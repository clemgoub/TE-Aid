"""Reader for Stockholm seed alignments (Dfam format).

A seed carries everything the default sheet needs and more: the copies, their
genomic loci, their alignment to the family consensus, and the consensus itself
on the ``#=GC RF`` line. Nothing is re-derived by search.

Format notes that matter here:

- One file may hold many records, each ``# STOCKHOLM 1.0`` … ``//``. Pick one by
  ``#=GF ID``.
- Alignments may be *interleaved* across several blocks; rows are accumulated by
  sequence name rather than assumed to arrive in one piece.
- Dfam uses ``.`` as the gap character; ``-`` is accepted too, since other
  producers emit it and the distinction carries no meaning for us.
- ``#=GF TP`` is Dfam's full classification path, semicolon separated
  (``Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;…``),
  not a short RepeatMasker-style ``LTR/Gypsy`` label.

Sequence identifiers are Smitten format, in either of two shapes::

    GCA_951799975.1:OX637595.1:15848-16090_+     assembly, sequence, span, strand
    OY720097.1:14692470-14693460_+               sequence, span, strand

Both appear in practice — RepeatModeler writes the short form. Coordinates are
**1-based fully closed** and are converted to the package's 0-based half-open
convention on read, as every other reader does.
"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

from ..records import Annotation, Copy, DivergenceKind

GAPS = frozenset(".-~")

# ...:<start>-<end>_<strand>, with an optional leading assembly accession.
_SMITTEN = re.compile(
    r"^(?:(?P<assembly>[^:]+):)?(?P<sequence>[^:]+):(?P<start>\d+)-(?P<end>\d+)_(?P<strand>[+-])$"
)


class MalformedSeedError(ValueError):
    """The file is not usable Stockholm.

    Raised rather than silently skipped: a caller processing many families needs
    to record *why* one failed, so the reason has to reach stderr intact.
    """


@dataclass(slots=True)
class SeedSequence:
    name: str
    aligned: str
    chrom: str | None = None
    start: int | None = None  # 0-based half-open
    end: int | None = None
    strand: str = "."
    assembly: str | None = None

    @property
    def has_locus(self) -> bool:
        return self.chrom is not None


@dataclass(slots=True)
class Seed:
    """One Stockholm record."""

    identifier: str | None = None
    features: dict[str, str] = field(default_factory=dict)  # #=GF
    columns: dict[str, str] = field(default_factory=dict)  # #=GC
    sequences: list[SeedSequence] = field(default_factory=list)
    source_path: str | None = None

    @property
    def expected_class(self) -> str | None:
        """``#=GF TP`` — the externally supplied classification, if present.

        This is the only claim on the sheet TE-Aid did not derive itself, and it
        is displayed as a claim to be checked, never as a conclusion.
        """
        return self.features.get("TP")

    @property
    def description(self) -> str | None:
        return self.features.get("DE")

    @property
    def reference_line(self) -> str | None:
        return self.columns.get("RF")

    @property
    def declared_count(self) -> int | None:
        raw = self.features.get("SQ")
        try:
            return int(raw) if raw is not None else None
        except ValueError:
            return None

    def __len__(self) -> int:
        return len(self.sequences)

    @property
    def alignment_width(self) -> int:
        return max((len(s.aligned) for s in self.sequences), default=0)

    def consensus(self) -> str:
        """The consensus, taken from ``#=GC RF`` with gap columns removed.

        Falls back to a per-column majority over the sequences when the seed has
        no RF line, so a seed from a producer that omits it still yields a sheet.
        """
        rf = self.reference_line
        if rf:
            return "".join(ch for ch in rf if ch not in GAPS).upper()
        return self._majority_consensus()

    def _majority_consensus(self) -> str:
        import collections

        out = []
        for index in range(self.alignment_width):
            column = [
                s.aligned[index].upper()
                for s in self.sequences
                if index < len(s.aligned) and s.aligned[index] not in GAPS
            ]
            if len(column) * 2 > len(self.sequences):  # occupied in a majority
                out.append(collections.Counter(column).most_common(1)[0][0])
        return "".join(out)

    def consensus_columns(self) -> list[int]:
        """Alignment column indices that correspond to a consensus position.

        Insert columns — where the reference has a gap — are not consensus
        positions; a copy with an insertion must not shift every coordinate
        after it.
        """
        rf = self.reference_line
        if rf:
            return [i for i, ch in enumerate(rf) if ch not in GAPS]
        width = self.alignment_width
        occupied = []
        for index in range(width):
            filled = sum(
                1
                for s in self.sequences
                if index < len(s.aligned) and s.aligned[index] not in GAPS
            )
            if filled * 2 > len(self.sequences):
                occupied.append(index)
        return occupied

    def depth(self) -> list[int]:
        """Per-consensus-position count of sequences aligned there."""
        return [
            sum(
                1
                for s in self.sequences
                if column < len(s.aligned) and s.aligned[column] not in GAPS
            )
            for column in self.consensus_columns()
        ]

    def mismatches(self) -> list[int]:
        """Per-consensus-position count of sequences differing from the reference.

        Depth alone can flatter a seed: a column with twelve sequences that
        disagree with each other is not twelve sequences of support. Dfam's own
        browser colours its coverage bars by allele fraction for exactly this
        reason, and this is the quantity behind that.

        Positions where the reference is ambiguous (``N``) are not counted as
        disagreement, since there is nothing definite to disagree with.
        """
        reference = self.reference_line
        columns = self.consensus_columns()
        if not reference:
            return [0] * len(columns)

        counts = []
        for column in columns:
            base = reference[column].upper()
            if base in GAPS or base == "N":
                counts.append(0)
                continue
            counts.append(
                sum(
                    1
                    for s in self.sequences
                    if column < len(s.aligned)
                    and s.aligned[column] not in GAPS
                    and s.aligned[column].upper() != base
                )
            )
        return counts

    def aligned_blocks(self) -> list[tuple[str, list[tuple[int, int]]]]:
        """Each sequence's aligned stretches, in consensus coordinates.

        One entry per sequence: its name, and the half-open runs of consensus
        positions where it actually has a base. Internal deletions appear as the
        gaps *between* runs, which is what makes a pileup show where a copy is
        truncated and where it is interrupted — two very different things that a
        single start-to-end bar conflates.

        Sorted by start position, then by span, so a pileup reads top-left to
        bottom-right the way a genome browser lays reads out.
        """
        columns = self.consensus_columns()
        out: list[tuple[str, list[tuple[int, int]]]] = []

        for sequence in self.sequences:
            runs: list[tuple[int, int]] = []
            start: int | None = None
            for position, column in enumerate(columns):
                filled = (
                    column < len(sequence.aligned)
                    and sequence.aligned[column] not in GAPS
                )
                if filled and start is None:
                    start = position
                elif not filled and start is not None:
                    runs.append((start, position))
                    start = None
            if start is not None:
                runs.append((start, len(columns)))
            if runs:
                out.append((sequence.name, runs))

        out.sort(key=lambda item: (item[1][0][0], -(item[1][-1][1] - item[1][0][0])))
        return out

    def to_annotation(self) -> Annotation:
        """Copies, consensus coordinates and divergence, all read from the seed."""
        columns = self.consensus_columns()
        index_of = {column: position for position, column in enumerate(columns)}
        annotation = Annotation(source_path=self.source_path, source_format="stockholm")
        reference = self.reference_line

        for line_number, sequence in enumerate(self.sequences, start=1):
            occupied = [
                c for c in columns if c < len(sequence.aligned) and sequence.aligned[c] not in GAPS
            ]
            if not occupied:
                continue
            consensus_start = index_of[occupied[0]]
            consensus_end = index_of[occupied[-1]] + 1

            divergence = None
            kind = DivergenceKind.NONE
            if reference:
                compared = mismatched = 0
                for column in occupied:
                    reference_base = reference[column].upper()
                    if reference_base in GAPS or reference_base == "N":
                        continue
                    compared += 1
                    if sequence.aligned[column].upper() != reference_base:
                        mismatched += 1
                if compared:
                    # Divergence from the seed's own consensus: the same
                    # quantity RepeatMasker's perc_div reports, computed here
                    # from the alignment rather than taken on trust.
                    divergence = 100.0 * mismatched / compared
                    kind = DivergenceKind.CONSENSUS

            annotation.copies.append(
                Copy(
                    chrom=sequence.chrom if sequence.has_locus else sequence.name,
                    start=sequence.start if sequence.start is not None else 0,
                    end=sequence.end
                    if sequence.end is not None
                    else len(occupied),
                    strand=sequence.strand,
                    family=self.identifier or "seed",
                    class_label=self.expected_class,
                    divergence=divergence,
                    divergence_kind=kind,
                    consensus_start=consensus_start,
                    consensus_end=consensus_end,
                    consensus_left=len(columns) - consensus_end,
                    source_line=line_number,
                )
            )
        return annotation


def parse_smitten(name: str) -> tuple[str | None, str, int, int, str] | None:
    """Split a Smitten identifier into (assembly, sequence, start, end, strand).

    Returns 0-based half-open coordinates, or None when the name is not Smitten
    — a seed may legitimately carry plain sequence names, and that is not an
    error, just an absence of loci.
    """
    match = _SMITTEN.match(name.strip())
    if match is None:
        return None
    start, end = int(match.group("start")), int(match.group("end"))
    if start > end:
        start, end = end, start
    return (
        match.group("assembly"),
        match.group("sequence"),
        start - 1,
        end,
        match.group("strand"),
    )


def _open(path: Path) -> io.TextIOBase:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def read_all(path: str | Path) -> list[Seed]:
    """Every record in a Stockholm file, in file order."""
    path = Path(path)
    seeds: list[Seed] = []
    current: Seed | None = None
    rows: dict[str, list[str]] = {}
    order: list[str] = []
    saw_header = False

    def close() -> None:
        nonlocal current, rows, order
        if current is None:
            return
        current.sequences = [
            _make_sequence(name, "".join(rows[name])) for name in order
        ]
        seeds.append(current)
        current, rows, order = None, {}, []

    with _open(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n").rstrip("\r")
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("# STOCKHOLM"):
                close()
                saw_header = True
                current = Seed(source_path=str(path))
                continue

            if current is None:
                # Tolerate leading junk before the first header, but a file with
                # no header at all is not Stockholm.
                continue

            if stripped == "//":
                close()
                continue

            if stripped.startswith("#=GF"):
                parts = stripped.split(None, 2)
                if len(parts) >= 3:
                    tag, value = parts[1], parts[2].strip()
                    # A repeated tag is continuation, not replacement.
                    current.features[tag] = (
                        f"{current.features[tag]} {value}" if tag in current.features else value
                    )
                    if tag == "ID":
                        current.identifier = current.features["ID"]
                continue

            if stripped.startswith("#=GC"):
                parts = stripped.split(None, 2)
                if len(parts) >= 3:
                    tag, value = parts[1], parts[2].strip()
                    current.columns[tag] = current.columns.get(tag, "") + value
                continue

            if stripped.startswith("#"):
                continue  # #=GS, #=GR and free comments are not used here

            parts = stripped.split(None, 1)
            if len(parts) != 2:
                continue
            name, chunk = parts[0], parts[1].strip()
            if name not in rows:
                rows[name] = []
                order.append(name)
            rows[name].append(chunk)

    close()

    if not saw_header:
        raise MalformedSeedError(f"{path}: no '# STOCKHOLM' header found")
    if not seeds:
        raise MalformedSeedError(f"{path}: no records found")
    return seeds


def _make_sequence(name: str, aligned: str) -> SeedSequence:
    parsed = parse_smitten(name)
    if parsed is None:
        return SeedSequence(name=name, aligned=aligned)
    assembly, sequence, start, end, strand = parsed
    return SeedSequence(
        name=name,
        aligned=aligned,
        chrom=sequence,
        start=start,
        end=end,
        strand=strand,
        assembly=assembly,
    )


def read(path: str | Path, family: str | None = None) -> Seed:
    """One record: the named family, or the only record when there is one."""
    seeds = read_all(path)
    if family is None:
        if len(seeds) != 1:
            raise MalformedSeedError(
                f"{path} holds {len(seeds)} records; name one with --family"
            )
        return seeds[0]

    wanted = family.split("#", 1)[0].strip().casefold()
    for seed in seeds:
        if seed.identifier and seed.identifier.split("#", 1)[0].strip().casefold() == wanted:
            return seed
    raise MalformedSeedError(f"{path}: no record with #=GF ID {family!r}")


def identifiers(path: str | Path) -> list[str]:
    """Every ``#=GF ID`` in the file, for error messages and listings."""
    return [s.identifier for s in read_all(path) if s.identifier]
