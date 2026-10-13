"""The shared record struct every annotation dialect is read into.

Coordinate convention
---------------------
Every coordinate stored on a :class:`Copy` is **0-based half-open**, genomic and
consensus alike, regardless of the dialect it was read from. RepeatMasker
``.out``/``.gff3`` (1-based fully closed) and Smitten identifiers in Stockholm
seeds (1-based fully closed) are converted at read time; BED16 is already
0-based half-open and passes through.

This is deliberate: mixing conventions silently is the single easiest way to
produce off-by-one errors that survive review because the plots still look
plausible. Convert once, at the edge. Use :meth:`Copy.genomic_1based` and
:meth:`Copy.consensus_1based` when reporting to a human or writing a format that
expects closed intervals.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterator


class DivergenceKind(str, Enum):
    """What a divergence value actually measures.

    ``perc_div`` means different things depending on which tool produced it, and
    the three meanings are not interchangeable on a plot axis:

    - :attr:`CONSENSUS` — divergence from a library consensus. A proxy for the
      age of the insertion; this is the one panel 1 is designed around.
    - :attr:`ARRAY_HOMOGENEITY` — how uniform a tandem array is. Says nothing
      about age; a young satellite can be wildly heterogeneous. No reader infers
      it; a BED16 caller sets it with ``--divergence-kind``.
    - :attr:`BLASTN_IDENTITY` — ``100 - pident`` of one blastn HSP, from the
      genome route. An identity over one local alignment whose length counts
      gaps, where RepeatMasker reports substitutions against the consensus and
      puts indels in separate columns, so the same copy can read older or
      younger than in a RepeatMasker annotation. v1 labelled it "divergence".
    - :attr:`NONE` — the producer emits no divergence at all.
    """

    CONSENSUS = "consensus"
    ARRAY_HOMOGENEITY = "array_homogeneity"
    BLASTN_IDENTITY = "blastn_identity"
    NONE = "none"

    @property
    def axis_label(self) -> str:
        return {
            DivergenceKind.CONSENSUS: "divergence from consensus (%)",
            DivergenceKind.ARRAY_HOMOGENEITY: "array heterogeneity (%)",
            DivergenceKind.BLASTN_IDENTITY: "100 − blastn identity (%)",
            DivergenceKind.NONE: "divergence not reported",
        }[self]


@dataclass(slots=True)
class Copy:
    """One annotated genomic copy of one TE family."""

    chrom: str
    start: int  # 0-based half-open, genomic
    end: int
    # '+', '-', or '.' for strandless features. RepeatMasker 'C' is normalised
    # to '-' at read time. '.' is not hypothetical: tandem and satellite callers
    # legitimately report no orientation.
    strand: str
    family: str
    class_label: str | None = None
    divergence: float | None = None
    divergence_kind: DivergenceKind = DivergenceKind.NONE
    # Position of this copy within the family consensus, 0-based half-open.
    consensus_start: int | None = None
    consensus_end: int | None = None
    # Bases of consensus remaining past consensus_end; lets us infer consensus
    # length from any single copy without reading the FASTA.
    consensus_left: int | None = None
    # Set on RepeatMasker rows flagged '*': a lower-scoring hit overlapping a
    # better one. Kept rather than dropped so callers can decide.
    is_overlapping_lower_score: bool = False
    # Groups fragments of one interrupted insertion (BED16 column 16; the
    # RepeatMasker 'ID' column). Copies sharing a hit_id are pieces of a single
    # element split by an indel or a nested insertion, not independent copies —
    # counting them separately is exactly the fragment inflation that makes
    # search-driven copy discovery misleading.
    hit_id: str | None = None
    source_line: int | None = None
    # Alignment score, higher is better. Only the genome route sets it (blastn
    # bitscore), for -r to choose between hits on the same stretch of genome.
    score: float | None = None

    def __post_init__(self) -> None:
        if self.end < self.start:
            raise ValueError(
                f"genomic end {self.end} precedes start {self.start} "
                f"({self.chrom}, {self.family})"
            )
        if (
            self.consensus_start is not None
            and self.consensus_end is not None
            and self.consensus_end < self.consensus_start
        ):
            raise ValueError(
                f"consensus end {self.consensus_end} precedes start "
                f"{self.consensus_start} ({self.family})"
            )

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def consensus_length_estimate(self) -> int | None:
        """Total consensus length implied by this copy, if it can be inferred."""
        if self.consensus_end is None or self.consensus_left is None:
            return None
        return self.consensus_end + self.consensus_left

    @property
    def has_divergence(self) -> bool:
        """True only for a real value of a kind that is meaningful to plot.

        Absent divergence is never zero. A record with no divergence is left out
        of panel 1's y-axis entirely rather than being drawn on the zero line,
        where it would masquerade as a pristine, very recent insertion.
        """
        return self.divergence is not None and self.divergence_kind is not DivergenceKind.NONE

    def genomic_1based(self) -> tuple[int, int]:
        """(start, end) as 1-based fully closed, for display and report output."""
        return self.start + 1, self.end

    def consensus_1based(self) -> tuple[int, int] | None:
        if self.consensus_start is None or self.consensus_end is None:
            return None
        return self.consensus_start + 1, self.consensus_end


class JunctionKind(str, Enum):
    """What lies between two fragments of one insertion, relative to consensus.

    - :attr:`INSERTION` — the genome carries sequence the consensus does not,
      most often another element nested in this copy. It uses no consensus,
      so on a consensus axis the two fragments abut.
    - :attr:`DELETION` — the copy lacks consensus sequence that its two
      neighbouring fragments bracket. On a consensus axis it is a real gap.
    - :attr:`FLUSH` — the two differences cancel to within
      :data:`JUNCTION_TOLERANCE`, as when an alignment was split for its score.
    """

    INSERTION = "insertion"
    DELETION = "deletion"
    FLUSH = "flush"


# Base pairs by which the genomic and consensus gaps may differ before a join
# counts as an insertion or a deletion. RepeatMasker fragments routinely
# overlap or leave a few unaligned bases at a join, and those are not indels.
# On the GenomeArk goby (11,629 classifiable joins) this gives 68.9% insertions
# (median 243 bp), 20.2% deletions (median 218 bp) and 10.9% flush, with the
# same split on either strand.
JUNCTION_TOLERANCE = 50


@dataclass(slots=True, frozen=True)
class Junction:
    """The join between two fragments that are adjacent in the genome.

    ``consensus_at`` and ``consensus_to`` bound the consensus interval the join
    skips (a deletion) or the point where it sits (an insertion or a flush
    join, where the two are equal). ``size`` is the indel length in bp, 0 when
    flush. A join can skip some consensus and also carry extra genome; it is
    classified by the net difference, and both gaps are kept here.
    """

    consensus_at: int
    consensus_to: int
    kind: JunctionKind
    size: int
    genomic_gap: int
    consensus_gap: int


def _junction(prev: Copy, nxt: Copy, strand: str) -> Junction | None:
    """Classify the join between two genomically adjacent fragments.

    On the minus strand the element runs backwards along the genome, so the
    fragment later in the genome sits *earlier* in the consensus. Ignoring that
    turns every minus-strand join into a spurious insertion.
    """
    if None in (prev.consensus_start, prev.consensus_end,
                nxt.consensus_start, nxt.consensus_end) or strand not in "+-":
        return None
    genomic_gap = nxt.start - prev.end
    if strand == "+":
        at, to = prev.consensus_end, nxt.consensus_start
    else:
        at, to = nxt.consensus_end, prev.consensus_start
    consensus_gap = to - at
    delta = consensus_gap - genomic_gap
    if delta > JUNCTION_TOLERANCE:
        kind, size = JunctionKind.DELETION, delta
    elif delta < -JUNCTION_TOLERANCE:
        kind, size = JunctionKind.INSERTION, -delta
    else:
        kind, size = JunctionKind.FLUSH, 0
    if kind is not JunctionKind.DELETION or to < at:
        to = at  # a point on the consensus axis, not an interval
    return Junction(at, to, kind, size, genomic_gap, consensus_gap)


def _merge(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(blocks):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


@dataclass(slots=True)
class Insertion:
    """One genomic insertion: the annotation rows RepeatMasker split it into.

    Panels 1 and 2 count insertions, not rows. A copy interrupted by a nested
    element arrives as two rows sharing one ID, and counting both is the copy
    inflation v2 exists to fix. An ungrouped row is an insertion of one
    fragment, so every consumer handles one type.
    """

    fragments: list[Copy]  # sorted by genomic start

    @property
    def chrom(self) -> str:
        return self.fragments[0].chrom

    @property
    def family(self) -> str:
        return self.fragments[0].family

    @property
    def class_label(self) -> str | None:
        return self.fragments[0].class_label

    @property
    def hit_id(self) -> str | None:
        return self.fragments[0].hit_id

    @property
    def start(self) -> int:
        return min(f.start for f in self.fragments)

    @property
    def end(self) -> int:
        return max(f.end for f in self.fragments)

    @property
    def strand(self) -> str:
        strands = {f.strand for f in self.fragments}
        return strands.pop() if len(strands) == 1 else "."

    @property
    def n_fragments(self) -> int:
        return len(self.fragments)

    @property
    def is_fragmented(self) -> bool:
        return len(self.fragments) > 1

    def genomic_1based(self) -> tuple[int, int]:
        return self.start + 1, self.end

    @property
    def consensus_blocks(self) -> list[tuple[int, int]]:
        """Consensus intervals the fragments align to, merged, ascending."""
        return _merge([
            (f.consensus_start, f.consensus_end) for f in self.fragments
            if f.consensus_start is not None and f.consensus_end is not None
        ])

    @property
    def consensus_covered(self) -> int:
        """Consensus bases aligned by any fragment: the union, not the span.

        A 5' and a 3' fragment of a 5 kb element span nearly all of it while
        covering a few hundred bases, and full-length has to measure the
        second. Deletions inside a single alignment are not visible here,
        because a row gives only its two endpoints.
        """
        return sum(end - start for start, end in self.consensus_blocks)

    @property
    def consensus_start(self) -> int | None:
        blocks = self.consensus_blocks
        return blocks[0][0] if blocks else None

    @property
    def consensus_end(self) -> int | None:
        blocks = self.consensus_blocks
        return blocks[-1][1] if blocks else None

    @property
    def divergence_kind(self) -> DivergenceKind:
        kind, _ = self._divergence()
        return kind

    @property
    def divergence(self) -> float | None:
        _, value = self._divergence()
        return value

    @property
    def has_divergence(self) -> bool:
        kind, value = self._divergence()
        return value is not None and kind is not DivergenceKind.NONE

    def _divergence(self) -> tuple[DivergenceKind, float | None]:
        """Mean divergence, weighted by each fragment's consensus extent.

        This is RepeatMasker's own rule when it merges two annotations
        (ProcessRepeats, where the weight is SubjEnd - SubjStart + 1, the
        consensus length of each piece). Fragments with no divergence do not
        vote; if the fragments disagree on what divergence means, the kind
        carrying the most consensus wins and only its fragments are averaged.
        """
        weights: dict[DivergenceKind, list[tuple[float, int]]] = {}
        for f in self.fragments:
            if not f.has_divergence:
                continue
            extent = (
                f.consensus_end - f.consensus_start
                if f.consensus_start is not None and f.consensus_end is not None
                else f.length
            )
            weights.setdefault(f.divergence_kind, []).append((f.divergence, max(extent, 1)))
        if not weights:
            return self.fragments[0].divergence_kind, None
        kind = max(weights, key=lambda k: sum(w for _, w in weights[k]))
        pairs = weights[kind]
        total = sum(w for _, w in pairs)
        return kind, sum(d * w for d, w in pairs) / total

    @property
    def junctions(self) -> list[Junction]:
        """One per pair of genomically adjacent fragments, when classifiable.

        A pair is skipped when either fragment lacks consensus coordinates or
        the insertion's fragments disagree on strand.
        """
        strand = self.strand
        found = []
        for prev, nxt in zip(self.fragments, self.fragments[1:]):
            junction = _junction(prev, nxt, strand)
            if junction is not None:
                found.append(junction)
        return found


@dataclass(slots=True)
class Annotation:
    """All copies read from one annotation file, plus where they came from."""

    copies: list[Copy] = field(default_factory=list)
    source_path: str | None = None
    source_format: str | None = None
    # Rows the reader could not parse: (line number, reason). Surfaced in the
    # report rather than raised, so one malformed line does not lose a genome.
    skipped: list[tuple[int, str]] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.copies)

    def __iter__(self) -> Iterator[Copy]:
        return iter(self.copies)

    def families(self) -> list[str]:
        return sorted({c.family for c in self.copies})

    def for_family(self, name: str, *, exact: bool = True) -> "Annotation":
        """Subset to one family.

        With ``exact=False``, matches case-insensitively and ignores a
        ``#class/family`` suffix on either side, since a consensus FASTA header
        (``rnd-1_family-257#Unknown``) and an annotation row
        (``rnd-1_family-257``) routinely disagree about it.
        """
        if exact:
            picked = [c for c in self.copies if c.family == name]
        else:
            want = _bare_name(name)
            picked = [c for c in self.copies if _bare_name(c.family) == want]
        return Annotation(
            copies=picked,
            source_path=self.source_path,
            source_format=self.source_format,
            skipped=self.skipped,
        )

    def divergence_kinds(self) -> set[DivergenceKind]:
        return {c.divergence_kind for c in self.copies}

    def fragment_groups(self) -> list[list[Copy]]:
        """Copies grouped into insertions, joining fragments that share a hit_id.

        Copies without a hit_id are each their own group. Order within a group
        follows genomic position, so callers can measure the span of a
        reassembled element.
        """
        groups: dict[str, list[Copy]] = {}
        singles: list[list[Copy]] = []
        for copy in self.copies:
            if copy.hit_id is None:
                singles.append([copy])
            else:
                groups.setdefault(f"{copy.chrom}\t{copy.hit_id}", []).append(copy)
        joined = [sorted(g, key=lambda c: c.start) for g in groups.values()]
        return joined + singles

    def insertions(self, *, group: bool = True) -> list[Insertion]:
        """The copies as insertions, ordered by locus.

        With ``group=False`` every row is an insertion of one fragment, which
        reproduces row counting without changing the type consumers see.
        """
        groups = self.fragment_groups() if group else [[c] for c in self.copies]
        found = [Insertion(sorted(g, key=lambda c: c.start)) for g in groups]
        found.sort(key=lambda i: (i.chrom, i.start))
        return found

    def consensus_length(self) -> int | None:
        """Best estimate of consensus length from the copies themselves.

        Takes the maximum implied length: a copy spanning further into the
        consensus than its neighbours is more informative than a truncated one.
        Returns None when no copy carries consensus coordinates.
        """
        estimates = [
            e for e in (c.consensus_length_estimate for c in self.copies) if e is not None
        ]
        return max(estimates) if estimates else None


def _bare_name(name: str) -> str:
    return name.split("#", 1)[0].strip().casefold()
