"""The v1 input route: copies found by ``blastn`` of the consensus against a genome.

v1 took a consensus and a genome and nothing else, and this module keeps that
working (``-q/--query`` with ``-g/--genome``, or ``--blastn GENOME``). It is the
leakiest route v2 offers, and the sheet says so. Megablast misses diverged
copies. Each HSP is one copy, because a search has no ID column to join the
pieces of an interrupted insertion. An annotation the user already trusts gives
a better sheet from the same genome.

What it copies from v1, deliberately:

- the database is built beside the genome (``makeblastdb -in G -out G``), so an
  existing v1 database is reused as it is;
- ``blastn`` runs with its default task (megablast) at v1's E-value, 10e-8.

What it does not copy: v1 plotted ``100 - pident`` as "divergence to consensus".
That is an uncorrected identity from one local alignment, not RepeatMasker's
Kimura-corrected divergence, so it is stored as
:attr:`DivergenceKind.BLASTN_IDENTITY` and panel 1's axis names it for what it is.

Nor is ``-r/--remove-redundant`` v1's rule; see :func:`remove_redundant`.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

from .records import Annotation, Copy, DivergenceKind

# v1's default, written the way v1 wrote it (TE-Aid shell script, EVALUE).
DEFAULT_EVALUE = 10e-8

_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# The most genome -r bridges when it joins two pieces. RepeatMasker's own joins
# in the GenomeArk goby bridge a median of 144 bp, and 99% of them bridge under
# 4,863 bp.
JOIN_MAX_GAP = 5000
# Two hits on one strand compete for the same stretch of genome, and two pieces
# are too alike along the consensus to be parts of one copy, when they share at
# least this fraction of the shorter one: the rule homology.best_per_region uses.
OVERLAP = 0.5
_BIN = 10_000


class BlastnNotFound(RuntimeError):
    pass


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise BlastnNotFound(
            f"{tool} not found on PATH; the genome route needs NCBI BLAST+ "
            "(https://blast.ncbi.nlm.nih.gov/)"
        )
    return path


def _is_database(prefix: Path) -> bool:
    if shutil.which("blastdbcmd") is None:
        return any(prefix.with_name(prefix.name + ext).exists() for ext in (".nin", ".nal"))
    result = subprocess.run(["blastdbcmd", "-db", str(prefix), "-info"],
                            capture_output=True, text=True, check=False)
    return result.returncode == 0


def ensure_database(genome: Path, *, log=print) -> Path:
    """The blastn database for ``genome``, built on first use.

    Built beside the genome, as v1 did, so a v1 database is reused. When that
    directory is not writable it goes to the cache instead, keyed on the
    genome's path, size and modification time.
    """
    from .proteins import cache_root

    if _is_database(genome):
        return genome
    _require("makeblastdb")
    stat = genome.stat()
    key = hashlib.sha256(f"{genome.resolve()}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    fallback = cache_root() / "genomes" / key.hexdigest()[:16] / genome.name
    if _is_database(fallback):
        return fallback
    reason = "no writable place for it"
    for prefix in (genome, fallback):
        try:
            prefix.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            reason = str(exc)
            continue
        log(f"building a blastn database for {genome.name} at {prefix.parent} (once per genome)")
        # No -parse_seqids, as in v1: it rejects some headers v1 accepted.
        result = subprocess.run(
            ["makeblastdb", "-in", str(genome), "-out", str(prefix), "-dbtype", "nucl"],
            capture_output=True, text=True, check=False,
        )
        if result.returncode == 0 and _is_database(prefix):
            return prefix
        reason = result.stderr.strip()[:300] or f"exit {result.returncode}"
    raise RuntimeError(f"makeblastdb could not build a database for {genome}: {reason}")


def search(sequence: str, name: str, database: Path, *, evalue: float = DEFAULT_EVALUE,
           family: str | None = None, class_label: str | None = None) -> Annotation:
    """Genomic copies of one consensus, one per HSP, as an :class:`Annotation`."""
    blastn = _require("blastn")
    with tempfile.TemporaryDirectory(prefix="teaid-genome-") as tmp:
        query = Path(tmp) / "consensus.fa"
        query.write_text(f">{name}\n{sequence}\n")
        result = subprocess.run(
            [blastn, "-query", str(query), "-db", str(database), "-evalue", str(evalue),
             "-outfmt", _OUTFMT],
            capture_output=True, text=True, check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"blastn failed: {result.stderr.strip()[:300]}")
    annotation = parse(result.stdout, family=family or name, class_label=class_label,
                       consensus_length=len(sequence))
    annotation.source_path = str(database)
    return annotation


def parse(text: str, *, family: str, class_label: str | None,
          consensus_length: int) -> Annotation:
    """``-outfmt 6`` rows as copies. Coordinates become 0-based half-open."""
    annotation = Annotation(source_format="blastn")
    for lineno, line in enumerate(text.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) < 12:
            if line.strip():
                annotation.skipped.append((lineno, f"expected 12 fields, got {len(fields)}"))
            continue
        try:
            pident = float(fields[2])
            q_start, q_end = int(fields[6]), int(fields[7])
            s_start, s_end = int(fields[8]), int(fields[9])
            bitscore = float(fields[11])
        except ValueError as exc:
            annotation.skipped.append((lineno, str(exc)))
            continue
        strand = "+" if s_start <= s_end else "-"
        annotation.copies.append(Copy(
            chrom=fields[1],
            start=min(s_start, s_end) - 1,
            end=max(s_start, s_end),
            strand=strand,
            family=family,
            class_label=class_label,
            divergence=100.0 - pident,
            divergence_kind=DivergenceKind.BLASTN_IDENTITY,
            consensus_start=min(q_start, q_end) - 1,
            consensus_end=max(q_start, q_end),
            consensus_left=max(0, consensus_length - max(q_start, q_end)),
            source_line=lineno,
            score=bitscore,
        ))
    return annotation


@dataclass(frozen=True)
class Reduction:
    """What ``-r`` did to one family's hits, for the sheet's note."""

    hits: int
    placed: int
    copies: int
    ties: int


def _shared(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """Overlap of two intervals as a fraction of the shorter one."""
    shorter = min(a_end - a_start, b_end - b_start)
    if shorter <= 0:
        return 0.0
    return max(0, min(a_end, b_end) - max(a_start, b_start)) / shorter


def _follows(prev: Copy, nxt: Copy, strand: str) -> bool:
    """Whether ``nxt``, later in the genome, continues ``prev`` along the consensus."""
    if None in (prev.consensus_start, prev.consensus_end, nxt.consensus_start, nxt.consensus_end):
        return False
    if strand == "-":
        # On the minus strand the piece later in the genome sits earlier in the
        # consensus.
        prev, nxt = nxt, prev
    if nxt.consensus_start <= prev.consensus_start or nxt.consensus_end <= prev.consensus_end:
        return False
    return _shared(prev.consensus_start, prev.consensus_end,
                   nxt.consensus_start, nxt.consensus_end) < OVERLAP


def _bins(copy: Copy) -> range:
    return range(copy.start // _BIN, max(copy.end - 1, copy.start) // _BIN + 1)


def remove_redundant(annotation: Annotation, *, max_gap: int = JOIN_MAX_GAP
                     ) -> tuple[Annotation, Reduction]:
    """``-r``: one placement per stretch of genome, and interrupted copies joined.

    1. Of the hits on one strand that share at least half of the shorter one's
       locus, keep the highest-scoring. On a tie keep the lower consensus start,
       and count a tie when the hit it beat sat elsewhere on the consensus: a
       solo LTR matches both LTRs of the consensus equally well.
    2. Join pieces that follow each other along the consensus, on the same
       sequence and strand, with at most ``max_gap`` bp of genome between them
       and less than half of the shorter piece's consensus in common. Joined
       pieces share a ``hit_id``, so they become one insertion with its joins
       classified, as RepeatMasker's split rows do.

    This replaces v1's ``reduce.cpp`` rule, which dropped a hit only when
    another hit's locus contained it. Under that rule, which of two identical
    loci survived depended on blastn's output order, and the pieces of an
    interrupted copy were never joined. On a synthetic LTR element with three
    real insertions, blastn gives 9 hits, v1's rule 5 and this one 3.
    """
    order = sorted(annotation.copies,
                   key=lambda c: (-(c.score or 0.0), c.consensus_start or 0, c.chrom, c.start))
    occupied: dict[tuple[str, str, int], list[Copy]] = {}
    placed: list[Copy] = []
    tied: set[int] = set()
    for copy in order:
        rival = next(
            (other for b in _bins(copy)
             for other in occupied.get((copy.chrom, copy.strand, b), ())
             if _shared(other.start, other.end, copy.start, copy.end) >= OVERLAP),
            None,
        )
        if rival is None:
            placed.append(copy)
            for b in _bins(copy):
                occupied.setdefault((copy.chrom, copy.strand, b), []).append(copy)
        elif rival.score == copy.score and rival.consensus_start != copy.consensus_start:
            tied.add(id(rival))

    by_strand: dict[tuple[str, str], list[Copy]] = {}
    for copy in placed:
        by_strand.setdefault((copy.chrom, copy.strand), []).append(copy)
    kept: list[Copy] = []
    n_chains = n_copies = 0
    for (_, strand), pieces in by_strand.items():
        pieces.sort(key=lambda c: (c.start, c.end))
        chains: list[list[Copy]] = []
        open_chains: list[list[Copy]] = []
        for piece in pieces:
            still_open = []
            for chain in open_chains:
                (still_open if piece.start - chain[-1].end <= max_gap else chains).append(chain)
            open_chains = still_open
            fits = [c for c in open_chains if _follows(c[-1], piece, strand)]
            if fits:
                max(fits, key=lambda c: c[-1].end).append(piece)
            else:
                open_chains.append([piece])
        chains += open_chains
        n_copies += len(chains)
        for chain in chains:
            if len(chain) == 1:
                kept.append(chain[0])
                continue
            n_chains += 1
            kept += [replace(c, hit_id=f"r{n_chains}") for c in chain]
    kept.sort(key=lambda c: (c.chrom, c.start))
    reduced = Annotation(copies=kept, source_path=annotation.source_path,
                         source_format=annotation.source_format, skipped=annotation.skipped)
    return reduced, Reduction(hits=len(annotation.copies), placed=len(placed),
                              copies=n_copies, ties=len(tied))
