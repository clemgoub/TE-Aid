"""Derived quantities for the panels: pileup, full-length calls, self-blast."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .records import Annotation, Copy


def coverage(annotation: Annotation, consensus_length: int) -> np.ndarray:
    """Per-position count of copies overlapping each consensus base.

    Returns a length-``consensus_length`` integer array; index *i* is the number
    of copies covering 0-based consensus position *i*.

    v1 built a dense ``n_copies x consensus_length`` matrix of zeros and summed
    the columns, which is O(n x L) memory — roughly 2 GB for 50k copies on a 5 kb
    consensus, and the largest family in a fish genome reaches that easily. A
    difference array gives the identical result in O(n + L).
    """
    if consensus_length <= 0:
        raise ValueError(f"consensus_length must be positive, got {consensus_length}")

    diff = np.zeros(consensus_length + 1, dtype=np.int64)
    for copy in annotation:
        if copy.consensus_start is None or copy.consensus_end is None:
            continue
        start = max(0, copy.consensus_start)
        end = min(consensus_length, copy.consensus_end)
        if end > start:
            diff[start] += 1
            diff[end] -= 1
    return np.cumsum(diff[:-1])


def full_length(
    annotation: Annotation, consensus_length: int, threshold: float = 0.9
) -> list[Copy]:
    """Copies spanning at least ``threshold`` of the consensus.

    Measured on the consensus interval, so a copy interrupted by a genomic
    insertion is judged on how much of the *family* it represents rather than
    how much genome it occupies.

    v1 computed the span as ``abs(qend - qstart)`` on 1-based closed
    coordinates, one short of the true length; a copy covering a consensus
    exactly end to end scored ``L-1`` and could miss a threshold of 1.0. Working
    in half-open coordinates removes the off-by-one rather than reproducing it.
    """
    cutoff = threshold * consensus_length
    return [
        c
        for c in annotation
        if c.consensus_start is not None
        and c.consensus_end is not None
        and (c.consensus_end - c.consensus_start) >= cutoff
    ]


@dataclass(slots=True)
class SelfHit:
    """One HSP of the consensus against itself."""

    q_start: int  # 0-based half-open, consensus coordinates
    q_end: int
    s_start: int
    s_end: int
    identity: float
    evalue: float
    bitscore: float

    @property
    def is_reverse(self) -> bool:
        return self.s_end < self.s_start

    @property
    def is_trivial_diagonal(self) -> bool:
        """The full-length self match every self-comparison necessarily finds."""
        return (
            not self.is_reverse
            and self.q_start == self.s_start
            and self.q_end == self.s_end
        )


class BlastNotFound(RuntimeError):
    pass


def self_blast(
    sequence: str,
    name: str = "consensus",
    *,
    word_size: int = 7,
    evalue: float = 1e-3,
    dust: bool = False,
) -> list[SelfHit]:
    """blastn of a consensus against itself, for the dot-plot and terminal repeats.

    Defaults follow the v2 decision in ``docs/BRIEF_v2.md`` §5.2:
    ``-word_size 7 -dust no -evalue 1e-3``, which is more sensitive than v1's
    ``-word_size 11`` with dust left on. Dust matters: it masks low-complexity
    sequence, and TE termini are frequently AT-rich, so leaving it enabled hides
    exactly the terminal repeats this panel exists to show.

    ``-task blastn`` is mandatory here — the default megablast task rejects word
    sizes below 16, so omitting it turns a sensitive search into an error.
    """
    if shutil.which("blastn") is None:
        raise BlastNotFound(
            "blastn not found on PATH; install NCBI BLAST+ to draw the dot-plot"
        )

    with tempfile.TemporaryDirectory(prefix="teaid-selfblast-") as tmp:
        fasta = Path(tmp) / "query.fa"
        fasta.write_text(f">{name}\n{sequence}\n")

        command = [
            "blastn",
            "-task", "blastn",
            "-query", str(fasta),
            "-subject", str(fasta),  # -subject avoids building a blastdb
            "-word_size", str(word_size),
            "-evalue", str(evalue),
            "-dust", "yes" if dust else "no",
            "-outfmt", "6 qstart qend sstart send pident evalue bitscore",
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=False)

    if result.returncode != 0:
        raise RuntimeError(f"blastn failed: {result.stderr.strip()}")

    hits: list[SelfHit] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        q_start, q_end, s_start, s_end, pident, ev, bits = line.split("\t")
        hits.append(
            SelfHit(
                # blast reports 1-based closed on both axes; the subject side
                # descends on reverse hits, so only the start is shifted and the
                # descending orientation is preserved for drawing.
                q_start=int(q_start) - 1,
                q_end=int(q_end),
                s_start=int(s_start) - 1,
                s_end=int(s_end),
                identity=float(pident),
                evalue=float(ev),
                bitscore=float(bits),
            )
        )
    return hits


def terminal_repeats(
    hits: list[SelfHit], consensus_length: int, *, min_span: float = 0.5
) -> dict[str, list[SelfHit]]:
    """Split self-hits into LTR- and TIR-like candidates.

    A terminal repeat is a pair of copies of the same sequence at opposite ends
    of the element: same-strand for an LTR, opposite-strand for a TIR. The test
    used here is that the two arms do not overlap and that the region they
    bracket covers at least ``min_span`` of the consensus — if a repeated pair
    brackets most of the element, its arms are near the termini by construction.

    Anchoring instead on "the first arm starts within x% of position 0" was
    tried and is too strict on real data: RepeatModeler consensuses frequently
    carry extra sequence beyond the element's true 5' boundary, which pushes a
    genuine LTR pair inward and hides it.

    blastn reports every off-diagonal repeat twice, once per direction, so the
    reciprocal duplicate is collapsed here; a two-LTR element yields one
    candidate, not two.

    These are suggestions for a curator to judge, never assertions: an internal
    segmental duplication produces the same signature as an LTR.
    """
    minimum = min_span * consensus_length
    candidates: list[SelfHit] = []

    for hit in hits:
        if hit.is_trivial_diagonal:
            continue
        s_lo, s_hi = sorted((hit.s_start, hit.s_end))
        # Overlapping arms are one region matching itself, not a repeat pair.
        if min(hit.q_end, s_hi) > max(hit.q_start, s_lo):
            continue
        span = max(hit.q_end, s_hi) - min(hit.q_start, s_lo)
        if span < minimum:
            continue
        candidates.append(hit)

    # blastn reports one terminal repeat many times over: once per direction,
    # and again for each slightly different extension of the same alignment.
    # Left alone that turns a single pair of LTRs into a dozen near-identical
    # candidates and a dozen lanes in the structure panel. Keep the
    # best-scoring representative of each group of mutually overlapping
    # candidates, comparing arm to arm.
    ltr: list[SelfHit] = []
    tir: list[SelfHit] = []
    for hit in sorted(candidates, key=lambda h: h.bitscore, reverse=True):
        kept = tir if hit.is_reverse else ltr
        if any(_same_repeat(hit, other) for other in kept):
            continue
        kept.append(hit)

    ltr.sort(key=lambda h: min(h.q_start, h.s_start, h.s_end))
    tir.sort(key=lambda h: min(h.q_start, h.s_start, h.s_end))
    return {"LTR": ltr, "TIR": tir}


def _arms(hit: SelfHit) -> tuple[tuple[int, int], tuple[int, int]]:
    """The hit's two arms as ascending intervals, ordered left to right."""
    query = (hit.q_start, hit.q_end)
    subject = tuple(sorted((hit.s_start, hit.s_end)))
    return tuple(sorted((query, subject)))  # type: ignore[return-value]


def _same_repeat(a: SelfHit, b: SelfHit) -> bool:
    """True when two hits describe the same terminal repeat.

    Compares left arm to left arm and right arm to right arm, so two genuinely
    different repeats that happen to share one endpoint are not merged.
    """
    a_left, a_right = _arms(a)
    b_left, b_right = _arms(b)
    return _overlaps(a_left, b_left) and _overlaps(a_right, b_right)


def _overlaps(a: tuple[int, int], b: tuple[int, int]) -> bool:
    return min(a[1], b[1]) > max(a[0], b[0])
