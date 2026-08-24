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
      about age; a young satellite can be wildly heterogeneous.
    - :attr:`NONE` — the producer emits no divergence at all.
    """

    CONSENSUS = "consensus"
    ARRAY_HOMOGENEITY = "array_homogeneity"
    NONE = "none"

    @property
    def axis_label(self) -> str:
        return {
            DivergenceKind.CONSENSUS: "divergence from consensus (%)",
            DivergenceKind.ARRAY_HOMOGENEITY: "array heterogeneity (%)",
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
