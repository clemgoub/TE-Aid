"""Row layout for panel 5 — ORFs with their domains inlaid on them.

The "Stack" design: one row per ORF, its protein hits drawn as arrows on top of
it; a hit with no open frame in its own register gets a bare row instead. That
silhouette — a compact block of framed rows, then a run of bare ones — is the
picture the tool exists to produce, and it is not producible by an
ORF-finder-then-align search at all.

**Housing tests reading frame, not merely overlap.** This is the load-bearing
rule. Coordinate overlap alone will occasionally seat a frameshifted domain
inside a stop-to-stop block belonging to a *different* register, and the panel
then asserts "intact coding domain" for the exact finding it was built to
surface. Frame is one arithmetic line and it removes that whole class of lie.

It also makes the geometry honest for free: a frameshifted alignment changes
register by definition, so it can never be same-frame-housed by two ORFs, and
"a notched bar sitting inside a rectangle" becomes impossible to draw.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# A hit must share this much of itself with an ORF to be housed by it.
HOUSING_OVERLAP = 0.5

# Tick lines that fit in one sub-lane's height. A row's tick names every domain
# on it, so a row with six domains needs three sub-lanes of vertical room even
# when none of them overlap — sizing by sub-lanes alone lets the rosters of
# adjacent rows collide, which is unreadable and was the first thing to break.
NAMES_PER_SUBLANE = 2


def reading_frame(start: int, end: int, strand: str, consensus_length: int) -> int:
    """Reading frame 1-3 for a feature, counted from the strand it reads on.

    ``start``/``end`` are 0-based half-open and already normalised ascending —
    BATH reports reverse hits with descending coordinates, and
    ``homology._parse_tblout`` collapses that before anything reaches here.

    A reverse-strand feature is translated from the other end, so its frame is
    counted from the distance to the consensus end, not from ``start``. Getting
    this sign wrong silently misaligns every hit against every ORF, which is
    worse than having no frame logic at all — hence the tests.
    """
    if strand == "-":
        return ((consensus_length - end) % 3) + 1
    return (start % 3) + 1


@dataclass(slots=True)
class Row:
    """One row of panel 5: an ORF and whatever domains sit in its register."""

    orf: object | None = None
    hits: list = field(default_factory=list)
    # Sub-lane per hit, for the rare case where two housed hits overlap in x.
    lanes: list[int] = field(default_factory=list)

    @property
    def start(self) -> int:
        candidates = [h.start for h in self.hits]
        if self.orf is not None:
            candidates.append(self.orf.start)
        return min(candidates) if candidates else 0

    @property
    def height(self) -> int:
        """Sub-lanes needed, by both the marks *and* the tick roster.

        Whichever needs more room wins: overlapping domains need their own
        sub-lanes, and a long roster needs the lines to print it.
        """
        packed = (max(self.lanes) + 1) if self.lanes else 1
        names = -(-max(1, len(self.hits)) // NAMES_PER_SUBLANE)  # ceil
        return max(1, packed, names)

    @property
    def is_housed(self) -> bool:
        return self.orf is not None


def _houses(orf, hit, consensus_length: int) -> bool:
    if orf.strand != hit.strand:
        return False
    overlap = min(orf.end, hit.end) - max(orf.start, hit.start)
    span = hit.end - hit.start
    if span <= 0 or overlap <= 0 or overlap / span < HOUSING_OVERLAP:
        return False
    return reading_frame(orf.start, orf.end, orf.strand, consensus_length) == reading_frame(
        hit.start, hit.end, hit.strand, consensus_length
    )


def _pack(hits: list) -> list[int]:
    """Sub-lane per hit, so two that overlap in x are never drawn on top of each other."""
    lanes: list[int] = []
    occupied: list[int] = []  # rightmost x used, per sub-lane
    for hit in sorted(hits, key=lambda h: h.start):
        for lane, edge in enumerate(occupied):
            if hit.start >= edge:
                lanes.append(lane)
                occupied[lane] = hit.end
                break
        else:
            lanes.append(len(occupied))
            occupied.append(hit.end)
    return lanes


def build(orfs: list, hits: list, consensus_length: int) -> list[Row]:
    """Lay panel 5 out: housed rows first by position, then the bare ones.

    Rows are ordered by the consensus start of their leftmost feature, with no
    above/below zoning. A zone rule would import a "demoted" valence and would
    conflate four different reasons a hit is unhoused — opposite strand, another
    register, an ORF below the length floor, or getorf not having run at all.
    """
    remaining = list(hits)
    rows: list[Row] = []

    for orf in orfs:
        housed = [h for h in remaining if _houses(orf, h, consensus_length)]
        for hit in housed:
            remaining.remove(hit)
        ordered = sorted(housed, key=lambda h: h.start)
        rows.append(Row(orf=orf, hits=ordered, lanes=_pack(ordered)))

    for hit in remaining:
        rows.append(Row(orf=None, hits=[hit], lanes=[0]))

    rows.sort(key=lambda r: r.start)
    return rows
