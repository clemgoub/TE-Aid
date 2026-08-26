"""Assembly of the evidence sheet.

The layout is v1's 2x2 quadrant grid, deliberately preserved — reading the whole
sheet at a glance, with the dot-plot square, is what the tool is used for::

    +---------------------------+---------------------------+
    | 1  copies vs divergence   | 2  consensus coverage     |
    +---------------------------+---------------------------+
    |                           | 4  structure              |
    | 3  self dot-plot (square) +---------------------------+
    |                           | 5  homology evidence      |
    +---------------------------+---------------------------+

Panels 4 and 5 are the split of v1's single crowded structure quadrant and share
one consensus x-axis so features line up vertically. Until there is homology
evidence to show, panel 4 takes the whole bottom-right quadrant rather than
leaving a visible hole in the grid.

Every quadrant is square and every panel carries the *same* x-range, so a
consensus position lands at the same screen x in all four. The dot-plot holds a
true 1:1 data aspect via ``scaleanchor`` with ``constrain="domain"`` on both of
its axes: Plotly shrinks the plotting box to fit the aspect rather than widening
the data range, which would have made the dot-plot silently wider than its
neighbours.

Segments are drawn as a *single* trace per series, with ``None`` separating each
segment, rather than one trace per copy. A family with 10,000 copies otherwise
produces 10,000 Plotly traces, which no browser renders comfortably.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from . import annotation_rows, classcheck, tetypes
from .analysis import SelfHit
from .orfs import ORF
from .records import Annotation, Copy, DivergenceKind
from .theme import FONT_FAMILY, Theme

# Terminal-repeat lanes drawn per type before the list is cut. The structure
# quadrant is a quarter of the sheet; past a handful of lanes the arrows are too
# thin to read, and the dot-plot already shows every self-match.
MAX_REPEAT_LANES = 4

# Fraction of the consensus added beyond each end of every x-axis.
#
# Without it the axis stops exactly at 0 and at the consensus length, and any
# feature drawn at a terminus is clipped: an arrowhead is a fixed-size marker
# centred on its coordinate, so half of it falls outside the plot. Autoscaling
# the panel fixes the clipping but gives that panel its own range, which is
# worse — the four quadrants are meant to be read against one shared consensus
# scale. Building the margin into the range keeps every panel identical *and*
# uncropped, on load and after a reset alike.
#
# 2.5% clears a 12px arrowhead comfortably at any usable figure width, and sits
# just under the 4% padding R applies by default, which is what v1's panels had.
X_PAD_FRACTION = 0.025

# Sheet geometry. The height is derived, not chosen: it is whatever makes each
# quadrant square.
#
# The dot-plot needs a square plotting box to hold a 1:1 data aspect, and it
# needs the same pixel width as the panel above it or a feature at a given
# consensus position lands at a different screen x in the two panels, which is
# exactly the vertical comparison the grid exists for. Both hold only when the
# quadrant itself is square, so the figure is sized to make it so — the same
# reason v1 used a 12x12 inch page with a 2x2 grid of 6x6 inch panels.
_MARGIN_L, _MARGIN_R, _MARGIN_T, _MARGIN_B = 68, 26, 96, 128
_H_SPACING, _V_SPACING = 0.10, 0.085

SHEET_WIDTH = 1180
_PLOT_WIDTH = SHEET_WIDTH - _MARGIN_L - _MARGIN_R
_CELL = _PLOT_WIDTH * (1 - _H_SPACING) / 2


def sheet_height(cell_rows: float = 2) -> int:
    """Figure height that keeps every cell square for a given number of rows."""
    return round(cell_rows * _CELL / (1 - _V_SPACING) + _MARGIN_T + _MARGIN_B)


SHEET_HEIGHT = sheet_height(2)

# Distance below the plotting area, in pixels, for the legend and the notes
# line. Both must stay inside the bottom margin.
_LEGEND_OFFSET_PX = 62
_NOTE_OFFSET_PX = 104
# Height of one wrapped legend row; the notes line is pushed down by this much
# per extra row so a long legend cannot land on top of it.
_LEGEND_ROW_PX = 22
# Extra top margin for the provenance line under the subtitle.
_SOURCE_LINE_PX = 17
# Average glyph width of a 12px subplot title, and the gap before its '?'.
_TITLE_CHAR_PX = 6.6
_HELP_GAP_PX = 7


def _legend_rows(fig: go.Figure) -> int:
    """How many rows the horizontal legend will wrap onto.

    Plotly wraps a horizontal legend silently, and the extra row lands on top of
    the notes line beneath it. Estimating the wrap from the labels' own width is
    cheaper and more reliable than reserving space for a worst case that usually
    is not there.
    """
    labels = [t.name for t in fig.data if t.name and t.showlegend is not False]
    if not labels:
        return 0
    # ~6.2px per character at 11px, plus the swatch and padding per entry.
    width = sum(len(label) * 6.2 + 46 for label in labels)
    return max(1, int(width / _PLOT_WIDTH) + 1)


def _plot_height(cell_rows: float) -> float:
    return sheet_height(cell_rows) - _MARGIN_T - _MARGIN_B

# Dfam requires a seed alignment to hold at least this many sequences. Drawn as
# a rule on panel 6 so a curator sees at a glance whether the seed clears it and
# where it runs thin.
DFAM_SEQUENCE_FLOOR = 3

# Height of the appended seed-QC row, relative to one main quadrant row. Full
# height because panel 6 carries a pileup as well as a coverage band.
SEED_ROW_HEIGHT = 1.0

# Share of panel 6's cell given to the coverage band; the rest holds the pileup.
_COVERAGE_BASE = 0.52

# Pileup lanes drawn before the list is cut. Past this the lanes are thinner
# than a pixel and the coverage band above already carries the depth.
MAX_PILEUP_LANES = 40

# Protein hits drawn in panel 5 before the list is cut. Hits arrive sorted by
# E-value, so the cut keeps the strongest.
MAX_HOMOLOGY_LANES = 14

# Panel 5 row metrics, in axis units. A sub-lane is one domain's height; a row
# is as many sub-lanes as its busiest overlap needs.
_SUBLANE = 1.0
_ROW_HALF = 0.42
_ROW_GAP = 0.55

# Status colour for the one threshold on the sheet that is a pass/fail
# requirement rather than a measurement. It always ships with the label naming
# the floor, so colour alone carries nothing.
#
# It is the same hex as theme.orf_reverse (v1's reverse-strand ORF red), and on
# a --seed-qc sheet both are visible: reverse ORFs in panel 5, the Dfam floor in
# panel 6. Two meanings for one red. They never share a panel, and both readings
# are labelled, so this is recorded rather than resolved -- see the open items
# in docs/BRIEF_v2.md §6.
STATUS_BELOW_FLOOR = "#d03b3b"

# What each panel is, and the trap it exists to avoid. Shown on hovering the '?'
# beside a panel title. These are the things a curator would otherwise have to
# take on trust or read the docs for — especially where two panels look alike.
PANEL_HELP = {
    1: (
        "One horizontal line per annotated copy, drawn across the stretch of "
        "consensus it matches, at its divergence.<br><br>"
        "Copies whose source reports <b>no</b> divergence are left out entirely "
        "rather than drawn at zero — zero here would read as a pristine, very "
        "recent insertion, the opposite of 'unknown'. Any omitted are counted in "
        "the axis label.<br><br>"
        "Highlighted copies span at least the full-length threshold of the "
        "consensus."
    ),
    2: (
        "How many copies <b>span</b> each consensus position: everything between "
        "a copy's first and last aligned base, internal deletions included."
        "<br><br>"
        "This is not the same quantity as panel 6, which counts only sequences "
        "that contribute an actual base. Where the two differ, the difference is "
        "internal deletion."
    ),
    3: (
        "The consensus aligned against itself, at a true 1:1 aspect — so an "
        "off-diagonal repeat reads as parallel to the main diagonal.<br><br>"
        "Same-strand matches near both termini suggest LTRs; opposite-strand "
        "matches suggest TIRs. An internal segmental duplication produces the "
        "same signature, so these are suggestions for you to judge, never calls."
    ),
    4: (
        "Terminal-repeat candidates, on the consensus axis.<br><br>"
        "Arrowheads carry orientation: a <b>direct</b> pair points the same way "
        "(→ … →), an <b>inverted</b> pair points inward (→ … ←).<br><br>"
        "Lanes are labelled by what was measured — direct or inverted — not by a "
        "class. On an internally repetitive consensus the same signature is a "
        "tandem unit, not a terminal repeat."
    ),
    5: (
        "ORFs as rectangles — <b>black outline forward, red reverse</b>, as v1 "
        "drew them — with domain hits as arrows on the ORF sharing their "
        "reading frame. A domain in another frame gets its own bare row.<br><br>"
        "A label like <b>2fs</b> or <b>1⊗</b> counts frameshifts and in-frame "
        "stops in that hit: a domain that was once coding and has decayed. "
        "Finding a broken domain and finding nothing are different results, and "
        "an ORF-finder-then-align search cannot tell them apart.<br><br>"
        "A hollow arrow marked <b>=</b> is a named element from blastp, not a "
        "domain model. Hits are evidence, never a classification."
    ),
    6: (
        "The seed alignment, in the shape of Dfam's own seed track.<br><br>"
        "<b>Above:</b> coverage split into sequences that match the consensus and "
        "sequences that differ. Depth alone flatters a seed — a column with 43 "
        "sequences of which 25 disagree is not 43 sequences of support.<br><br>"
        "<b>Below:</b> one lane per seed sequence, drawn as its aligned runs, so "
        "an internal deletion is a gap rather than being smoothed over.<br><br>"
        "The rule marks Dfam's 3-sequence minimum; stretches falling short are "
        "shaded. This counts <b>bases</b>, where panel 2 counts spans."
    ),
    8: (
        "The <code>#=GF TP</code> classification carried by the seed.<br><br>"
        "This label arrives <b>with the seed</b> — TE-Aid did not derive it and "
        "never asserts a classification of its own. The deliverable is whether "
        "the evidence on this sheet disagrees with it, which only means something "
        "because the rest of the sheet has no opinion."
    ),
}


@dataclass(slots=True)
class SheetData:
    """Everything the panels need, already computed."""

    family: str
    consensus_length: int
    copies: list[Copy]
    full_length: list[Copy]
    coverage: np.ndarray
    self_hits: list[SelfHit] = field(default_factory=list)
    terminal_repeats: dict[str, list[SelfHit]] = field(default_factory=dict)
    orfs: list[ORF] = field(default_factory=list)
    orf_min_size: int = 400
    class_label: str | None = None
    full_length_threshold: float = 0.9
    source_format: str | None = None
    # Where the evidence came from, one entry per input, shown under the title.
    # A sheet outlives the shell that produced it, and a curator looking at one
    # weeks later needs to know whether it was built from an annotation, a seed
    # or both, and which files. It also surfaces the *detected* annotation
    # format, which is otherwise chosen silently.
    sources: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    # Panel 5 content; empty until the BATH protein row and Dfam nucleotide row
    # land (work order steps 6 and 8).
    homology: list = field(default_factory=list)
    # Every protein hit found, including those collapsed by best_per_region and
    # those past the lane cap. Panel 5 draws a readable subset; this is what the
    # sheet's expandable table lists, so nothing found is only in a log file.
    homology_all: list = field(default_factory=list)
    # Seed-QC (--seed-qc): per-consensus-position alignment depth from the
    # Stockholm seed, and the externally supplied expected classification.
    seed_depth: np.ndarray | None = None
    # Sequences *spanning* each position (internal deletions included), the same
    # quantity panel 2 plots. Carried so panel 6 can draw it faintly behind the
    # base-level depth and make the difference between the two visible.
    seed_span: np.ndarray | None = None
    # Per-position count of sequences differing from the reference, and each
    # sequence's aligned runs, for the coverage split and the pileup.
    seed_mismatches: np.ndarray | None = None
    seed_blocks: list | None = None
    expected_class: str | None = None
    seed_sequence_count: int | None = None
    # Result of comparing the seed's TP against the protein hits. None when the
    # comparison has not been run -- which is different from "no disagreement".
    class_check: object | None = None

    @property
    def has_seed_qc(self) -> bool:
        return self.seed_depth is not None or self.expected_class is not None

    @property
    def n_copies(self) -> int:
        return len(self.copies)

    @property
    def has_homology(self) -> bool:
        return bool(self.homology)

    @property
    def divergence_kind(self) -> DivergenceKind:
        kinds = {c.divergence_kind for c in self.copies if c.has_divergence}
        if not kinds:
            return DivergenceKind.NONE
        if len(kinds) == 1:
            return next(iter(kinds))
        return DivergenceKind.CONSENSUS


def _segments(starts, ends, ys_start, ys_end):
    """Interleave segment endpoints with None breaks for a single Plotly trace."""
    xs: list[float | None] = []
    ys: list[float | None] = []
    for x0, x1, y0, y1 in zip(starts, ends, ys_start, ys_end):
        xs.extend((x0, x1, None))
        ys.extend((y0, y1, None))
    return xs, ys


def build_figure(data: SheetData, theme: Theme) -> go.Figure:
    """Compose the 2x2 evidence sheet."""
    # Titles are appended in the same order make_subplots walks the grid --
    # row-major over the non-None specs -- so they must be built alongside the
    # specs below, not assembled up front. Supplying one title short shifts
    # every later panel's title onto the wrong panel, silently.
    titles = [
        "1 · Annotated copies vs divergence",
        "2 · Consensus coverage",
        "3 · Self dot-plot",
        "4 · Structure",
    ]
    # The grid is assembled from blocks so the default sheet is never disturbed
    # by an optional one. The main block is always the 2x2 quadrants; the
    # bottom-right quadrant splits in two when there is homology evidence, and
    # --seed-qc appends a further row of full-width-cell panels beneath.
    specs: list[list] = [[{}, {}]]
    heights: list[float] = [1.0]
    # Panel 5 now carries the ORFs as well as the homology, so it takes the
    # larger share; panel 4 is down to terminal-repeat candidates alone.
    if data.has_homology or data.orfs:
        specs += [[{"rowspan": 2}, {}], [None, {}]]
        heights += [0.42, 0.58]
        structure_row, homology_row = 2, 3
        titles += ["5 · ORFs and protein homology"]
    else:
        specs += [[{}, {}]]
        heights += [1.0]
        structure_row, homology_row = 2, None

    # The seed-QC block is appended at half height: it holds a depth profile and
    # a text panel, neither of which needs the square cell the dot-plot does.
    # Keeping it short leaves the main 2x2 untouched and the sheet compact.
    seed_row = None
    if data.has_seed_qc:
        seed_row = len(specs) + 1
        specs += [[{}, {}]]
        heights += [SEED_ROW_HEIGHT]
        titles += ["6 · Seed depth", "8 · Expected class"]

    cell_rows = 2 + (SEED_ROW_HEIGHT if data.has_seed_qc else 0)
    # Normalise so each *cell* row gets an equal share, whatever the split above.
    total = sum(heights)
    fig = make_subplots(
        rows=len(specs),
        cols=2,
        specs=specs,
        row_heights=[h / total for h in heights],
        column_widths=[0.5, 0.5],
        horizontal_spacing=_H_SPACING,
        vertical_spacing=_V_SPACING / cell_rows * 2,
        subplot_titles=titles,
    )

    _panel_hits(fig, data, theme, row=1, col=1)
    _panel_coverage(fig, data, theme, row=1, col=2)
    _panel_dotplot(fig, data, theme, row=2, col=1)
    _panel_structure(fig, data, theme, row=structure_row, col=2)
    if homology_row is not None:
        _panel_homology(fig, data, theme, row=homology_row, col=2)
    if seed_row is not None:
        _panel_seed_depth(fig, data, theme, row=seed_row, col=1)
        _panel_expected_class(fig, data, theme, row=seed_row, col=2)

    _layout(
        fig,
        data,
        theme,
        structure_row=structure_row,
        homology_row=homology_row,
        seed_row=seed_row,
        cell_rows=cell_rows,
    )
    return fig


def _panel_hits(fig: go.Figure, data: SheetData, theme: Theme, row: int, col: int) -> None:
    """Panel 1: one horizontal segment per copy, at its divergence.

    Copies with no divergence are omitted rather than drawn at zero: a zero on
    this axis reads as a pristine, very recent insertion, which is the opposite
    of 'we do not know'. The count omitted is reported in the axis title.
    """
    plotted = [
        c
        for c in data.copies
        if c.has_divergence and c.consensus_start is not None and c.consensus_end is not None
    ]
    full_ids = {id(c) for c in data.full_length}
    fragments = [c for c in plotted if id(c) not in full_ids]
    full = [c for c in plotted if id(c) in full_ids]
    omitted = data.n_copies - len(plotted)

    for subset, colour, label, width in (
        (fragments, theme.base, "fragment", 1.2),
        (full, theme.highlight, f"full length ≥{data.full_length_threshold:.0%}", 2.0),
    ):
        if not subset:
            continue
        xs, ys = _segments(
            [c.consensus_start for c in subset],
            [c.consensus_end for c in subset],
            [c.divergence for c in subset],
            [c.divergence for c in subset],
        )
        customdata = []
        for c in subset:
            g0, g1 = c.genomic_1based()
            customdata.extend([[c.chrom, g0, g1, c.strand, c.divergence]] * 3)
        fig.add_trace(
            go.Scattergl(
                x=xs,
                y=ys,
                mode="lines",
                name=f"{label} (n={len(subset)})",
                line=dict(color=colour, width=width),
                # Semi-transparent so overplotted copies read as density, as in
                # v1; kept high enough that the legend swatch stays legible on
                # the dark surface.
                opacity=0.65 if label == "fragment" else 0.9,
                customdata=customdata,
                hovertemplate=(
                    "%{customdata[0]}:%{customdata[1]:,}-%{customdata[2]:,} "
                    "(%{customdata[3]})<br>"
                    "consensus %{x:,.0f} bp<br>divergence %{customdata[4]:.1f}%"
                    "<extra></extra>"
                ),
            ),
            row=row,
            col=col,
        )

    title = data.divergence_kind.axis_label
    if omitted:
        title += f"<br><span style='font-size:10px'>{omitted:,}/{data.n_copies:,} omitted: none reported</span>"
    fig.update_yaxes(title_text=title, row=row, col=col)
    fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)

    if not plotted:
        _empty_note(fig, "no copy carries a divergence value", theme, row, col)


def _panel_coverage(fig: go.Figure, data: SheetData, theme: Theme, row: int, col: int) -> None:
    """Panel 2: how many copies overlap each consensus position.

    v1 labelled this axis 'coverage (bp)' while plotting a count of BLAST HSPs.
    The quantity is a number of copies, and the label now says so.
    """
    positions = np.arange(data.consensus_length)
    fig.add_trace(
        go.Scattergl(
            x=positions,
            y=data.coverage,
            mode="lines",
            name="copies overlapping",
            showlegend=False,
            line=dict(color=theme.base, width=2),
            fill="tozeroy",
            fillcolor=_alpha(theme.base, 0.18),
            hovertemplate="consensus %{x:,.0f} bp<br>%{y:,.0f} copies<extra></extra>",
        ),
        row=row,
        col=col,
    )
    # 'Spanning', not 'covering'. A copy spans every position between its first
    # and last aligned base, internal deletions included. Panel 6 counts only
    # positions where a sequence actually contributes a base, so the two axes
    # are deliberately named for different quantities — see _panel_seed_depth.
    fig.update_yaxes(title_text="copies spanning", rangemode="tozero", row=row, col=col)
    fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)


def _panel_dotplot(fig: go.Figure, data: SheetData, theme: Theme, row: int, col: int) -> None:
    """Panel 3: the consensus against itself, at a true 1:1 aspect.

    Same-strand off-diagonal matches near both termini suggest LTRs;
    opposite-strand matches suggest TIRs. Both are suggestions for the curator,
    never assertions — a segmental duplication produces the same signature.
    """
    if not data.self_hits:
        _empty_note(fig, "self dot-plot unavailable (blastn not run)", theme, row, col)
        fig.update_yaxes(title_text="consensus (bp)", row=row, col=col)
        fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)
        return

    direct = [h for h in data.self_hits if not h.is_reverse]
    inverted = [h for h in data.self_hits if h.is_reverse]

    for subset, colour, label in (
        (direct, theme.base, "direct (same strand)"),
        (inverted, theme.inverted, "inverted (opposite strand)"),
    ):
        if not subset:
            continue
        xs, ys = _segments(
            [h.q_start for h in subset],
            [h.q_end for h in subset],
            [h.s_start for h in subset],
            [h.s_end for h in subset],
        )
        customdata = []
        for h in subset:
            customdata.extend([[h.identity, h.evalue]] * 3)
        fig.add_trace(
            go.Scattergl(
                x=xs,
                y=ys,
                mode="lines",
                name=f"{label} (n={len(subset)})",
                line=dict(color=colour, width=2),
                customdata=customdata,
                hovertemplate=(
                    "query %{x:,.0f} bp<br>subject %{y:,.0f} bp<br>"
                    "identity %{customdata[0]:.1f}%<br>E %{customdata[1]:.1g}"
                    "<extra></extra>"
                ),
            ),
            row=row,
            col=col,
        )

    fig.update_yaxes(title_text="consensus (bp)", row=row, col=col)
    fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)


def _panel_structure(fig: go.Figure, data: SheetData, theme: Theme, row: int, col: int) -> None:
    """Panel 4: terminal-repeat candidates and ORFs on the consensus axis.

    Colour carries the same meaning as in the dot-plot — blue for direct
    repeats, aqua for inverted — so the two panels read together. ORF strand is
    shown by arrowhead direction rather than colour, keeping direction as a
    shape channel and leaving colour free to mean 'coding feature'.
    """
    y = 0.0
    tickvals: list[float] = []
    ticktext: list[str] = []
    drew_anything = False

    ltr = data.terminal_repeats.get("LTR", [])
    tir = data.terminal_repeats.get("TIR", [])

    # v1's genuinely good idea, kept: draw both arms of a terminal repeat as
    # arrows and let the arrowheads carry the orientation. Direct repeats point
    # the same way (-> ... ->), inverted repeats point at each other
    # (-> ... <-), so LTR and TIR are distinguishable by shape alone rather
    # than by colour.
    #
    # Three v1 bugs are not reproduced. It placed each pair at y = i, the row
    # index of the *unfiltered* self-blast table, so a consensus with 50 hits
    # and 3 drawn pairs scattered 3 arrows over 50 empty lanes; lanes here are
    # packed consecutively. It coloured by rainbow(n)[i], which carries no
    # meaning, shifts with the hit count and is not colourblind-safe; colour
    # here is the LTR/TIR distinction and matches the dot-plot. And it drew
    # every self-similarity under the heading 'structure', including purely
    # internal repeats; only terminal-repeat candidates reach this panel.
    seen_types: set[str] = set()
    # The lane names what was measured — a direct or inverted repeat pair — and
    # the legend offers the interpretation. Calling a lane 'LTR' would assert a
    # classification, and on an internally repetitive consensus the same
    # signature is a tandem unit, not a long terminal repeat.
    for all_hits, colour, label, lane_label in (
        (ltr, theme.base, "direct pair (LTR-like)", "direct"),
        (tir, theme.inverted, "inverted pair (TIR-like)", "inverted"),
    ):
        # An internally repetitive consensus can yield many genuine repeat
        # pairs. Show the strongest few rather than letting one family fill the
        # quadrant, and say so when the list is cut.
        hits = sorted(all_hits, key=lambda h: h.bitscore, reverse=True)[:MAX_REPEAT_LANES]
        hits.sort(key=lambda h: min(h.q_start, h.s_start, h.s_end))
        hidden = len(all_hits) - len(hits)
        if hidden:
            data.notes.append(
                f"{hidden} further {lane_label} repeat pair"
                f"{'s' if hidden != 1 else ''} not drawn "
                f"(showing the {MAX_REPEAT_LANES} strongest)"
            )
        for hit in hits:
            arm_hover = (
                f"{label}<br>arms %{{x:,.0f}} bp<br>identity {hit.identity:.1f}%"
            )
            # Arrow direction is decided from the arms' positions, not from
            # blastn's coordinate order. blastn reports a repeat pair in
            # whichever direction it happened to find it, so an inverted pair
            # came out pointing inward or outward depending on which reciprocal
            # survived deduplication. Inverted repeats always point inward here:
            # that is the palindrome the curator is looking for, and it makes
            # the two types distinguishable at a glance even in grayscale.
            left_arm, right_arm = _ordered_arms(hit)
            inverted = hit.is_reverse
            for index, (tail, head) in enumerate(
                (
                    (left_arm[0], left_arm[1]),
                    (right_arm[1], right_arm[0]) if inverted else (right_arm[0], right_arm[1]),
                )
            ):
                _arrow(
                    fig, tail, head, y, colour, arm_hover, row, col,
                    show_legend=index == 0 and label not in seen_types,
                    legend_name=f"{label} (n={len(all_hits)})",
                )
            seen_types.add(label)
            tickvals.append(y)
            ticktext.append(lane_label)
            y -= 1.0
            drew_anything = True

    # ORFs have moved to panel 5, where they can be drawn together with the
    # protein hits that sit in them. Panel 4 now holds only what the self
    # dot-plot above it implies, so the two read as one statement.

    # A rail across the full consensus, so how much of the element the repeats
    # bracket is visible rather than inferred from tick labels.
    if drew_anything:
        fig.add_trace(
            go.Scatter(
                x=[0, data.consensus_length],
                y=[y + 0.25, y + 0.25],
                mode="lines",
                line=dict(color=theme.axis, width=1),
                showlegend=False,
                hovertemplate=f"consensus 0-{data.consensus_length:,} bp<extra></extra>",
            ),
            row=row,
            col=col,
        )
        tickvals.append(y + 0.25)
        ticktext.append(f"{data.consensus_length:,} bp")
        y -= 0.6

    if not drew_anything:
        _empty_note(fig, "no terminal repeat candidates", theme, row, col)

    fig.update_yaxes(
        tickvals=tickvals,
        ticktext=ticktext,
        range=[y - 0.4, 0.8],
        showgrid=False,
        zeroline=False,
        tickfont=dict(size=10, color=theme.muted),
        row=row,
        col=col,
    )
    fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)


def _panel_seed_depth(fig: go.Figure, data: SheetData, theme: Theme, row: int, col: int) -> None:
    """Panel 6: the seed alignment, laid out as Dfam's own seed track.

    Two stacked parts sharing the consensus axis, the way a genome browser shows
    an alignment — Dfam's family browser uses igv.js for exactly this:

    * **Coverage**, on top, split into sequences that match the consensus and
      sequences that differ. Depth alone flatters a seed: a column with 43
      sequences of which 25 disagree is not 43 sequences of support, and only
      the split shows it. Dfam colours its coverage bars by allele fraction for
      the same reason.
    * **The pileup**, below: one lane per sequence, drawn as its aligned runs so
      an internal deletion appears as a gap in the lane rather than being
      smoothed over. Truncation and interruption look different, which is the
      point.

    Dfam's floor of three sequences is a rule across the coverage part, with the
    stretches that fall short shaded.

    **This is not panel 2 restated.** Panel 2 counts copies that *span* a
    position, internal deletions included; this counts sequences that contribute
    a *base*. Across the 409 GenomeArk seeds for one assembly, 362 differ, by as
    much as 37 sequences at a single position — so the spanning curve is drawn
    here too, faintly, and the gap between the two is the deletion structure.
    """
    depth = data.seed_depth
    if depth is None or not len(depth):
        _empty_note(fig, "no seed alignment", theme, row, col)
        fig.update_yaxes(title_text="seed alignment", row=row, col=col)
        fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)
        return

    positions = np.arange(len(depth))
    mismatch = data.seed_mismatches
    if mismatch is None or len(mismatch) != len(depth):
        mismatch = np.zeros_like(depth)
    agreeing = np.maximum(depth - mismatch, 0)

    ceiling = max(
        int(depth.max()),
        int(data.seed_span.max()) if data.seed_span is not None and len(data.seed_span) else 0,
        DFAM_SEQUENCE_FLOOR,
    ) + 1

    # The cell is split: coverage in the upper band, pileup lanes beneath. Both
    # are drawn on one normalised axis so they share the consensus scale exactly.
    def to_y(value: float) -> float:
        return _COVERAGE_BASE + (value / ceiling) * (1.0 - _COVERAGE_BASE)

    # An explicit baseline at the foot of the coverage band. The fills stack
    # onto it rather than onto y=0, which is the floor of the whole cell -- a
    # 'tozeroy' fill here floods the pileup underneath.
    fig.add_trace(
        go.Scatter(
            x=positions, y=np.full(len(positions), _COVERAGE_BASE),
            mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip",
        ),
        row=row, col=col,
    )
    # Matching sequences, then the disagreeing remainder stacked on top, so the
    # total height is the depth and the coloured cap is the disagreement.
    fig.add_trace(
        go.Scatter(
            x=positions, y=[to_y(v) for v in agreeing],
            mode="lines", name="matches consensus",
            line=dict(color=theme.base, width=1.5),
            fill="tonexty", fillcolor=_alpha(theme.base, 0.30),
            hovertemplate="consensus %{x:,.0f} bp<extra>matching</extra>",
        ),
        row=row, col=col,
    )
    if mismatch.any():
        fig.add_trace(
            go.Scatter(
                x=positions, y=[to_y(v) for v in depth],
                mode="lines", name="differs from consensus",
                line=dict(color=theme.highlight, width=1.5),
                fill="tonexty", fillcolor=_alpha(theme.highlight, 0.45),
                customdata=np.stack([depth, mismatch], axis=-1),
                hovertemplate=(
                    "consensus %{x:,.0f} bp<br>%{customdata[0]:,.0f} aligned, "
                    "%{customdata[1]:,.0f} differing<extra></extra>"
                ),
            ),
            row=row, col=col,
        )

    span = data.seed_span
    if span is not None and len(span) == len(depth) and not np.array_equal(span, depth):
        fig.add_trace(
            go.Scatter(
                x=positions, y=[to_y(v) for v in span],
                mode="lines", name="sequences spanning",
                line=dict(color=theme.muted, width=1.2, dash="dot"),
                hovertemplate="consensus %{x:,.0f} bp<br>%{y}<extra>spanning</extra>",
            ),
            row=row, col=col,
        )
        deleted = int((span - depth).sum())
        if deleted:
            data.notes.append(
                f"{deleted:,} sequence-positions fall inside internal deletions "
                f"(the gap between spanning and aligned)"
            )

    floor_y = to_y(DFAM_SEQUENCE_FLOOR)
    fig.add_trace(
        go.Scatter(
            x=[-len(depth), len(depth) * 2], y=[floor_y, floor_y],
            mode="lines", name=f"Dfam floor ({DFAM_SEQUENCE_FLOOR} sequences)",
            line=dict(color=STATUS_BELOW_FLOOR, width=1.5, dash="dash"),
            hoverinfo="skip",
        ),
        row=row, col=col,
    )

    # The pileup: one lane per sequence, drawn as its aligned runs.
    blocks = data.seed_blocks or []
    shown = blocks[:MAX_PILEUP_LANES]
    if len(blocks) > len(shown):
        data.notes.append(
            f"pileup shows {len(shown)} of {len(blocks):,} seed sequences"
        )
    if shown:
        lane_top = _COVERAGE_BASE - 0.06
        step = lane_top / (len(shown) + 1)
        xs: list[float | None] = []
        ys: list[float | None] = []
        for index, (_, runs) in enumerate(shown):
            y = lane_top - index * step
            for run_start, run_end in runs:
                xs.extend((run_start, run_end, None))
                ys.extend((y, y, None))
        fig.add_trace(
            go.Scattergl(
                x=xs, y=ys, mode="lines",
                name=f"seed sequences (n={len(blocks)})",
                line=dict(color=theme.base, width=max(1.0, min(4.0, 160 / len(shown)))),
                opacity=0.75,
                hovertemplate="consensus %{x:,.0f} bp<extra>seed sequence</extra>",
            ),
            row=row, col=col,
        )

    for start, end in _runs_below(depth, DFAM_SEQUENCE_FLOOR):
        fig.add_vrect(
            x0=start, x1=end, fillcolor=STATUS_BELOW_FLOOR, opacity=0.16,
            line_width=0, layer="below", row=row, col=col,
        )

    thin = int((depth < DFAM_SEQUENCE_FLOOR).sum())
    if thin:
        data.notes.append(
            f"seed is below the Dfam {DFAM_SEQUENCE_FLOOR}-sequence floor at "
            f"{thin:,} of {len(depth):,} consensus positions"
        )

    ticks = _depth_ticks(ceiling)
    fig.update_yaxes(
        title_text="sequences · pileup",
        range=[0, 1.0],
        tickvals=[to_y(v) for v in ticks],
        ticktext=[str(v) for v in ticks],
        showgrid=False,
        zeroline=False,
        row=row, col=col,
    )
    fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)


def _depth_ticks(ceiling: int) -> list[int]:
    """A few round depth values to label the coverage band with."""
    if ceiling <= 6:
        return list(range(0, ceiling + 1))
    step = max(1, round(ceiling / 4))
    return [v for v in range(0, ceiling + 1, step)]


def _panel_expected_class(
    fig: go.Figure, data: SheetData, theme: Theme, row: int, col: int
) -> None:
    """Panel 8: the seed's ``#=GF TP``, and whether the evidence contradicts it.

    The label is *not* TE-Aid's opinion — it arrives with the seed. The
    deliverable is the flag, and the flag only means something because the rest
    of the sheet has no opinion of its own (§1). With no homology evidence drawn
    yet, the comparison is stated as pending rather than guessed at.
    """
    lines: list[str] = []
    if data.expected_class:
        # Dfam's TP is a full semicolon-separated path; show it broken up rather
        # than as one unreadable run.
        path = [p.strip() for p in data.expected_class.split(";") if p.strip()]
        lines.append("<b>#=GF TP</b> (supplied with the seed)")
        lines.append("<br>".join(f"&#8195;{'└ ' if i else ''}{p}" for i, p in enumerate(path)))
    else:
        lines.append(
            f"<span style='color:{theme.muted}'>the seed carries no "
            f"<b>#=GF TP</b>; nothing to check the evidence against</span>"
        )

    if data.seed_sequence_count is not None:
        qualifies = data.seed_sequence_count >= DFAM_SEQUENCE_FLOOR
        colour = theme.text_secondary if qualifies else STATUS_BELOW_FLOOR
        verdict = "meets" if qualifies else "below"
        lines.append(
            f"<br><span style='color:{colour}'>{data.seed_sequence_count} sequence"
            f"{'s' if data.seed_sequence_count != 1 else ''} — {verdict} the Dfam "
            f"minimum of {DFAM_SEQUENCE_FLOOR}</span>"
        )

    if data.expected_class:
        check = data.class_check
        if check is None:
            lines.append(
                f"<br><span style='color:{theme.muted}'>agreement with homology "
                f"evidence: pending (panel 5)</span>"
            )
        elif check.disagrees:
            lines.append(
                f"<br><span style='color:{STATUS_BELOW_FLOOR}'><b>⚠ disagrees with "
                f"the evidence</b><br>{_wrap_help(check.detail)}</span>"
            )
        else:
            lines.append(
                f"<br><span style='color:{theme.muted}'>no disagreement: "
                f"{check.detail}</span>"
            )

    fig.add_trace(
        go.Scatter(
            x=[0], y=[0], mode="markers",
            marker=dict(opacity=0, size=1), showlegend=False, hoverinfo="skip",
        ),
        row=row,
        col=col,
    )
    fig.add_annotation(
        text="<br>".join(lines),
        xref="x domain",
        yref="y domain",
        x=0.02,
        y=0.94,
        xanchor="left",
        yanchor="top",
        align="left",
        showarrow=False,
        font=dict(size=11, color=theme.text_secondary, family=FONT_FAMILY),
        row=row,
        col=col,
    )
    fig.update_xaxes(
        showticklabels=False, showgrid=False, zeroline=False, title_text=None,
        row=row, col=col,
    )
    fig.update_yaxes(
        showticklabels=False, showgrid=False, zeroline=False, row=row, col=col
    )


def homology_lanes(hits: list, limit: int = MAX_HOMOLOGY_LANES) -> list:
    """The hits panel 5 has room to draw, chosen by strength and laid out by position.

    One function so the panel and the sheet's hit table cannot disagree about
    which hits were drawn — they did, briefly, and the table then claimed a hit
    was on the panel when it was not.

    Slicing the position-ordered list instead keeps the leftmost hits and drops
    the strongest, while the note claims the opposite. Strength is ranked within
    a source, since a blastp bitscore and a profile-HMM bit score are not
    comparable, and taken round-robin so one source cannot crowd out the other.
    """
    if len(hits) <= limit:
        return list(hits)

    by_source: dict[str, list] = {}
    for hit in hits:
        by_source.setdefault(getattr(hit, "source", "pfam"), []).append(hit)
    for group in by_source.values():
        group.sort(key=lambda h: (h.evalue, -h.score))

    shown: list = []
    index = 0
    while len(shown) < limit and any(index < len(g) for g in by_source.values()):
        for group in by_source.values():
            if index < len(group):
                shown.append(group[index])
                if len(shown) == limit:
                    break
        index += 1
    shown.sort(key=lambda h: h.start)
    return shown


def _panel_homology(fig: go.Figure, data: SheetData, theme: Theme, row: int, col: int) -> None:
    """Panel 5: ORFs with their protein domains inlaid on them.

    One row per ORF, drawn as a rectangle outlined by strand — black forward,
    red reverse, as v1 did — with the domains that sit in that ORF's *register*
    drawn as arrows on top of it. A domain with no open frame in its own
    register gets a bare row over a dotted ground rule instead.

    Domain colour is v1's TE-class scheme (green LTR, blue LINE, salmon DNA
    transposon), so a curator reads type at a glance the way they always have.
    Colour is never the only channel: every row names its contents in the tick
    gutter, which is the one text space on the sheet that cannot collide.

    The silhouette is the point. A compact block of framed rows above a run of
    bare, notched rows is a decayed element, and no ORF-finder-then-align search
    can draw it — the bare rows are precisely the domains such a search cannot
    see. The sheet still asserts nothing: what it states is "no ORF >= the
    length floor in this register here" and "the alignment carries N
    frameshifts".
    """
    hits = homology_lanes(list(data.homology))
    if len(data.homology) > len(hits):
        hidden = len(data.homology) - len(hits)
        data.notes.append(
            f"{hidden} further protein hit{'s' if hidden != 1 else ''} not drawn "
            f"(showing the {MAX_HOMOLOGY_LANES} strongest); all are listed under the sheet"
        )

    rows = annotation_rows.build(list(data.orfs), hits, data.consensus_length)
    if not rows:
        _empty_note(fig, "no ORFs and no protein homology", theme, row, col)
        fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=row, col=col)
        fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)
        return

    orders = classcheck.domain_orders()
    tickvals: list[float] = []
    ticktext: list[str] = []
    seen_classes: dict[str, str] = {}
    seen_strands: set[str] = set()
    y = 0.0

    for entry in rows:
        height = entry.height
        # A row is as tall as its deepest sub-lane, so a busy ORF never draws
        # two domains on top of each other.
        centre = y - (height - 1) * _SUBLANE / 2

        if entry.orf is not None:
            orf = entry.orf
            stroke = theme.orf_forward if orf.strand == "+" else theme.orf_reverse
            half = (height - 1) * _SUBLANE / 2 + _ROW_HALF
            fig.add_trace(
                go.Scatter(
                    x=[orf.start, orf.end, orf.end, orf.start, orf.start],
                    y=[centre - half, centre - half, centre + half, centre + half, centre - half],
                    mode="lines",
                    line=dict(color=stroke, width=1.5),
                    fill="toself",
                    fillcolor=_alpha(stroke, 0.10),
                    name=f"ORF {orf.strand} strand",
                    showlegend=orf.strand not in seen_strands,
                    legendgroup=f"orf{orf.strand}",
                    hovertemplate=(
                        f"ORF {orf.coordinates_1based()[0]:,}-{orf.coordinates_1based()[1]:,} "
                        f"({orf.strand})<br>{orf.length_aa:,} aa"
                        f"<extra></extra>"
                    ),
                ),
                row=row,
                col=col,
            )
            seen_strands.add(orf.strand)
        else:
            # A positive mark for "the frame is not open along this row", in
            # theme.axis rather than theme.grid, which is too faint to read.
            fig.add_trace(
                go.Scatter(
                    x=[0, data.consensus_length],
                    y=[centre, centre],
                    mode="lines",
                    line=dict(color=theme.axis, width=1, dash="dot"),
                    showlegend=False,
                    hovertemplate=(
                        f"no ORF ≥ {data.orf_min_size} bp in this reading frame"
                        f"<extra></extra>"
                    ),
                ),
                row=row,
                col=col,
            )

        # Spread the arrows over the sub-lanes actually *used*, not over the row's
        # nominal height. The height is also driven by the tick roster
        # (NAMES_PER_SUBLANE), so a row with six domains on one lane is three
        # units tall -- sizing the spread by that pushed the single row of arrows
        # to the bottom of the row while its label stayed centred.
        used = (max(entry.lanes) + 1) if entry.lanes else 1
        roster: list[tuple[float, int, str]] = []
        for hit, lane in zip(entry.hits, entry.lanes):
            hy = centre + (lane - (used - 1) / 2) * _SUBLANE
            key = tetypes.classify_hit(hit, orders)
            tint = tetypes.colour(key)
            disrupted = getattr(hit, "is_disrupted", False)
            derived = getattr(hit, "source", "pfam") == "repeatpeps-blastp"
            label = getattr(hit, "display_name", hit.query)

            _hit_bar(fig, hit, hy, tint, theme, row, col, hollow=derived)
            if disrupted:
                # A notch cut out of the bar, in the surface colour, at a fixed
                # pixel size: it reads at every bar width from 2 px to 150 px,
                # where a dash pattern needs ~25 px of shaft and most domains
                # are narrower than that.
                fig.add_trace(
                    go.Scatter(
                        x=[(hit.start + hit.end) / 2],
                        y=[hy],
                        mode="markers",
                        marker=dict(symbol="line-ns", size=9,
                                    line=dict(width=2.5, color=theme.surface)),
                        showlegend=False,
                        hoverinfo="skip",
                    ),
                    row=row,
                    col=col,
                )
            if key not in seen_classes:
                seen_classes[key] = tint

            mark = ""
            if disrupted:
                bits = []
                if hit.frameshifts:
                    bits.append(f"{hit.frameshifts}fs")
                if hit.stop_codons:
                    bits.append(f"{hit.stop_codons}⊗")
                mark = " " + " ".join(bits)
            roster.append((hy, hit.start, ("= " if derived else "") + label + mark))

        # The roster must read top-to-bottom in the order the arrows actually
        # sit, or it names the wrong bar. Sub-lane 0 draws at the *bottom* of the
        # row (hy rises with the lane index), while the hits arrive sorted by
        # start -- so listing them as they come inverted the pairing on every
        # multi-lane row. Sort by descending y, then left-to-right among arrows
        # sharing a lane.
        roster.sort(key=lambda item: (-item[0], item[1]))
        names = [label for _, _, label in roster]

        if names:
            ticktext.append("<br>".join(names))
        elif entry.orf is not None:
            ticktext.append(f"{entry.orf.length_aa:,} aa {entry.orf.strand}")
        else:
            ticktext.append("")
        tickvals.append(centre)
        y -= height * _SUBLANE + _ROW_GAP

    # One legend entry per TE class actually present, so the key stays as short
    # as the sheet allows.
    for key, tint in seen_classes.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None], mode="lines",
                line=dict(color=tint, width=7),
                name=key, hoverinfo="skip",
            ),
            row=row, col=col,
        )

    fig.update_yaxes(
        tickvals=tickvals,
        ticktext=ticktext,
        range=[y - _ROW_GAP, _ROW_HALF + _ROW_GAP],
        showgrid=False,
        zeroline=False,
        tickfont=dict(size=9, color=theme.muted),
        row=row,
        col=col,
    )
    fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)


def _hit_bar(fig, hit, y, tint, theme, row, col, *, hollow: bool) -> None:
    """One domain: a bar with a scaled arrowhead at its 3' end.

    ``hollow`` marks a tier-3 hit, whose consensus coordinates are computed
    *from* the ORF it came out of, so its containment is arithmetic rather than
    a finding. Drawing it solid would turn a tautology into evidence.
    """
    tail, head = (hit.start, hit.end) if hit.strand == "+" else (hit.end, hit.start)
    fig.add_trace(
        go.Scatter(
            x=[tail, head],
            y=[y, y],
            mode="lines",
            line=dict(color=tint, width=2 if hollow else 7),
            opacity=0.95,
            showlegend=False,
            hovertemplate=(
                f"<b>{getattr(hit, 'display_name', hit.query)}</b> "
                f"{hit.query_accession}<br>"
                f"consensus {hit.start + 1:,}-{hit.end:,} ({hit.strand})<br>"
                f"E {hit.evalue:.1g} · score {hit.score:.0f} · {hit.identity:.0f}% id<br>"
                f"model coverage {hit.coverage:.0%}"
                + (f"<br>{hit.frameshifts} frameshift(s), {hit.stop_codons} stop(s)"
                   if getattr(hit, "is_disrupted", False) else "")
                + "<extra></extra>"
            ),
        ),
        row=row,
        col=col,
    )
    fig.add_trace(
        go.Scatter(
            x=[head],
            y=[y],
            mode="markers",
            marker=dict(
                symbol="triangle-right" if head >= tail else "triangle-left",
                size=8,
                color=theme.surface if hollow else tint,
                line=dict(width=1.5, color=tint),
            ),
            showlegend=False,
            hoverinfo="skip",
        ),
        row=row,
        col=col,
    )


def _ordered_arms(hit: SelfHit) -> tuple[tuple[int, int], tuple[int, int]]:
    """The hit's two arms as ascending intervals, ordered left to right."""
    query = (hit.q_start, hit.q_end)
    subject = tuple(sorted((hit.s_start, hit.s_end)))
    left, right = sorted((query, subject))
    return left, right  # type: ignore[return-value]


def _arrow(
    fig: go.Figure,
    x_tail: float,
    x_head: float,
    y: float,
    colour: str,
    hover: str,
    row: int,
    col: int,
    *,
    show_legend: bool = False,
    legend_name: str | None = None,
) -> None:
    """One horizontal arrow from ``x_tail`` to ``x_head`` on lane ``y``.

    The arrowhead is placed at ``x_head`` and points in the direction of travel,
    so a feature's orientation is readable from its shape rather than from a
    colour key. ``x_head < x_tail`` is normal and means the feature runs right
    to left.
    """
    fig.add_trace(
        go.Scatter(
            x=[x_tail, x_head],
            y=[y, y],
            mode="lines",
            line=dict(color=colour, width=8),
            opacity=0.85,
            name=legend_name or "",
            showlegend=show_legend,
            legendgroup=legend_name or None,
            hovertemplate=f"{hover}<extra></extra>",
        ),
        row=row,
        col=col,
    )
    fig.add_trace(
        go.Scatter(
            x=[x_head],
            y=[y],
            mode="markers",
            marker=dict(
                symbol="triangle-right" if x_head >= x_tail else "triangle-left",
                size=12,
                color=colour,
            ),
            opacity=0.85,
            showlegend=False,
            legendgroup=legend_name or None,
            hoverinfo="skip",
        ),
        row=row,
        col=col,
    )


def _empty_note(fig: go.Figure, text: str, theme: Theme, row: int, col: int) -> None:
    """Explain why a quadrant is empty, keeping the quadrant itself drawn.

    Plotly only renders a subplot's axes if some trace references them, so an
    empty panel would otherwise vanish — taking the 2x2 grid with it and
    stranding this annotation over a neighbouring panel. An invisible anchor
    point keeps the axes, and therefore the grid, intact.
    """
    fig.add_trace(
        go.Scatter(
            x=[0],
            y=[0],
            mode="markers",
            marker=dict(opacity=0, size=1),
            showlegend=False,
            hoverinfo="skip",
        ),
        row=row,
        col=col,
    )
    fig.add_annotation(
        text=text,
        xref="x domain",
        yref="y domain",
        x=0.5,
        y=0.5,
        showarrow=False,
        font=dict(color=theme.muted, size=11, family=FONT_FAMILY),
        row=row,
        col=col,
    )


def _layout(
    fig: go.Figure,
    data: SheetData,
    theme: Theme,
    *,
    structure_row: int,
    homology_row: int | None,
    seed_row: int | None = None,
    cell_rows: int = 2,
) -> None:
    # A wrapped legend needs the space it takes, or it lands on the notes line.
    _extra_bottom = max(0, _legend_rows(fig) - 1) * _LEGEND_ROW_PX

    subtitle_bits = [f"{data.consensus_length:,} bp", f"{data.n_copies:,} copies"]
    if data.full_length:
        subtitle_bits.append(
            f"{len(data.full_length):,} full length ≥{data.full_length_threshold:.0%}"
        )
    if data.orfs:
        subtitle_bits.append(f"{len(data.orfs)} ORF{'s' if len(data.orfs) != 1 else ''}")
    if data.class_label:
        subtitle_bits.append(data.class_label)

    title_text = (
        f"<b>{data.family}</b><br>"
        f"<span style='font-size:12px;color:{theme.text_secondary}'>"
        f"{' &nbsp;·&nbsp; '.join(subtitle_bits)}</span>"
    )
    if data.sources:
        title_text += (
            f"<br><span style='font-size:11px;color:{theme.muted}'>"
            f"from &nbsp;{' &nbsp;+&nbsp; '.join(data.sources)}</span>"
        )
    # The provenance line needs its own room in the top margin.
    top_margin = _MARGIN_T + (_SOURCE_LINE_PX if data.sources else 0)

    fig.update_layout(
        title=dict(
            text=title_text,
            x=0.012,
            xanchor="left",
            yref="container",
            y=0.985,
            yanchor="top",
            font=dict(size=18, color=theme.text_primary, family=FONT_FAMILY),
        ),
        width=SHEET_WIDTH,
        height=sheet_height(cell_rows) + _extra_bottom + (top_margin - _MARGIN_T),
        paper_bgcolor=theme.paper,
        plot_bgcolor=theme.surface,
        font=dict(family=FONT_FAMILY, color=theme.text_secondary, size=11),
        # The legend sits below the grid, centred. In a 2x2 the top-right
        # subplot title occupies the same band as a top-anchored legend, so
        # anywhere at the top collides with it by construction.
        legend=dict(
            orientation="h",
            yanchor="top",
            # Offsets below the plot are expressed in pixels and converted, so
            # they stay put as the sheet grows taller with optional blocks. A
            # fixed paper fraction scales with plot height and pushes the note
            # past the bottom margin on a tall sheet.
            y=-_LEGEND_OFFSET_PX / _plot_height(cell_rows),
            xanchor="center",
            x=0.5,
            bgcolor="rgba(0,0,0,0)",
            font=dict(color=theme.text_secondary, size=11),
        ),
        margin=dict(
            l=_MARGIN_L, r=_MARGIN_R, t=top_margin, b=_MARGIN_B + _extra_bottom
        ),
        hovermode="closest",
        dragmode="pan",
    )

    fig.update_xaxes(
        showgrid=True,
        gridcolor=theme.grid,
        gridwidth=1,
        zeroline=False,
        linecolor=theme.axis,
        tickcolor=theme.axis,
        tickfont=dict(color=theme.muted, size=10),
        title_font=dict(color=theme.text_secondary, size=11),
        range=_padded_span(data.consensus_length),
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor=theme.grid,
        gridwidth=1,
        zeroline=False,
        linecolor=theme.axis,
        tickcolor=theme.axis,
        tickfont=dict(color=theme.muted, size=10),
        title_font=dict(color=theme.text_secondary, size=11),
    )

    # The dot-plot is the one panel with a meaningful data aspect: both axes are
    # consensus base pairs, so a 1:1 ratio is what makes an off-diagonal repeat
    # read as parallel to the main diagonal. constrain='domain' shrinks the
    # plotting box to honour it instead of widening the range.
    dotplot_axis = fig.get_subplot(2, 1)
    span = _padded_span(data.consensus_length)
    fig.update_yaxes(
        range=span,
        scaleanchor=dotplot_axis.xaxis.anchor.replace("y", "x"),
        scaleratio=1,
        constrain="domain",
        row=2,
        col=1,
    )
    # Both axes of the dot-plot carry the same padded span, and both satisfy the
    # 1:1 constraint by shrinking their domain. That keeps the plotting box
    # square while leaving the *ranges* untouched, so the dot-plot's x-axis
    # still matches the other three exactly. Letting the range float instead —
    # Plotly's default — makes the dot-plot silently wider than its neighbours
    # and destroys the alignment the grid exists for.
    fig.update_xaxes(range=span, constrain="domain", row=2, col=1)

    # Panels 4 and 5 share one consensus x-axis so features line up vertically.
    if homology_row is not None:
        structure_axis = fig.get_subplot(structure_row, 2)
        fig.update_xaxes(
            matches=structure_axis.xaxis.anchor.replace("y", "x"),
            row=homology_row,
            col=2,
        )
        fig.update_xaxes(title_text=None, showticklabels=False, row=structure_row, col=2)

    # Anything the sheet chose not to draw is stated on the sheet. A silent
    # truncation reads as 'this is everything', which is the one thing an
    # evidence sheet must never imply.
    notes = list(dict.fromkeys(data.notes))
    if notes:
        fig.add_annotation(
            text=" · ".join(notes),
            xref="paper",
            yref="paper",
            x=0,
            y=-(_NOTE_OFFSET_PX + _extra_bottom) / _plot_height(cell_rows),
            xanchor="left",
            yanchor="top",
            showarrow=False,
            font=dict(size=10, color=theme.muted, family=FONT_FAMILY),
        )

    titles = []
    for annotation in fig.layout.annotations:
        number = _panel_number(annotation.text)
        if number is None:
            continue
        annotation.update(
            font=dict(size=12, color=theme.text_primary, family=FONT_FAMILY),
            x=annotation.x,
            xanchor="center",
        )
        titles.append((number, annotation))

    for number, annotation in titles:
        help_text = PANEL_HELP.get(number)
        if help_text:
            _add_help_marker(fig, annotation, help_text, theme)


# A **raw** string. The body is JavaScript and CSS, so a backslash in it is
# meant for the browser, not for Python: as a normal string, the `'\t'` and
# `'\n'` in the TSV export below became a real tab and a real newline, and a
# literal newline inside a JS string literal is a SyntaxError that kills the
# whole <script> block -- taking the reset button, which sits in the same block
# and had nothing wrong with it, down with it. Keep the `r` prefix.
_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  :root {{ color-scheme: {scheme}; }}
  body {{
    margin: 0;
    background: {paper};
    color: {text};
    font-family: {font};
  }}
  .teaid-bar {{
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 10px 16px 0;
    max-width: 1180px;
    margin: 0 auto;
  }}
  .teaid-reset {{
    font: inherit;
    font-size: 12px;
    color: {text};
    background: {surface};
    border: 1px solid {axis};
    border-radius: 6px;
    padding: 5px 12px;
    cursor: pointer;
  }}
  .teaid-reset:hover {{ border-color: {accent}; color: {accent}; }}
  .teaid-reset:focus-visible {{ outline: 2px solid {accent}; outline-offset: 2px; }}
  .teaid-hint {{ font-size: 11px; color: {muted}; }}
  .teaid-plot {{ max-width: 1180px; margin: 0 auto; }}
  .teaid-hits {{
    max-width: 1180px; margin: 4px auto 40px; padding: 0 16px;
    font-size: 12px; color: {text};
  }}
  .teaid-hits > summary {{
    cursor: pointer; padding: 7px 0; color: {muted};
    border-top: 1px solid {axis}; list-style: revert;
  }}
  .teaid-hits > summary:hover {{ color: {accent}; }}
  .teaid-hits > summary:focus-visible {{ outline: 2px solid {accent}; outline-offset: 2px; }}
  .teaid-hits-actions {{ display: flex; gap: 8px; margin: 8px 0 10px; }}
  .teaid-hits-actions button {{
    font: inherit; font-size: 11px; color: {text}; background: {surface};
    border: 1px solid {axis}; border-radius: 5px; padding: 4px 11px; cursor: pointer;
  }}
  .teaid-hits-actions button:hover {{ border-color: {accent}; color: {accent}; }}
  .teaid-hits-scroll {{ overflow-x: auto; max-height: 420px; overflow-y: auto; }}
  .teaid-hits table {{ border-collapse: collapse; width: 100%; }}
  .teaid-hits th {{
    position: sticky; top: 0; background: {paper}; text-align: left;
    font-weight: 500; font-size: 10px; letter-spacing: .06em; text-transform: uppercase;
    color: {muted}; padding: 5px 10px 5px 0; border-bottom: 1px solid {axis};
  }}
  .teaid-hits td {{
    padding: 4px 10px 4px 0; border-bottom: 1px solid {axis}; color: {muted};
    white-space: nowrap; font-variant-numeric: tabular-nums;
  }}
  .teaid-hits td.name {{ color: {text}; }}
  .teaid-hits tr.collapsed td {{ opacity: .62; }}
  .teaid-hits .flag {{ color: {alarm}; }}
</style>
</head>
<body>
<div class="teaid-bar">
  <button class="teaid-reset" id="teaid-reset" type="button">Reset view</button>
  <span class="teaid-hint">drag to pan · scroll to zoom · double-click a panel to autoscale it</span>
</div>
<div class="teaid-plot">{plot}</div>
{hits_table}
<script>
  // Restore every axis to the range the sheet was built with, so one click
  // re-centres all four quadrants on the full consensus. Plotly's own "reset
  // axes" autoscales each panel to its data instead, which leaves the panels
  // on different x-ranges and breaks the vertical alignment between them.
  (function () {{
    var RESET = {reset};
    var gd = document.getElementById({div_id!r});
    var button = document.getElementById('teaid-reset');
    if (!gd || !button) return;
    button.addEventListener('click', function () {{
      Plotly.relayout(gd, RESET);
    }});
  }})();

  // Protein-hit table: export what the panel could not draw.
  (function () {{
    var rows = window.__teaidHits;
    if (!rows || !rows.length) return;
    var head = Object.keys(rows[0]);
    var tsv = [head.join('\t')].concat(
      rows.map(function (r) {{ return head.map(function (k) {{ return r[k]; }}).join('\t'); }})
    ).join('\n');

    var save = document.getElementById('teaid-hits-download');
    if (save) save.addEventListener('click', function () {{
      // A blob URL rather than a data: URI — a large table exceeds what some
      // browsers accept in a URL.
      var url = URL.createObjectURL(new Blob([tsv], {{type: 'text/tab-separated-values'}}));
      var a = document.createElement('a');
      a.href = url; a.download = {tsv_name!r};
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(function () {{ URL.revokeObjectURL(url); }}, 1000);
    }});

    var copy = document.getElementById('teaid-hits-copy');
    if (copy) copy.addEventListener('click', function () {{
      // Clipboard works where a download is blocked, e.g. a sandboxed viewer.
      navigator.clipboard.writeText(tsv).then(
        function () {{ copy.textContent = 'Copied'; setTimeout(function () {{ copy.textContent = 'Copy TSV'; }}, 1600); }},
        function () {{ copy.textContent = 'Copy failed — select the table'; }}
      );
    }});
  }})();
</script>
</body>
</html>
"""


def _hits_table(data: SheetData, theme: Theme, tsv_name: str) -> tuple[str, str]:
    """The expandable table of every protein hit, and its JSON for export.

    Panel 5 draws a readable subset — competitors collapsed, lanes capped — so
    without this the rest exists only in a note saying how many were dropped.
    A curator checking a marginal call needs the actual numbers.
    """
    import html as _html
    import json as _json

    hits = list(data.homology_all or data.homology)
    if not hits:
        return "", "null"

    drawn = {id(h) for h in homology_lanes(list(data.homology))}
    rows = []
    for hit in sorted(hits, key=lambda h: (h.start, h.evalue)):
        start, end = hit.coordinates_1based()
        rows.append({
            "name": getattr(hit, "display_name", hit.query),
            "class": getattr(hit, "source_class", None) or "",
            "accession": hit.query_accession,
            "start": start,
            "end": end,
            "strand": hit.strand,
            "evalue": f"{hit.evalue:.2g}",
            "score": f"{hit.score:.1f}",
            "identity_pct": f"{hit.identity:.1f}",
            "model_coverage_pct": f"{hit.coverage * 100:.0f}",
            "frameshifts": hit.frameshifts,
            "stop_codons": hit.stop_codons,
            "source": getattr(hit, "source", "pfam"),
            "drawn": "yes" if id(hit) in drawn else "no",
        })

    headers = ["name", "class", "accession", "start", "end", "strand", "evalue",
               "score", "identity_pct", "model_coverage_pct", "frameshifts",
               "stop_codons", "source", "drawn"]
    labels = ["hit", "class", "accession", "start", "end", "str", "E-value",
              "score", "% id", "% model", "shifts", "stops", "tier", "drawn"]

    body = []
    for row in rows:
        classes = [] if row["drawn"] == "yes" else ["collapsed"]
        cells = []
        for key in headers:
            value = _html.escape(str(row[key]))
            css = "name" if key == "name" else ""
            if key in ("frameshifts", "stop_codons") and row[key]:
                css = (css + " flag").strip()
            cells.append(f'<td class="{css}">{value}</td>' if css else f"<td>{value}</td>")
        body.append(f'<tr class="{" ".join(classes)}">{"".join(cells)}</tr>')

    n_drawn = sum(1 for r in rows if r["drawn"] == "yes")
    hidden = len(rows) - n_drawn
    summary = (
        f"{len(rows)} protein hit{'s' if len(rows) != 1 else ''}"
        + (f" — {n_drawn} drawn, {hidden} collapsed or capped" if hidden else " — all drawn")
    )

    markup = (
        '<details class="teaid-hits">'
        f"<summary>{_html.escape(summary)}</summary>"
        '<div class="teaid-hits-actions">'
        '<button type="button" id="teaid-hits-download">Download TSV</button>'
        '<button type="button" id="teaid-hits-copy">Copy TSV</button>'
        "</div>"
        '<div class="teaid-hits-scroll"><table><thead><tr>'
        + "".join(f"<th>{_html.escape(l)}</th>" for l in labels)
        + "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></div></details>"
    )
    return markup, _json.dumps(rows)


def write_html(fig: go.Figure, path, data: SheetData, theme: Theme) -> None:
    """Write the interactive sheet, with a reset control the modebar lacks."""
    import json

    div_id = "teaid-sheet"
    plot_div = fig.to_html(
        include_plotlyjs="cdn",
        full_html=False,
        div_id=div_id,
        default_width="100%",
        config={
            "scrollZoom": True,
            "displaylogo": False,
            "responsive": True,
            "toImageButtonOptions": {
                "filename": f"{data.family}.teaid",
                "format": "png",
                "scale": 2,
            },
        },
    )

    reset: dict[str, object] = {}
    for name in fig.layout:
        if not (name.startswith("xaxis") or name.startswith("yaxis")):
            continue
        axis = fig.layout[name]
        # The layout attribute path is the full axis name ('xaxis2.range').
        # Traces reference axes by the short form ('x2'), which is not a valid
        # relayout key -- using it fails silently, leaving the button inert.
        if axis.range is not None:
            reset[f"{name}.range"] = list(axis.range)
        else:
            # Panels whose y-scale depends on the data (divergence, coverage)
            # go back to autoscaling rather than to a range invented here.
            reset[f"{name}.autorange"] = True

    tsv_name = f"{data.family}.protein_hits.tsv"
    hits_markup, hits_json = _hits_table(data, theme, tsv_name)
    if hits_markup:
        plot_div += f"\n<script>window.__teaidHits = {hits_json};</script>"

    path = str(path)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(
            _HTML_TEMPLATE.format(
                title=f"{data.family} — TE-Aid",
                scheme="dark" if theme.is_dark else "light",
                paper=theme.paper,
                surface=theme.surface,
                text=theme.text_primary,
                muted=theme.muted,
                axis=theme.axis,
                accent=theme.base,
                font=FONT_FAMILY,
                plot=plot_div,
                reset=json.dumps(reset),
                div_id=div_id,
                hits_table=hits_markup,
                tsv_name=tsv_name,
                alarm=STATUS_BELOW_FLOOR,
            )
        )


def _runs_below(values: np.ndarray, floor: int) -> list[tuple[int, int]]:
    """Half-open [start, end) spans where ``values`` sits under ``floor``."""
    below = values < floor
    if not below.any():
        return []
    # Difference of the padded mask marks every rising and falling edge.
    edges = np.diff(np.concatenate(([0], below.view(np.int8), [0])))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    return list(zip(starts.tolist(), ends.tolist()))


HELP_MARKER = "?"

# Visible characters per line in a help tooltip. Plotly hover labels do not wrap
# on their own, so an unwrapped paragraph becomes one line as wide as the page
# and covers the panels it is meant to explain.
_HELP_WRAP = 62


def _wrap_help(text: str) -> str:
    """Wrap help text to a readable column width, respecting its markup.

    Explicit ``<br>`` breaks are kept as paragraph boundaries. Line length is
    measured on visible characters only, so a ``<b>`` tag does not push the
    wrap; and breaks are only ever inserted between words, never inside a tag.
    """
    import re as _re

    visible = lambda s: len(_re.sub(r"<[^>]+>", "", s))
    out_paragraphs = []
    for paragraph in text.split("<br><br>"):
        line, lines = "", []
        for word in paragraph.split(" "):
            candidate = f"{line} {word}".strip()
            if line and visible(candidate) > _HELP_WRAP:
                lines.append(line)
                line = word
            else:
                line = candidate
        if line:
            lines.append(line)
        out_paragraphs.append("<br>".join(lines))
    return "<br><br>".join(out_paragraphs)


def _panel_number(text: str | None) -> int | None:
    """The leading panel number of a subplot title, e.g. '3 · Self dot-plot'."""
    if not text:
        return None
    head = text.split("·", 1)[0].strip()
    return int(head) if head.isdigit() else None


def _add_help_marker(fig: go.Figure, title, help_text: str, theme: Theme) -> None:
    """Put a faded '?' after a panel title that explains the panel on hover.

    A separate annotation rather than part of the title, so the hover target is
    the marker itself and the title stays clean. Its position is estimated from
    the title's own width: subplot titles are centred, so the marker sits half a
    title-width to the right of that centre. Paper x spans the plotting area, so
    pixels convert by dividing by its width.
    """
    half_title_px = len(title.text) * _TITLE_CHAR_PX / 2
    offset = (half_title_px + _HELP_GAP_PX) / _PLOT_WIDTH

    fig.add_annotation(
        text=HELP_MARKER,
        x=title.x + offset,
        y=title.y,
        xref=title.xref,
        yref=title.yref,
        xanchor="left",
        yanchor=title.yanchor or "bottom",
        showarrow=False,
        font=dict(size=11, color=theme.muted, family=FONT_FAMILY),
        hovertext=_wrap_help(help_text),
        hoverlabel=dict(
            bgcolor=theme.surface,
            bordercolor=theme.axis,
            font=dict(color=theme.text_primary, size=11, family=FONT_FAMILY),
        ),
        captureevents=True,
    )


def strip_help_markers(fig: go.Figure) -> go.Figure:
    """Remove the '?' markers, for static export.

    They are an interactive affordance: in a PDF or PNG the tooltip cannot open,
    so the marker is a question mark with no answer — and these files are the
    ones that end up in papers.
    """
    fig.layout.annotations = tuple(
        a for a in fig.layout.annotations if a.text != HELP_MARKER
    )
    return fig


def _padded_span(consensus_length: int) -> list[float]:
    """The shared x-range: the consensus plus a margin at each end."""
    pad = consensus_length * X_PAD_FRACTION
    return [-pad, consensus_length + pad]


def _alpha(hex_colour: str, alpha: float) -> str:
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"
