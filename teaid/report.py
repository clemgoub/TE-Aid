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

The dot-plot holds a true 1:1 data aspect via ``scaleanchor`` with
``constrain="domain"``: Plotly shrinks the plotting area to fit the aspect
instead of widening the data range, so the square survives any figure size.

Segments are drawn as a *single* trace per series, with ``None`` separating each
segment, rather than one trace per copy. A family with 10,000 copies otherwise
produces 10,000 Plotly traces, which no browser renders comfortably.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .analysis import SelfHit
from .orfs import ORF
from .records import Annotation, Copy, DivergenceKind
from .theme import FONT_FAMILY, Theme


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
    notes: list[str] = field(default_factory=list)
    # Panel 5 content; empty until the BATH protein row and Dfam nucleotide row
    # land (work order steps 6 and 8).
    homology: list = field(default_factory=list)

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
    titles = [
        "1 · Annotated copies vs divergence",
        "2 · Consensus coverage",
        "3 · Self dot-plot",
        "4 · Structure",
    ]
    if data.has_homology:
        fig = make_subplots(
            rows=3,
            cols=2,
            specs=[[{}, {}], [{"rowspan": 2}, {}], [None, {}]],
            row_heights=[0.5, 0.27, 0.23],
            column_widths=[0.5, 0.5],
            horizontal_spacing=0.10,
            vertical_spacing=0.085,
            subplot_titles=titles + ["5 · Homology evidence"],
        )
        structure_row, homology_row = 2, 3
    else:
        fig = make_subplots(
            rows=2,
            cols=2,
            specs=[[{}, {}], [{}, {}]],
            row_heights=[0.5, 0.5],
            column_widths=[0.5, 0.5],
            horizontal_spacing=0.10,
            vertical_spacing=0.085,
            subplot_titles=titles,
        )
        structure_row, homology_row = 2, None

    _panel_hits(fig, data, theme, row=1, col=1)
    _panel_coverage(fig, data, theme, row=1, col=2)
    _panel_dotplot(fig, data, theme, row=2, col=1)
    _panel_structure(fig, data, theme, row=structure_row, col=2)
    if homology_row is not None:
        _panel_homology(fig, data, theme, row=homology_row, col=2)

    _layout(fig, data, theme, structure_row=structure_row, homology_row=homology_row)
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
    fig.update_yaxes(title_text="copies overlapping", rangemode="tozero", row=row, col=col)
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
    lanes: list[tuple[str, str]] = []  # (tick label, ...)
    y = 0.0
    tickvals: list[float] = []
    ticktext: list[str] = []
    drew_anything = False

    ltr = data.terminal_repeats.get("LTR", [])
    tir = data.terminal_repeats.get("TIR", [])

    if ltr or tir:
        for hits, colour, label in (
            (ltr, theme.base, "LTR-like"),
            (tir, theme.inverted, "TIR-like"),
        ):
            if not hits:
                continue
            xs: list[float | None] = []
            ys: list[float | None] = []
            for hit in hits:
                s_lo, s_hi = sorted((hit.s_start, hit.s_end))
                # Both arms of the terminal repeat, on one lane.
                xs.extend((hit.q_start, hit.q_end, None, s_lo, s_hi, None))
                ys.extend((y, y, None, y, y, None))
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="lines",
                    name=f"{label} (n={len(hits)})",
                    line=dict(color=colour, width=9),
                    opacity=0.85,
                    hovertemplate=f"{label} arm<br>%{{x:,.0f}} bp<extra></extra>",
                ),
                row=row,
                col=col,
            )
        tickvals.append(y)
        ticktext.append("terminal<br>repeats")
        y -= 1.0
        drew_anything = True

    for orf in data.orfs:
        tip = orf.end if orf.strand == "+" else orf.start
        fig.add_trace(
            go.Scatter(
                x=[orf.start, orf.end],
                y=[y, y],
                mode="lines",
                line=dict(color=theme.highlight, width=9),
                opacity=0.85,
                showlegend=False,
                hovertemplate=(
                    f"ORF {orf.coordinates_1based()[0]:,}-{orf.coordinates_1based()[1]:,} "
                    f"({orf.strand})<br>{orf.length_aa:,} aa<extra></extra>"
                ),
            ),
            row=row,
            col=col,
        )
        fig.add_trace(
            go.Scatter(
                x=[tip],
                y=[y],
                mode="markers",
                marker=dict(
                    symbol="triangle-right" if orf.strand == "+" else "triangle-left",
                    size=11,
                    color=theme.highlight,
                ),
                showlegend=False,
                hoverinfo="skip",
            ),
            row=row,
            col=col,
        )
        tickvals.append(y)
        ticktext.append(f"{orf.length_aa:,} aa {orf.strand}")
        y -= 1.0
        drew_anything = True

    if data.orfs:
        # One legend entry for the ORF track as a whole; the per-ORF traces stay
        # out of the legend so it does not fill with one row per frame.
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines",
                line=dict(color=theme.highlight, width=9),
                name=f"ORF ≥{data.orf_min_size} bp (n={len(data.orfs)})",
            ),
            row=row,
            col=col,
        )

    if not drew_anything:
        _empty_note(
            fig,
            f"no terminal repeats and no ORF ≥ {data.orf_min_size} bp",
            theme,
            row,
            col,
        )

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


def _panel_homology(fig: go.Figure, data: SheetData, theme: Theme, row: int, col: int) -> None:
    """Panel 5: best protein and nucleotide hits (work order steps 6 and 8)."""
    _empty_note(fig, "homology evidence pending", theme, row, col)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=row, col=col)
    fig.update_xaxes(title_text="consensus (bp)", row=row, col=col)


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
) -> None:
    subtitle_bits = [f"{data.consensus_length:,} bp", f"{data.n_copies:,} copies"]
    if data.full_length:
        subtitle_bits.append(
            f"{len(data.full_length):,} full length ≥{data.full_length_threshold:.0%}"
        )
    if data.orfs:
        subtitle_bits.append(f"{len(data.orfs)} ORF{'s' if len(data.orfs) != 1 else ''}")
    if data.class_label:
        subtitle_bits.append(data.class_label)

    fig.update_layout(
        title=dict(
            text=(
                f"<b>{data.family}</b><br>"
                f"<span style='font-size:12px;color:{theme.text_secondary}'>"
                f"{' &nbsp;·&nbsp; '.join(subtitle_bits)}</span>"
            ),
            x=0.012,
            xanchor="left",
            yref="container",
            y=0.985,
            yanchor="top",
            font=dict(size=18, color=theme.text_primary, family=FONT_FAMILY),
        ),
        width=1180,
        height=1180,
        paper_bgcolor=theme.paper,
        plot_bgcolor=theme.surface,
        font=dict(family=FONT_FAMILY, color=theme.text_secondary, size=11),
        # The legend sits below the grid, centred. In a 2x2 the top-right
        # subplot title occupies the same band as a top-anchored legend, so
        # anywhere at the top collides with it by construction.
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.055,
            xanchor="center",
            x=0.5,
            bgcolor="rgba(0,0,0,0)",
            font=dict(color=theme.text_secondary, size=11),
        ),
        margin=dict(l=68, r=26, t=96, b=104),
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
        range=[0, data.consensus_length],
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
    fig.update_yaxes(
        range=[0, data.consensus_length],
        scaleanchor=dotplot_axis.xaxis.anchor.replace("y", "x"),
        scaleratio=1,
        constrain="domain",
        row=2,
        col=1,
    )

    # Panels 4 and 5 share one consensus x-axis so features line up vertically.
    if homology_row is not None:
        structure_axis = fig.get_subplot(structure_row, 2)
        fig.update_xaxes(
            matches=structure_axis.xaxis.anchor.replace("y", "x"),
            row=homology_row,
            col=2,
        )
        fig.update_xaxes(title_text=None, showticklabels=False, row=structure_row, col=2)

    for annotation in fig.layout.annotations:
        if annotation.text in {
            "1 · Annotated copies vs divergence",
            "2 · Consensus coverage",
            "3 · Self dot-plot",
            "4 · Structure",
            "5 · Homology evidence",
        }:
            annotation.update(
                font=dict(size=12, color=theme.text_primary, family=FONT_FAMILY),
                x=annotation.x,
                xanchor="center",
            )


def _alpha(hex_colour: str, alpha: float) -> str:
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"
