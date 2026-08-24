"""Assembly of the evidence sheet: panels 1-3 over a shared consensus axis.

Panel 1 — annotated copies vs divergence
Panel 2 — consensus coverage pileup
Panel 3 — self dot-plot

All three share one consensus-coordinate x-axis so features line up vertically;
panels 4 and 5 will slot into the same stack.

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
from .records import Annotation, Copy, DivergenceKind
from .theme import FONT_FAMILY, Theme


@dataclass(slots=True)
class SheetData:
    """Everything panels 1-3 need, already computed."""

    family: str
    consensus_length: int
    copies: list[Copy]
    full_length: list[Copy]
    coverage: np.ndarray
    self_hits: list[SelfHit] = field(default_factory=list)
    class_label: str | None = None
    full_length_threshold: float = 0.9
    source_format: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def n_copies(self) -> int:
        return len(self.copies)

    @property
    def divergence_kind(self) -> DivergenceKind:
        kinds = {c.divergence_kind for c in self.copies if c.has_divergence}
        if not kinds:
            return DivergenceKind.NONE
        if len(kinds) == 1:
            return next(iter(kinds))
        # Mixed producers in one file: refuse to average two different
        # quantities onto one axis and say so instead.
        return DivergenceKind.CONSENSUS


def _segments(
    starts, ends, ys_start, ys_end
) -> tuple[list[float | None], list[float | None]]:
    """Interleave segment endpoints with None breaks for a single Plotly trace."""
    xs: list[float | None] = []
    ys: list[float | None] = []
    for x0, x1, y0, y1 in zip(starts, ends, ys_start, ys_end):
        xs.extend((x0, x1, None))
        ys.extend((y0, y1, None))
    return xs, ys


def build_figure(data: SheetData, theme: Theme) -> go.Figure:
    """Compose the three-panel sheet."""
    fig = make_subplots(
        rows=3,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.07,
        row_heights=[0.32, 0.20, 0.48],
        subplot_titles=(
            "Annotated copies vs divergence",
            "Consensus coverage",
            "Self dot-plot",
        ),
    )

    _panel_hits(fig, data, theme, row=1)
    _panel_coverage(fig, data, theme, row=2)
    _panel_dotplot(fig, data, theme, row=3)

    _layout(fig, data, theme)
    return fig


def _panel_hits(fig: go.Figure, data: SheetData, theme: Theme, row: int) -> None:
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
        (full, theme.highlight, f"full length (≥{data.full_length_threshold:.0%})", 2.0),
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
                legendgroup=label,
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
            col=1,
        )

    kind = data.divergence_kind
    title = kind.axis_label
    if omitted:
        title += f"<br><span style='font-size:11px'>{omitted:,} of {data.n_copies:,} copies omitted: no divergence reported</span>"
    fig.update_yaxes(title_text=title, row=row, col=1)

    if not plotted:
        fig.add_annotation(
            text="no copy carries a divergence value",
            xref=f"x{row if row > 1 else ''} domain",
            yref=f"y{row if row > 1 else ''} domain",
            x=0.5,
            y=0.5,
            showarrow=False,
            font=dict(color=theme.text_secondary, size=13, family=FONT_FAMILY),
        )


def _panel_coverage(fig: go.Figure, data: SheetData, theme: Theme, row: int) -> None:
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
        col=1,
    )
    fig.update_yaxes(title_text="copies overlapping", rangemode="tozero", row=row, col=1)


def _panel_dotplot(fig: go.Figure, data: SheetData, theme: Theme, row: int) -> None:
    """Panel 3: the consensus against itself.

    Same-strand off-diagonal matches near both termini suggest LTRs;
    opposite-strand matches suggest TIRs. Both are suggestions for the curator,
    never assertions - a segmental duplication produces the same signature.
    """
    if not data.self_hits:
        fig.add_annotation(
            text="self dot-plot unavailable (blastn not run)",
            xref=f"x{row} domain",
            yref=f"y{row} domain",
            x=0.5,
            y=0.5,
            showarrow=False,
            font=dict(color=theme.text_secondary, size=13, family=FONT_FAMILY),
        )
        fig.update_yaxes(title_text="consensus (bp)", row=row, col=1)
        return

    direct = [h for h in data.self_hits if not h.is_reverse]
    inverted = [h for h in data.self_hits if h.is_reverse]

    for subset, colour, label in (
        (direct, theme.base, "same strand (direct)"),
        (inverted, theme.inverted, "opposite strand (inverted)"),
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
            col=1,
        )

    # No scaleanchor here. A square dot-plot would force its own x-range, and
    # because all three panels share one consensus axis that range would
    # propagate upward and desynchronise panels 1 and 2. The shared axis is the
    # point of the sheet - features must line up vertically - so the dot-plot
    # gives up its square aspect rather than the stack giving up alignment.
    fig.update_yaxes(
        title_text="consensus (bp)",
        range=[0, data.consensus_length],
        row=row,
        col=1,
    )


def _layout(fig: go.Figure, data: SheetData, theme: Theme) -> None:
    subtitle_bits = [f"{data.consensus_length:,} bp", f"{data.n_copies:,} annotated copies"]
    if data.full_length:
        subtitle_bits.append(
            f"{len(data.full_length):,} full length "
            f"(≥{data.full_length_threshold:.0%} of consensus)"
        )
    if data.class_label:
        subtitle_bits.append(data.class_label)

    fig.update_layout(
        title=dict(
            text=(
                f"<b>{data.family}</b><br>"
                f"<span style='font-size:13px;color:{theme.text_secondary}'>"
                f"{' &nbsp;·&nbsp; '.join(subtitle_bits)}</span>"
            ),
            x=0.012,
            xanchor="left",
            yref="container",
            y=0.985,
            yanchor="top",
            font=dict(size=19, color=theme.text_primary, family=FONT_FAMILY),
        ),
        height=1000,
        paper_bgcolor=theme.paper,
        plot_bgcolor=theme.surface,
        font=dict(family=FONT_FAMILY, color=theme.text_secondary, size=12),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.012,
            xanchor="right",
            x=1,
            bgcolor="rgba(0,0,0,0)",
            font=dict(color=theme.text_secondary),
        ),
        # The top margin holds two stacked rows of chrome: the two-line title,
        # then the horizontal legend beneath it. Sized for both so neither the
        # legend nor the first subplot title overlaps the subtitle.
        margin=dict(l=78, r=32, t=126, b=58),
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
        tickfont=dict(color=theme.muted),
        range=[0, data.consensus_length],
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor=theme.grid,
        gridwidth=1,
        zeroline=False,
        linecolor=theme.axis,
        tickcolor=theme.axis,
        tickfont=dict(color=theme.muted),
        title_font=dict(color=theme.text_secondary, size=12),
    )
    fig.update_xaxes(title_text="consensus position (bp)", row=3, col=1)

    for annotation in fig.layout.annotations[:3]:
        annotation.update(
            font=dict(size=13, color=theme.text_primary, family=FONT_FAMILY),
            x=0,
            xanchor="left",
        )


def _alpha(hex_colour: str, alpha: float) -> str:
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"
