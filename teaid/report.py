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

# Terminal-repeat lanes drawn per type before the list is cut. The structure
# quadrant is a quarter of the sheet; past a handful of lanes the arrows are too
# thin to read, and the dot-plot already shows every self-match.
MAX_REPEAT_LANES = 4


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
            # blastn always reports the query interval ascending, so the query
            # arm points right; the subject arm follows its own coordinate
            # order, which is descending exactly when the match is inverted.
            arm_hover = (
                f"{label}<br>arms %{{x:,.0f}} bp<br>identity {hit.identity:.1f}%"
            )
            _arrow(
                fig, hit.q_start, hit.q_end, y, colour, arm_hover, row, col,
                show_legend=label not in seen_types,
                legend_name=f"{label} (n={len(all_hits)})",
            )
            seen_types.add(label)
            _arrow(fig, hit.s_start, hit.s_end, y, colour, arm_hover, row, col,
                   legend_name=f"{label} (n={len(all_hits)})")
            tickvals.append(y)
            ticktext.append(lane_label)
            y -= 1.0
            drew_anything = True

    for index, orf in enumerate(data.orfs):
        # Drawn tail-to-head so the arrowhead sits at the ORF's 3' end,
        # pointing the way it is translated.
        tail, head = (orf.start, orf.end) if orf.strand == "+" else (orf.end, orf.start)
        _arrow(
            fig,
            tail,
            head,
            y,
            theme.highlight,
            f"ORF {orf.coordinates_1based()[0]:,}-{orf.coordinates_1based()[1]:,} "
            f"({orf.strand})<br>{orf.length_aa:,} aa",
            row,
            col,
            show_legend=index == 0,
            legend_name=f"ORF ≥{data.orf_min_size} bp (n={len(data.orfs)})",
        )
        tickvals.append(y)
        ticktext.append(f"{orf.length_aa:,} aa {orf.strand}")
        y -= 1.0
        drew_anything = True

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
        margin=dict(l=68, r=26, t=96, b=128),
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
            y=-0.098,
            xanchor="left",
            yanchor="top",
            showarrow=False,
            font=dict(size=10, color=theme.muted, family=FONT_FAMILY),
        )

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


_HTML_TEMPLATE = """<!doctype html>
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
</style>
</head>
<body>
<div class="teaid-bar">
  <button class="teaid-reset" id="teaid-reset" type="button">Reset view</button>
  <span class="teaid-hint">drag to pan · scroll to zoom · double-click a panel to autoscale it</span>
</div>
<div class="teaid-plot">{plot}</div>
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
</script>
</body>
</html>
"""


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
            )
        )


def _alpha(hex_colour: str, alpha: float) -> str:
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"
