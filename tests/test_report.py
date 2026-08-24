"""Tests for the sheet: the quadrant grid, the reset control, and static export."""

from __future__ import annotations

import json
import re

import numpy as np
import pytest

from teaid import report, theme
from teaid.analysis import SelfHit
from teaid.orfs import ORF
from teaid.records import Copy, DivergenceKind


def make_data(**overrides) -> report.SheetData:
    copies = [
        Copy(
            chrom="chr1",
            start=0,
            end=900,
            strand="+",
            family="fam",
            divergence=5.0,
            divergence_kind=DivergenceKind.CONSENSUS,
            consensus_start=0,
            consensus_end=900,
        )
    ]
    base = dict(
        family="fam",
        consensus_length=1000,
        copies=copies,
        full_length=copies,
        coverage=np.ones(1000, dtype=np.int64),
    )
    base.update(overrides)
    return report.SheetData(**base)


class TestQuadrantGrid:
    def test_four_quadrants_without_homology(self):
        """The 2x2 grid is the layout; panel 4 takes the whole bottom-right
        quadrant until there is homology evidence to split it with."""
        fig = report.build_figure(make_data(), theme.LIGHT)
        axes = [k for k in fig.layout if re.fullmatch(r"xaxis\d*", k)]
        assert len(axes) == 4

    def test_five_panels_once_homology_exists(self):
        fig = report.build_figure(make_data(homology=["placeholder"]), theme.LIGHT)
        axes = [k for k in fig.layout if re.fullmatch(r"xaxis\d*", k)]
        assert len(axes) == 5

    def test_dotplot_axes_are_locked_to_a_square_data_aspect(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        y = fig.layout.yaxis3
        assert y.scaleanchor == "x3"
        assert y.scaleratio == 1
        # 'domain' shrinks the plotting box to honour the aspect; the default
        # would widen the data range instead, and because the quadrants share a
        # consensus scale that would desynchronise them.
        assert y.constrain == "domain"
        assert tuple(y.range) == (0, 1000)

    def test_every_panel_spans_the_full_consensus(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        for name in ("xaxis", "xaxis2", "xaxis3", "xaxis4"):
            assert tuple(fig.layout[name].range) == (0, 1000), name

    def test_empty_panels_keep_their_axes(self):
        """Plotly renders a subplot's axes only if a trace references them, so an
        empty quadrant would otherwise disappear and take the grid with it."""
        fig = report.build_figure(make_data(), theme.LIGHT)  # no self-hits, no ORFs
        referenced = {t.xaxis or "x" for t in fig.data}
        assert {"x3", "x4"} <= referenced


class TestResetControl:
    def _reset_payload(self, tmp_path, data) -> dict:
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        match = re.search(r"var RESET = (\{.*?\});", path.read_text(), re.S)
        assert match, "reset payload missing from the page"
        return json.loads(match.group(1))

    def test_keys_are_full_layout_paths(self, tmp_path):
        """Plotly.relayout takes 'xaxis2.range'. The short form 'x2.range' is how
        a *trace* references an axis; as a relayout key it is silently ignored,
        which leaves the button inert while every static check still passes."""
        payload = self._reset_payload(tmp_path, make_data())
        assert payload
        for key in payload:
            assert re.fullmatch(r"[xy]axis\d*\.(range|autorange)", key), key

    def test_every_axis_is_covered(self, tmp_path):
        payload = self._reset_payload(tmp_path, make_data())
        axes = {k.split(".")[0] for k in payload}
        assert axes == {
            "xaxis", "xaxis2", "xaxis3", "xaxis4",
            "yaxis", "yaxis2", "yaxis3", "yaxis4",
        }

    def test_x_axes_reset_to_the_consensus_span(self, tmp_path):
        payload = self._reset_payload(tmp_path, make_data())
        for name in ("xaxis", "xaxis2", "xaxis3", "xaxis4"):
            assert payload[f"{name}.range"] == [0, 1000]

    def test_data_dependent_y_axes_autoscale_rather_than_guess(self, tmp_path):
        payload = self._reset_payload(tmp_path, make_data())
        assert payload["yaxis.autorange"] is True  # divergence
        assert payload["yaxis2.autorange"] is True  # coverage
        assert payload["yaxis3.range"] == [0, 1000]  # dot-plot, square

    def test_page_wires_the_button_to_the_plot_div(self, tmp_path):
        data = make_data()
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        html = path.read_text()
        assert 'id="teaid-reset"' in html
        assert 'id="teaid-sheet"' in html
        assert "Plotly.relayout" in html


class TestStructurePanel:
    def test_terminal_repeats_draw_both_arms_as_arrows(self):
        """v1's idea, kept: direct repeats point the same way, inverted repeats
        point at each other, so orientation is readable from shape alone."""
        direct = SelfHit(100, 200, 800, 900, 98.0, 1e-30, 400.0)
        fig = report.build_figure(
            make_data(terminal_repeats={"LTR": [direct], "TIR": []}), theme.LIGHT
        )
        arrowheads = [
            t for t in fig.data
            if t.mode == "markers" and str(t.marker.symbol).startswith("triangle")
        ]
        assert len(arrowheads) == 2
        assert {str(t.marker.symbol) for t in arrowheads} == {"triangle-right"}

    def test_inverted_repeat_arrows_point_at_each_other(self):
        inverted = SelfHit(100, 200, 900, 800, 98.0, 1e-30, 400.0)
        fig = report.build_figure(
            make_data(terminal_repeats={"LTR": [], "TIR": [inverted]}), theme.LIGHT
        )
        symbols = [
            str(t.marker.symbol) for t in fig.data
            if t.mode == "markers" and str(t.marker.symbol).startswith("triangle")
        ]
        assert sorted(symbols) == ["triangle-left", "triangle-right"]

    def test_reverse_strand_orf_points_left(self):
        fig = report.build_figure(
            make_data(orfs=[ORF(100, 700, "-", "M" * 200)]), theme.LIGHT
        )
        symbols = [
            str(t.marker.symbol) for t in fig.data
            if t.mode == "markers" and str(t.marker.symbol).startswith("triangle")
        ]
        assert symbols == ["triangle-left"]

    def test_lane_count_is_capped_and_the_cut_is_stated(self):
        """A silent truncation reads as 'this is everything'."""
        hits = [
            SelfHit(i * 10, i * 10 + 50, 800, 900 - i, 98.0, 1e-30, 500.0 - i)
            for i in range(report.MAX_REPEAT_LANES + 3)
        ]
        data = make_data(terminal_repeats={"LTR": hits, "TIR": []})
        fig = report.build_figure(data, theme.LIGHT)
        # Each lane draws two arms, so count distinct lane positions.
        lanes = {
            t.y[0] for t in fig.data if t.mode == "lines" and t.xaxis == "x4"
        }
        assert len(lanes) == report.MAX_REPEAT_LANES
        assert any("not drawn" in note for note in data.notes)
        assert any("not drawn" in (a.text or "") for a in fig.layout.annotations)

    def test_lane_labels_name_the_observation_not_a_classification(self):
        """'LTR' on a lane would assert a class; on a repetitive consensus the
        same signature is a tandem unit."""
        fig = report.build_figure(
            make_data(
                terminal_repeats={
                    "LTR": [SelfHit(0, 100, 900, 1000, 99.0, 0.0, 500.0)],
                    "TIR": [SelfHit(0, 100, 1000, 900, 99.0, 0.0, 500.0)],
                }
            ),
            theme.LIGHT,
        )
        assert set(fig.layout.yaxis4.ticktext) == {"direct", "inverted"}


class TestStaticExport:
    @pytest.mark.parametrize("fmt", ["png", "pdf", "svg"])
    def test_writes_a_non_empty_file(self, tmp_path, fmt):
        pytest.importorskip("kaleido")
        fig = report.build_figure(make_data(), theme.LIGHT)
        path = tmp_path / f"sheet.{fmt}"
        fig.write_image(path, scale=2)
        assert path.stat().st_size > 5_000

    def test_pdf_is_a_real_pdf(self, tmp_path):
        pytest.importorskip("kaleido")
        fig = report.build_figure(make_data(), theme.LIGHT)
        path = tmp_path / "sheet.pdf"
        fig.write_image(path)
        assert path.read_bytes().startswith(b"%PDF")


class TestThemes:
    @pytest.mark.parametrize("name", ["light", "dark"])
    def test_both_themes_paint_an_explicit_background(self, name):
        chosen = theme.get(name)
        fig = report.build_figure(make_data(), chosen)
        assert fig.layout.paper_bgcolor == chosen.paper
        assert fig.layout.plot_bgcolor == chosen.surface

    def test_unknown_theme_is_rejected(self):
        with pytest.raises(ValueError):
            theme.get("solarized")
