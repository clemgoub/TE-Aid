"""Tests for the sheet: the quadrant grid, the reset control, and static export."""

from __future__ import annotations

import json
import re

import numpy as np
import pytest

from teaid import report, theme
from teaid.analysis import SelfHit
from teaid.homology import ProteinHit
from teaid.orfs import ORF
from teaid.records import Copy, DivergenceKind


def hit(**overrides) -> ProteinHit:
    base = dict(
        query="RVT_1", query_accession="PF00078.33", start=100, end=700, strand="+",
        evalue=1e-40, score=127.0, identity=27.0, hmm_from=1, hmm_to=200,
        hmm_len=200, frameshifts=0, stop_codons=0,
    )
    base.update(overrides)
    return ProteinHit(**base)


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
        fig = report.build_figure(make_data(homology=[hit()]), theme.LIGHT)
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
        # Both axes shrink their domain rather than floating their range, so the
        # box is square while the x-range still matches the other panels.
        assert fig.layout.xaxis3.constrain == "domain"
        assert tuple(y.range) == tuple(report._padded_span(1000))

    def test_every_panel_shares_one_identical_x_range(self):
        """The grid is only comparable if the four quadrants agree exactly. The
        range carries a margin so a feature drawn at a terminus is not clipped:
        without it the built range ends at the consensus length, an arrowhead
        there is half outside the plot, and the only cure is autoscaling that
        one panel — which is what breaks the shared scale."""
        fig = report.build_figure(make_data(), theme.LIGHT)
        expected = tuple(report._padded_span(1000))
        for name in ("xaxis", "xaxis2", "xaxis3", "xaxis4"):
            assert tuple(fig.layout[name].range) == expected, name
        assert expected[0] < 0 and expected[1] > 1000

    def test_sheet_height_is_derived_so_quadrants_are_square(self):
        """A square quadrant is what lets the dot-plot be square *and* the same
        pixel width as the panel above it. With a non-square quadrant you get one
        or the other: either the dot-plot's box narrows (breaking the vertical
        alignment between panel 1 and panel 3) or its data range widens
        (breaking the shared consensus scale)."""
        cell_w = (report.SHEET_WIDTH - report._MARGIN_L - report._MARGIN_R) * (
            1 - report._H_SPACING
        ) / 2
        cell_h = (report.SHEET_HEIGHT - report._MARGIN_T - report._MARGIN_B) * (
            1 - report._V_SPACING
        ) / 2
        assert abs(cell_w - cell_h) < 1.0

        fig = report.build_figure(make_data(), theme.LIGHT)
        assert fig.layout.width == report.SHEET_WIDTH
        assert fig.layout.height == report.SHEET_HEIGHT

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

    def test_x_axes_reset_to_the_shared_consensus_span(self, tmp_path):
        payload = self._reset_payload(tmp_path, make_data())
        expected = report._padded_span(1000)
        for name in ("xaxis", "xaxis2", "xaxis3", "xaxis4"):
            assert payload[f"{name}.range"] == expected

    def test_data_dependent_y_axes_autoscale_rather_than_guess(self, tmp_path):
        payload = self._reset_payload(tmp_path, make_data())
        assert payload["yaxis.autorange"] is True  # divergence
        assert payload["yaxis2.autorange"] is True  # coverage
        assert payload["yaxis3.range"] == report._padded_span(1000)  # square

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

    @staticmethod
    def _arrowheads(fig) -> list[tuple[float, str]]:
        """(tip position, symbol) for every arrowhead in the structure panel."""
        return sorted(
            (t.x[0], str(t.marker.symbol))
            for t in fig.data
            if t.mode == "markers" and str(t.marker.symbol).startswith("triangle")
        )

    def test_inverted_repeat_arms_point_inward(self):
        """The left arm points right and the right arm points left, so an
        inverted pair reads as the palindrome it is."""
        inverted = SelfHit(100, 200, 900, 800, 98.0, 1e-30, 400.0)
        fig = report.build_figure(
            make_data(terminal_repeats={"LTR": [], "TIR": [inverted]}), theme.LIGHT
        )
        assert self._arrowheads(fig) == [(200, "triangle-right"), (800, "triangle-left")]

    def test_inward_orientation_survives_blastn_reporting_either_reciprocal(self):
        """blastn reports a pair in whichever direction it found it, and
        deduplication keeps an arbitrary one of the two. Drawing straight from
        its coordinate order therefore pointed some inverted pairs outward."""
        forward = SelfHit(100, 200, 900, 800, 98.0, 1e-30, 400.0)
        reciprocal = SelfHit(800, 900, 200, 100, 98.0, 1e-30, 400.0)
        drawn = [
            self._arrowheads(
                report.build_figure(
                    make_data(terminal_repeats={"LTR": [], "TIR": [hit]}), theme.LIGHT
                )
            )
            for hit in (forward, reciprocal)
        ]
        assert drawn[0] == drawn[1]
        assert drawn[0] == [(200, "triangle-right"), (800, "triangle-left")]

    def test_direct_repeat_arms_both_point_the_same_way(self):
        direct = SelfHit(100, 200, 800, 900, 98.0, 1e-30, 400.0)
        fig = report.build_figure(
            make_data(terminal_repeats={"LTR": [direct], "TIR": []}), theme.LIGHT
        )
        assert self._arrowheads(fig) == [(200, "triangle-right"), (900, "triangle-right")]

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


class TestProvenance:
    """A sheet outlives the shell that produced it."""

    def test_sources_are_shown_under_the_title(self):
        data = make_data(sources=["seed <b>a.stk</b>", "annotation <b>b.out</b> (rmout)"])
        fig = report.build_figure(data, theme.LIGHT)
        assert "a.stk" in fig.layout.title.text
        assert "b.out" in fig.layout.title.text

    def test_no_provenance_line_when_there_is_nothing_to_say(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        assert "from" not in fig.layout.title.text

    def test_the_line_gets_its_own_room(self):
        """Otherwise it lands on the first row of panel titles."""
        plain = report.build_figure(make_data(), theme.LIGHT)
        sourced = report.build_figure(make_data(sources=["seed <b>a.stk</b>"]), theme.LIGHT)
        assert sourced.layout.margin.t > plain.layout.margin.t
        assert sourced.layout.height > plain.layout.height


class TestPanelHelp:
    @staticmethod
    def _markers(fig):
        return [a for a in fig.layout.annotations if a.text == report.HELP_MARKER]

    def test_every_drawn_panel_gets_a_marker(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        assert len(self._markers(fig)) == 4  # panels 1-4

    def test_seed_qc_panels_get_markers_too(self):
        data = make_data(seed_depth=np.ones(50, dtype=np.int64), seed_sequence_count=1)
        fig = report.build_figure(data, theme.LIGHT)
        assert len(self._markers(fig)) == 6  # panels 1-4, 6, 8

    def test_markers_carry_the_help_text_and_capture_hover(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        for marker in self._markers(fig):
            assert marker.captureevents is True
            assert len(marker.hovertext) > 50

    def test_marker_sits_to_the_right_of_its_title(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        titles = {
            report._panel_number(a.text): a
            for a in fig.layout.annotations
            if report._panel_number(a.text)
        }
        for marker in self._markers(fig):
            nearest = min(titles.values(), key=lambda t: abs(t.x - marker.x))
            assert marker.x > nearest.x

    def test_help_text_is_wrapped(self):
        """Plotly hover labels do not wrap; an unwrapped paragraph covers the
        panels it is meant to explain."""
        for text in report.PANEL_HELP.values():
            wrapped = report._wrap_help(text)
            for line in wrapped.replace("<br><br>", "<br>").split("<br>"):
                assert len(re.sub(r"<[^>]+>", "", line)) <= report._HELP_WRAP + 20

    def test_wrapping_never_breaks_inside_a_tag(self):
        wrapped = report._wrap_help(report.PANEL_HELP[2])
        for line in wrapped.split("<br>"):
            assert line.count("<") == line.count(">")

    def test_markers_are_stripped_for_static_export(self):
        """A question mark in a PDF has no answer, and PDFs go into papers."""
        fig = report.build_figure(make_data(), theme.LIGHT)
        assert self._markers(fig)
        report.strip_help_markers(fig)
        assert not self._markers(fig)
        # The panel titles themselves survive.
        assert any(report._panel_number(a.text) for a in fig.layout.annotations)
