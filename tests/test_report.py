"""Tests for the sheet: the quadrant grid, the reset control, and static export."""

from __future__ import annotations

import base64
import json
import re

import numpy as np
import pytest

from teaid import report, theme
from teaid.analysis import SelfHit
from teaid.homology import ProteinHit
from teaid.orfs import ORF
from teaid.records import Copy, DivergenceKind, Insertion


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
    insertions = [Insertion([c]) for c in copies]
    base = dict(
        family="fam",
        consensus_length=1000,
        insertions=insertions,
        full_length=insertions,
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


class TestTheLogo:
    def test_the_sheet_carries_the_logo_inline(self, tmp_path):
        """Users move a sheet around as one file, so the logo is embedded as a
        data URI. A relative link would break once the file leaves its folder."""
        data = make_data()
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        page = path.read_text()
        match = re.search(r'<img class="teaid-logo" src="data:image/png;base64,([A-Za-z0-9+/=]+)"', page)
        assert match, "logo missing from the toolbar"
        assert base64.b64decode(match.group(1))[:8] == b"\x89PNG\r\n\x1a\n"


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

    def test_orfs_are_no_longer_in_this_panel(self):
        """Panel 4 holds only what the dot-plot above it implies; ORFs moved to
        panel 5, where they can be drawn with the domains that sit in them."""
        fig = report.build_figure(
            make_data(orfs=[ORF(100, 700, "-", "M" * 200)]), theme.LIGHT
        )
        assert "aa" not in " ".join(fig.layout.yaxis4.ticktext or [])

    def test_lane_count_is_capped_and_the_cut_is_stated(self):
        """A silent truncation reads as 'this is everything'."""
        hits = [
            SelfHit(i * 10, i * 10 + 50, 800, 900 - i, 98.0, 1e-30, 500.0 - i)
            for i in range(report.MAX_REPEAT_LANES + 3)
        ]
        data = make_data(terminal_repeats={"LTR": hits, "TIR": []})
        fig = report.build_figure(data, theme.LIGHT)
        # Each lane draws two arms; the consensus-span rail is not a lane.
        lanes = {
            t.y[0] for t in fig.data
            if t.mode == "lines" and t.xaxis == "x4" and float(t.y[0]) == int(float(t.y[0]))
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
        assert {"direct", "inverted"} <= set(fig.layout.yaxis4.ticktext)


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
        # The cards paint the surface; the plot background is clear so they show.
        assert fig.layout.plot_bgcolor == "rgba(0,0,0,0)"
        cards = [s for s in fig.layout.shapes if s.name == "quadrant"]
        assert cards and all(s.fillcolor == chosen.surface for s in cards)

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


class TestNotesLine:
    """The notes say what the sheet left out, so the page must not cut one off."""

    def test_a_long_line_wraps_and_the_margin_grows_to_hold_it(self):
        words = " ".join(f"word{i}" for i in range(120))
        short = report.build_figure(make_data(notes=["one short note"]), theme.LIGHT)
        wrapped = report.build_figure(make_data(notes=[words]), theme.LIGHT)
        (text,) = [a.text for a in wrapped.layout.annotations if a.text and "word0" in a.text]
        lines = text.split("<br>")
        assert len(lines) > 1
        assert all(len(line) * report._NOTE_CHAR_PX <= report._PLOT_WIDTH for line in lines)
        assert " ".join(lines) == words
        extra = (len(lines) - 1) * report._NOTE_LINE_PX
        assert wrapped.layout.margin.b == short.layout.margin.b + extra
        assert wrapped.layout.height == short.layout.height + extra


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


class TestPanelTitleOrdering:
    """make_subplots walks the grid row-major over non-None specs, so a title
    list one entry short shifts every later title onto the wrong panel — and
    does it silently."""

    @staticmethod
    def _titles(fig):
        return [a.text for a in fig.layout.annotations if report._panel_number(a.text)]

    def test_default_sheet(self):
        assert self._titles(report.build_figure(make_data(), theme.LIGHT)) == [
            "1 · Annotated copies vs divergence", "2 · Consensus coverage",
            "3 · Self dot-plot", "4 · Structure",
        ]

    def test_with_homology(self):
        fig = report.build_figure(make_data(homology=[hit()]), theme.LIGHT)
        assert self._titles(fig)[-1] == "5 · ORFs and protein homology"

    def test_with_seed_qc(self):
        data = make_data(seed_depth=np.ones(50, dtype=np.int64), seed_sequence_count=1)
        assert self._titles(report.build_figure(data, theme.LIGHT))[-2:] == [
            "6 · Seed depth", "8 · Expected class",
        ]

    def test_with_both_every_panel_keeps_its_own_title(self):
        data = make_data(
            homology=[hit()], seed_depth=np.ones(50, dtype=np.int64), seed_sequence_count=1
        )
        fig = report.build_figure(data, theme.LIGHT)
        assert self._titles(fig) == [
            "1 · Annotated copies vs divergence", "2 · Consensus coverage",
            "3 · Self dot-plot", "4 · Structure", "5 · ORFs and protein homology",
            "6 · Seed depth", "8 · Expected class",
        ]

    def test_a_title_for_every_drawn_panel(self):
        data = make_data(
            homology=[hit()], seed_depth=np.ones(50, dtype=np.int64), seed_sequence_count=1
        )
        fig = report.build_figure(data, theme.LIGHT)
        axes = [k for k in fig.layout if re.fullmatch(r"xaxis\d*", k)]
        assert len(self._titles(fig)) == len(axes) == 7


def _domain_names(fig) -> list[str]:
    """The names panel 5 writes above its arrows, in drawing order."""
    for trace in fig.data:
        if trace.meta and trace.meta.get("teaidDomainName"):
            return list(trace.text)
    return []


class TestHomologyPanel:
    def test_repeatpeps_names_are_shortened_for_the_panel(self):
        """A tier-2 model is named after the whole FASTA header; the class after
        the '#' is redundant with the panel and runs far past the arrow."""
        fig = report.build_figure(
            make_data(homology=[hit(query="ACROBAT1_tnp#DNA/PiggyBac")]), theme.LIGHT
        )
        assert "ACROBAT1_tnp" in _domain_names(fig)[0]
        assert "#DNA/PiggyBac" not in _domain_names(fig)[0]

    def test_a_disrupted_hit_is_counted_on_the_sheet_not_only_in_hover(self):
        """The counts, not a bare mark: 2 frameshifts and 4 is a different
        finding from 1 and 0, and a static export has no hover."""
        fig = report.build_figure(
            make_data(homology=[hit(frameshifts=2, stop_codons=1)]), theme.LIGHT
        )
        name = _domain_names(fig)[0]
        assert "2fs" in name and "1⊗" in name

    def test_intact_hits_carry_no_mark(self):
        fig = report.build_figure(make_data(homology=[hit()]), theme.LIGHT)
        name = _domain_names(fig)[0]
        assert "fs" not in name and "⊗" not in name

    def test_a_row_with_domains_keeps_its_tick_gutter_empty(self):
        """Long model names ran out of the panel's card from the tick gutter,
        so the names sit above their arrows and only a bare ORF row has a
        gutter label."""
        orf = ORF(0, 900, "+", "M" * 300)
        fig = report.build_figure(
            make_data(orfs=[orf, ORF(1000, 1600, "-", "M" * 200)],
                      homology=[hit(start=100, end=400)]), theme.LIGHT
        )
        ticks = list(fig.layout.yaxis5.ticktext)
        assert "" in ticks
        assert any(t.endswith("aa -") for t in ticks)

    def test_lane_count_is_capped_and_stated(self):
        hits = [hit(query=f"D{i}", start=i * 400, end=i * 400 + 200, score=100.0 - i)
                for i in range(report.MAX_HOMOLOGY_LANES + 4)]
        data = make_data(homology=hits)
        fig = report.build_figure(data, theme.LIGHT)
        assert len(fig.layout.yaxis5.ticktext) == report.MAX_HOMOLOGY_LANES
        assert any("not drawn" in n for n in data.notes)


class TestHomologyLaneSelection:
    """One selection function, so the panel and the hit table cannot disagree
    about which hits were drawn."""

    @staticmethod
    def _hits(n, source="pfam", **over):
        return [hit(query=f"D{i}", start=i * 300, end=i * 300 + 200,
                    evalue=10.0 ** -(n - i), score=float(i), source=source, **over)
                for i in range(n)]

    def test_everything_is_drawn_when_it_fits(self):
        hits = self._hits(5)
        assert report.homology_lanes(hits) == hits

    def test_the_cap_keeps_the_strongest_not_the_leftmost(self):
        """Slicing the position-ordered list would keep the leftmost and drop
        the strongest, while the note claims the opposite."""
        hits = self._hits(report.MAX_HOMOLOGY_LANES + 6)
        shown = report.homology_lanes(hits)
        assert len(shown) == report.MAX_HOMOLOGY_LANES
        strongest = min(hits, key=lambda h: h.evalue)
        assert strongest in shown

    def test_output_is_ordered_along_the_consensus(self):
        shown = report.homology_lanes(self._hits(report.MAX_HOMOLOGY_LANES + 4))
        assert [h.start for h in shown] == sorted(h.start for h in shown)

    def test_neither_source_crowds_the_other_out(self):
        """A blastp bitscore dwarfs a profile-HMM bit score, so a shared ranking
        lets tier 3 take every lane."""
        pfam = self._hits(20, source="pfam")
        blast = [hit(query=f"E{i}", start=i * 300 + 50, end=i * 300 + 250,
                     evalue=1e-99, score=2000.0, source="repeatpeps-blastp")
                 for i in range(20)]
        shown = report.homology_lanes(pfam + blast)
        sources = {getattr(h, "source", "pfam") for h in shown}
        assert sources == {"pfam", "repeatpeps-blastp"}

    def test_tiers_one_and_two_are_ranked_together(self):
        """Tier-1 and tier-2 hits both have profile-HMM bit scores, so
        homology_lanes draws the stronger tier-2 hits instead of giving half the
        lanes to tier 1 by round-robin."""
        n = report.MAX_HOMOLOGY_LANES
        pfam = self._hits(n, source="pfam")
        gaps = [hit(query=f"G{i}", start=i * 300 + 50, end=i * 300 + 250,
                    evalue=1e-50, score=300.0, source="repeatpeps-gaps")
                for i in range(n)]
        shown = report.homology_lanes(pfam + gaps)
        assert {h.source for h in shown} == {"repeatpeps-gaps"}

    def test_the_table_agrees_with_the_panel(self, tmp_path):
        hits = self._hits(report.MAX_HOMOLOGY_LANES + 5)
        data = make_data(homology=hits, homology_all=hits)
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        html = path.read_text()
        drawn_in_table = html.count('<tr class="">')
        assert drawn_in_table == len(fig.layout.yaxis5.ticktext)


class TestHitsTable:
    def test_lists_every_hit_including_collapsed_ones(self, tmp_path):
        drawn = hit(query="kept", start=0, end=300)
        collapsed = hit(query="dropped", start=10, end=290, score=1.0)
        data = make_data(homology=[drawn], homology_all=[drawn, collapsed])
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        html = path.read_text()
        assert "kept" in html and "dropped" in html
        assert "2 protein hits" in html
        assert 'class="collapsed"' in html

    def test_offers_export(self, tmp_path):
        data = make_data(homology=[hit()], homology_all=[hit()])
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        html = path.read_text()
        assert "teaid-hits-download" in html
        assert "teaid-hits-copy" in html
        assert "__teaidHits" in html
        assert "fam.protein_hits.tsv" in html

    def test_no_table_when_there_are_no_hits(self, tmp_path):
        data = make_data()
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        # The stylesheet always carries the class; the element must not exist.
        assert '<details class="teaid-hits">' not in path.read_text()

    def test_frameshifts_are_flagged_in_the_table(self, tmp_path):
        broken = hit(frameshifts=3, stop_codons=1)
        data = make_data(homology=[broken], homology_all=[broken])
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        assert "flag" in path.read_text()


class TestTheInlineScriptParses:
    """The page's one <script> block holds both the reset handler and the TSV
    export. A syntax error anywhere in it silently disables both, and nothing
    server-side notices -- which is exactly what happened: ``_HTML_TEMPLATE``
    was a non-raw Python string, so the JS escapes ``'\\t'`` and ``'\\n'`` in the
    TSV builder became a real tab and a real newline, and a literal newline
    inside a JS string literal is an unterminated-string SyntaxError.
    """

    def test_the_template_holds_no_raw_control_characters(self):
        """A literal tab or a stray carriage return in the template can only
        arrive by Python interpreting an escape meant for the browser."""
        assert "\t" not in report._HTML_TEMPLATE
        assert "\r" not in report._HTML_TEMPLATE

    def test_tsv_escapes_reach_the_browser_as_escapes(self, tmp_path):
        data = make_data(homology=[hit()], homology_all=[hit()])
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        text = path.read_text()
        assert r"join('\t')" in text, "tab separator must stay a JS escape"
        assert r"join('\n')" in text, "row separator must stay a JS escape"

    def test_the_tsv_builder_stays_on_its_own_lines(self, tmp_path):
        """The specific failure mode: the row separator turned into a real
        newline, splitting ``.join('\\n')`` across two source lines and leaving
        an unterminated string. Pin that the expression is intact on one line."""
        data = make_data(homology=[hit()], homology_all=[hit()])
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        joins = [ln for ln in path.read_text().splitlines() if ".join('" in ln]
        assert joins, "expected the TSV builder in the page"
        for line in joins:
            assert line.count("'") % 2 == 0, (
                f"quote left open, string spans a newline: {line.strip()[:90]!r}"
            )


class TestNamesSitOnTheirArrows:
    """Panel 5 writes each domain's name above its arrow. The name's anchor
    must share the arrow's y, or it labels the wrong bar on a multi-lane row.
    """

    @staticmethod
    def _names_and_positions(fig):
        named = {}
        for trace in fig.data:
            if trace.meta and trace.meta.get("teaidDomainName"):
                named = dict(zip(trace.text, trace.y))
        drawn = {}
        for trace in fig.data:
            tpl = getattr(trace, "hovertemplate", None) or ""
            for label in ("LOWER", "UPPER"):
                if f"<b>{label}</b>" in tpl and label not in drawn:
                    drawn[label] = trace.y[0]
        return named, drawn

    def _two_lane_figure(self):
        # Both hits share the ORF's reading frame (frame 1) and overlap each
        # other, so the packer must give them separate sub-lanes.
        orf = ORF(0, 900, "+", "M" * 300)
        lower = hit(query="LOWER", start=0, end=600)
        upper = hit(query="UPPER", start=300, end=900)
        data = make_data(orfs=[orf], homology=[lower, upper])
        return report.build_figure(data, theme.LIGHT)

    def test_the_two_hits_really_do_land_on_different_lanes(self):
        """Guards the fixture: if they shared a lane the test would pass for
        the wrong reason."""
        _, drawn = self._names_and_positions(self._two_lane_figure())
        assert set(drawn) == {"LOWER", "UPPER"}
        assert drawn["LOWER"] != drawn["UPPER"]

    def test_each_name_is_anchored_on_its_own_arrow(self):
        named, drawn = self._names_and_positions(self._two_lane_figure())
        assert set(named) == {"LOWER", "UPPER"}
        for label, y in drawn.items():
            assert named[label] == pytest.approx(y), f"{label} is named at another y"

    def test_a_single_lane_row_centres_its_arrow_on_the_row(self):
        """Row height also grows with the name count, so sizing the arrow
        spread by height dropped a lone arrow below the row's centre."""
        orf = ORF(0, 6000, "+", "M" * 2000)
        hits = [hit(query=f"D{i}", start=i * 900, end=i * 900 + 300)
                for i in range(6)]
        data = make_data(consensus_length=8000, orfs=[orf], homology=hits,
                         coverage=np.ones(8000, dtype=np.int64))
        fig = report.build_figure(data, theme.LIGHT)
        ys = {t.y[0] for t in fig.data
              if "<b>D" in (getattr(t, "hovertemplate", None) or "")}
        assert len(ys) == 1, "six non-overlapping hits belong on one lane"
        centres = list(fig.layout.yaxis5.tickvals)
        assert len(centres) == 1, "expected one row"
        assert ys.pop() == pytest.approx(centres[0]), (
            "the arrow must sit on the row's centre"
        )


def fragmented_data(**overrides):
    """One copy split into three rows: an insertion join, then a deletion join."""
    rows = [
        Copy(chrom="chr1", start=0, end=300, strand="+", family="fam", divergence=5.0,
             divergence_kind=DivergenceKind.CONSENSUS, consensus_start=0,
             consensus_end=300, hit_id="7"),
        Copy(chrom="chr1", start=1300, end=1600, strand="+", family="fam", divergence=5.0,
             divergence_kind=DivergenceKind.CONSENSUS, consensus_start=300,
             consensus_end=600, hit_id="7"),
        Copy(chrom="chr1", start=1600, end=1800, strand="+", family="fam", divergence=5.0,
             divergence_kind=DivergenceKind.CONSENSUS, consensus_start=800,
             consensus_end=1000, hit_id="7"),
    ]
    insertions = [Insertion(rows)]
    return make_data(insertions=insertions, full_length=[], **overrides)


class TestInsertionsOnTheSheet:
    def test_panel_1_counts_copies_not_rows(self):
        fig = report.build_figure(fragmented_data(), theme.LIGHT)
        names = [t.name for t in fig.data if t.name]
        assert "partial (n=1)" in names

    def test_the_subtitle_names_both_counts(self):
        fig = report.build_figure(fragmented_data(), theme.LIGHT)
        title = fig.layout.title.text
        assert "1 copy " in title and "3 annotation rows" in title

    def test_both_junction_kinds_are_drawn_and_legended(self):
        fig = report.build_figure(fragmented_data(), theme.LIGHT)
        names = [t.name for t in fig.data if t.name]
        assert "insertion in the copy (n=1)" in names
        assert "consensus missing from the copy (n=1)" in names
        gap = [t for t in fig.data if t.name and t.name.startswith("consensus missing")][0]
        assert gap.line.dash == "dash"
        assert gap.line.color == theme.LIGHT.muted, "grey, not the copy's colour"
        assert gap.line.width < report.COPY_WIDTH, "thinner than the copy it belongs to"
        assert not [t for t in fig.data if t.mode == "text" and "]" in list(t.text or [])]

    def test_a_full_length_copy_shares_the_gap_trace(self):
        data = fragmented_data()
        data = make_data(insertions=data.insertions, full_length=list(data.insertions))
        fig = report.build_figure(data, theme.LIGHT)
        gaps = [t for t in fig.data if t.legendgroup == "missing" and t.mode == "lines"]
        assert len(gaps) == 1 and gaps[0].line.color == theme.LIGHT.muted

    def test_panel_1_lines_carry_their_resting_width(self):
        fig = report.build_figure(fragmented_data(), theme.LIGHT)
        lines = [t for t in fig.data if t.meta and t.meta.get("teaidWidth")]
        assert sorted(t.meta["teaidWidth"] for t in lines) == [report.GAP_WIDTH, report.COPY_WIDTH]
        for trace in lines:
            assert trace.line.width == trace.meta["teaidWidth"]
        assert "line.width" in report._HTML_TEMPLATE and "teaidWidth" in report._HTML_TEMPLATE

    def test_a_gap_that_is_not_a_deletion_is_still_drawn(self):
        """A join that skips 40 bp of consensus but carries 1,000 bp of extra
        genome is an insertion by net difference; the 40 bp it lacks still shows."""
        rows = [
            Copy(chrom="chr1", start=0, end=300, strand="+", family="fam", divergence=5.0,
                 divergence_kind=DivergenceKind.CONSENSUS, consensus_start=0,
                 consensus_end=300, hit_id="9"),
            Copy(chrom="chr1", start=1300, end=1560, strand="+", family="fam", divergence=5.0,
                 divergence_kind=DivergenceKind.CONSENSUS, consensus_start=340,
                 consensus_end=600, hit_id="9"),
        ]
        fig = report.build_figure(make_data(insertions=[Insertion(rows)], full_length=[]),
                                  theme.LIGHT)
        gap = [t for t in fig.data if t.name and t.name.startswith("consensus missing")][0]
        assert list(gap.x[:2]) == [300, 340]
        assert "insertion in the copy (n=1)" in [t.name for t in fig.data if t.name]

    def test_insertion_sizes_start_hidden_and_gap_lengths_stay_in_the_hover(self):
        """At rest the sizes are hover-only: printed at a join they landed on
        neighbouring copies as often as on their own (rnd-1_family-21,
        2026-11-04). The page script shows the insertion sizes above their
        ticks once the view is zoomed in past INSERTION_LABEL_ZOOM."""
        fig = report.build_figure(fragmented_data(), theme.LIGHT)
        tick = [t for t in fig.data if t.name and t.name.startswith("insertion in the copy")][0]
        assert tick.mode == "markers" and tick.marker.symbol == "line-ns-open"
        assert list(tick.x) == [300]
        labels = [t for t in fig.data if t.mode == "text"]
        assert len(labels) == 1
        label = labels[0]
        assert list(label.text) == ["1,000 bp"] and list(label.x) == [300]
        assert label.visible is False and label.textposition == "top center"
        assert label.meta == {"teaidInsertionLabel": True}
        assert "teaidInsertionLabel" in report._HTML_TEMPLATE
        assert "LABEL_ZOOM" in report._HTML_TEMPLATE


def _legend_of(fig, trace):
    return fig.layout[trace.legend or "legend"]


class TestQuadrantCards:
    def _cards(self, fig):
        return [s for s in fig.layout.shapes if s.name == "quadrant"]

    def test_one_card_per_quadrant(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        cards = self._cards(fig)
        assert len(cards) == 4
        for card in cards:
            assert (card.type, card.xref, card.yref, card.layer) == (
                "path", "paper", "paper", "below")
            assert card.fillcolor == theme.LIGHT.surface
            assert (card.line.color, card.line.width) == (theme.LIGHT.grid, 1)

    def test_corners_are_rounded_by_the_card_radius(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        card = self._cards(fig)[0]
        x0, x1, y0, y1 = report.card_bounds(card)
        first = card.path.split()[1]  # the path starts one radius in from the corner
        assert float(first.split(",")[0]) == pytest.approx(
            x0 + report._CARD_RADIUS_PX / report._PLOT_WIDTH, abs=1e-6)
        assert card.path.count("Q") == 4

    def test_cards_in_a_row_end_at_the_same_depth(self):
        """Panel 1 has a legend and panel 2 none; their cards still align."""
        fig = report.build_figure(make_data(), theme.LIGHT)
        tops = {}
        for card in self._cards(fig):
            x0, x1, y0, y1 = report.card_bounds(card)
            tops.setdefault(round(y1, 6), set()).add(round(y0, 6))
        assert len(tops) == 2 and all(len(bottoms) == 1 for bottoms in tops.values())

    def test_cards_go_under_a_panel_s_own_shading(self):
        """The seed-depth shading is a 'below' shape too; a card appended after
        it would cover it, so the cards go first."""
        fig = report.build_figure(make_data(), theme.LIGHT)
        names = [s.name for s in fig.layout.shapes]
        assert names[:4] == ["quadrant"] * 4

    def test_a_card_holds_its_panel_title_and_legend(self):
        fig = report.build_figure(make_data(), theme.LIGHT)
        xaxis, yaxis = fig.layout.xaxis, fig.layout.yaxis
        bounds = [report.card_bounds(c) for c in self._cards(fig)]
        x0, x1, y0, y1 = min(bounds, key=lambda b: (b[0], -b[3]))
        plot_height = report._plot_height(2)
        assert x0 < xaxis.domain[0] and x1 > xaxis.domain[1]
        assert y1 > yaxis.domain[1] + report._TITLE_PX / plot_height
        legend_bottoms = [
            fig.layout[name].y - report._LEGEND_ROW_PX / plot_height
            for name in fig.layout if name.startswith("legend")
            and abs(fig.layout[name].x - sum(xaxis.domain) / 2) < 1e-9
        ]
        assert legend_bottoms and y0 < min(legend_bottoms)

    def test_panels_4_and_5_share_one_card(self):
        data = fragmented_data(
            terminal_repeats={"LTR": [SelfHit(100, 200, 800, 900, 98.0, 1e-30, 400.0)],
                              "TIR": []},
            homology=[hit()],
        )
        fig = report.build_figure(data, theme.LIGHT)
        cards = self._cards(fig)
        assert len(cards) == 4
        structure = fig.get_subplot(2, 2)
        homology = fig.get_subplot(3, 2)
        bounds = [report.card_bounds(c) for c in cards]
        right = [b for b in bounds if b[1] > 0.9]
        shared = [b for b in right if b[3] > structure.yaxis.domain[1]
                  and b[2] < homology.yaxis.domain[0]]
        assert len(shared) == 1
        # Panel 4 shows no tick labels, so it grows down towards panel 5.
        gap_px = (structure.yaxis.domain[0] - homology.yaxis.domain[1]) * report._plot_height(2)
        assert abs(gap_px - report._STACK_GAP_PX) < 0.5

    def test_the_toolbar_stays_at_the_top_of_the_page(self, tmp_path):
        data = make_data()
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        page = path.read_text()
        bar = re.search(r"\.teaid-bar \{(.*?)\}", page, re.S)
        assert bar and "position: sticky" in bar.group(1)

    def test_the_page_gives_the_cards_a_shadow(self, tmp_path):
        """Plotly shapes have no shadow, so the page adds one in CSS to the
        shape layer below the traces, where the cards are drawn."""
        data = make_data()
        fig = report.build_figure(data, theme.LIGHT)
        path = tmp_path / "sheet.html"
        report.write_html(fig, path, data, theme.LIGHT)
        assert re.search(r"\.layer-below \.shapelayer path \{\s*filter: drop-shadow",
                         path.read_text())


class TestLegendsSitUnderTheirPanels:
    def _sheet(self):
        data = fragmented_data(
            terminal_repeats={"LTR": [SelfHit(100, 200, 800, 900, 98.0, 1e-30, 400.0)],
                              "TIR": []},
            homology=[hit()],
        )
        return report.build_figure(data, theme.LIGHT)

    def test_each_panel_legend_is_centred_below_its_own_panel(self):
        fig = self._sheet()
        for trace in fig.data:
            if not (trace.name and trace.showlegend is not False):
                continue
            legend = _legend_of(fig, trace)
            xaxis = fig.layout[(trace.xaxis or "x").replace("x", "xaxis", 1)]
            yaxis = fig.layout[(trace.yaxis or "y").replace("y", "yaxis", 1)]
            assert legend.x == pytest.approx(sum(xaxis.domain) / 2), trace.name
            assert legend.y < yaxis.domain[0], trace.name
            assert legend.xanchor == "center" and legend.yanchor == "top"

    def test_panel_1_entries_wrap_into_rows_no_wider_than_the_panel(self):
        """Three entries at ~5.8 px a character do not fit a 489 px panel on one
        row; Plotly would wrap them at the full plot width, across panel 2."""
        fig = self._sheet()
        panel_1 = [t for t in fig.data if (t.xaxis or "x") == "x"
                   and t.name and t.showlegend is not False]
        assert len(panel_1) == 3
        ids = list(dict.fromkeys(t.legend or "legend" for t in panel_1))
        assert len(ids) == 2
        first, second = (fig.layout[i] for i in ids)
        plot_height = report._plot_height(2)
        assert first.y - second.y == pytest.approx(report._LEGEND_ROW_PX / plot_height)
        assert first.y == pytest.approx(
            fig.layout.yaxis.domain[0] - report._LEGEND_OFFSET_PX / plot_height)

    def test_panels_do_not_share_a_legend(self):
        fig = self._sheet()
        by_panel = {}
        for trace in fig.data:
            if trace.name and trace.showlegend is not False:
                by_panel.setdefault(trace.xaxis or "x", set()).add(trace.legend or "legend")
        assert len(by_panel) >= 3
        seen = [i for ids in by_panel.values() for i in ids]
        assert len(seen) == len(set(seen))

    def test_a_group_member_follows_the_entry_that_shows_it(self):
        """The arrowhead of a terminal-repeat arm is drawn without a legend
        entry but in the arm's legend group; it must sit in the arm's legend
        or a click on the entry would hide the arm and not the head."""
        fig = self._sheet()
        arms = [t for t in fig.data if t.legendgroup and t.legendgroup.startswith("direct pair")]
        assert len(arms) >= 2 and len({t.legend or "legend" for t in arms}) == 1

    def test_the_structure_legend_sits_close_under_a_bare_axis(self):
        """Panel 4 shows no tick labels when panel 5 is beneath it, so its
        legend starts right under the axis line rather than 46 px down, where
        it would sit midway between the two panels."""
        fig = self._sheet()
        arm = [t for t in fig.data if t.name and t.name.startswith("direct pair")][0]
        xaxis = fig.layout[(arm.xaxis or "x").replace("x", "xaxis", 1)]
        yaxis = fig.layout[(arm.yaxis or "y").replace("y", "yaxis", 1)]
        assert xaxis.showticklabels is False
        assert _legend_of(fig, arm).y == pytest.approx(
            yaxis.domain[0] - report._LEGEND_OFFSET_BARE_PX / report._plot_height(2))

    def test_the_notes_clear_the_deepest_bottom_legend(self):
        """Six class keys and two ORF strands under panel 5 are more than one
        row holds; the notes line and the bottom margin move down by the
        extra rows."""
        orders = ["PF00078", "PF00077", "PF00517", "PF02994", "PF01541", "PF00589",
                  "PF01609", "PF01359"]
        many = [hit(query=f"dom{i}", query_accession=f"{acc}.1", start=100 + 60 * i,
                    end=150 + 60 * i) for i, acc in enumerate(orders)]
        orfs = [ORF(0, 600, "+", "M" * 200), ORF(600, 900, "-", "M" * 100)]
        data = fragmented_data(homology=many, orfs=orfs)
        data.notes.append("a note")
        fig = report.build_figure(data, theme.LIGHT)
        plot_height = report._plot_height(2)
        panel_5 = [t for t in fig.data if t.xaxis == "x5" and t.name and t.showlegend is not False]
        rows = {t.legend for t in panel_5}
        assert len(rows) >= 2
        deepest = min(fig.layout[k].y for k in rows)
        note = [a for a in fig.layout.annotations if a.text == "a note"][0]
        assert note.y < deepest - report._LEGEND_ROW_PX / plot_height * 0.5
        assert fig.layout.margin.b >= report._MARGIN_B + (len(rows) - 1) * report._LEGEND_ROW_PX


class TestResetIsNotAliased:
    def test_the_button_hands_plotly_a_copy_of_the_reset_ranges(self):
        """Plotly adopts the arrays passed to relayout as its live axis ranges and
        edits them in place on a scroll-zoom. Passing RESET itself worked once,
        then restored the last zoom: seen in a browser on ltr-1_family-26, where
        panel 1 reset to [-93, 10921] instead of the padded [-290, 11883]."""
        assert "Plotly.relayout(gd, RESET)" not in report._HTML_TEMPLATE
        assert "JSON.parse(JSON.stringify(RESET))" in report._HTML_TEMPLATE


class TestThePageScriptRuns:
    def test_the_page_script_is_valid_javascript(self, tmp_path):
        """Static checks missed a SyntaxError here once. Parse it for real."""
        import shutil
        import subprocess

        node = shutil.which("node")
        if node is None:
            pytest.skip("node not installed")
        data = fragmented_data(homology=[hit()], homology_all=[hit()])
        fig = report.build_figure(data, theme.LIGHT)
        page = tmp_path / "sheet.html"
        report.write_html(fig, page, data, theme.LIGHT)
        text = page.read_text()
        # The sheet's own block is the last one; plotly's come before it.
        script = text[text.rindex("<script>") + len("<script>"):text.rindex("</script>")]
        source = tmp_path / "sheet.js"
        source.write_text(script)
        result = subprocess.run([node, "--check", str(source)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr


class TestStaticExportIsBounded:
    def test_a_normal_export_still_works(self, tmp_path):
        pytest.importorskip("kaleido")
        path = tmp_path / "sheet.svg"
        report.export_static(report.build_figure(make_data(), theme.LIGHT), path)
        assert path.stat().st_size > 5_000

    def test_a_stalled_export_is_killed_at_the_deadline(self, tmp_path, monkeypatch):
        """kaleido's headless Chrome stalled once in a 50-family batch and held
        the whole run. A stand-in child that never finishes must be killed."""
        import time

        monkeypatch.setattr(report, "_EXPORT", "import time\ntime.sleep(60)\n")
        started = time.perf_counter()
        with pytest.raises(TimeoutError):
            report.export_static(report.build_figure(make_data(), theme.LIGHT),
                                 tmp_path / "sheet.png", timeout=1.0)
        assert time.perf_counter() - started < 10
