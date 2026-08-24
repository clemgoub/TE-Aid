"""Tests for the Stockholm seed reader and the seed-QC panels."""

from __future__ import annotations

import gzip

import numpy as np
import pytest

from teaid import report, theme
from teaid.readers import stockholm
from teaid.records import DivergenceKind

# Two records, the second interleaved across blocks and using the short (2-part)
# Smitten identifier RepeatModeler writes. Gap character is Dfam's '.'.
SEED = """\
# STOCKHOLM 1.0
#=GF ID    fam-one
#=GF DE    a test family
#=GF TP    Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition
#=GF SQ    3
#=GC RF    ACGTACGTAC
OY720097.1:101-110_+    ACGTACGTAC
OY720097.1:201-210_-    ACGTACGTAT
OY720097.1:301-306_+    ACGTAC....
//
# STOCKHOLM 1.0
#=GF ID    fam-two
#=GF SQ    2
#=GC RF    ACGTA
seqA:1-10_+    ACGTA
seqB:1-10_+    ACGTA
#=GC RF    CCCCC
seqA:1-10_+    CCCCC
seqB:1-10_+    CCCGC
//
"""


@pytest.fixture
def seed_file(tmp_path):
    path = tmp_path / "families.stk"
    path.write_text(SEED)
    return path


class TestParsing:
    def test_reads_every_record(self, seed_file):
        assert [s.identifier for s in stockholm.read_all(seed_file)] == ["fam-one", "fam-two"]

    def test_selects_a_record_by_id(self, seed_file):
        assert stockholm.read(seed_file, "fam-two").identifier == "fam-two"

    def test_unknown_id_is_a_malformed_seed_error(self, seed_file):
        with pytest.raises(stockholm.MalformedSeedError):
            stockholm.read(seed_file, "absent")

    def test_multi_record_file_needs_a_family(self, seed_file):
        with pytest.raises(stockholm.MalformedSeedError):
            stockholm.read(seed_file)

    def test_interleaved_blocks_are_joined(self, seed_file):
        """Standard Stockholm may split an alignment across blocks; rows are
        accumulated by name rather than assumed to arrive in one piece."""
        seed = stockholm.read(seed_file, "fam-two")
        assert seed.reference_line == "ACGTACCCCC"
        assert [s.aligned for s in seed.sequences] == ["ACGTACCCCC", "ACGTACCCGC"]

    def test_a_file_without_a_header_is_rejected(self, tmp_path):
        path = tmp_path / "not.stk"
        path.write_text("ACGT\nACGT\n")
        with pytest.raises(stockholm.MalformedSeedError):
            stockholm.read_all(path)

    def test_gzip_is_transparent(self, tmp_path):
        path = tmp_path / "families.stk.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(SEED)
        assert len(stockholm.read_all(path)) == 2

    def test_features_and_declared_count(self, seed_file):
        seed = stockholm.read(seed_file, "fam-one")
        assert seed.declared_count == 3
        assert seed.description == "a test family"
        assert seed.expected_class.startswith("Interspersed_Repeat;")


class TestSmittenIdentifiers:
    def test_short_form_without_an_assembly(self):
        assert stockholm.parse_smitten("OY720097.1:101-110_+") == (
            None, "OY720097.1", 100, 110, "+",
        )

    def test_long_form_with_an_assembly(self):
        assert stockholm.parse_smitten("GCA_951799975.1:OX637595.1:15848-16090_+") == (
            "GCA_951799975.1", "OX637595.1", 15847, 16090, "+",
        )

    def test_coordinates_convert_from_1_based_closed(self):
        """Smitten is 1-based fully closed; everything here is 0-based half-open,
        and mixing the two silently is the documented trap (§7)."""
        _, _, start, end, _ = stockholm.parse_smitten("chr1:1-10_+")
        assert (start, end) == (0, 10)
        assert end - start == 10

    def test_a_plain_name_is_not_an_error(self):
        """A seed may legitimately carry non-Smitten names; that is an absence
        of loci, not a malformed file."""
        assert stockholm.parse_smitten("just_a_name") is None

    def test_non_smitten_sequences_still_parse(self, tmp_path):
        path = tmp_path / "plain.stk"
        path.write_text("# STOCKHOLM 1.0\n#=GF ID x\n#=GC RF ACGT\nplain_name ACGT\n//\n")
        seed = stockholm.read(path, "x")
        assert seed.sequences[0].has_locus is False
        assert len(seed.to_annotation()) == 1


class TestDerivedQuantities:
    def test_consensus_strips_reference_gaps(self, tmp_path):
        """Insert columns are not consensus positions."""
        path = tmp_path / "ins.stk"
        path.write_text(
            "# STOCKHOLM 1.0\n#=GF ID x\n#=GC RF AC.GT\n"
            "chr1:1-5_+ ACTGT\nchr1:11-15_+ AC.GT\n//\n"
        )
        seed = stockholm.read(path, "x")
        assert seed.consensus() == "ACGT"
        assert seed.consensus_columns() == [0, 1, 3, 4]

    def test_an_insertion_does_not_shift_later_coordinates(self, tmp_path):
        path = tmp_path / "ins.stk"
        path.write_text(
            "# STOCKHOLM 1.0\n#=GF ID x\n#=GC RF AC.GT\nchr1:1-5_+ ACTGT\n//\n"
        )
        copy = stockholm.read(path, "x").to_annotation().copies[0]
        assert copy.consensus_1based() == (1, 4)

    def test_depth_counts_sequences_per_consensus_position(self, seed_file):
        seed = stockholm.read(seed_file, "fam-one")
        # Third sequence covers only the first six columns.
        assert seed.depth() == [3, 3, 3, 3, 3, 3, 2, 2, 2, 2]

    def test_divergence_is_computed_against_the_reference(self, seed_file):
        copies = stockholm.read(seed_file, "fam-one").to_annotation().copies
        assert copies[0].divergence == pytest.approx(0.0)  # identical to RF
        assert copies[1].divergence == pytest.approx(10.0)  # one mismatch in ten
        assert copies[0].divergence_kind is DivergenceKind.CONSENSUS

    def test_loci_and_strand_come_from_the_identifiers(self, seed_file):
        copies = stockholm.read(seed_file, "fam-one").to_annotation().copies
        assert copies[0].chrom == "OY720097.1"
        assert copies[0].genomic_1based() == (101, 110)
        assert copies[1].strand == "-"

    def test_consensus_falls_back_to_a_majority_without_an_rf_line(self, tmp_path):
        path = tmp_path / "norf.stk"
        path.write_text(
            "# STOCKHOLM 1.0\n#=GF ID x\nchr1:1-4_+ ACGT\nchr1:11-14_+ ACGA\n//\n"
        )
        assert stockholm.read(path, "x").consensus() == "ACGT"


class TestSeedQcPanels:
    @staticmethod
    def _sheet(depth, expected_class=None, count=None):
        from teaid.records import Copy

        copies = [
            Copy(
                chrom="c", start=0, end=900, strand="+", family="f",
                divergence=5.0, divergence_kind=DivergenceKind.CONSENSUS,
                consensus_start=0, consensus_end=900,
            )
        ]
        return report.SheetData(
            family="f",
            consensus_length=len(depth) if depth is not None else 1000,
            copies=copies,
            full_length=copies,
            coverage=np.ones(1000, dtype=np.int64),
            seed_depth=depth,
            expected_class=expected_class,
            seed_sequence_count=count,
        )

    def test_seed_qc_is_off_unless_asked_for(self):
        assert self._sheet(None).has_seed_qc is False
        fig = report.build_figure(self._sheet(None), theme.LIGHT)
        assert len([k for k in fig.layout if k.startswith("xaxis")]) == 4

    def test_seed_qc_appends_a_row_without_disturbing_the_quadrants(self):
        fig = report.build_figure(self._sheet(np.array([5] * 1000)), theme.LIGHT)
        assert len([k for k in fig.layout if k.startswith("xaxis")]) == 6
        # The main 2x2 keeps its shared, padded consensus range.
        expected = tuple(report._padded_span(1000))
        for name in ("xaxis", "xaxis2", "xaxis3", "xaxis4"):
            assert tuple(fig.layout[name].range) == expected

    def test_thin_stretches_are_shaded_and_nothing_else_is(self):
        """A full-width sub-floor band muddies every seed including good ones;
        shading only the shortfall answers 'where is it thin' directly."""
        depth = np.array([5] * 400 + [1] * 200 + [5] * 400)
        fig = report.build_figure(self._sheet(depth), theme.LIGHT)
        shapes = fig.layout.shapes or []
        assert [(s.x0, s.x1) for s in shapes] == [(400, 600)]

    def test_a_seed_clearing_the_floor_is_not_shaded(self):
        fig = report.build_figure(self._sheet(np.array([9] * 1000)), theme.LIGHT)
        assert not (fig.layout.shapes or [])

    def test_shortfall_is_stated_on_the_sheet(self):
        data = self._sheet(np.array([2] * 300 + [8] * 700))
        report.build_figure(data, theme.LIGHT)
        assert any("below the Dfam" in note for note in data.notes)

    def test_expected_class_is_shown_as_a_supplied_claim(self):
        fig = report.build_figure(
            self._sheet(np.array([5] * 100), "Interspersed_Repeat;Foo", 5), theme.LIGHT
        )
        text = " ".join(a.text or "" for a in fig.layout.annotations)
        assert "#=GF TP" in text
        assert "supplied with the seed" in text
        # The comparison is stated as pending, never guessed at, until there is
        # homology evidence to check against.
        assert "pending" in text

    def test_a_seed_without_tp_says_so_rather_than_inventing_one(self):
        fig = report.build_figure(self._sheet(np.array([5] * 100), None, 5), theme.LIGHT)
        text = " ".join(a.text or "" for a in fig.layout.annotations)
        assert "no <b>#=GF TP</b>" in text

    def test_sequence_count_is_judged_against_the_dfam_floor(self):
        below = report.build_figure(self._sheet(np.array([2] * 50), None, 2), theme.LIGHT)
        assert "below the Dfam minimum" in " ".join(
            a.text or "" for a in below.layout.annotations
        )
        meets = report.build_figure(self._sheet(np.array([5] * 50), None, 5), theme.LIGHT)
        assert "meets the Dfam minimum" in " ".join(
            a.text or "" for a in meets.layout.annotations
        )


class TestRunsBelow:
    def test_finds_each_contiguous_shortfall(self):
        values = np.array([5, 1, 1, 5, 5, 0, 5])
        assert report._runs_below(values, 3) == [(1, 3), (5, 6)]

    def test_returns_nothing_when_all_values_clear_the_floor(self):
        assert report._runs_below(np.array([9, 9, 9]), 3) == []

    def test_handles_a_shortfall_at_both_edges(self):
        assert report._runs_below(np.array([0, 9, 0]), 3) == [(0, 1), (2, 3)]
