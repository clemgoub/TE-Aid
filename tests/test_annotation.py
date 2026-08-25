"""Tests for panel 5's row layout and v1's TE-class colour scheme."""

from __future__ import annotations

import pytest

from teaid import annotation_rows, orfs as orfs_mod, tetypes
from teaid.annotation_rows import reading_frame
from teaid.homology import ProteinHit
from teaid.orfs import ORF


def hit(start, end, strand="+", **over) -> ProteinHit:
    base = dict(
        query="D", query_accession="PF00078", start=start, end=end, strand=strand,
        evalue=1e-20, score=100.0, identity=30.0, hmm_from=1, hmm_to=100,
        hmm_len=100, frameshifts=0, stop_codons=0,
    )
    base.update(over)
    return ProteinHit(**base)


class TestReadingFrame:
    """A sign error here silently misaligns every hit against every ORF, which
    is worse than having no frame logic at all."""

    @pytest.mark.parametrize("start,expected", [(0, 1), (1, 2), (2, 3), (3, 1), (99, 1)])
    def test_forward_frame_counts_from_the_start(self, start, expected):
        assert reading_frame(start, start + 30, "+", 1000) == expected

    @pytest.mark.parametrize("end,expected", [(1000, 1), (999, 2), (998, 3), (997, 1)])
    def test_reverse_frame_counts_from_the_consensus_end(self, end, expected):
        """A reverse feature is translated from the other end, so its register
        is set by the distance to the end, not by its start."""
        assert reading_frame(end - 30, end, "-", 1000) == expected

    def test_the_two_strands_are_computed_differently(self):
        assert reading_frame(0, 30, "+", 1000) != reading_frame(0, 30, "-", 1000)

    def test_frame_is_stable_under_length_of_the_feature(self):
        assert reading_frame(9, 39, "+", 1000) == reading_frame(9, 300, "+", 1000)


class TestHousing:
    def test_a_hit_in_the_same_frame_is_housed(self):
        orf = ORF(300, 900, "+", "M" * 200)
        h = hit(300, 600)
        rows = annotation_rows.build([orf], [h], 2000)
        assert len(rows) == 1 and rows[0].hits == [h]

    def test_a_hit_in_another_frame_is_not_housed(self):
        """The load-bearing rule. Overlap alone would seat a frameshifted domain
        in a stop-to-stop block of a different register, and the panel would
        then assert 'intact coding domain' for the finding it exists to show."""
        orf = ORF(300, 900, "+", "M" * 200)
        h = hit(301, 601)  # same span, one base over — different frame
        rows = annotation_rows.build([orf], [h], 2000)
        assert len(rows) == 2
        assert any(r.orf is None and r.hits == [h] for r in rows)

    def test_the_opposite_strand_is_not_housed(self):
        orf = ORF(300, 900, "+", "M" * 200)
        h = hit(300, 600, strand="-")
        rows = annotation_rows.build([orf], [h], 2000)
        assert any(r.orf is None for r in rows)

    def test_a_hit_barely_touching_an_orf_is_not_housed(self):
        orf = ORF(300, 900, "+", "M" * 200)
        h = hit(negative := 0, 330)          # only 30 of 330 bp inside
        rows = annotation_rows.build([orf], [h], 2000)
        assert any(r.orf is None and r.hits == [h] for r in rows)

    def test_an_orf_with_no_hits_still_gets_a_row(self):
        rows = annotation_rows.build([ORF(300, 900, "+", "M" * 200)], [], 2000)
        assert len(rows) == 1 and rows[0].orf is not None and rows[0].hits == []

    def test_a_hit_is_housed_by_only_one_orf(self):
        a = ORF(300, 900, "+", "M" * 200)
        b = ORF(300, 1200, "+", "M" * 300)
        rows = annotation_rows.build([a, b], [hit(300, 600)], 2000)
        assert sum(len(r.hits) for r in rows) == 1


class TestRowLayout:
    def test_rows_are_ordered_along_the_consensus(self):
        rows = annotation_rows.build(
            [ORF(2000, 2600, "+", "M" * 200)], [hit(0, 300)], 4000
        )
        assert rows[0].start < rows[1].start

    def test_overlapping_hits_get_their_own_sublanes(self):
        orf = ORF(0, 1200, "+", "M" * 400)
        a, b = hit(0, 600), hit(300, 900)
        rows = annotation_rows.build([orf], [a, b], 2000)
        assert sorted(rows[0].lanes) == [0, 1]

    def test_non_overlapping_hits_share_one_sublane(self):
        orf = ORF(0, 1200, "+", "M" * 400)
        rows = annotation_rows.build([orf], [hit(0, 300), hit(600, 900)], 2000)
        assert set(rows[0].lanes) == {0}

    def test_row_height_makes_room_for_the_tick_roster(self):
        """Six domains need three sub-lanes of room even when none overlap:
        sizing by packing alone lets adjacent rosters collide."""
        orf = ORF(0, 6000, "+", "M" * 2000)
        hits = [hit(i * 900, i * 900 + 300) for i in range(6)]
        rows = annotation_rows.build([orf], hits, 8000)
        assert set(rows[0].lanes) == {0}
        assert rows[0].height == 3

    def test_a_bare_row_is_one_unit_tall(self):
        rows = annotation_rows.build([], [hit(0, 300)], 2000)
        assert rows[0].height == 1 and rows[0].orf is None


class TestV1Palette:
    """v1's exact hex values, because curators have read TE-Aid sheets with them
    for years: green is LTR, blue is LINE, salmon is a DNA transposon."""

    @pytest.mark.parametrize("label,key,hexcode", [
        ("LTR/Gypsy", "LTR", "#00cc44"),
        ("LINE/L1", "LINE", "#3399ff"),
        ("SINE/tRNA", "SINE", "#800080"),
        ("DNA/hAT-Ac", "TIR", "#ff6666"),
        ("RC/Helitron", "RC", "#ff6600"),
    ])
    def test_the_main_classes_keep_v1_colours(self, label, key, hexcode):
        assert tetypes.class_of_repeatmasker_label(label) == key
        assert tetypes.colour(key) == hexcode

    @pytest.mark.parametrize("label,key", [
        ("DNA/Maverick", "MAV"),
        ("DNA/Crypton", "CRY"),
        ("DNA/Cryp", "CRY"),
        ("LINE/Penelope", "PLE"),
        ("LTR/DIRS", "DIRS"),
    ])
    def test_v1_relabelling_beats_the_general_prefix(self, label, key):
        """v1 rewrote these before collapsing DNA -> TIR and they must still
        win, or Maverick reads as an ordinary DNA transposon."""
        assert tetypes.class_of_repeatmasker_label(label) == key

    def test_an_unrecognised_label_is_unknown_not_a_guess(self):
        assert tetypes.class_of_repeatmasker_label("Something/Else") == "Unknown"
        assert tetypes.class_of_repeatmasker_label(None) == "Unknown"

    def test_simple_repeat_is_a_real_colour(self):
        """v1 wrote '#8686ac' into a table whose other values omitted the '#',
        so it rendered as '##8686ac' and fell back to a default."""
        assert tetypes.colour("Simple_repeat") == "#8686ac"
        assert all(c.startswith("#") and len(c) == 7 for c in tetypes.CLASS_COLOURS.values())

    def test_order_agnostic_domains_are_not_given_an_order(self):
        """A bare reverse transcriptase is shared by every Class I order, so
        colouring it LTR or LINE would assert something the hit cannot support."""
        assert tetypes.class_of_order("Class I") == "Unknown"
        assert tetypes.class_of_order("domesticated") == "Unknown"

    def test_specific_orders_do_get_their_colour(self):
        assert tetypes.class_of_order("TIR/hAT") == "TIR"
        assert tetypes.class_of_order("RC/Helitron") == "RC"
        assert tetypes.class_of_order("LTR/ERV") == "LTR"

    def test_a_tier3_hit_is_classified_from_its_own_label(self):
        h = hit(0, 300, query="Gypsy-1_DR_pol#LTR/Gypsy", query_accession="-")
        assert tetypes.classify_hit(h, {}) == "LTR"

    def test_a_pfam_hit_is_classified_from_the_curated_table(self):
        h = hit(0, 300, query="Dimer_Tnp_hAT", query_accession="PF05699.1")
        assert tetypes.classify_hit(h, {"PF05699": "TIR/hAT"}) == "TIR"


class TestOrfInvocation:
    def test_getorf_is_pinned_to_stop_to_stop_blocks(self):
        """The 'domains in one frame cannot overlap' guarantee, and 'the
        rectangle's edge is where the frame closes', both depend on getorf's
        default -find 0 (translate between stops). Passing a different -find
        would invalidate both, silently."""
        command = orfs_mod.build_command("q.fa", "o.fa", 400, True)
        assert "-find" not in command
        assert command[:1] == ["getorf"]
        assert "-minsize" in command and "400" in command

    def test_reverse_orfs_are_opt_out_not_opt_in(self):
        assert "-reverse" in orfs_mod.build_command("q", "o", 400, True)
        assert "-noreverse" in orfs_mod.build_command("q", "o", 400, False)
