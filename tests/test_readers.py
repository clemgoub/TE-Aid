"""Reader tests, centred on the traps each format sets.

Fixtures are written inline rather than committed as data files: every row here
exists to pin one documented behaviour, and keeping the row beside its assertion
is what makes the test readable.
"""

from __future__ import annotations

import gzip

import pytest

from teaid import readers
from teaid.readers import bed16, gff3, rmout
from teaid.records import DivergenceKind

# Real rows from GenomeArk GCA_963082875.1. The third and fourth are 'C' rows,
# where the consensus columns are written '(left) end begin'. Family
# rnd-1_family-257 is 1755 bp, so row 4 (308..400 with 1355 left) is the
# arithmetic check: 400 + 1355 == 1755.
RMOUT = """\
   SW   perc perc perc  query       position in query              matching          repeat                position in repeat
score   div. del. ins.  sequence    begin    end          (left)   repeat            class/family      begin   end    (left)     ID

 6759    7.0  0.1  0.0  OY720097.1         4     8148 (23708455) + (ACCCTG)n         Simple_repeat           1   8154     (0)     1
14672    5.0  1.8  0.4  OY720097.1      8303    10606 (23705997) C rnd-1_family-257  Unknown               (0)   1755       1     3
  368   20.4  0.0  1.1  OY720097.1     10674    10767 (23705836) C rnd-1_family-257  Unknown            (1355)    400     308     4
 4110    9.0  1.2  1.9  OY720097.1     11531    12323 (23704280) C rnd-1_family-21   Unknown               (0)   1302     608     8 *
"""


@pytest.fixture
def rmout_file(tmp_path):
    path = tmp_path / "sample.fa.out"
    path.write_text(RMOUT)
    return path


class TestRepeatMaskerOut:
    def test_skips_headers_and_reads_every_data_row(self, rmout_file):
        annotation = rmout.read(rmout_file)
        assert len(annotation) == 4
        assert annotation.skipped == []

    def test_plus_strand_consensus_columns_are_begin_end_left(self, rmout_file):
        simple = rmout.read(rmout_file).copies[0]
        assert (simple.consensus_start, simple.consensus_end) == (0, 8154)
        assert simple.consensus_left == 0

    def test_c_strand_consensus_columns_are_reversed(self, rmout_file):
        """The central trap: on 'C' rows the order is (left) end begin."""
        copy = rmout.read(rmout_file).copies[2]
        # Row reads '(1355) 400 308' -> begin 308, end 400, left 1355.
        assert copy.consensus_1based() == (308, 400)
        assert copy.consensus_left == 1355
        # Reading positionally instead would give begin=1355, end=400 and a
        # transposed, negative-length interval.
        assert copy.consensus_end > copy.consensus_start

    def test_c_strand_becomes_minus(self, rmout_file):
        assert [c.strand for c in rmout.read(rmout_file)] == ["+", "-", "-", "-"]

    def test_consensus_length_inferred_from_end_plus_left(self, rmout_file):
        family = rmout.read(rmout_file).for_family("rnd-1_family-257")
        # 1755 + 0 from the full-length copy, 400 + 1355 from the fragment.
        assert {c.consensus_length_estimate for c in family} == {1755}
        assert family.consensus_length() == 1755

    def test_genomic_coordinates_convert_to_half_open(self, rmout_file):
        copy = rmout.read(rmout_file).copies[0]
        assert (copy.start, copy.end) == (3, 8148)  # from 1-based 4..8148
        assert copy.genomic_1based() == (4, 8148)
        assert copy.length == 8145

    def test_trailing_star_marks_lower_scoring_overlap(self, rmout_file):
        flags = [c.is_overlapping_lower_score for c in rmout.read(rmout_file)]
        assert flags == [False, False, False, True]


ROW = "  368   20.4  0.0  1.1  chr1   10674  10767 (23705836) C fam  Unknown  (1355)  400  308"


class TestRepeatMaskerFragments:
    """The ID column is what joins the pieces of one interrupted insertion."""

    @staticmethod
    def parse(tmp_path, *rows):
        path = tmp_path / "x.fa.out"
        path.write_text("\n".join(rows) + "\n")
        return rmout.read(path)

    def test_reads_the_id_column(self, tmp_path):
        copy = self.parse(tmp_path, ROW + "  7").copies[0]
        assert copy.hit_id == "7"

    def test_rows_sharing_an_id_group_into_one_insertion(self, tmp_path):
        annotation = self.parse(tmp_path, ROW + "  7", ROW + "  7", ROW + "  8")
        assert sorted(len(g) for g in annotation.fragment_groups()) == [1, 2]

    def test_the_star_flag_is_not_read_as_an_id(self, tmp_path):
        copy = self.parse(tmp_path, ROW + "  7 *").copies[0]
        assert copy.hit_id == "7" and copy.is_overlapping_lower_score

    def test_a_no_id_file_still_parses(self, tmp_path):
        """RepeatMasker -no_id drops the column; every such row used to be skipped."""
        annotation = self.parse(tmp_path, ROW, ROW)
        assert len(annotation) == 2 and annotation.skipped == []
        assert all(c.hit_id is None for c in annotation)

    def test_a_no_id_row_with_a_star_has_no_id(self, tmp_path):
        copy = self.parse(tmp_path, ROW + "  *").copies[0]
        assert copy.hit_id is None and copy.is_overlapping_lower_score

    def test_every_row_reads_divergence_from_consensus(self, rmout_file):
        """Decided 2026-10-13: Simple_repeat rows included; array homogeneity is VGP's."""
        copies = rmout.read(rmout_file).copies
        assert {c.divergence_kind for c in copies} == {DivergenceKind.CONSENSUS}

    def test_unparsable_row_is_skipped_not_raised(self, tmp_path):
        path = tmp_path / "bad.out"
        path.write_text(RMOUT + " 100    1.0  0.0  0.0  chr1  1  10 (0) ? fam  Unknown  1 10 (0) 9\n")
        annotation = rmout.read(path)
        assert len(annotation) == 4
        assert len(annotation.skipped) == 1
        assert "strand" in annotation.skipped[0][1]


BED16_HEADER = (
    "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tSW_score\tperc_div\t"
    "perc_del\tperc_ins\tquery_left\trepeat_class_family\trepeat_start\t"
    "repeat_end\trepeat_left\thit_id\n"
)
BED16_ROWS = (
    "OX637595.1\t15848\t16090\t(ACTACT)n\t171\t+\t171\t9.4\t0.0\t0.0\t"
    "76485642\tSimple_repeat\t1\t242\t0\t1\n"
    "OX637595.1\t21044\t21395\tL2-3_DR\t421\t-\t421\t22.1\t1.4\t0.7\t"
    "76480337\tLINE/L2\t2891\t3247\t1204\t2\n"
    # A structural caller: NA everywhere it cannot measure, strand '.'.
    "OX637595.1\t48210\t53887\tLTR_retro_1\t0\t.\tNA\tNA\tNA\tNA\tNA\t"
    "LTR/Gypsy\tNA\tNA\tNA\tTE_struc_1\n"
)


@pytest.fixture
def bed16_file(tmp_path):
    path = tmp_path / "sample.bed"
    path.write_text(BED16_HEADER + BED16_ROWS)
    return path


class TestBed16:
    def test_header_is_skipped_and_rows_read(self, bed16_file):
        annotation = bed16.read(bed16_file)
        assert len(annotation) == 3
        assert annotation.skipped == []

    def test_genomic_coordinates_pass_through_unshifted(self, bed16_file):
        """BED16 is already 0-based half-open, unlike .out."""
        copy = bed16.read(bed16_file).copies[0]
        assert (copy.start, copy.end) == (15848, 16090)

    def test_na_becomes_none_never_zero(self, bed16_file):
        """The load-bearing distinction: NA means 'cannot report', 0 is a measurement."""
        structural = bed16.read(bed16_file).copies[2]
        assert structural.divergence is None
        assert structural.divergence_kind is DivergenceKind.NONE
        assert structural.has_divergence is False
        assert structural.consensus_start is None
        assert structural.consensus_left is None

    def test_strandless_features_are_legal(self, bed16_file):
        assert bed16.read(bed16_file).copies[2].strand == "."

    def test_every_row_reads_divergence_from_consensus(self, bed16_file):
        copies = bed16.read(bed16_file).copies
        assert copies[0].divergence_kind is DivergenceKind.CONSENSUS  # Simple_repeat
        assert copies[1].divergence_kind is DivergenceKind.CONSENSUS  # LINE/L2

    def test_divergence_kind_can_be_forced(self, bed16_file):
        annotation = bed16.read(bed16_file, divergence_kind=DivergenceKind.ARRAY_HOMOGENEITY)
        assert annotation.copies[1].divergence_kind is DivergenceKind.ARRAY_HOMOGENEITY
        # A record with no divergence stays NONE regardless of the override.
        assert annotation.copies[2].divergence_kind is DivergenceKind.NONE

    def test_hit_id_groups_interrupted_fragments(self, tmp_path):
        path = tmp_path / "frag.bed"
        rows = (
            "chr1\t100\t200\tfamA\t0\t+\tNA\t5.0\tNA\tNA\tNA\tLINE/L2\t1\t100\t400\tcopy7\n"
            "chr1\t300\t600\tfamA\t0\t+\tNA\t5.0\tNA\tNA\tNA\tLINE/L2\t101\t400\t100\tcopy7\n"
            "chr1\t900\t950\tfamA\t0\t+\tNA\t5.0\tNA\tNA\tNA\tLINE/L2\t1\t50\t450\tcopy8\n"
        )
        path.write_text(BED16_HEADER + rows)
        groups = bed16.read(path).fragment_groups()
        assert sorted(len(g) for g in groups) == [1, 2]

    def test_a_star_in_column_16_is_not_an_id(self, tmp_path):
        """A -no_id .out converted by rmout2bed puts RepeatMasker's '*' flag there."""
        path = tmp_path / "star.bed"
        rows = (
            "chr1\t100\t200\tfamA\t0\t+\tNA\t5.0\tNA\tNA\tNA\tLINE/L2\t1\t100\t400\t*\n"
            "chr1\t900\t950\tfamA\t0\t+\tNA\t5.0\tNA\tNA\tNA\tLINE/L2\t1\t50\t450\t*\n"
        )
        path.write_text(BED16_HEADER + rows)
        annotation = bed16.read(path)
        assert [c.hit_id for c in annotation.copies] == [None, None]
        assert len(annotation.fragment_groups()) == 2

    def test_repeatmasker_c_strand_reads_as_minus(self, tmp_path):
        path = tmp_path / "c.bed"
        path.write_text(BED16_HEADER + "chr1\t100\t200\tfamA\t0\tC\tNA\t5.0\tNA\tNA\tNA"
                                       "\tLINE/L2\t1\t100\t400\tcopy1\n")
        assert bed16.read(path).copies[0].strand == "-"

    def test_gzip_is_transparent(self, tmp_path):
        path = tmp_path / "sample.bed.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(BED16_HEADER + BED16_ROWS)
        assert len(bed16.read(path)) == 3

    def test_short_row_is_skipped(self, tmp_path):
        path = tmp_path / "short.bed"
        path.write_text(BED16_HEADER + "chr1\t1\t2\tfam\n")
        annotation = bed16.read(path)
        assert len(annotation) == 0
        assert len(annotation.skipped) == 1


GFF3 = """\
##gff-version 3
OY720097.1\tRepeatMasker\tdispersed_repeat\t4\t8148\t6759\t+\t.\tTarget=(ACCCTG)n 1 8154;Class=Simple_repeat
OY720097.1\tRepeatMasker\tsimilarity\t8303\t10606\t14672\t-\t.\tTarget "Motif:rnd-1_family-257" 1 1755;Divergence=5.0;Class=Unknown
"""


class TestGff3:
    def test_both_target_dialects_parse(self, tmp_path):
        path = tmp_path / "sample.gff3"
        path.write_text(GFF3)
        annotation = gff3.read(path)
        assert [c.family for c in annotation] == ["(ACCCTG)n", "rnd-1_family-257"]
        assert annotation.copies[1].consensus_1based() == (1, 1755)

    def test_coordinates_convert_to_half_open(self, tmp_path):
        path = tmp_path / "sample.gff3"
        path.write_text(GFF3)
        assert gff3.read(path).copies[0].start == 3

    def test_score_column_is_never_read_as_divergence(self, tmp_path):
        """Column 6 is a Smith-Waterman score; putting it on a divergence axis
        would plot an alignment score as a percentage."""
        path = tmp_path / "sample.gff3"
        path.write_text(GFF3)
        no_div = gff3.read(path).copies[0]
        assert no_div.divergence is None
        assert no_div.divergence_kind is DivergenceKind.NONE

    def test_row_without_target_is_skipped(self, tmp_path):
        path = tmp_path / "no_target.gff3"
        path.write_text("##gff-version 3\nchr1\tx\tmatch\t1\t10\t.\t+\t.\tID=a\n")
        annotation = gff3.read(path)
        assert len(annotation) == 0
        assert "Target" in annotation.skipped[0][1]


class TestDetection:
    def test_detects_each_format_by_content(self, tmp_path, rmout_file, bed16_file):
        gff_path = tmp_path / "x.gff3"
        gff_path.write_text(GFF3)
        assert readers.detect(rmout_file) == "rmout"
        assert readers.detect(bed16_file) == "bed16"
        assert readers.detect(gff_path) == "gff3"

    def test_detection_ignores_a_misleading_extension(self, tmp_path):
        """The hub ships .bed files converted from .out; extensions lie."""
        path = tmp_path / "actually_rmout.bed"
        path.write_text(RMOUT)
        assert readers.detect(path) == "rmout"

    def test_unrecognisable_file_raises(self, tmp_path):
        path = tmp_path / "prose.txt"
        path.write_text("this is not an annotation\n")
        with pytest.raises(readers.FormatDetectionError):
            readers.detect(path)

    def test_read_dispatches_on_detected_format(self, rmout_file):
        assert readers.read(rmout_file).source_format == "rmout"
