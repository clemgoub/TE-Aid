"""Tests for the derived quantities behind panels 1-3."""

from __future__ import annotations

import numpy as np
import pytest

from teaid import analysis
from teaid.records import Annotation, Copy, DivergenceKind
from teaid.sequences import read_fasta


def copy(cons_start, cons_end, *, divergence=5.0, family="fam", chrom="chr1"):
    return Copy(
        chrom=chrom,
        start=0,
        end=cons_end - cons_start,
        strand="+",
        family=family,
        divergence=divergence,
        divergence_kind=DivergenceKind.CONSENSUS,
        consensus_start=cons_start,
        consensus_end=cons_end,
    )


class TestCoverage:
    def test_counts_copies_overlapping_each_position(self):
        annotation = Annotation(copies=[copy(0, 10), copy(5, 15), copy(5, 8)])
        result = analysis.coverage(annotation, 20)
        assert result[0] == 1  # only the first copy
        assert result[5] == 3  # all three
        assert result[9] == 2  # first and second
        assert result[14] == 1  # second only
        assert result[15] == 0  # past every copy

    def test_matches_a_naive_dense_implementation(self):
        """Same answer as v1's dense matrix, in O(n + L) instead of O(n x L)."""
        rng = np.random.default_rng(20260823)
        length = 500
        copies = []
        for _ in range(200):
            start = int(rng.integers(0, length - 1))
            end = int(rng.integers(start + 1, length + 1))
            copies.append(copy(start, end))
        annotation = Annotation(copies=copies)

        naive = np.zeros((len(copies), length), dtype=np.int64)
        for row, c in enumerate(copies):
            naive[row, c.consensus_start : c.consensus_end] = 1

        assert np.array_equal(analysis.coverage(annotation, length), naive.sum(axis=0))

    def test_copies_without_consensus_coordinates_are_ignored(self):
        blind = Copy(chrom="chr1", start=0, end=100, strand="+", family="fam")
        annotation = Annotation(copies=[copy(0, 10), blind])
        assert analysis.coverage(annotation, 20)[0] == 1

    def test_clamps_copies_extending_past_the_consensus(self):
        annotation = Annotation(copies=[copy(0, 50)])
        result = analysis.coverage(annotation, 10)
        assert len(result) == 10
        assert result.max() == 1

    def test_rejects_a_nonpositive_length(self):
        with pytest.raises(ValueError):
            analysis.coverage(Annotation(), 0)


class TestFullLength:
    def test_selects_copies_spanning_the_threshold(self):
        # Spans of 100, 50, 90 and 89 against a 90 bp cutoff: the 90 is included
        # (the comparison is >=), the 89 is not.
        annotation = Annotation(
            copies=[copy(0, 100), copy(0, 50), copy(10, 100), copy(10, 99)]
        )
        selected = analysis.full_length(annotation, 100, 0.9)
        assert [(c.consensus_start, c.consensus_end) for c in selected] == [
            (0, 100),
            (10, 100),
        ]

    def test_a_copy_spanning_the_whole_consensus_passes_a_threshold_of_one(self):
        """v1 measured the span as abs(qend - qstart) on closed coordinates, so a
        copy covering the consensus end to end scored L-1 and missed."""
        annotation = Annotation(copies=[copy(0, 100)])
        assert len(analysis.full_length(annotation, 100, 1.0)) == 1

    def test_ignores_copies_without_consensus_coordinates(self):
        blind = Copy(chrom="chr1", start=0, end=100, strand="+", family="fam")
        assert analysis.full_length(Annotation(copies=[blind]), 100, 0.9) == []


class TestSelfBlast:
    @staticmethod
    @pytest.fixture(scope="class")
    def ltr_element():
        """A synthetic LTR retrotransposon: a 120 bp repeat at both termini."""
        rng = np.random.default_rng(7)
        bases = np.array(list("ACGT"))
        ltr = "".join(rng.choice(bases, 120))
        interior = "".join(rng.choice(bases, 600))
        return ltr + interior + ltr

    @staticmethod
    @pytest.fixture(scope="class")
    def tir_element():
        """A synthetic TIR element: terminal repeat, reverse-complemented at the 3' end."""
        rng = np.random.default_rng(11)
        bases = np.array(list("ACGT"))
        tir = "".join(rng.choice(bases, 120))
        interior = "".join(rng.choice(bases, 600))
        rc = tir[::-1].translate(str.maketrans("ACGT", "TGCA"))
        return tir + interior + rc

    def test_finds_the_trivial_diagonal(self, ltr_element):
        try:
            hits = analysis.self_blast(ltr_element)
        except analysis.BlastNotFound:
            pytest.skip("blastn not installed")
        assert any(h.is_trivial_diagonal for h in hits)

    def test_detects_direct_terminal_repeats_as_ltr_candidates(self, ltr_element):
        try:
            hits = analysis.self_blast(ltr_element)
        except analysis.BlastNotFound:
            pytest.skip("blastn not installed")
        terminal = analysis.terminal_repeats(hits, len(ltr_element))
        assert terminal["LTR"], "a 120 bp direct terminal repeat should be found"
        assert not terminal["TIR"]

    def test_detects_inverted_terminal_repeats_as_tir_candidates(self, tir_element):
        try:
            hits = analysis.self_blast(tir_element)
        except analysis.BlastNotFound:
            pytest.skip("blastn not installed")
        terminal = analysis.terminal_repeats(hits, len(tir_element))
        assert terminal["TIR"], "a 120 bp inverted terminal repeat should be found"
        assert not terminal["LTR"]

    def test_reverse_hits_keep_their_descending_orientation(self, tir_element):
        try:
            hits = analysis.self_blast(tir_element)
        except analysis.BlastNotFound:
            pytest.skip("blastn not installed")
        for hit in hits:
            if hit.is_reverse:
                assert hit.s_end < hit.s_start


class TestSequences:
    def test_matches_across_a_class_suffix(self, tmp_path):
        """A library header says 'fam#Unknown'; the annotation says 'fam'."""
        path = tmp_path / "lib.fa"
        path.write_text(">rnd-1_family-257#Unknown\nACGTACGTAC\n>other#LTR/Gypsy\nTTTT\n")
        library = read_fasta(path)
        assert library.get("rnd-1_family-257").bare_name == "rnd-1_family-257"
        assert library.get("rnd-1_family-257#Unknown") is not None
        assert library.get("RND-1_FAMILY-257") is not None

    def test_reads_class_label_and_length(self, tmp_path):
        path = tmp_path / "lib.fa"
        path.write_text(">fam#LTR/Gypsy\nACGT\nACGT\n")
        entry = read_fasta(path).get("fam")
        assert entry.class_label == "LTR/Gypsy"
        assert len(entry) == 8

    def test_suggests_near_misses(self, tmp_path):
        path = tmp_path / "lib.fa"
        path.write_text(">rnd-1_family-257\nACGT\n>rnd-1_family-258\nACGT\n")
        assert "rnd-1_family-257" in read_fasta(path).suggest("rnd-1_family-25")


class TestRecords:
    def test_absent_divergence_is_never_plottable(self):
        blind = Copy(chrom="c", start=0, end=10, strand="+", family="f")
        assert blind.has_divergence is False

    def test_zero_divergence_is_a_real_value(self):
        perfect = Copy(
            chrom="c",
            start=0,
            end=10,
            strand="+",
            family="f",
            divergence=0.0,
            divergence_kind=DivergenceKind.CONSENSUS,
        )
        assert perfect.has_divergence is True

    def test_transposed_interval_is_rejected_at_construction(self):
        with pytest.raises(ValueError):
            Copy(chrom="c", start=100, end=10, strand="+", family="f")

    def test_family_subsetting_ignores_class_suffix_when_inexact(self):
        annotation = Annotation(copies=[copy(0, 10, family="famA#LTR/Gypsy"), copy(0, 10, family="famB")])
        assert len(annotation.for_family("famA", exact=False)) == 1
        assert len(annotation.for_family("famA", exact=True)) == 0

    def test_divergence_kind_axis_labels_differ(self):
        labels = {k.axis_label for k in DivergenceKind}
        assert len(labels) == 3
