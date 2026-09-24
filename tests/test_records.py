"""Insertions: the unit panels 1 and 2 count, built from fragments sharing an ID."""

from __future__ import annotations

import pytest

from teaid.records import (Annotation, Copy, DivergenceKind, Insertion, JunctionKind,
                           JUNCTION_TOLERANCE)


def frag(start, end, cs, ce, *, strand="+", div=10.0, hit_id="1", chrom="chr1",
         kind=DivergenceKind.CONSENSUS):
    return Copy(chrom=chrom, start=start, end=end, strand=strand, family="fam",
                divergence=div, divergence_kind=kind, consensus_start=cs,
                consensus_end=ce, consensus_left=5000 - ce, hit_id=hit_id)


class TestInsertion:
    def test_a_lone_row_is_an_insertion_of_one_fragment(self):
        one = Insertion([frag(100, 400, 0, 300)])
        assert one.n_fragments == 1 and not one.is_fragmented
        assert one.consensus_blocks == [(0, 300)] and one.junctions == []

    def test_overlapping_fragments_merge_into_one_block(self):
        both = Insertion([frag(0, 300, 0, 300), frag(1000, 1300, 250, 550)])
        assert both.consensus_blocks == [(0, 550)]
        assert both.consensus_covered == 550

    def test_divergence_is_weighted_by_consensus_extent_as_in_processrepeats(self):
        # 900 bp at 5% and 100 bp at 25%: weighted 7.0, not the plain mean 15.0
        # and not the longest fragment's 5.0.
        both = Insertion([frag(0, 900, 0, 900, div=5.0), frag(2000, 2100, 900, 1000, div=25.0)])
        assert both.divergence == pytest.approx(7.0)

    def test_a_fragment_without_divergence_does_not_vote(self):
        both = Insertion([frag(0, 100, 0, 100, div=8.0),
                          frag(500, 900, 100, 500, div=None, kind=DivergenceKind.NONE)])
        assert both.divergence == pytest.approx(8.0)
        assert both.has_divergence

    def test_span_and_union_differ_for_a_split_insertion(self):
        """Pins the full-length decision: a 5' and a 3' piece of a 5 kb element
        span almost all of it and cover 400 bp."""
        ends = Insertion([frag(0, 200, 0, 200), frag(900, 1100, 4800, 5000)])
        assert ends.consensus_end - ends.consensus_start == 5000
        assert ends.consensus_covered == 400


class TestJunctions:
    def test_extra_genomic_sequence_is_an_insertion(self):
        # Consensus continues 300 -> 300 while the genome skips 1,000 bp.
        nested = Insertion([frag(0, 300, 0, 300), frag(1300, 1600, 300, 600)])
        (j,) = nested.junctions
        assert j.kind is JunctionKind.INSERTION and j.size == 1000
        assert j.consensus_at == j.consensus_to == 300

    def test_missing_consensus_is_a_deletion(self):
        # The genome continues while the consensus skips 300 -> 800.
        deleted = Insertion([frag(0, 300, 0, 300), frag(300, 500, 800, 1000)])
        (j,) = deleted.junctions
        assert j.kind is JunctionKind.DELETION and j.size == 500
        assert (j.consensus_at, j.consensus_to) == (300, 800)

    def test_a_join_within_tolerance_is_flush(self):
        split = Insertion([frag(0, 300, 0, 300),
                           frag(320, 600, 300 + JUNCTION_TOLERANCE, 580 + JUNCTION_TOLERANCE)])
        (j,) = split.junctions
        assert j.kind is JunctionKind.FLUSH and j.size == 0

    def test_a_minus_strand_element_runs_backwards_through_the_consensus(self):
        """The later fragment in the genome sits earlier in the consensus. Read
        in genomic order without that, this deletion looks like an insertion."""
        minus = Insertion([frag(0, 200, 800, 1000, strand="-"),
                           frag(200, 500, 0, 300, strand="-")])
        (j,) = minus.junctions
        assert j.kind is JunctionKind.DELETION and j.size == 500
        assert (j.consensus_at, j.consensus_to) == (300, 800)

    def test_fragments_that_disagree_on_strand_give_no_junction(self):
        mixed = Insertion([frag(0, 300, 0, 300), frag(1300, 1600, 300, 600, strand="-")])
        assert mixed.strand == "." and mixed.junctions == []


class TestAnnotationInsertions:
    def test_rows_sharing_an_id_become_one_insertion(self):
        annotation = Annotation(copies=[
            frag(0, 300, 0, 300, hit_id="7"), frag(1300, 1600, 300, 600, hit_id="7"),
            frag(5000, 5300, 0, 300, hit_id="8"),
        ])
        assert [i.n_fragments for i in annotation.insertions()] == [2, 1]

    def test_group_false_counts_every_row(self):
        annotation = Annotation(copies=[
            frag(0, 300, 0, 300, hit_id="7"), frag(1300, 1600, 300, 600, hit_id="7"),
        ])
        assert [i.n_fragments for i in annotation.insertions(group=False)] == [1, 1]

    def test_insertions_are_ordered_by_locus(self):
        annotation = Annotation(copies=[
            frag(9000, 9300, 0, 300, hit_id="2", chrom="chr2"),
            frag(5000, 5300, 0, 300, hit_id="3"), frag(100, 400, 0, 300, hit_id="4"),
        ])
        assert [(i.chrom, i.start) for i in annotation.insertions()] == [
            ("chr1", 100), ("chr1", 5000), ("chr2", 9000)]
