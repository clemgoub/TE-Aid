"""Tests for the protein-homology layer: BATH output, hit selection, class check."""

from __future__ import annotations

import pytest

from teaid import classcheck, homology, proteins
from teaid.homology import ProteinHit

# Real bathsearch --tblout output: PF00078 (RVT_1) against ltr-1_family-26.
# Row 2 is on the reverse strand, which bathsearch reports with descending
# alignment coordinates rather than a strand column.
TBLOUT = """\
# hit ID  target name         accession  query name           accession   hmm len  hmm from    hmm to   seq len  ali from    ali to    E-value  score  bias   PID  shifts  stops  description of target
#------- ------------------- ---------- -------------------- ---------- --------- --------- --------- --------- --------- ---------  --------- ------ ----- ----- ------- ------ ---------------------
       1 ltr-1_family-26      -          RVT_1                PF00078.33      200         1       200      11593      1350      1826   7.7e-40  127.0   0.0 27.00       0      0 -
       2 ltr-1_family-26      -          RVT_1                PF00078.33      200        57       179      11593      7611      7243     3e-08   23.8   0.0 17.69       2      1 -
#
# Program:         bathsearch
"""


def hit(**over) -> ProteinHit:
    base = dict(
        query="RVT_1", query_accession="PF00078.33", start=0, end=300, strand="+",
        evalue=1e-20, score=100.0, identity=30.0, hmm_from=1, hmm_to=200,
        hmm_len=200, frameshifts=0, stop_codons=0,
    )
    base.update(over)
    return ProteinHit(**base)


class TestTbloutParsing:
    def test_reads_every_data_row(self):
        assert len(homology._parse_tblout(TBLOUT)) == 2

    def test_coordinates_convert_to_half_open(self):
        first = homology._parse_tblout(TBLOUT)[0]
        assert (first.start, first.end) == (1349, 1826)
        assert first.coordinates_1based() == (1350, 1826)

    def test_descending_coordinates_mean_the_reverse_strand(self):
        """bathsearch has no strand column; direction is in the coordinate order."""
        second = homology._parse_tblout(TBLOUT)[1]
        assert second.strand == "-"
        assert second.start < second.end, "the interval is normalised ascending"
        assert (second.start, second.end) == (7242, 7611)

    def test_frameshifts_and_stops_are_carried(self):
        """The whole reason the search is BATH: a disrupted domain is reported
        *and* marked, where an ORF-finder-then-align search sees nothing."""
        intact, disrupted = homology._parse_tblout(TBLOUT)
        assert intact.frameshifts == 0 and intact.stop_codons == 0
        assert intact.is_disrupted is False
        assert disrupted.frameshifts == 2 and disrupted.stop_codons == 1
        assert disrupted.is_disrupted is True

    def test_hits_are_sorted_by_evalue(self):
        parsed = homology._parse_tblout(TBLOUT)
        assert parsed[0].evalue < parsed[1].evalue

    def test_model_coverage(self):
        first = homology._parse_tblout(TBLOUT)[0]
        assert first.coverage == pytest.approx(1.0)

    def test_comments_and_short_rows_are_skipped(self):
        assert homology._parse_tblout("# only a comment\n\nbroken row\n") == []


class TestBestPerRegion:
    def test_keeps_the_strongest_of_competing_hits(self):
        """RepeatClassifier's rule, borrowed for evidence selection only."""
        strong = hit(query="strong", score=200.0, start=0, end=300)
        weak = hit(query="weak", score=50.0, start=50, end=280)
        kept = homology.best_per_region([weak, strong])
        assert [h.query for h in kept] == ["strong"]

    def test_keeps_hits_that_do_not_overlap(self):
        a = hit(query="a", start=0, end=200)
        b = hit(query="b", start=400, end=600)
        assert len(homology.best_per_region([a, b])) == 2

    def test_slight_overlap_is_not_competition(self):
        a = hit(query="a", start=0, end=300, score=200.0)
        b = hit(query="b", start=280, end=600, score=100.0)
        assert len(homology.best_per_region([a, b])) == 2

    def test_result_is_ordered_along_the_consensus(self):
        kept = homology.best_per_region(
            [hit(query="late", start=800, end=900), hit(query="early", start=0, end=100)]
        )
        assert [h.query for h in kept] == ["early", "late"]

    def test_no_library_means_no_hits_rather_than_an_error(self):
        assert homology.search("ACGT", None) == []


class TestClassCheck:
    def test_class_i_seed_with_class_i_evidence_does_not_flag(self):
        tp = "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon"
        result = classcheck.check(tp, [hit(query_accession="PF00078")])
        assert result.expected == "I"
        assert result.disagrees is False

    def test_class_ii_seed_with_only_class_i_evidence_flags(self):
        """The example from the brief: TP says a DNA transposon, the only
        protein hit is a reverse transcriptase."""
        tp = "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT"
        result = classcheck.check(tp, [hit(query_accession="PF00078")])
        assert result.expected == "II"
        assert result.disagrees is True
        assert "Class II" in result.detail and "Class I" in result.detail

    def test_no_tp_never_flags(self):
        assert classcheck.check(None, [hit(query_accession="PF00078")]).disagrees is False

    def test_no_hits_never_flags(self):
        tp = "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition"
        result = classcheck.check(tp, [])
        assert result.disagrees is False
        assert "no protein hit" in result.detail

    def test_an_unclassifiable_tp_never_flags(self):
        assert classcheck.check("Interspersed_Repeat;Unknown",
                                [hit(query_accession="PF00078")]).disagrees is False

    def test_domesticated_domains_never_drive_a_contradiction(self):
        """A domesticated gag says TE ancestry, not how a live element moves."""
        tp = "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition"
        # PF03732 is the PEG10 capsid domain, te_order 'domesticated'.
        assert classcheck.check(tp, [hit(query_accession="PF03732")]).disagrees is False

    def test_mixed_evidence_including_the_expected_class_does_not_flag(self):
        tp = "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition"
        hits = [hit(query_accession="PF00078"), hit(query_accession="PF00872")]
        assert classcheck.check(tp, hits).disagrees is False

    def test_the_class_map_covers_every_order_it_should(self):
        orders = set(classcheck.domain_orders().values())
        unmapped = orders - set(classcheck._ORDER_CLASS)
        # DIRS/Crypton spans both worlds and 'domesticated' says nothing about
        # transposition, so both are deliberately unmapped.
        assert unmapped == {"DIRS/Crypton", "domesticated"}


class TestCuratedTable:
    def test_ships_with_the_package(self):
        accessions = proteins.curated_accessions()
        assert len(accessions) > 100
        assert all(a.startswith("PF") for a, _ in accessions)

    def test_names_are_present_for_every_accession(self):
        assert all(name for _, name in proteins.curated_accessions())

    def test_no_duplicate_accessions(self):
        accessions = [a for a, _ in proteins.curated_accessions()]
        assert len(accessions) == len(set(accessions))

    def test_the_pfam_blind_superfamilies_are_documented(self):
        """Each one is why RepeatPeps cannot simply be dropped."""
        assert set(proteins.PFAM_BLIND_SUPERFAMILIES) == {
            "DNA/PiggyBac", "DNA/Maverick", "DNA/Cryp", "LINE/Penelope"
        }
        assert all(reason for reason in proteins.PFAM_BLIND_SUPERFAMILIES.values())

    def test_extracts_only_the_pfam_blind_entries(self, tmp_path):
        source = tmp_path / "peps.lib"
        source.write_text(
            ">a_pol#LTR/Gypsy desc\nMKV\n"
            ">b_tnp#DNA/PiggyBac desc\nMKW\n"
            ">c_pol#LINE/Penelope desc\nMKY\n"
            ">d_tnp#DNA/hAT-Ac desc\nMKZ\n"
        )
        dest = tmp_path / "subset.fa"
        assert proteins.extract_pfam_blind(source, dest) == 2
        text = dest.read_text()
        assert "PiggyBac" in text and "Penelope" in text
        assert "Gypsy" not in text and "hAT-Ac" not in text

    def test_signature_changes_with_the_table(self, tmp_path, monkeypatch):
        """An edited domain table must rebuild rather than reuse a stale library."""
        before = proteins._signature(False, None)
        fake = tmp_path / "te_domains.tsv"
        fake.write_text("pfam_acc\tname\n PF99999\tmade up\n")
        monkeypatch.setattr(proteins, "TE_DOMAINS", fake)
        assert proteins._signature(False, None) != before

    def test_deep_and_default_get_different_caches(self):
        assert proteins._signature(True, None) != proteins._signature(False, None)
