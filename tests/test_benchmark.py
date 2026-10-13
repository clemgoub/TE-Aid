"""tools/benchmark.py: the parts of the protein benchmark that need no search tools."""

from __future__ import annotations

import importlib.util
import random
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location("benchmark", ROOT / "tools" / "benchmark.py")
bm = importlib.util.module_from_spec(_spec)
sys.modules["benchmark"] = bm  # dataclasses look their module up there
_spec.loader.exec_module(bm)


def hit(family="f", label="LTR/Gypsy", start=0, end=300, evalue=1e-20, score=50.0,
        coverage=0.9, model="rve", model_class=None, accession="PF00665", source="pfam"):
    return bm.Hit(family, label, model, model_class, accession, source, start, end,
                  evalue, score, coverage)


class TestShuffle:
    def test_keeps_length_composition_and_every_dinucleotide(self):
        rng = random.Random(1)
        seq = "".join(rng.choice("ACGT") for _ in range(2000))
        out = bm.dinucleotide_shuffle(seq, random.Random(7))
        pairs = lambda s: Counter(s[i:i + 2] for i in range(len(s) - 1))
        assert len(out) == len(seq) and out != seq
        assert pairs(out) == pairs(seq)
        assert out[0] == seq[0] and out[-1] == seq[-1]

    def test_the_seed_fixes_the_result(self):
        seq = "ACGTTGCAACGGTACCATGA" * 20
        assert bm.dinucleotide_shuffle(seq, random.Random(3)) == bm.dinucleotide_shuffle(seq, random.Random(3))


class TestWilson:
    def test_all_found_has_the_closed_form_lower_bound(self):
        """For k = n the lower bound is n / (n + z²)."""
        lo, hi = bm.wilson(20, 20)
        assert lo == pytest.approx(20 / (20 + 1.959964 ** 2))
        assert hi == pytest.approx(1.0)

    def test_no_regions_bounds_nothing(self):
        assert bm.wilson(0, 0) == (0.0, 1.0)


class TestCompare:
    def test_a_region_found_in_two_pieces_reports_its_union(self):
        ref = {"f": [hit(start=0, end=1000)]}
        arm = {"f": [hit(start=0, end=300, score=40), hit(start=500, end=900, score=30)]}
        (row,) = bm.compare(ref, arm)
        assert row["found"] and row["pieces"] == 2
        assert row["span"] == pytest.approx(0.7)

    def test_a_hit_overlapping_less_than_half_of_the_shorter_does_not_count(self):
        ref = {"f": [hit(start=0, end=1000)]}
        arm = {"f": [hit(start=900, end=1300)]}  # 100 shared of a 400 bp hit
        (row,) = bm.compare(ref, arm)
        assert not row["found"] and row["span"] == 0

    def test_hits_above_the_e_value_are_not_regions(self):
        ref = {"f": [hit(evalue=0.5)]}
        assert bm.compare(ref, ref) == []


class TestSpecificity:
    def test_the_rule_picks_the_floor_that_removes_a_short_contradicting_hit(self):
        hits = {
            "a": [hit(family="a", start=0, end=300, coverage=0.9),
                  hit(family="a", start=1000, end=1100, coverage=0.25, model="Polinton-1_int",
                      model_class="DNA/Maverick", accession="", score=30)],
            "b": [hit(family="b", start=0, end=300, coverage=0.8)],
        }
        result = bm.specificity(hits, dict(bm.DRAFT_THRESHOLDS))
        rows = {r["floor"]: r for r in result["grid"] if r["evalue"] == 1e-3}
        assert (rows[0.0]["agree"], rows[0.0]["contradict"]) == (2, 1)
        assert rows[0.3]["contradict"] == 0 and rows[0.3]["kept"] == 1.0
        assert result["chosen"]["0.001"] == 0.3

    def test_no_floor_meets_both_when_agreeing_hits_go_first(self):
        hits = {"a": [hit(family="a", start=0, end=300, coverage=0.2),
                      hit(family="a", start=1000, end=1300, coverage=0.95, model="Polinton-1_int",
                          model_class="DNA/Maverick", accession="")]}
        assert bm.specificity(hits, dict(bm.DRAFT_THRESHOLDS))["chosen"]["0.001"] is None


class TestCriteria:
    @pytest.mark.parametrize("label, group", [
        ("DNA/PiggyBac", "piggyBac"), ("DNA/Maverick", "Maverick"), ("DNA/Crypton-A", "Crypton"),
        ("LINE/Penelope", "Penelope"), ("PLE", None), ("LTR/Gypsy", None),
    ])
    def test_pfam_blind_groups(self, label, group):
        assert bm.pfam_blind_group(label) == group

    def test_a3_cannot_pass_with_nothing_to_test(self):
        regions = [{"label": "LTR/Gypsy", "found": True, "span": 1.0}]
        result = bm.criteria(regions, {"ratio": 3.0}, dict(bm.DRAFT_THRESHOLDS))
        assert result["A1"]["found"] == 1 and not result["A3"]["pass"]
        assert result["A4"]["pass"]


def test_a_batched_search_is_rescaled_to_each_consensus_alone():
    """E scales with database size: one consensus of L bp is L / 1e6 of a 1 Mb batch."""
    from dataclasses import dataclass

    @dataclass(frozen=True)
    class H:
        target: str
        evalue: float

    found = [H("short", 2000.0), H("short", 1e5), H("long", 0.5)]
    out = bm.rescale(found, {"short": 100, "long": 2000})
    assert [h.evalue for h in out["short"]] == [pytest.approx(0.2)]
    assert [h.evalue for h in out["long"]] == [pytest.approx(0.001)]


class TestGuard:
    def test_results_are_refused_until_the_protocol_is_signed(self, monkeypatch):
        monkeypatch.setattr(bm, "protocol_state", lambda: {
            "commit": None, "committed_at": None, "modified": True, "signed": False})
        with pytest.raises(SystemExit, match="refusing to write under bench/results/"):
            bm.guard(bm.RESULTS / "run1")

    def test_other_directories_are_allowed(self, tmp_path, monkeypatch):
        monkeypatch.setattr(bm, "protocol_state", lambda: {
            "commit": None, "committed_at": None, "modified": True, "signed": False})
        bm.guard(tmp_path / "scratch")

    def test_a_signed_committed_protocol_lets_results_through(self, monkeypatch):
        monkeypatch.setattr(bm, "protocol_state", lambda: {
            "commit": "abc", "committed_at": 1, "modified": False, "signed": True})
        bm.guard(bm.RESULTS / "run1")


def test_hits_round_trip_with_their_class_label(tmp_path):
    arm = tmp_path / "default"
    arm.mkdir()
    (arm / "hits.tsv").write_text(
        "\t".join(bm.HIT_FIELDS) + "\n"
        "f\tLTR/Gypsy\tPolinton-1_int\tDNA/Maverick\t\tpfam\t10\t310\t+\t1e-05\t30.0\t1\t50\t200\t0.25\t0\t0\n")
    (arm / "timing.tsv").write_text("family\tseconds\thits\nf\t1.0\t1\ng\t1.0\t0\n")
    hits = bm.read_hits(arm / "hits.tsv")
    assert set(hits) == {"f", "g"} and hits["g"] == []
    (h,) = hits["f"]
    assert h.source_class == "DNA/Maverick" and h.coverage == 0.25
    assert bm.classcheck.judge(h.label, h).verdict == "contradict"
