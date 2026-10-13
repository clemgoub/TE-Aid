"""Order-level judging for the protein benchmark (bench/PROTOCOL.md)."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from teaid import classcheck


def pfam(accession):
    return SimpleNamespace(query_accession=accession, source_class=None)


def repeatpeps(label):
    return SimpleNamespace(query_accession="", source_class=label)


ORDERS = {"PF00078": "Class I", "PF00665": "LTR", "PF03184": "TIR/Mutator",
          "PF14529": "RC/Helitron", "PF00000": "domesticated", "PF99999": "DIRS/Crypton"}


@pytest.mark.parametrize("label, expected", [
    ("LTR/Gypsy", ("LTR", "I")),
    ("LTR/DIRS", ("DIRS", "I")),
    ("LTR/Ngaro", ("DIRS", "I")),
    ("LINE/L1", ("LINE", "I")),
    ("LINE/Penelope", ("PLE", "I")),
    ("PLE", ("PLE", "I")),
    ("SINE/Alu", ("SINE", "I")),
    ("DNA/hAT-Ac", ("TIR", "II")),
    ("DNA/PiggyBac", ("TIR", "II")),
    ("DNA/Maverick", ("Maverick", "II")),
    ("DNA/Crypton-A", ("Crypton", "II")),
    ("RC/Helitron", ("Helitron", "II")),
    ("DNA", (None, "II")),
    ("Retroposon/SVA", (None, "I")),
    ("Unknown", (None, None)),
    ("Satellite/centromeric", (None, None)),
    ("tRNA", (None, None)),
    (None, (None, None)),
])
def test_labels_map_to_an_order_and_a_class(label, expected):
    assert classcheck.order_of_label(label) == expected


def test_a_maverick_protein_on_an_ltr_family_contradicts_at_order_level():
    """The case a curator reported: Class II only because of where Maverick is filed."""
    call = classcheck.judge("LTR/Gypsy", repeatpeps("DNA/Maverick"), ORDERS)
    assert (call.verdict, call.level) == ("contradict", "order")


def test_an_integrase_on_an_ltr_family_agrees():
    call = classcheck.judge("LTR/Copia", pfam("PF00665.20"), ORDERS)
    assert (call.verdict, call.level, call.hit_order) == ("agree", "order", "LTR")


def test_a_bare_reverse_transcriptase_is_judged_at_class_level():
    agree = classcheck.judge("LINE/L2", pfam("PF00078"), ORDERS)
    contradict = classcheck.judge("DNA/hAT", pfam("PF00078"), ORDERS)
    assert (agree.verdict, agree.level) == ("agree", "class")
    assert (contradict.verdict, contradict.level) == ("contradict", "class")


def test_a_family_labelled_only_by_class_is_judged_at_class_level():
    call = classcheck.judge("DNA", pfam("PF03184"), ORDERS)
    assert (call.verdict, call.level) == ("agree", "class")


def test_any_tir_superfamily_domain_counts_as_tir():
    assert classcheck.judge("DNA/MULE-MuDR", pfam("PF03184"), ORDERS).verdict == "agree"


def test_helitron_domains_agree_with_rc_families():
    assert classcheck.judge("RC/Helitron", pfam("PF14529"), ORDERS).verdict == "agree"


@pytest.mark.parametrize("label, hit", [
    ("Unknown", pfam("PF00665")),
    ("LTR/Gypsy", pfam("PF00000")),  # domesticated
    ("LTR/Gypsy", pfam("PF99999")),  # DIRS/Crypton recombinase names neither
    ("LTR/Gypsy", repeatpeps(None)),
])
def test_nothing_is_judged_without_a_class_on_both_sides(label, hit):
    assert classcheck.judge(label, hit, ORDERS).verdict is None


def test_the_shipped_table_maps_every_te_order_it_uses():
    """A te_order the judge cannot read would silently leave its hits unjudged."""
    readable = {"domesticated", "DIRS/Crypton"}
    for accession, te_order in classcheck.domain_orders().items():
        order, klass = classcheck.order_of_hit(pfam(accession))
        assert klass is not None or te_order in readable, (accession, te_order)
