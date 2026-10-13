"""Compare a seed's declared class against the homology evidence on the sheet.

This is the one place TE-Aid *evaluates* rather than renders, and it is only
legitimate because both sides come from somewhere else: ``#=GF TP`` arrives with
the seed, and the domain-to-order mapping comes from the curated table. The
sheet still has no opinion of its own — it reports that two external claims
disagree (``docs/BRIEF_v2.md`` §1, §4 panel 8).

**The flag is the deliverable, not the label**, so this is deliberately
conservative: it fires only on a Class I / Class II contradiction, where the
seed says the element transposes via an RNA intermediate and every protein hit
says DNA, or the reverse. Finer disagreements — hAT versus Mutator, Copia versus
Gypsy — are left to the curator, because at that resolution a single domain hit
is not strong enough to contradict a curated label, and a flag that cries wolf
is worse than no flag.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

DATA = Path(__file__).parent / "data"
TE_DOMAINS = DATA / "te_domains.tsv"

# Tokens in a Dfam #=GF TP path that pin the transposition class.
_CLASS_I_TOKENS = ("class_i_retrotransposition", "retrotransposon", "retroposon", "ltr", "line", "sine")
_CLASS_II_TOKENS = ("class_ii_dna_transposition", "transposase", "helicase", "tir", "helitron")

# te_order values in the curated table, mapped to a transposition class.
_ORDER_CLASS = {
    "Class I": "I", "LTR": "I", "LTR/ERV": "I", "LINE": "I", "PLE": "I",
    "TIR": "II", "TIR/Tc1-mariner": "II", "TIR/hAT": "II", "TIR/Mutator": "II",
    "TIR/CACTA": "II", "TIR/PIF-Harbinger": "II", "TIR/P-element": "II",
    "TIR/KDZ": "II", "RC/Helitron": "II", "RC": "II",
    # Deliberately unmapped: DIRS/Crypton use a tyrosine recombinase and span
    # both worlds, and 'domesticated' says nothing about how a live element
    # would move. Neither should ever drive a contradiction.
}

# Tier-2 hits come from RepeatPeps rather than Pfam, so they have no Pfam
# accession to look up — but their RepeatMasker class label states the class
# directly. Without this the four superfamilies tier 2 exists to cover
# (piggyBac, Maverick, Crypton, Penelope) could never reach the check at all.
_REPEATPEPS_CLASS = {
    "DNA": "II",
    "RC": "II",
    "LTR": "I",
    "LINE": "I",
    "SINE": "I",
}


@dataclass(slots=True)
class ClassCheck:
    expected: str | None  # the transposition class the seed's TP implies
    observed: set[str]  # classes implied by the protein hits
    disagrees: bool
    detail: str


def domain_orders() -> dict[str, str]:
    """Accession -> te_order, from the curated table."""
    if not TE_DOMAINS.exists():
        return {}
    out = {}
    for line in TE_DOMAINS.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 3:
            out[parts[0].split(".")[0]] = parts[2]
    return out


def class_of_tp(tp: str | None) -> str | None:
    """The transposition class a Dfam ``#=GF TP`` path implies, if any."""
    if not tp:
        return None
    tokens = [t.strip().casefold() for t in tp.split(";")]
    for token in tokens:
        if any(k in token for k in _CLASS_II_TOKENS):
            return "II"
    for token in tokens:
        if any(k in token for k in _CLASS_I_TOKENS):
            return "I"
    return None


def _class_of_repeatpeps(label: str | None) -> str | None:
    """Transposition class implied by a RepeatMasker class label, e.g. 'DNA/hAT-Ac'.

    Crypton is deliberately excluded: it is filed under ``DNA`` but integrates
    with a tyrosine recombinase, so it should not be evidence for or against
    either class.
    """
    if not label:
        return None
    head, _, rest = label.partition("/")
    if head == "DNA" and rest.startswith("Cryp"):
        return None
    return _REPEATPEPS_CLASS.get(head)


def check(tp: str | None, hits) -> ClassCheck:
    """Flag a Class I / Class II contradiction between the seed and the evidence."""
    expected = class_of_tp(tp)
    orders = domain_orders()

    observed: set[str] = set()
    names: dict[str, list[str]] = {}
    for hit in hits:
        accession = (getattr(hit, "query_accession", "") or "").split(".")[0]
        klass = _ORDER_CLASS.get(orders.get(accession) or "")
        if klass is None:
            klass = _class_of_repeatpeps(getattr(hit, "source_class", None))
        if klass:
            observed.add(klass)
            names.setdefault(klass, []).append(
                getattr(hit, "display_name", None) or getattr(hit, "query", accession)
            )

    if expected is None:
        return ClassCheck(None, observed, False,
                          "the seed's TP does not pin a transposition class")
    if not observed:
        return ClassCheck(expected, observed, False,
                          "no protein hit maps to a transposition class")

    if expected in observed:
        return ClassCheck(expected, observed, False,
                          f"the evidence includes Class {expected} domains, as the seed says")

    contrary = sorted(observed)[0]
    examples = ", ".join(dict.fromkeys(names.get(contrary, [])))[:80]
    return ClassCheck(
        expected, observed, True,
        f"the seed declares Class {expected}, but every domain that maps to a class "
        f"is Class {contrary} ({examples})",
    )


# --- Order level, for the protein benchmark (bench/PROTOCOL.md) -------------
#
# The sheet's flag above judges a family as a whole, at class level. The
# benchmark judges each hit, and at order level where both sides name an order:
# a LINE reverse transcriptase on an LTR family, or a TIR transposase on a
# Helitron, is wrong at order level but not at class level. Nothing here
# reaches the sheet.

# Each order's transposition class.
ORDER_CLASS = {
    "LTR": "I", "DIRS": "I", "LINE": "I", "PLE": "I", "SINE": "I",
    "TIR": "II", "Helitron": "II", "Maverick": "II", "Crypton": "II",
}

# te_order values in the curated domain table that name an order. "Class I"
# names only a class; "domesticated" and the shared tyrosine recombinase of
# DIRS and Crypton name neither, so hits from them are not judged.
_TE_ORDER_ORDER = {
    "LTR": "LTR", "LTR/ERV": "LTR", "LINE": "LINE", "PLE": "PLE",
    "RC": "Helitron", "RC/Helitron": "Helitron",
}
_TE_ORDER_CLASS_ONLY = {"Class I": "I"}


def order_of_label(label: str | None) -> tuple[str | None, str | None]:
    """(order, class) of a RepeatMasker or Dfam label such as ``LTR/Gypsy``.

    Either can be None: ``DNA`` alone gives a class but no order, and
    ``Unknown``, satellites and RNA genes give neither.
    """
    if not label:
        return None, None
    head, _, sub = label.partition("/")
    sub = sub.casefold()
    if head == "LTR":
        order = "DIRS" if sub.startswith(("dirs", "ngaro", "viper")) else "LTR"
    elif head == "LINE":
        order = "PLE" if sub.startswith("penelope") else "LINE"
    elif head == "PLE":
        order = "PLE"
    elif head == "SINE":
        order = "SINE"
    elif head == "DNA":
        if sub.startswith("maverick") or sub.startswith("polinton"):
            order = "Maverick"
        elif sub.startswith("crypton"):
            order = "Crypton"
        elif sub.startswith("helitron"):
            order = "Helitron"
        else:
            order = "TIR" if sub else None
    elif head == "RC":
        order = "Helitron"
    elif head == "Retroposon":
        return None, "I"
    else:
        return None, None
    if order is None:
        return None, "II"
    return order, ORDER_CLASS[order]


def order_of_hit(hit, orders: dict[str, str] | None = None) -> tuple[str | None, str | None]:
    """(order, class) a protein hit stands for.

    A Pfam domain takes its te_order from the curated table; a RepeatPeps
    protein takes the label after the ``#`` in its name.
    """
    if orders is None:
        orders = domain_orders()
    accession = (getattr(hit, "query_accession", "") or "").split(".")[0]
    te_order = orders.get(accession)
    if te_order is not None:
        if te_order in _TE_ORDER_CLASS_ONLY:
            return None, _TE_ORDER_CLASS_ONLY[te_order]
        if te_order.startswith("TIR"):
            return "TIR", "II"
        order = _TE_ORDER_ORDER.get(te_order)
        return (order, ORDER_CLASS[order]) if order else (None, None)
    return order_of_label(getattr(hit, "source_class", None))


@dataclass(frozen=True, slots=True)
class OrderCall:
    """Whether one hit agrees with its family, and at which level that was judged."""

    verdict: str | None  # "agree", "contradict", or None when not judged
    level: str | None  # "order" or "class"
    family_order: str | None
    hit_order: str | None


def judge(family_label: str | None, hit, orders: dict[str, str] | None = None) -> OrderCall:
    """Judge a hit against its family's label, at order level where both sides allow.

    A hit that names only a class (a bare reverse transcriptase), or a family
    whose label names only a class, is judged at class level. A family of
    unknown class, and a hit that names neither, are not judged.
    """
    family_order, family_class = order_of_label(family_label)
    hit_order, hit_class = order_of_hit(hit, orders)
    if family_class is None or hit_class is None:
        return OrderCall(None, None, family_order, hit_order)
    if family_order and hit_order:
        verdict = "agree" if family_order == hit_order else "contradict"
        return OrderCall(verdict, "order", family_order, hit_order)
    verdict = "agree" if family_class == hit_class else "contradict"
    return OrderCall(verdict, "class", family_order, hit_order)
