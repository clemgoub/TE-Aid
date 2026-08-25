"""v1's TE-class colour scheme, carried forward.

These are the exact hex values v1 used for protein hits (``TE-Aid`` shell script,
the ``awk -v LINEcol=...`` block), kept because curators have been reading TE-Aid
sheets with them for years: green means LTR, blue means LINE, salmon means a DNA
transposon. Changing them would cost more in re-learning than any palette
improvement could return.

v1's own relabelling is preserved too, and the order matters — the specific
cases must be tested before the general ones they sit inside::

    DNA/Maverick   -> MAV      LINE/Penelope -> PLE
    DNA/Crypton    -> CRY      LTR/DIRS      -> DIRS
    DNA/*          -> TIR

Two things this scheme is not: it is not colourblind-validated (v1's palette
predates that concern here), and colour is never the only channel — every hit
carries its name in the tick gutter, so the scheme is a fast index rather than
the sole carrier of identity. One v1 bug is not carried over: ``Simple_repeat``
was defined as ``#8686ac`` in a table whose other values omitted the ``#``, so it
rendered as ``##8686ac`` and silently fell back to a default.
"""

from __future__ import annotations

# v1's palette, normalised to real CSS colours.
CLASS_COLOURS: dict[str, str] = {
    "LINE": "#3399ff",
    "SINE": "#800080",
    "TIR": "#ff6666",
    "LTR": "#00cc44",
    "RC": "#ff6600",
    "PLE": "#b2edba",
    "DIRS": "#fce7bd",
    "CRY": "#8f1800",
    "MAV": "#669999",
    "Low_complexity": "#d1d1e0",
    "Satellite": "#ff99ff",
    "Simple_repeat": "#8686ac",
    "Unknown": "#c9c9c9",
}

UNKNOWN = CLASS_COLOURS["Unknown"]

# te_order values in teaid/data/te_domains.tsv -> v1 class key.
#
# 'Class I' stays Unknown deliberately: a bare reverse transcriptase domain is
# shared by every Class I order, so colouring it LINE or LTR would assert an
# order the evidence does not support. 'domesticated' likewise says TE ancestry,
# not which element.
_ORDER_TO_CLASS = {
    "LTR": "LTR",
    "LTR/ERV": "LTR",
    "LINE": "LINE",
    "SINE": "SINE",
    "PLE": "PLE",
    "DIRS/Crypton": "DIRS",
    "TIR": "TIR",
    "TIR/Tc1-mariner": "TIR",
    "TIR/hAT": "TIR",
    "TIR/Mutator": "TIR",
    "TIR/CACTA": "TIR",
    "TIR/PIF-Harbinger": "TIR",
    "TIR/P-element": "TIR",
    "TIR/KDZ": "TIR",
    "RC/Helitron": "RC",
    "RC": "RC",
    "Class I": "Unknown",
    "domesticated": "Unknown",
}


def class_of_repeatmasker_label(label: str | None) -> str:
    """v1's class key for a RepeatMasker label such as ``LTR/Gypsy``.

    Order is v1's: the special cases are recognised before the general prefix
    they are nested inside, so ``DNA/Maverick`` is MAV rather than TIR and
    ``LINE/Penelope`` is PLE rather than LINE.
    """
    if not label:
        return "Unknown"
    text = label.strip()
    lowered = text.casefold()

    if "maverick" in lowered or "polinton" in lowered:
        return "MAV"
    if "crypton" in lowered or lowered.startswith("dna/cryp"):
        return "CRY"
    if "penelope" in lowered:
        return "PLE"
    if "dirs" in lowered:
        return "DIRS"

    head = text.split("/", 1)[0]
    if head == "DNA":
        return "TIR"
    if head in CLASS_COLOURS:
        return head
    for key in ("Low_complexity", "Satellite", "Simple_repeat"):
        if lowered.startswith(key.casefold()):
            return key
    return "Unknown"


def class_of_order(te_order: str | None) -> str:
    """v1's class key for a ``te_order`` from the curated Pfam table."""
    return _ORDER_TO_CLASS.get(te_order or "", "Unknown")


def colour(class_key: str) -> str:
    return CLASS_COLOURS.get(class_key, UNKNOWN)


def classify_hit(hit, orders: dict[str, str] | None = None) -> str:
    """v1 class key for one protein hit, whichever tier produced it.

    A tier-3 hit carries its RepeatMasker label directly. A tier-1/2 hit is a
    Pfam domain, so its order comes from the curated table.
    """
    label = getattr(hit, "source_class", None)
    if label:
        return class_of_repeatmasker_label(label)
    if orders:
        accession = (getattr(hit, "query_accession", "") or "").split(".")[0]
        order = orders.get(accession)
        if order:
            return class_of_order(order)
    return "Unknown"
