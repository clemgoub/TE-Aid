#!/usr/bin/env python3
"""Build ``teaid/data/te_domains.tsv`` — the curated Pfam TE-domain set.

Run with a Pfam metadata cache produced by ``tools/fetch_pfam.py``. Every
accession's *name* is taken from InterPro, never from this file, so a typo or a
dead accession fails the build instead of shipping.

Why this table exists at all, and why it is not a clan selection: see
``docs/BRIEF_v2.md`` §5.4. Selecting by clan fails in both directions at once,
and the fetched clan memberships quantify it:

    CL0219 RNase_H          121 members, ~75 of them host enzymes
                            (RNase H, exonucleases, Argonaute/Piwi, Cas9, RuvC,
                            DNA pol subunits, the spliceosome's PRP8 …)
    CL0027 RdRP              15 members, 13 of them *viral RNA-dependent RNA
                            polymerases* — only PF00078 and PF07727 are RTs
    CL0169 Rep               24 members, 22 of them plasmid/virus replication
                            proteins — only the Helitron and Replitron domains
                            are transposon-derived
    CL0523 GAG-polyprotein   16 members, ~6 of them domesticated host genes
                            (PEG10, Arc, PNMA, RTL1, LDOC)

and in the other direction the hAT dimerisation domain and the Helitron-associated
DUF4216 sit in no clan at all.

The inclusion test is one question: **would a hit here, on a TE consensus, be
evidence of transposable-element origin?** Host enzymes fail it — a hit to
cellular RNase H manufactures TE evidence where there is none. Domains from
domesticated, TE-derived host genes *pass* it, but are labelled ``domesticated``
because a hit means "TE-derived", which may be a live element or a host gene
built from one — a distinction the curator needs to make, and which the sheet
must not make for them.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "teaid" / "data" / "te_domains.tsv"

# ---------------------------------------------------------------------------
# The curated set: accession -> (te_order, why_included)
#
# te_order follows Wicker et al. where it applies. 'Class I' / 'Class II' mean
# the domain is shared across an order's families. 'domesticated' marks host
# genes of TE origin. Rationales say what the domain *is*, not what it is called.
# ---------------------------------------------------------------------------
CURATED: dict[str, tuple[str, str]] = {}


def add(order: str, why: str, *accessions: str) -> None:
    for acc in accessions:
        if acc in CURATED:
            raise SystemExit(f"{acc} listed twice")
        CURATED[acc] = (order, why)


# --- Class I: the reverse transcriptase core ------------------------------
add("Class I", "reverse transcriptase catalytic core; the defining enzyme of retrotransposition",
    "PF00078", "PF07727", "PF13456")
add("Class I", "flanking domains of the RT fold, retained in diverged elements where the catalytic core has drifted",
    "PF13655", "PF06817", "PF06815", "PF13966")
add("Class I", "RNase H domain *of a reverse transcriptase* — the retroelement-specific model, unlike the cellular RNase H of PF00075",
    "PF17917", "PF17919", "PF05380")

# --- Class I / LTR: integrase, protease, gag ------------------------------
add("LTR", "DDE integrase that inserts the cDNA; the Class I counterpart of a DNA transposase",
    "PF00665", "PF13333", "PF13683", "PF24764")
add("LTR", "integrase accessory domains: zinc-binding and DNA-binding modules flanking the catalytic core",
    "PF02022", "PF17921", "PF00552", "PF18103", "PF18697")
add("LTR", "aspartyl protease that cleaves the pol polyprotein; retroelement-specific, unlike the eukaryotic aspartyl proteases of PF00026",
    "PF00077", "PF08284", "PF12382", "PF12384", "PF13975")
add("LTR", "gag capsid protein — the structural half of a retroelement, and often the only ORF left in a decayed copy",
    "PF01021", "PF14223", "PF14244", "PF17241", "PF19259", "PF13976")
add("LTR", "gag-associated domain of the copia/gypsy lineages, named as a DUF but confined to retroelement gag",
    "PF03564", "PF13961", "PF23055", "PF23309")
add("LTR/ERV", "retroviral structural proteins: capsid, matrix, nucleocapsid and envelope of endogenous retroviruses",
    "PF00607", "PF19317", "PF30021", "PF02337", "PF02813", "PF08705", "PF08723",
    "PF03276", "PF00517", "PF25597")
add("LTR", "lineage-specific polyprotein segments of named LTR families",
    "PF04195", "PF24614", "PF29702", "PF29703", "PF29829", "PF20167", "PF10599")

# --- Class I / LINE -------------------------------------------------------
add("LINE", "L1 ORF1 RNA-binding domain — the nucleic-acid chaperone of the L1 machinery",
    "PF02994", "PF29741")
add("LINE", "apurinic endonuclease that nicks the target site for target-primed reverse transcription. "
            "Lower specificity: the fold is shared with cellular AP endonucleases, so treat a lone hit as weak",
    "PF03372", "PF14529")
add("LINE", "retrotransposon-associated domain of the hot-spot protein families",
    "PF20445")

# --- Class I / PLE and DIRS ----------------------------------------------
add("PLE", "GIY-YIG endonuclease of Penelope-like elements. Lower specificity: the fold also occurs in "
           "cellular UvrC-type nucleases, so corroborate with an RT hit",
    "PF01541")
add("DIRS/Crypton", "tyrosine recombinase used for integration instead of a DDE enzyme. Lower specificity: "
                    "shared with phage and bacterial integrases, so corroborate",
    "PF00589")

# --- Class II / TIR: the DDE transposases ---------------------------------
add("TIR", "DDE transposase catalytic core; the cut-and-paste enzyme of Class II elements",
    "PF01609", "PF13843", "PF13701", "PF13737", "PF13751", "PF13586", "PF13612",
    "PF13610", "PF03184", "PF13358", "PF13359", "PF13546")
add("TIR", "domains flanking a DDE core: zinc fingers, HTH and associated modules that persist when the catalytic core is degraded",
    "PF13808", "PF13842", "PF13613", "PF13963", "PF14706", "PF14319", "PF12760")
add("TIR/Tc1-mariner", "Tc1/mariner transposase and its DNA-binding HTH; the superfamily behind Sleeping Beauty and Mos1",
    "PF01359", "PF17906", "PF25787", "PF11427", "PF21517", "PF01498", "PF03221", "PF04236")
add("TIR/hAT", "hAT superfamily: the C-terminal dimerisation region and the RNase-H-like catalytic fold. "
               "The dimerisation domain sits in no Pfam clan, so clan-based selection misses it entirely",
    "PF05699", "PF14372", "PF10683")
add("TIR/Mutator", "MULE/Mutator transposase and its accessory domains",
    "PF00872", "PF10551", "PF20700", "PF03108", "PF18221")
add("TIR/CACTA", "CACTA (En/Spm) transposase and the plant transposon proteins of that superfamily",
    "PF03017", "PF04827")
add("TIR/PIF-Harbinger", "Harbinger-associated domain, carried by the second ORF of the superfamily",
    "PF04937")
add("TIR/P-element", "P-element transposase and the THAP9 domesticated copy of it",
    "PF12017", "PF12596", "PF22824")
add("TIR/KDZ", "Kyakuja-Dileera-Zisupton transposases and the CxC cysteine clusters that accompany them",
    "PF18758", "PF18759", "PF18717", "PF18718", "PF18721", "PF18802", "PF18803",
    "PF18804", "PF18866")
add("TIR", "transposase families described only as domains of unknown function, but confined to elements",
    "PF25273", "PF21787", "PF21789", "PF21804")

# --- Class II / rolling circle -------------------------------------------
add("RC/Helitron", "Helitron helicase-like domain and the Helitron-associated DUF4216, which is in no clan. "
                   "Rolling-circle transposition leaves no DDE core, so these are the only catalytic evidence",
    "PF14214", "PF13952")
add("RC", "HUH endonuclease of Replitrons and the rolling-circle initiator fold they share",
    "PF21859", "PF18106")

# --- Domesticated: TE-derived host genes ---------------------------------
add("domesticated", "host gene built from a transposase. A hit means TE ancestry, which may be a live element "
                    "*or* a domesticated gene — the curator decides which, and the sheet must not",
    "PF26100", "PF14291", "PF27039", "PF27041", "PF27046", "PF27073")
add("domesticated", "host gene built from a retroelement gag. Same caveat: TE ancestry, not necessarily a live element",
    "PF03732", "PF30901", "PF16297", "PF14893", "PF18162", "PF21395", "PF29000", "PF29013")

# ---------------------------------------------------------------------------
# Deliberate exclusions worth recording, so the reasoning survives review.
# ---------------------------------------------------------------------------
EXCLUDED: dict[str, str] = {
    "PF00075": "cellular RNase H. Retroelements do carry an RNase H, but this model matches the host enzyme "
               "too; the RT-specific models PF17917/PF17919 are included instead",
    "PF00098": "zinc knuckle. Present in retroelement gag, but far more abundant in host RNA-binding proteins",
    "PF00136": "DNA polymerase family B. Mavericks/Polintons carry one, but so does every eukaryotic replisome",
    "PF00026": "eukaryotic aspartyl protease — the host counterpart of the retroviral protease already included",
    "PF01693": "RNase H1 N-terminal domain — host",
    "PF08615": "RNase H2 non-catalytic subunit — host",
    "PF31260": "RNase H2 subunit C — host",
    "PF11474": "telomerase reverse transcriptase. Telomerase is a host RT of retroelement ancestry, but a hit "
               "would flag the host telomerase machinery, not an element",
    "PF17984": "telomerase RT thumb domain — as above",
    "PF21399": "telomerase RT C-terminal extension — as above",
    "PF23119": "telomerase RT CTE domain — as above",
    "PF02171": "Piwi domain — host silencing machinery, and the pathway that *represses* TEs",
    "PF13017": "piRNA pathway germ-plasm component — host silencing machinery",
    "PF22474": "Argonaute middle domain — host",
    "PF02075": "RuvC Holliday-junction resolvase — host recombination",
    "PF22702": "Cas9 RuvC domain — prokaryotic immunity",
    "PF22126": "C2c1 CRISPR-Cas endonuclease — prokaryotic immunity",
    "PF12134": "PRP8 domain IV — the spliceosome",
    "PF06333": "Mediator complex subunit 13 — transcription machinery",
    "PF00929": "exonuclease — host",
    "PF01612": "3'-5' exonuclease — host",
    "PF01351": "ribonuclease HII — host",
    "PF03104": "DNA polymerase B exonuclease domain — host",
    "PF05188": "MutS domain II — host mismatch repair",
    "PF04493": "endonuclease V — host",
    "PF06550": "presenilin aspartyl protease — host",
    "PF16641": "CLIP1 zinc knuckle — host",
    "PF21890": "Lin-28A zinc knuckle — host",
    "PF03455": "dDENN domain — host",
}

# Whole categories left out of the default set, with the reason.
EXCLUDED_CATEGORIES = {
    "prokaryotic IS and phage": (
        "IS transposases, phage integrases and conjugative-transposon proteins (PF01527, PF01710, PF01797, "
        "PF02371, PF13005-PF13007, the PF02899/PF12482/PF12834 phage-integrase family, and the "
        "conjugative Tra/Tcp/Tnp proteins). They are genuine transposases, but not of the elements TE-Aid "
        "curates. Including them would let bacterial contamination in a eukaryotic assembly read as a TE hit — "
        "arguably useful, but a different question, and one the curator should opt into rather than meet by surprise."
    ),
    "PD-(D/E)XK nuclease superfamily": (
        "a large fold shared by restriction enzymes, phage nucleases and host repair proteins "
        "(PF04411, PF07788, PF08011, PF11645 and ~20 more). Almost none of it is transposon-specific."
    ),
    "viral RNA-dependent RNA polymerases": (
        "the other 13 members of CL0027. They share the RdRP fold with reverse transcriptase but do not make "
        "DNA, and a hit means an RNA virus, not a retroelement."
    ),
}


def load_names(cache: Path, verified: Path) -> dict[str, str]:
    """Authoritative Pfam names, from InterPro — never from this file."""
    names: dict[str, str] = {}
    if cache.exists():
        blob = json.loads(cache.read_text())
        for group in (*blob.get("clans", {}).values(), *blob.get("search", {}).values()):
            for row in group:
                names[row["acc"]] = row["name"]
    if verified.exists():
        for acc, row in json.loads(verified.read_text()).items():
            if "error" not in row and row.get("name"):
                names[acc] = row["name"]
    return names


def provenance(acc: str, cache: Path) -> str:
    """Which clan or search first surfaced this accession."""
    if not cache.exists():
        return "manual"
    blob = json.loads(cache.read_text())
    for clan, rows in blob.get("clans", {}).items():
        if any(r["acc"] == acc for r in rows):
            return f"clan:{clan}"
    for term, rows in blob.get("search", {}).items():
        if any(r["acc"] == acc for r in rows):
            return f"search:{term}"
    return "manual"


def main() -> int:
    scratch = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    cache, verified = scratch / "pfam_cache.json", scratch / "pfam_verified.json"
    names = load_names(cache, verified)

    missing = sorted(a for a in CURATED if a not in names)
    if missing:
        print(f"ERROR: no verified InterPro name for {len(missing)}: {missing}", file=sys.stderr)
        return 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    order_rank = {o: i for i, o in enumerate(
        ["Class I", "LTR", "LTR/ERV", "LINE", "PLE", "DIRS/Crypton", "TIR",
         "TIR/Tc1-mariner", "TIR/hAT", "TIR/Mutator", "TIR/CACTA",
         "TIR/PIF-Harbinger", "TIR/P-element", "TIR/KDZ", "RC/Helitron", "RC",
         "domesticated"])}
    rows = sorted(
        CURATED.items(),
        key=lambda kv: (order_rank.get(kv[1][0], 99), kv[0]),
    )

    with OUT.open("w") as handle:
        handle.write("pfam_acc\tname\tte_order\tprovenance\twhy_included\n")
        for acc, (order, why) in rows:
            handle.write(f"{acc}\t{names[acc]}\t{order}\t{provenance(acc, cache)}\t{why}\n")

    print(f"wrote {OUT} — {len(rows)} domains across {len({o for o, _ in CURATED.values()})} orders")
    print(f"recorded {len(EXCLUDED)} explicit exclusions and "
          f"{len(EXCLUDED_CATEGORIES)} excluded categories")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
