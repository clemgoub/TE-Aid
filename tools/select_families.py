#!/usr/bin/env python3
"""Pick a stratified sample of families from a RepeatModeler2 library.

Two uses, one selector:

- **Example sheets.** A pool spanning every class label, the full range of
  consensus length and the full range of copy number, for eyeballing the sheet
  design against real data rather than against one convenient family.
- **The §5.5 benchmark set.** That section requires a *recorded* family list —
  "a set re-derived per run is not a benchmark". This script is deterministic:
  no randomness, no dict-ordering dependence, so the same library always yields
  the same list, and the list can be regenerated from the two input files rather
  than pasted around.

Selection runs in three passes, in priority order, so the guaranteed coverage is
taken before the filler:

1. families named in ``KNOWN`` — each exercises a behaviour already documented
   in ``docs/BRIEF_v2.md`` §7b, so a sample without them misses known edge cases;
2. one representative of every distinct RepeatMasker class label, so no type the
   library contains goes unrepresented, however rare;
3. fill to ``--count`` by walking size bands within each broad TE class,
   preferring high-copy families first and stepping deeper on each pass.

Families with **no annotated copy** are excluded: the ``--annot`` route exits 5
(``no-evidence``) for them, so they would be 50 failures rather than 50 examples.

``--count`` is a **floor, not a cap.** Passes 1 and 2 are guaranteed coverage and
run to completion regardless — asking for 20 families out of a library holding 33
distinct class labels cannot both honour the count and keep every type
represented, and coverage is the property worth keeping. On the dev goby the two
coverage passes alone yield 40, so ``--count`` only bites above that. The printed
row count is the truth; do not assume it equals ``--count``.

Usage::

    tools/select_families.py --annot X.fa.out --consensus X-families.fa [--count 50]
"""
from __future__ import annotations

import argparse
from collections import defaultdict

from teaid import readers, tetypes
from teaid.sequences import read_fasta

# Behaviours worth keeping in any sample. Documented in docs/BRIEF_v2.md §7b;
# names are library-specific, so a different assembly simply skips them.
KNOWN = {
    "ltr-1_family-65": "direct + inverted terminal repeat",
    "ltr-1_family-26": "LTR pair pushed inward; 4 ORFs",
    "ltr-1_family-22": "tandem gag + pol ORFs",
    "rnd-1_family-137": "515 genomic copies vs 44 in seed; heavy deletion",
    "rnd-1_family-30": "internally repetitive; repeat-lane cap",
    "ltr-1_family-11": "2 seed sequences - below the Dfam floor",
    "rnd-1_family-117": "10,851 copies - performance case",
    "ltr-1_family-51": "RT/integrase/rve all frameshifted, no covering ORF",
}

BANDS = [(0, 300), (300, 1000), (1000, 3000), (3000, 7000), (7000, 10**9)]


def select(annot: str, consensus: str, count: int) -> list[dict]:
    lib = read_fasta(consensus)
    ann = readers.read(annot)

    copies: dict[str, int] = defaultdict(int)
    for c in ann.copies:
        copies[c.family.split("#", 1)[0]] += 1

    fam = {}
    for entry in lib:
        label = entry.class_label or "NONE"
        fam[entry.bare_name] = dict(
            name=entry.bare_name,
            label=label,
            klass=tetypes.class_of_repeatmasker_label(label),
            length=len(entry.sequence),
            copies=copies.get(entry.bare_name, 0),
        )

    eligible = {k: v for k, v in fam.items() if v["copies"] > 0}
    chosen: list[str] = []
    why: dict[str, str] = {}

    def take(name: str, reason: str) -> None:
        if name in eligible and name not in why:
            chosen.append(name)
            why[name] = reason

    for name, reason in KNOWN.items():
        take(name, reason)

    by_label: dict[str, list[dict]] = defaultdict(list)
    for v in eligible.values():
        by_label[v["label"]].append(v)
    for label in sorted(by_label):
        # Highest copy number is the most typical example of its label.
        best = max(by_label[label], key=lambda v: (v["copies"], v["name"]))
        take(best["name"], f"only/typical {label}")

    by_class: dict[str, list[dict]] = defaultdict(list)
    for v in eligible.values():
        by_class[v["klass"]].append(v)

    depth = 0
    while len(chosen) < count and depth < 12:
        for klass in sorted(by_class):
            if len(chosen) >= count:
                break
            for lo, hi in BANDS:
                band = [v for v in by_class[klass]
                        if lo <= v["length"] < hi and v["name"] not in why]
                if not band:
                    continue
                band.sort(key=lambda v: (-v["copies"], v["name"]))
                pick = band[min(depth, len(band) - 1)]
                edge = hi if hi < 10**8 else ""
                take(pick["name"], f"{klass}, {lo}-{edge}bp band")
                if len(chosen) >= count:
                    break
        depth += 1

    chosen.sort(key=lambda n: (fam[n]["klass"], fam[n]["length"], n))
    return [dict(fam[n], why=why[n]) for n in chosen]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--annot", required=True, help="RepeatMasker .out")
    parser.add_argument("--consensus", required=True, help="family FASTA")
    parser.add_argument(
        "--count", type=int, default=50,
        help="minimum families to emit; the class-label coverage pass may "
             "exceed it (40 on the dev goby) and is never trimmed to fit",
    )
    args = parser.parse_args()

    for row in select(args.annot, args.consensus, args.count):
        print("\t".join(str(row[k]) for k in
                        ("name", "label", "klass", "length", "copies", "why")))


if __name__ == "__main__":
    main()
