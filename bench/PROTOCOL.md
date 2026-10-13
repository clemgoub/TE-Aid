# TE-Aid protein benchmark: protocol

Status: **signed** by Clement Goubert, 2026-10-13. Nothing in `bench/results/`
may exist before this file is committed on its own, and `tools/benchmark.py
verify` checks that the commit that signs it predates every result.

We run the benchmark to answer two questions, in this order:

1. **Specificity.** Are tier 1–2 hits too permissive? A curator reported short
   hits from the wrong order that survive wherever no stronger hit overlaps them.
   On the goby example sheets, 91 hits from Polinton integrase models sit on three
   LTR families, spanning 20–59% of their model, at E-values from 3e-14 to 5e-6.
   We use the answer to set `--min-model-coverage` and to decide whether the
   E-value default moves.
2. **Sensitivity and cost.** Do tiers 1–2 find what the full `--deep` library finds,
   at a cost worth the tiering? We use the answer to choose the default protein
   path (§5.1a).

## Inputs, pinned

| Input | Pin |
| --- | --- |
| Curated families, *D. melanogaster* | `dev-data/dfam/dmel.*`: Dfam 4.0 (May 2026), 390 families, FASTA sha256 `c9cd7af0…86dc` |
| Curated families, human | `dev-data/dfam/human.*`: Dfam 4.0 (May 2026), 1,403 families, sha256 in its provenance file |
| Curated Maverick families, for A3 | `dev-data/dfam/maverick.*`: Dfam 4.0 (May 2026), 15 families from every clade (`--type DNA --subtype Maverick`), FASTA sha256 `787a9588…86fc` |
| Uncurated set, cost only | the 50 goby families of `dev-data/examples/selection.tsv`, committed as `bench/families.tsv` with the inputs' sha256 |
| Protein library, tiers 1–2 | cache signature `07532cf00f976128`, stamp `pfam:130,repeatpeps-gaps`, `pfam.hmm.count` 130/130 |
| `--deep` library | built once from the same `RepeatPeps.lib`; its signature is recorded when built |
| Tools | versions of `bathsearch`, `hmmsearch`, `blastp`, `getorf` recorded in each run's `manifest.json` |

The two curated sets carry each family's classification, the only ground truth
we have. We use it to judge **agreement**, not correctness: a hit agrees when its
domain's order is compatible with the family's curated class.

## Arms

| Arm | Search | Question it answers |
| --- | --- | --- |
| `default` | `bathsearch --fs`, tiers 1–2 | the incumbent |
| `nofs` | `bathsearch` without `--fs`, tiers 1–2 | what frameshift-awareness adds, with nothing else changed |
| `pfam-orf` | `hmmsearch` of the 130 tier-1 models on `getorf` peptides | whether an ORF-first search misses what BATH finds, or only fragments it |
| `tier3` | `blastp` of `getorf` peptides against RepeatPeps | v1's arm |
| `deep` | `bathsearch --fs`, all of RepeatPeps as single-sequence pHMMs | the superset tiers 1–2 are measured against |

Every arm runs at E ≤ 1 and is filtered afterwards, so one run yields every
threshold. The primary threshold is the shipped 1e-3. `--cpus 2` throughout.

## Rules for question 1 (specificity)

- A hit **contradicts** its family when its domain's order and the family's curated
  class disagree at order level (for example a Maverick integrase on an LTR/Gypsy
  family), and **agrees** when they match. Hits from order-agnostic domains (a bare
  reverse transcriptase, te_order `Class I`) are judged at class level only.
  Families of unknown class, satellites and RNA genes are left out.
- For each E-value threshold (1e-3, 1e-5, 1e-10) and model-coverage floor (0 to
  0.9 in steps of 0.1), record the contradiction rate among hits that would be
  drawn, and the share of agreeing hits kept.
- **Decision rule:** the shipped floor is the smallest one at which the
  contradiction rate among drawn hits is at most 2.5% on both species, while at
  least 94% of agreeing drawn hits are kept. If no floor meets both, the floor
  stays at 0 and the E-value default is revisited instead.

## Rules for question 2 (sensitivity and cost)

Two hits "find the same domain" when they overlap by at least half of the shorter
one, the rule `homology.best_per_region` already uses. Every match also records
the fraction of the `deep` region it spans, so a domain found in fragments is not
scored the same as one found whole.

| Criterion | Passes when |
| --- | --- |
| A1 recall | the lower 95% (Wilson) bound on the share of `deep` regions that `default` also finds is ≥ 0.85 |
| A2 span | over those regions, the median span is ≥ 0.80 and at least 60% span ≥ 0.90 |
| A3 Pfam-blind superfamilies | `default` misses none of the `deep` regions on piggyBac, Maverick, Crypton or Penelope families |
| A4 cost | `deep` costs at least 2× `default` per family; otherwise the tiering is not worth keeping |

A3 has few families to work with. The two curated sets hold 7 piggyBac (1 fly,
6 human), 9 Crypton and 1 Penelope family, and no Maverick. The 15 curated
Maverick families come from a separate fetch across every clade, run through the
same arms as a third set.

Verdict: tiers 1–2 become the validated default only if A1, A2, A3 and A4 all pass.
`nofs` against `default`, and `pfam-orf` against `default`, are reported beside it:
they measure what frameshift-awareness adds, and neither gates the verdict.

## Settled at signing

Clement Goubert chose the thresholds with a companion page showing each rule
drawn and the Question 1 rule tried on a pilot of the 50 goby example families
(uncurated labels, one species, not a benchmark set):
https://claude.ai/artifact/W5N3BfdYHxUB7YzmuizNVx. These points, which the rules
above do not state, were settled at signing:

1. "Agreeing drawn hits kept" is measured against floor 0 at the same E-value,
   as on the companion page.
2. The contradiction rate pools hits across families.
3. The ceiling applies to the raw rate, with no minimum count of drawn hits.
4. A3 applies to Penelope although the curated sets hold one Penelope family.
5. The negative control uses five shuffled sequences, as below. On the goby
   pilot, 6 of 50 shuffled sequences carried a hit at E ≤ 1e-3 and none at 1e-5.

`tools/benchmark.py` implements `run`, `shuffle`, `compare`, `report` and
`verify`. It refuses to write under `bench/results/` while this file is
uncommitted, modified or unsigned.

## Checks before any number is trusted

- Each arm matched against itself must score recall 1.000 and span 1.000.
- Five dinucleotide-shuffled consensus sequences, seed recorded, must give no hits
  at 1e-3 in any arm.
- The `default` arm's hits must equal the CLI's own for three families.
- Two runs of the same arm must give identical hits.

## Where the verdict goes

`bench/results/<run>/verdict.md` holds every criterion, its threshold and its
observed value. §5.1a of the brief gets one paragraph with the decision; the
shipped defaults change in their own commit, citing the run.

## Signature

```
Signed: Clement Goubert, 2026-10-13

- Question 1 decision rule: the shipped floor is the smallest one at which the
  contradiction rate among drawn hits is at most 2.5% on both species, while
  at least 94% of agreeing drawn hits are kept.
- A1 recall: the lower 95% (Wilson) bound on the share of deep regions that
  default also finds is >= 0.85.
- A2 span: the median span is >= 0.80 and at least 60% of regions
  span >= 0.90.
- A3 Pfam-blind superfamilies: default misses none of the deep regions on
  piggyBac, Maverick, Crypton or Penelope families.
- A4 cost: deep costs at least 2x default per family.
```
