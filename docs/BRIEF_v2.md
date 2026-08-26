# TE-Aid v2 — working brief

Self-contained brief for a fresh Claude Code session. Everything needed is
here; no prior conversation required. Canonical home of this brief:
`docs/BRIEF_v2.md` on the `v2` branch (any copy elsewhere is stale).

Repo: `https://github.com/clemgoub/TE-Aid` (maintainer: Clément Goubert).
Work happens on the `v2` branch. `main` is v1 (shell + R), published
in Mobile DNA 2022 (doi:10.1186/s13100-021-00259-7) and preserved unchanged
on the `v1-legacy` branch.

---

## 1. What TE-Aid is, and what v2 changes

TE-Aid produces a **single evidence sheet for one TE family**, to support
manual curation. v1 takes a consensus FASTA + a reference genome, runs
`blastn` (consensus vs genome), a self-`blastn`, EMBOSS `getorf`, and `blastp`
of ORF peptides vs the RepeatMasker TE-protein library, and draws four panels:
genomic hits vs divergence · consensus coverage pileup · self dot-plot ·
structure (TIR/LTR suggestions, ORFs, protein hits).

**v2 makes four changes:**

1. **Annotation-driven instead of search-driven.** v1 rediscovers copies with
   `blastn`, which is leaky — it misses diverged copies and splits elements at
   indels. v2's default input is an annotation the user already trusts
   (RepeatMasker `.out`, `.gff3`, or BED16), so copies, coordinates,
   divergence and strand are read, not re-derived. No `makeblastdb`, no
   genome-wide search: seconds per family.
2. **Python rewrite with an interactive HTML sheet** (static PDF/PNG retained
   for papers). v1's R/shell implementation is not extended; it is preserved
   on a `v1-legacy` branch.
3. **Stockholm (Dfam seed) input** — copies, coordinates, the alignment *and*
   the consensus all come from the seed (§7). Combines with `--annot`.
4. **Stronger homology evidence**: frameshift-aware translated pHMM search
   (BATH) for protein hits, and optional nucleotide hits against a taxonomic
   slice of Dfam.

**Hard scope boundary: TE-Aid never asserts a classification.** It renders
evidence — best protein hits, best nucleotide hits, structure — and leaves the
verdict to the curator. This is deliberate and load-bearing: in seed-QC mode
the sheet flags disagreement between an *externally supplied* expected class
and the evidence shown, which only means something if the sheet itself has no
opinion. Do not add a "predicted class" line.

---

## 2. Settled decisions (do not re-litigate)

| Topic | Decision |
|---|---|
| Standalone input priority | 1) `--annot` (RM `.out`/`.gff3`/BED16) · 2) `--blastn` (legacy) · 3) `--stk` |
| Pipeline mode | `--pipeline` reverses to Stk-first |
| Language / output | Python 3; interactive HTML sheet + static export |
| Plotting stack | Plotly, with `kaleido` for static export. The HTML sheet loads plotly.js from a CDN, so it needs network to *view*; static exports are self-contained |
| v1 compatibility | **Not** maintained in v2. `v1-legacy` cut from `main` (done); README points at it. |
| EMBOSS | **Kept** — `getorf` for the ORF track (§5.3). `dotmatcher` retired (self-blastn dot plot replaces it, §5.2) |
| Seed-QC panels | Flag-gated `--seed-qc`; off by default. Appended as a row *below* the 2×2 grid, so the default sheet has the **same layout** whatever the input was — only the provenance line under the title differs |
| Classification | Never asserted (§1) |
| Protein library | **Three tiers (decided 2026-08-25, see §5.1a).** Default = curated Pfam pHMMs ∪ pHMMs for only the superfamilies Pfam cannot model — 233 MB searched. Fallback = `blastp` vs full RepeatPeps, as v1. `--deep` = full RepeatPeps as single-sequence pHMMs (~6.4 GB), built on demand and cached; it **replaces** tiers 1–2 and turns off tier 3. The original "Pfam ∪ all of RepeatPeps as pHMMs" is **superseded**: it is a 6.4 GB model file |
| Nucleotide library | Dfam FamDB only requirement; user-supplied merged famdb works via the same path |
| Benchmark set | ~20 families from `GCA_963082875.1` (the dev-data goby) stratified by order, **plus** a Dfam-curated subset for recall — see §5.5 |
| Sheet layout | **v1's 2×2 quadrant grid is load-bearing and must be preserved** — the whole sheet readable at a glance, with the dot-plot at a true 1:1 aspect. This is a large part of why v1 is used. Do not replace it with a vertical stack. |
| Panel 4 / 5 | The bottom-right quadrant splits into stacked sub-panels (0.42 / 0.58) **by provenance**: panel 4 is what the *dot-plot* implies (terminal-repeat candidates only), panel 5 is what *sequence annotation* implies (ORFs with their domain hits). Revised 2026-08-25 — ORFs moved out of panel 4 |
| Package naming | Importable `teaid`; CLI entry point `teaid` **plus** a `TE-Aid` alias that prints a one-line notice pointing at the new name |
| Third-party branches (`Jiangzhao`, `simonorozcoarias`) | Left alone — independent active projects |
| RepeatAfterMe | Out of scope |

---

## 3. Input formats

One reader module, three annotation dialects, one shared record struct:
`chrom · start · end · strand · family · class_label · divergence ·
divergence_kind · consensus_start · consensus_end · consensus_left`.

- **RepeatMasker `.out`** — fixed-width (not tab-delimited; split on
  whitespace after skipping 3 header lines). Strand column is `+` or **`C`**
  (complement), not `-`. A leading `*` marks a lower-scoring overlapping hit.
  Consensus coordinates are columns `repeat_start/repeat_end/repeat_left`,
  and on `C` rows they are written in the order `left, end, start` — the
  value in parentheses is the one *past* the matching end. Normalize to
  `consensus_start <= consensus_end` at read time.
- **RepeatMasker `.gff3`** — attribute-keyed; consensus coordinates live in
  the `Target=` attribute when present.
- **BED16** — the VGP repeat-hub interchange format. Column definitions are in
  `VGP_TEbed/docs/INPUT_FORMAT.md` (repo confirmed at `~/Documents/VGP_TEbed`);
  do not guess the columns. 0-based half-open, unlike `.out`.
- **Stockholm** — Dfam seed; see §7 for the identifier convention.

**`divergence_kind` is not optional bookkeeping.** `perc_div` means different
things depending on the producer: divergence from a library consensus (an age
proxy) for homology tools, but *array homogeneity* for tandem-repeat finders,
and nothing at all for pure detectors. Carry
`consensus | array_homogeneity | none` per record and **never plot absent
divergence as zero** — leave the point out and say so in the axis label.

---

## 4. Panels

**Default sheet — v1's 2×2 quadrant grid, preserved:**

```
+---------------------------+---------------------------+
| 1  copies vs divergence   | 2  consensus coverage     |
+---------------------------+---------------------------+
|                           | 4  terminal repeats  0.42 |
| 3  self dot-plot (square) +---------------------------+
|                           | 5  ORFs + domains    0.58 |
+---------------------------+---------------------------+
        what the dot-plot says  |  what annotation says
```

1. Genomic hits vs divergence (from the annotation / seed / blastn).
2. Consensus coverage pileup — copies **spanning** each position.
3. Self dot-plot (see §5.2), at a **true 1:1 data aspect** — both axes are
   consensus base pairs, and the square is what makes an off-diagonal repeat
   read as parallel to the main diagonal.

**Three geometry rules make the grid comparable, and they interlock:**

- **Every quadrant is square.** The sheet's height is derived from its width to
  make it so, exactly as v1 used a 12×12 in page for a 2×2 grid of 6×6 in
  panels. Without it the dot-plot can be square *or* the same pixel width as the
  panel above it, never both.
- **Every panel carries the identical x-range**, so a consensus position lands
  at the same screen x in all four. The dot-plot satisfies its 1:1 aspect by
  shrinking its *domain* (`constrain="domain"` on both axes), never by widening
  its range — the default behaviour makes it silently wider than its neighbours.
- **That range is padded 2.5% past each end of the consensus.** A feature at a
  terminus is otherwise clipped, because an arrowhead is a fixed-size marker
  centred on its coordinate. Autoscaling the panel un-clips it but gives that
  panel its own range, which breaks the first two rules.
4. **Structure** — terminal-repeat candidates **only**: what the self dot-plot
   in panel 3 implies, so the two quadrants read as one statement. ORFs moved
   to panel 5 (decided 2026-08-25), where they can be drawn together with the
   protein hits that sit in them. A rail across the full consensus shows how
   much of the element the repeats bracket. Repeat pairs are drawn as **arrows**,
   carrying over v1's good idea: a repeat pair's two arms point the same way
   when direct (`→ … →`, LTR-like) and at each other when inverted
   (`→ … ←`, TIR-like), so orientation is readable from shape rather than from
   a colour key. **Inverted pairs always point inward**, decided from the arms'
   positions rather than from blastn's coordinate order — blastn reports a pair
   in whichever direction it found it, so drawing it verbatim made some inverted
   pairs point outward. Lanes are labelled `direct`/`inverted` — what was measured —
   with the LTR/TIR reading offered in the legend, because on an internally
   repetitive consensus the same signature is a tandem unit and naming the lane
   `LTR` would assert a class (§1).
5. **ORFs and protein homology** — the "Stack" layout, chosen 2026-08-25 from
   three candidates rendered at true panel size. One row per ORF, drawn as a
   rectangle outlined by strand (**black forward, red reverse, as v1 did**),
   with the protein hits that sit in that ORF's *reading frame* drawn as arrows
   on top of it. A hit with no open frame in its own register gets a bare row
   over a dotted ground rule. Domain colour is **v1's TE-class scheme** (green
   LTR, blue LINE, salmon DNA transposon) so type reads at a glance the way
   curators are used to; colour is never the only channel, since every row names
   its contents in the tick gutter. Also: best nucleotide hits (§5.1), pending.

   Two rules the layout depends on, both argued in full where the code lives:
   **a hit is housed by reading frame, not by coordinate overlap**
   (`annotation_rows.py`), and **`getorf` keeps its default `-find 0`**
   (`orfs.py`). Either one changed silently invalidates the panel — the first by
   letting it assert "intact coding domain" for a frameshifted hit, the second
   by making an ORF rectangle's edge stop meaning "where the frame closes".

   The silhouette is the deliverable: a compact block of framed rows above a run
   of bare, notched ones is a decayed element, and no ORF-finder-then-align
   search can draw it — the bare rows are precisely what such a search cannot
   see.

Panels 4 and 5 are the split of v1's single crowded structure quadrant. They
occupy that quadrant stacked, sharing one consensus x-axis so features line up
vertically; the quadrant grid itself is not disturbed. The quadrant splits as
soon as there is **either** an ORF track **or** homology evidence to draw; with
neither, panel 4 takes the whole quadrant rather than leaving a hole.

**The grid is the product, not a styling choice.** Everything visible at once,
in fixed positions, is how curators read these sheets; a vertical stack that
forces scrolling loses it. Any layout change must keep the four quadrants and
the square dot-plot.

**`--seed-qc` adds:**

6. **Seed depth**, shaped after [Dfam's own seed-alignment track](https://dfam.org/family/DF000001423/browser)
   (decided 2026-08-25): a coverage band split into sequences that **match** the
   consensus and sequences that **differ**, drawn over a pileup of the individual
   seed sequences as their aligned runs, so an internal deletion shows as a gap
   in a lane rather than being smoothed away. The **Dfam ≥3-sequences floor** is a
   horizontal rule across the band, and stretches falling short are shaded.

   Splitting by agreement is the point: a column with 43 aligned sequences of
   which 25 disagree is not 43 sequences of support, and a single depth curve
   hides that. Both series are counts of sequences, so the panel has **one**
   y-axis (`sequences · pileup`), not two.

   Two silent caps, because an unbounded pileup stops being readable:
   `MAX_PILEUP_LANES = 40` sequences drawn, and `MAX_HOMOLOGY_LANES = 14` for
   panel 5's hits. Both are noted on the sheet rather than applied invisibly.

   **Panel 6 is not panel 2 restated, and the difference is load-bearing.**
   Panel 2 counts copies that *span* a position — everything between a copy's
   first and last aligned base, internal deletions included. Panel 6 counts
   sequences that actually *contribute a base* there. A copy with a 200 bp
   internal deletion is present in panel 2 across its whole span and absent
   from panel 6 inside the deletion, so the gap between the two curves *is* the
   family's internal deletion structure.

   This is not a corner case: across the 409 GenomeArk seeds for
   `GCA_963082875.1`, **362 differ**, by as many as 37 sequences at a single
   position (`rnd-1_family-137` position 179: 44 copies span it, 21 of them
   have a deletion there, so the depth is 23).

   **Consequence for anyone gating seeds on "depth ≥ 3".** Decide which quantity
   you mean. A gate computed from consensus coordinates measures *spanning* and
   is optimistic: a family can pass it and still have fewer than three real bases
   at positions inside common deletions. Dfam's requirement concerns bases at a
   column, which is the base-level quantity panel 6 draws. (§7 records one real
   pipeline whose gate has exactly this issue.)
7. **Consensus vs contributing library entries** — the rebuilt consensus
   aligned against each source library entry: end extension means a truncation
   was fixed; mid-sequence disagreement is a chimera warning.
8. **Expected-class label** — read `#=GF TP` from the seed, display it, and
   **flag disagreement** with the evidence in panel 5 (e.g. `TP` says
   `DNA/CMC-EnSpm` but the only protein hit is a reverse transcriptase). The
   flag is the deliverable, not the label.

---

## 5. Evidence and methods

### 5.1 Homology evidence (panel 5)

Borrowed from RepeatClassifier's **evidence search** — explicitly *not* its
classification decision. If RepeatClassifier's best-hit selection rules are
reused for overlapping/competing hits, cite it as the source in the docs.

**Protein row — BATH** (`https://github.com/TravisWheelerLab/BATH`,
Wheeler Lab): frameshift-aware translated pHMM search. Use
`bathsearch --fs`. The reason this replaces `getorf` + `blastp` as the
*homology* evidence: frameshifted and pseudogenized ORFs are invisible to an
ORF-finder-then-align approach by construction, and those are precisely the
copies that make TE ORFs hard to annotate. (The ORF *track* itself stays
`getorf`-based — see §5.3.)

Query set = two parts:
- **Curated Pfam TE domain list** — see §5.4; this is a deliverable, not an
  existing file.
- **RepeatPeps** (ships with RepeatMasker, Creative Commons licensed) built
  into single-sequence pHMMs. Verified: `bathbuild` accepts *unaligned*
  sequence files and builds one pHMM per sequence —
  `bathbuild RepeatPeps.bhmm RepeatPeps.lib`. No alignment work needed.

### 5.1a The protein library: three tiers, and why

**Decided 2026-08-25 by the maintainer, after the measurements below. Do not
re-derive this; the numbers are here so the decision survives.**

§5.1's original plan — the curated Pfam set *plus all of RepeatPeps* as
single-sequence pHMMs — does not survive contact with the sizes:

- **RepeatPeps as one pHMM per protein is ~6.4 GB.** Measured, not estimated:
  `bathbuild` had written 718 MB at 2,022 of 18,011 sequences before the build
  was stopped. 16.2 M residues at ~396 bytes each.
- **Reducing redundancy first does not help.** cd-hit at 90% identity removes
  0.5% of RepeatPeps (18,011 → 17,916); at 80%, 3%. It is already a
  non-redundant curated library, so the size is intrinsic.
- **But RepeatPeps cannot simply be dropped**, because Pfam has *no model at
  all* for **piggyBac**, **Maverick/Polinton** and **Crypton**, and only a
  generic GIY-YIG for **Penelope**. RepeatPeps covers those four with 870
  proteins. The two halves are genuinely complementary, not redundant.

So the library is **three tiers**, with `--deep` as an opt-in fourth. The
measured sizes and the mechanics live in **`teaid/proteins.py`'s docstring**,
beside the code that builds them; what belongs here is only the decision:

- **Tiers 1–2 are always searched**, frameshift-aware, and cover the conserved
  catalytic domains — the curated Pfam set plus pHMMs for *only* the four
  superfamilies Pfam cannot model. 233 MB, ~10 s per family.
- **Tier 3 is `blastp` vs the whole RepeatPeps FASTA, as v1 did**, kept because
  it answers a different question — *which named element does this resemble*,
  rather than *which domain does it encode* — at negligible cost.
- **`--deep` replaces tiers 1–2** with all of RepeatPeps as single-sequence
  pHMMs (~6.4 GB), and turns tier 3 off as well. For anyone who wants
  frameshift-aware search across everything and has the disk.

**What is deliberately *not* claimed here.** This is a cost decision, not a
sensitivity one. Whether tiers 1–2 actually recover what tier 3 or `--deep`
would is exactly what the §5.5 benchmark measures, and the framework is being
built first so that benchmark has something to run. Do not present the tiering
as validated until it has been.

**Nucleotide row — `nhmmer` against a taxonomic slice of Dfam FamDB.**
- `--species <taxon>` extracts a famdb slice, cached under
  `~/.teaid/famdb_slices/<taxon>-<dfam_version>/`.
- **Documentation must be extensive here** (explicit maintainer request):
  how to find valid taxon names (`famdb.py names`), that FamDB rejects
  informal names like "mammal", what a slice contains (lineage-specific +
  ancestral families), cache invalidation when the Dfam version changes, and
  disk cost per slice.
- Without `--species`, **skip the nucleotide row** rather than searching all
  of Dfam. Speed is a feature.
- Licensing, for the record: Dfam is open and is the only required source
  (RepeatMasker ≥ 4.1.7 ships no libraries; you download FamDB partitions).
  Repbase itself remains closed, but much Repbase-derived content is already
  inside open Dfam from its early releases. If a user points at a locally
  merged "withRBRM" famdb it flows through the same code path — identical
  format, no special casing, nothing redistributed.

### 5.2 Dot-plot method, and the TSD correction

The maintainer asked what the best dot-plot method is for identifying LTRs,
TIRs, and TSDs. Three separate answers:

- **LTRs** (direct terminal repeats) and **TIRs** (inverted terminal repeats)
  *are* findable from a consensus self-comparison: LTRs appear as
  off-diagonal same-strand hits near both termini, TIRs as opposite-strand
  hits. Self-`blastn` with a reduced word size and dust disabled
  (`-word_size 7 -dust no -evalue 1e-3`) is the practical default — terminal
  repeats within a *consensus* are usually well conserved, since the consensus
  is already an average over copies. Keep this; retire EMBOSS `dotmatcher`
  (the extra output added little that self-blastn does not show).
  If sensitivity on degraded terminals proves insufficient, LAST/`lastz` is
  the upgrade path — evaluate only if a real case demands it.
- **TSDs cannot be found from a dot plot at all.** A target-site duplication
  is a short (typically 2–10 bp) direct repeat *in the flanking genomic
  sequence* on either side of an insertion — it is not part of the consensus,
  so no self-comparison can reveal it. TSD detection requires the **flanks of
  multiple genomic copies**: align the flanks of the copies and look for a
  short direct repeat immediately abutting the element boundary, consistently
  across copies. Two consequences: (a) the extraction step must keep flanks
  (±500 bp is a reasonable default), and (b) a TSD consensus computed this way
  fills Dfam's `#=GF TD` field, which any seed producer needs anyway. Implement TSD detection in the flank-analysis path, and label
  it clearly as flank-derived, not dot-plot-derived.
- **Caveat: do not assume a perfect end-to-end consensus.** Flank-derived TSD
  detection presumes the consensus termini are the element's true boundaries.
  Sometimes they are — but consensuses can be chimeric (a TE embedded in
  another TE's consensus), in which case the extracted "flanks" are the host
  element, not genomic sequence, and no consistent TSD will emerge. Report
  that absence honestly rather than forcing a call. When the boundaries *are*
  right the payoff is real: several Class II (DNA) superfamilies carry
  diagnostic TSDs (e.g. `TA` for Tc1/mariner, ~8 bp for hAT, 9–11 bp for
  Mutator), so report the detected TSD sequence and length as evidence the
  curator can weigh — never as a classification (§1).
### 5.3 ORF track

ORF finding stays EMBOSS `getorf`, exactly as in v1 (EMBOSS remains a
dependency; only `dotmatcher` is retired, §5.2). Note that BATH's translated
alignments already mark coding regions *including frameshifted ones*, so the
ORF track and the protein row are complementary: ORFs show open frames, BATH
shows homology regardless of frame integrity. Draw both — for now. The
benchmark (§5.5) should additionally measure whether the ORF track is
redundant once BATH's coding annotation is drawn; if it never shows anything
BATH does not, retiring `getorf` can be revisited then, with data.

### 5.4 `teaid/data/te_domains.tsv` — a deliverable, and a trap to avoid

The protein query set needs a TE-domain list. **Do not select by Pfam clan.**
Verified clan memberships show clan selection failing in both directions at
once:

| domain | Pfam | clan | clan size |
|---|---|---|---|
| Reverse transcriptase (RVT_1) | PF00078 | CL0027 RdRP | 15 |
| Integrase core | PF00665 | **CL0219 RNase_H** | **121** |
| RNase H | PF00075 | CL0219 RNase_H | 121 |
| DDE endonuclease | PF03184 / PF13358 | CL0219 RNase_H | 121 |
| Mutator transposase | PF00872 | CL0219 RNase_H | 121 |
| Copia gag | PF14223 | CL0523 GAG-polyprotein | 16 |
| Helitron helicase-like | PF14214 | CL0169 Rep | — |
| hAT C-term dimerisation | PF05699 | **no clan** | — |
| DUF4216 (Helitron-assoc.) | PF13952 | **no clan** | — |

- **Over-inclusion:** CL0219 has 121 members and mixes integrase, RNase H and
  the DDE nucleases with *host* enzymes — RNase H is a cellular enzyme.
  Taking the clan manufactures false TE evidence.
- **Under-inclusion:** hAT dimerisation and DUF4216 are TE-diagnostic and have
  **no clan at all**; clan selection misses them entirely.

So build an **accession-level curated table**: `pfam_acc · name · te_order ·
provenance · why_included`. Seed it from the clans above minus host-enzyme
members, add the clanless TE-specific domains, and cross-check against what
RepeatPeps already covers. **Draft it, then send the table to the maintainer
for review** — he specifically wants eyes on the CL0219 host-enzyme
exclusions. Keep it versioned in the repo; it is small and reviewable, and it
is what makes the protein panel defensible.

**The shipped TSV is generated — do not hand-edit it.** The chain lives in
`tools/`, run in order:

| Script | Does |
|---|---|
| `fetch_pfam.py` | pulls candidate domains and clan memberships from InterPro |
| `verify_pfam.py` | checks each accession resolves and reports clan sizes |
| `build_te_domains.py` | applies the inclusion/exclusion rules → `teaid/data/te_domains.tsv` |
| `build_review_page.py` | renders the maintainer's review artifact |

**The exclusions awaiting review live in `build_te_domains.py`**, as
`EXCLUDED` (29 accessions, each with a reason) and `EXCLUDED_CATEGORIES` (3) —
not in the shipped TSV, which records only what survived. Review the script,
not just the table.

`te_order` is a **closed vocabulary with two hard-coded consumers that must be
edited together**: `classcheck._ORDER_CLASS` (Class I/II — drives the `#=GF TP`
disagreement flag) and `tetypes._ORDER_TO_CLASS` (v1's colours). A value absent
from either fails *silently*: it simply never votes, or renders Unknown grey.

The 17 values currently in use — `Class I`, `DIRS/Crypton`, `LINE`, `LTR`,
`LTR/ERV`, `PLE`, `RC`, `RC/Helitron`, `TIR`, `TIR/CACTA`, `TIR/hAT`,
`TIR/KDZ`, `TIR/Mutator`, `TIR/P-element`, `TIR/PIF-Harbinger`,
`TIR/Tc1-mariner`, `domesticated` — are covered, with **three deliberate gaps**:

| Value | Class vote | Colour | Why |
|---|---|---|---|
| `Class I` | I | Unknown grey | implies a class but no single order — an RT is shared by every Class I order |
| `DIRS/Crypton` | *none* | DIRS | the class assignment is contested; the colour is safe, the vote is not |
| `domesticated` | *none* | Unknown grey | a host-domesticated domain is evidence of neither |

Adding a value means adding it to both dicts **and** deciding which of these
three shapes it takes. Do not let it default.

Editing the TSV changes `proteins._signature()` and therefore **invalidates the
whole cached protein library**, forcing a ~5-minute rebuild on the next run.

### 5.5 Benchmark: cost and relative sensitivity

There is no ground-truth TE annotation for the test genome, so the benchmark
measures cost and *relative* sensitivity, not accuracy:

- **Cost:** wall-time per family per arm → decides which arm is the default
  and which becomes `--deep`.
- **Relative sensitivity:** domains recovered *only* by BATH with `--fs`
  enabled; hand-inspect a handful to confirm they are real degraded ORFs.
- **Partial recall check that is not genome-specific:** Dfam **curated**
  families carry known classifications — run the arms over a curated subset
  for the only honest accuracy signal available.
- **ORF-track redundancy (§5.3):** compare the `getorf` ORF track against
  BATH's coding regions across the benchmark families — quantify how often
  the ORF track shows something BATH does not. This decides whether `getorf`
  stays long-term.

Arms: v1 (`getorf` + `blastp` vs RepeatPeps) · `hmmscan` vs the curated Pfam
set · **BATH `--fs`** vs (Pfam set + RepeatPeps pHMMs) · `nhmmer` vs a Dfam
slice · **`--deep`** vs all of RepeatPeps as pHMMs.

The `--deep` arm is what makes §5.1a's claim testable — without it the benchmark
cannot say whether tiers 1–2 recover what `--deep` would, which is the specific
thing §5.1a refuses to assert. Budget ~6.4 GB and a long one-off build for it;
run it last, and if it is dropped, record here that §5.1a's claim stays open.

**What is still undecided, and must be fixed before the first arm runs.** These
are the things that make results comparable; deciding them mid-run invalidates
everything already measured.

- **Family set.** ~20 families from `GCA_963082875.1` (the goby already in
  `dev-data/`), stratified by the `.out` class label so each of LTR / LINE /
  SINE / TIR / RC / Unknown is represented. Write the chosen list into this
  section — a set re-derived per run is not a benchmark.
- **Dfam curated subset.** Name the release, the partition, the family count and
  the local path. "Dfam curated families" is not yet a reproducible input.
- **One E-value across all arms.** The defaults differ — `homology.search` uses
  `1e-3`, `search_repeatpeps` uses `1e-5` — so an arm comparison at defaults
  measures the thresholds, not the methods.
- **What counts as "the same domain" found by two arms.** Reuse
  `best_per_region`'s 0.5 reciprocal-overlap rule rather than inventing a second
  one.
- **A pre-registered pass/fail rule, written before the numbers exist**, plus
  where the verdict is recorded (here, in §5.1a, or both).

**Two invocation traps.** `--proteins FILE` constructs `Library(repeatpeps=None)`
and so **silently disables tier 3** — an arm meant to include it must not use
that flag. And no CLI flag exposes `homology.search(frameshift_aware=…)` or
`cpus`, so the `--fs` on/off comparison and any threading have to drive the
Python API directly, not the CLI.

**One cost is already measured:** ~10 s per family for `bathsearch` against the
233 MB tier-1+2 library, dominated by loading the models rather than by the
search. BATH has no `hmmpress` equivalent, so there is no index step to amortise
it — which is why a whole-library run wants one library load reused across
families rather than more parallelism alone (§6).

---

## 6. Work order and status

Detail that a reader needs *while changing the code* lives in the module
docstrings, not here — `analysis.py`, `homology.py`, `proteins.py`,
`annotation_rows.py`, `report.py` and `readers/*` each explain their own traps.
This section is status and direction only.

| # | Step | State |
|---|---|---|
| 1 | `v1-legacy` cut, `v2` opened, README pointer | done |
| 2 | Package, CLI, annotation readers, panels 1–3 | done — v1 parity without blastn |
| 3 | Interactive HTML + PNG/PDF/SVG export, both themes | done |
| 4 | Panels 4–5 split; `getorf` ORF track; self-blastn dot-plot | done |
| 5 | `--stk`, `--pipeline`, `--seed-qc`, panels 6 and 8, fail-soft contract | done bar panel 7 |
| 6 | `te_domains.tsv`, BATH, three-tier library, panel 5, `TP` flag | done, table **awaiting review** |
| 7 | Benchmark (§5.5) → set the default protein path | next |
| 8 | Nucleotide row: famdb slice + `nhmmer` + cache + docs | |
| 9 | TSD detection in the flank path (§5.2) + `#=GF TD` | |
| 10 | `--blastn` legacy path, with the leakiness warning | |
| 11 | Tag `v2.0` | |

### Open items, smallest first

- **Cosmetic, panel 5.** Tick rosters still crowd where two long rows abut, and
  a row's roster lists names without mapping them to individual bars within the
  row.
- **Cosmetic, one red means two things.** `report.STATUS_BELOW_FLOOR` (the Dfam
  floor warning, panel 6) is the same `#d03b3b` as `theme.orf_reverse` (v1's
  reverse-strand ORF outline, panel 5). On a `--seed-qc` sheet both are visible.
  They never share a panel and both are labelled, so nothing is unreadable — but
  a reserved status colour that is also a routine data colour is a smell. Either
  move the status red, or accept it and delete the "reserved" framing.
- **Judgement call, panel 5 colour.** `RVT_1` and other order-agnostic domains
  render Unknown grey, because a bare reverse transcriptase is shared by every
  Class I order and colouring it LTR would assert what the hit cannot support.
  It does make the most important domains look uninformative. The alternative
  is to let such a domain inherit the colour when every element-level hit on the
  sheet agrees. Undecided.
- **Panel 1 does not yet use `Annotation.fragment_groups()`.** BED16 column 16
  (and the `.out` `ID` column) groups fragments of one interrupted insertion;
  counting them separately inflates copy number, which is the specific failure
  that motivates v2. The grouping is implemented and tested; panel 1 still
  counts rows.
- **Panel 7** (consensus vs contributing library entries) needs the source
  entries, which only a seed *producer* has. Design the hook generically — a
  FASTA of candidate source entries — rather than around any one pipeline.
- **`te_domains.tsv` review** is with the maintainer: 130 domains, 29 explicit
  exclusions, 3 excluded categories.

### Direction after step 7 (raised 2026-08-25, not yet scoped)

- **Whole-library runs and parallelism.** Today `teaid` is one family per
  invocation. A 400-family library at ~10 s of `bathsearch` each is over an
  hour serially, and the model load dominates that (§5.5). Wants: a batch mode,
  a worker pool, and probably one library load reused across families.

  **Everything a batch mode must hoist out of the per-family path** — all of it
  is currently repeated for every single family:

  | Repeated work | Where |
  |---|---|
  | full re-parse of the `.fa.out` (147k records) | `cli.py` |
  | `stockholm.read()` reads *all* 409 records to return one | `readers/stockholm.py` |
  | `makeblastdb` over RepeatPeps into a fresh temp dir | `homology.py` |
  | cold protein-library build, unlocked and racy | `proteins.build()` |

  Two constraints on the obvious shortcuts. **One multi-family `bathsearch`
  cannot be split back apart**: `ProteinHit` drops tblout's `target` column, so
  there is nothing to group hits by — adding that field is a prerequisite. And
  **`proteins.build()` takes no lock**, so a pool must warm the cache with one
  serial run before fanning out. `homology.search(cpus=)` exists but no flag
  reaches it.
- **CLI update.** The flag surface has grown organically across steps 2–6 and
  deserves a pass: grouping, defaults, and a `--version`-style summary of which
  external tools were found.
- **GUI.** Unscoped. The sheet is already a standalone HTML page, so the
  smallest useful version may be an index over many sheets rather than a new
  application. One correction to that idea: the sheet is **not** self-contained
  — `report.py` writes `include_plotlyjs="cdn"`, so it renders blank without
  network. A GUI meant for offline or air-gapped use has to switch that to an
  inlined bundle (~3 MB per sheet) or serve one shared copy.

## 7. Scope, input routes, and downstream integrators

### TE-Aid v2 is a general-purpose tool, and this branch keeps it that way

**Nothing in this repo is specific to any one pipeline or consortium.** TE-Aid
inspects one TE family from whatever evidence the user has, and the routes are
designed to compose. A downstream project that needs project-specific behaviour
**forks TE-Aid and modifies its fork** — it does not push its architecture back
into this branch. Keeping the tool general is what makes it worth forking.

Practically, when a change is proposed, ask: *would this make sense to someone
curating TEs who has never heard of the pipeline that asked for it?* If not, it
belongs in the fork.

### The four input routes

**The routes and their contract are specified in `docs/INTEGRATION.md` §2** —
that file is the one callers read, so it is the one that stays right. What
belongs here is only the design intent behind them:

Routes `--annot` (a), `--blastn` (b, unimplemented), `--stk` (c) and both
together (d) are **designed to compose, not to compete**. Route **d** is the
point: the seed says which copies were *chosen*, the annotation says which
*exist*, and the ratio between them ("the seed used 44 of 515 copies") is a QC
signal neither input produces alone. It is available whenever a RepeatModeler2
run supplies both its `.out` and its `.stk`, which is the common case.

Two seed use cases, both first-class, and they differ in what is *trustworthy*
rather than in what is drawn:

1. **Raw seed evaluation** — a `.stk` straight from a RepeatModeler2 run, judged
   before anyone invests curation effort in it. The sampling is RepeatModeler's.
2. **QC of a purpose-built seed** — a seed some pipeline assembled from a chosen
   copy set; the question is whether that choice was sound.

### What an integrator can rely on

Anything below is a contract: it will not change without a version bump. Full
detail, including worked invocations, is in **`docs/INTEGRATION.md`** — read
that file first if you are building a tool that calls or forks TE-Aid.

- **Invocation** is per family. `--pipeline` reverses the standalone input
  priority to seed-first, for callers whose primary artefact is a seed.
- **Fail-soft**: failures TE-Aid diagnoses carry a distinct exit code *and* a
  stable stderr slug. Argparse-level errors and uncaught exceptions do not —
  `INTEGRATION.md` §3 has the full table and the exit-0-is-not-complete caveat.
  TE-Aid does not validate seeds beyond what it needs to draw them.
- **Stockholm identifiers** are Smitten, in either the 4-part or 2-part shape,
  1-based fully closed on the way in. Specified in `INTEGRATION.md` §4.
- **Dfam seed conventions**: `.` as the gap character, `#=GF` fields (`DE`, `AU`,
  `TP`, `OC`, `SQ`) and a `#=GC RF` consensus line. Spec: `Dfam_Seeds.md` in
  `https://github.com/Dfam-consortium/dfam-curator`, which also ships `stk lint`
  — the natural machine acceptance gate for a seed *producer*. TE-Aid is the
  human-facing half of that gate, not a replacement for it.
- **`#=GF TP`** is Dfam's full semicolon-separated path
  (`Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;…`), not
  a short `LTR/Gypsy` label.

### Known downstream consumer (context only, not a requirement)

One planned consumer is a VGP "low hanging fruits" pipeline that screens
candidate families from a multi-tool repeat-annotation track hub, rebuilds each
consensus, and emits Stockholm seeds for a curator to approve or reject, with
TE-Aid as the curator-facing QC step. **It will fork this repo.** Recorded here
only so its findings are not lost, and explicitly *not* as design constraints:

- Its candidate gate **G2** ("MSA depth ≥3 over ≥99% of the consensus") is
  computed from consensus coordinates, so it measures *spanning* and is
  therefore optimistic — the trap set out under §4 panel 6. Recorded here as a
  real instance of it, not as a second explanation.
- Its deposition unit is a *cross-tool cluster*, not one program's family, so
  "contributing library entries" for panel 7 means every cluster member's
  sequence. That shape is a property of that pipeline, and panel 7's input hook
  should be designed generically (a FASTA of candidate source entries) rather
  than around it.

---

## 7b. Resuming work on this machine

Not in git, so not obvious from a clean checkout:

- **Virtualenv `.venv-teaid/`** at the repo root, with an editable install and
  `kaleido`, `pytest`, `pillow`. Use `./.venv-teaid/bin/teaid` and
  `./.venv-teaid/bin/python -m pytest -q`. Recreate with
  `python3 -m venv .venv-teaid && ./.venv-teaid/bin/pip install -e '.[static,dev]'`.
- **`dev-data/`** (gitignored) holds the development triple for GenomeArk
  assembly `GCA_963082875.1` — `.fa.out` (19 MB), `-families.fa`, and
  `-families.stk` (409 records) — plus a BED16 conversion made with VGP_TEbed's
  `scripts/rmout2bed.py`. Re-fetch from
  `https://genomeark.s3.amazonaws.com/downstream_analyses/repeats/systematic_annotations/RepeatModeler-v2.0.8/`
  under `RepeatMasker/`, `fasta/` and `stk/`. 482 assemblies are available, each
  with all three files; this one is the smallest complete triple.
- **`getorf` is MacPorts, at `/opt/local/bin/getorf`**, and that directory *is*
  on this machine's PATH already — no export needed. If it ever is missing,
  `orfs.py` says so explicitly rather than silently dropping the track.
- **BATH 2.0** is built at `~/Documents/BATH/opt/bin/` (`bathsearch`,
  `bathbuild`, `bathconvert`, `bathfetch`, `bathstat`) and is **not** on the
  default PATH — `export PATH=$PATH:$HOME/Documents/BATH/opt/bin` before any run
  that should draw the protein row, or it degrades to a warning and exit 0.
  Two traps if it is ever rebuilt: `easel` must be cloned separately and checked
  out on its own `BATH` branch, and `--prefix` must not be `…/BATH/install`,
  which collides with the repo's own `INSTALL` file on a case-insensitive
  filesystem.
- **The protein library** caches under `$TEAID_CACHE` (or `~/.teaid/proteins`).
  A prebuilt, verified-complete tier-1+2 library sits in
  `dev-data/protein-cache/` (gitignored; 233 MB searched, 473 MB on disk because
  it also keeps the rebuildable tier parts). **`export TEAID_CACHE=dev-data/protein-cache`**
  to avoid a rebuild, which re-fetches 130 Pfam models from InterPro and takes
  ~5 minutes. The test suite isolates itself (`tests/conftest.py`), so plain
  `pytest` never touches either cache.
- **`RepeatPeps.lib`** is auto-discovered. Discovery returns the **first** of
  `~/RepeatPeps.lib`, `~/Downloads/RepeatMasker/Libraries/RepeatPeps.lib`,
  `/usr/local/…`, `/opt/…`, `/opt/homebrew/…` — on this machine both of the
  first two exist and `~/RepeatPeps.lib` wins. It ships inside RepeatMasker and
  is no longer in that project's git. `proteins._signature()` keys on
  `name:size`, so any same-named copy of the same file keeps a prebuilt cache
  valid.
- **Verifying the HTML in a browser**: `file://` URLs are blocked by the Chrome
  tooling, so serve the directory (`python3 -m http.server`) and open
  `http://localhost:…`. Worth doing for anything interactive — two bugs got
  through every static check and were caught only this way (a reset button that
  silently did nothing, and hover tooltips rendering a page wide).

Useful families in the dev data, for eyeballing a change:

| Family | Why |
|---|---|
| `ltr-1_family-65` | both a direct and an inverted terminal repeat |
| `ltr-1_family-26` | LTR pair pushed inward by extra 5' sequence; 4 ORFs |
| `ltr-1_family-22` | tandem gag + pol ORFs (856 aa, 1,155 aa) |
| `rnd-1_family-137` | 515 genomic copies vs 44 in the seed; heavy internal deletion |
| `rnd-1_family-30` | internally repetitive; exercises the repeat-lane cap |
| `ltr-1_family-11` | 2 sequences — below the Dfam floor |
| `rnd-1_family-117` | 10,851 copies — the performance case |

---

## 8. Ask the maintainer before assuming

Everything else raised here has been answered and folded into §7b (tool
inventory and paths) and §5.4 (test data). One item is still open:

- **`te_domains.tsv` review** — 130 domains, plus the 29 explicit exclusions and
  3 excluded categories that live in `tools/build_te_domains.py` rather than in
  the shipped table (§5.4). The maintainer specifically wants eyes on the CL0219
  host-enzyme exclusions. This is the same item as the last bullet of §6's open
  list; it is recorded twice on purpose, because it is the only thing blocking
  on a human rather than on code.
