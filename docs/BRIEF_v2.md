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
3. **Stockholm (Dfam seed) input**, for use with a seed-building pipeline
   (see §7): copies, coordinates and the alignment all come from the seed.
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
| Plotting stack | Whatever looks best — Plotly is the recommended default (polished out of the box, `kaleido` for static export). Prototype panel 1 in Plotly before committing. |
| v1 compatibility | **Not** maintained in v2. `v1-legacy` cut from `main` (done); README points at it. |
| EMBOSS | **Kept** — `getorf` for the ORF track (§5.3). `dotmatcher` retired (self-blastn dot plot replaces it, §5.2) |
| Seed-QC panels | Flag-gated `--seed-qc`; off by default |
| Classification | Never asserted (§1) |
| Protein library | Curated Pfam accession list ∪ RepeatPeps as single-sequence pHMMs |
| Nucleotide library | Dfam FamDB only requirement; user-supplied merged famdb works via the same path |
| Benchmark set | ~20 goby families stratified by order **plus** a Dfam-curated subset for recall |
| Panel 4 | Split into stacked sub-panels: structure ∥ homology evidence |
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

**Default sheet:**

1. Genomic hits vs divergence (from the annotation / seed / blastn).
2. Consensus coverage pileup.
3. Self dot-plot (see §5.2).
4. **Structure** — TIR/LTR suggestions, ORFs.
5. **Homology evidence** — best protein hits, best nucleotide hits (§5.1).

Panels 4 and 5 are the split of v1's single crowded structure panel, sharing
one consensus-coordinate x-axis so features line up vertically.

**`--seed-qc` adds:**

6. **Seed depth** — per-consensus-position alignment depth from the Stockholm
   seed, with the **Dfam ≥3-sequences floor drawn as a horizontal rule**.
   Curators need to see instantly whether the seed meets the requirement and
   where it is thin.
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
  (the companion pipeline extracts ±500 bp), and (b) a TSD consensus computed
  this way fills Dfam's `#=GF TD` field, which the seed-building pipeline
  needs anyway. Implement TSD detection in the flank-analysis path, and label
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
slice.

---

## 6. Work order

1. ~~Cut `v1-legacy` from current `main`; README pointer. Open `v2`.~~ **Done.**
2. ~~Python package skeleton + CLI (`teaid`, with the `TE-Aid` alias notice) +
   annotation reader (`.out` / `.gff3` / BED16) + panels 1–3 from `--annot`.~~
   **Done — milestone reached: v1 parity with no blastn.** Notes:
   - Format is detected by *content*, not extension (the hub ships `.bed` files
     converted from `.out`).
   - The `.out` and BED16 readers were cross-validated: the same GenomeArk
     annotation routed through VGP_TEbed's independent `rmout2bed.py` converter
     yields byte-identical records from both readers (147,231 copies).
   - Coverage uses a difference array, O(n + L); v1's dense
     `n_copies × consensus_length` matrix reached ~2 GB on large families.
   - Full-length is measured on half-open consensus coordinates, removing v1's
     off-by-one (`abs(qend - qstart)` was one short of the true span).
   - BED16 column 16 (`hit_id`) groups fragments of one interrupted insertion;
     `Annotation.fragment_groups()` exposes this. It is the direct fix for the
     fragment inflation that motivates v2, and panel 1 should use it once
     panels 4–5 land.
3. Interactive HTML sheet + static export.
4. Panels 4–5 split; `getorf` ORF track; self-blastn dot-plot per §5.2.
5. `--stk` reader + `--seed-qc` panels 6–8 + `TP` mismatch flag.
   **This is the join point with the seed-building pipeline (§7).**
6. `te_domains.tsv` draft → maintainer review; RepeatPeps → pHMMs via
   `bathbuild`; BATH protein row.
7. Benchmark (§5.5) → set the default protein path.
8. Nucleotide row: famdb slice + `nhmmer` + cache + the documentation of §5.1.
9. TSD detection in the flank path (§5.2) + `#=GF TD` output.
10. `--blastn` legacy path (port, or shell out to v1) with the leakiness
    warning documented.
11. Tag `v2.0`.

Steps 1–4 and 6–8 are independent of the companion pipeline; step 5 is the
join.

---

## 7. Context: the companion seed-building pipeline

TE-Aid v2 will be pinned as a **git submodule** of a separate pipeline repo
(planned name `TEbed-seeds`) that screens candidate TE families from a
multi-tool VGP repeat-annotation track hub, rebuilds each family's consensus
from its genomic copies, and emits **Dfam-compliant Stockholm seed
alignments** for a curator to approve or reject. TE-Aid v2 is that pipeline's
QC step, invoked with `--pipeline --stk`.

What that means for interfaces here:

- **Stockholm sequence identifiers are Smitten format**:
  `GCA_951799975.1:OX637595.1:15848-16090_+` — assembly accession, sequence
  name, **1-based fully-closed** coordinates, strand. When reading a seed,
  parse these to recover genomic loci; when reporting coordinates, do not
  silently mix them with BED16's 0-based half-open convention.
- Dfam seeds use `.` as the gap character and carry required `#=GF` fields
  (`DE`, `AU`, `TP`, `OC`, `SQ`) plus a `#=GC RF` consensus line. Spec:
  `Dfam_Seeds.md` in `https://github.com/Dfam-consortium/dfam-curator`; that
  repo also ships `stk lint`, which the pipeline uses as its acceptance gate.
  TE-Aid does not need to validate seeds, but should fail gracefully and
  informatively on a malformed one.
- Fail-soft contract: if TE-Aid errors on a packet, the pipeline still queues
  that packet for the curator with a note. Exit codes and stderr should make
  the failure reason machine-readable.

---

## 8. Ask the maintainer before assuming

1. ~~Path to the VGP repeat-hub repo~~ **Answered:** `~/Documents/VGP_TEbed`
   (`docs/INPUT_FORMAT.md` confirmed present).
2. ~~Whether a BATH install already exists locally~~ **Answered by inspection:**
   no `bathsearch`/`bathbuild` on this machine, so step 6 must include building
   BATH from source. Present and usable: `blastn`/`blastp`/`makeblastdb`
   (Homebrew), EMBOSS `getorf` and `dotmatcher` (MacPorts, `/opt/local/bin`),
   `nhmmer` and `hmmscan` (Homebrew). Python 3.13.7.
3. ~~Test genome + consensus sequences for development~~ **Answered:** use
   RepeatMasker `.out` + `.fasta` libraries from GenomeArk systematic
   annotations:
   `https://genomeark.s3.amazonaws.com/index.html?prefix=downstream_analyses/repeats/systematic_annotations/RepeatModeler-v2.0.8/`
4. `te_domains.tsv` review, once drafted (§5.4). **Open.**
