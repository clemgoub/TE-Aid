# TE-Aid v2 [![publication](https://img.shields.io/badge/v1_publication:-Mobile_DNA-blue)](https://doi.org/10.1186/s13100-021-00259-7)

<img src=https://i.imgur.com/pxxR3Ec.png width="500">

> **🚧 This is the `v2` development branch** — a Python rewrite, still incomplete.
> For the published, stable version (`shell`+`R`, [Mobile DNA 2022](https://doi.org/10.1186/s13100-021-00259-7)),
> use the [`v1-legacy`](https://github.com/clemgoub/TE-Aid/tree/v1-legacy) branch,
> which also carries the full v1 documentation.

**TE-Aid** builds a single **evidence sheet for one TE family** to support manual
curation. It renders evidence and leaves the verdict to the curator — it never
asserts a classification.

The working brief for this rewrite, including every settled design decision, is
[`docs/BRIEF_v2.md`](docs/BRIEF_v2.md).

## What changes in v2

1. **Annotation-driven, not search-driven.** v1 rediscovered copies with
   `blastn`, which misses diverged copies and splits elements at indels. v2
   reads an annotation you already trust — RepeatMasker `.out`, GFF3, or BED16 —
   so coordinates, divergence and strand are read, not re-derived. No
   `makeblastdb`, no genome-wide search: seconds per family.
2. **Python, with an interactive HTML sheet** (static PDF/PNG retained for papers).
3. **Stockholm (Dfam seed) input** — a seed carries the copies, their loci, the
   alignment and the consensus in one file.
4. **Stronger homology evidence**: frameshift-aware translated pHMM search (BATH),
   and optional nucleotide hits against a taxonomic slice of Dfam.

## Status

| Step | State |
|---|---|
| Package skeleton, CLI, annotation readers (`.out` / GFF3 / BED16) | ✅ done |
| Panels 1–3 (copies vs divergence · coverage · self dot-plot) | ✅ done |
| Interactive HTML sheet + static export, light and dark | ✅ done |
| Panel 4 (structure: terminal-repeat candidates) | ✅ done |
| `--stk` Stockholm seed input, `--pipeline`, `--seed-qc` | ✅ done |
| Panel 5 (ORFs + protein homology, BATH) | ✅ done |
| `#=GF TP` disagreement flag | ✅ done |
| Benchmark (cost and relative sensitivity) | ⬜ **next** |
| Nucleotide row (Dfam slice + nhmmer) | ⬜ |
| TSD detection, `#=GF TD` output | ⬜ |
| `--blastn` legacy path | ⬜ |
| Panel 7 (consensus vs library entries) | ⬜ needs an agreed input hook |
| Whole-library runs, batch mode · CLI pass · GUI | ⬜ raised, not yet scoped |

The remaining steps are numbered and tracked in
[`docs/BRIEF_v2.md`](docs/BRIEF_v2.md) §6, which is the authority on order and
on what "done" means for each; this table is a summary of it.

The sheet keeps **v1's 2×2 quadrant layout** — copies vs divergence and coverage
on top, the self dot-plot (square, 1:1) and structure below — so the whole family
is readable at a glance. The bottom-right quadrant splits: panel 4 shows what the
dot-plot implies (terminal repeats), panel 5 shows what sequence annotation
implies (ORFs with their protein domains inlaid).

ORFs keep v1's strand convention — black outline forward, red reverse — and
protein domains keep **v1's TE-class colours**: green LTR, blue LINE, salmon DNA
transposon.

## Install

```bash
git clone -b v2 https://github.com/clemgoub/TE-Aid.git
cd TE-Aid
python3 -m venv .venv && source .venv/bin/activate
pip install -e .              # add '[static]' for PNG/PDF export via kaleido
```

Requires Python ≥ 3.10, plus external tools on `PATH`:

- `blastn` ([NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/)) for the self dot-plot
- `getorf` ([EMBOSS](https://emboss.sourceforge.net/)) for the ORF track
- `bathsearch`, `bathbuild`, `bathconvert` ([BATH](https://github.com/TravisWheelerLab/BATH))
  for the protein row

Any may be absent: the sheet still renders, and the affected quadrant says what
is missing rather than disappearing.

> **The v1 files are still in this branch and are not used by v2.** `TE-Aid`,
> `consensus2genome.R`, `blastndotplot.R`, `Run-c2g.R`, `reduce.cpp`,
> `loop_TE-Aid.sh`, `extractfasta.sh`, `getlength.sh`, `dev/` and `Example/` are
> the published shell+R implementation, kept here until the `v2.0` tag. In
> particular **do not use `TE_AID.yml`** — it is v1's conda environment, pinning
> R and EMBOSS with no Python and no BATH. Use the venv above.

### The protein library

Panel 5 searches a library in three tiers:

| Tier | What | Size |
|---|---|---|
| 1 | 130 curated Pfam TE domains ([`teaid/data/te_domains.tsv`](teaid/data/te_domains.tsv)) | 7.3 MB |
| 2 | pHMMs for the superfamilies Pfam cannot model — piggyBac, Maverick, Crypton, Penelope | 225 MB |
| 3 | `blastp` of ORF peptides against all of RepeatPeps, as v1 did | 17 MB |

**The first run builds it, and that takes a while:** roughly 130 sequential
requests to InterPro plus two `bathbuild` passes — about **5 minutes and ~240 MB**,
network permitting. Every later run reuses the cache and costs ~10 s per family.
The build is safe to interrupt: it writes through temporary files and renames, so
a cancelled build leaves nothing behind and simply restarts.

The cache lives under `$TEAID_CACHE`, or `~/.teaid/proteins/` if that is unset,
and is keyed by a hash of the domain table, so editing the table rebuilds rather
than silently reusing a stale library. Three flags cover the awkward cases:

- `--proteins FILE` — search a **prebuilt** library and skip the build entirely.
  The route for offline or read-only machines. Note it also disables tier 3.
- `--rebuild-proteins` — force a rebuild, if a cache is ever suspect.
- `--repeatpeps FILE` — point at `RepeatPeps.lib` explicitly.

Tiers 2 and 3 need **`RepeatPeps.lib`**, which is not separately downloadable: it
ships inside a [RepeatMasker](https://www.repeatmasker.org/) installation, under
`Libraries/`. TE-Aid looks in `~/`, `~/Downloads/RepeatMasker/Libraries/`,
`/usr/local/RepeatMasker/Libraries/`, `/opt/RepeatMasker/Libraries/` and
`/opt/homebrew/share/RepeatMasker/Libraries/`. **Without it you still get a
useful protein row** from tier 1 alone — you lose only the four superfamilies
Pfam cannot model, and the sheet says so.

`--deep` searches the whole of RepeatPeps as pHMMs, which is **~6.4 GB**. It
*replaces* tiers 1–2 rather than adding to them, and turns off tier 3 too, so a
deep run has no Pfam accessions and panel 5 loses its order-derived colours. See
[`docs/BRIEF_v2.md`](docs/BRIEF_v2.md) §5.1a for why it is not the default.

A hit labelled **`2fs`** or **`1⊗`** carries that many frameshifts or in-frame
stops: a domain that was once coding and has since been disrupted. Finding one is
a different result from finding nothing, and an ORF-finder-then-align search
cannot find it at all — which is the whole reason the search is BATH. A hollow
arrow labelled **`=`** is a tier-3 `blastp` match to a *named element* rather
than a domain model.

Panel 5 draws only as many hits as it can show legibly, choosing the strongest
within each tier. **Every hit found is listed in an expandable table under the
grid**, with the collapsed ones dimmed, and can be exported with *Download TSV*
or *Copy TSV* — so nothing the search found exists only in a note about how many
were dropped.

## Usage

```bash
teaid --annot genome.fa.out \
      --consensus families.fa \
      --family rnd-1_family-257 \
      --output sheets/
```

This writes `sheets/rnd-1_family-257.teaid.html`. Add `--static png` (or `pdf`,
`svg`) for a publication figure alongside it, and `--theme dark` for a dark sheet.

The sheet states its own provenance under the title — which inputs built it, and
the **detected** annotation format, so a misdetection is visible rather than
silent:

```
from  seed families.stk (Stockholm, 44 sequences)  +  annotation genome.fa.out (rmout)
```

In the HTML sheet you can drag to pan, scroll to zoom, and double-click a panel
to autoscale it. **Reset view** returns all four quadrants to the full consensus
span at once — unlike Plotly's built-in reset, which autoscales each panel
independently and so leaves them on different x-ranges. Each panel title carries
a faded **?**: hover it for what the panel shows and the trap it exists to avoid.
The markers are dropped from static exports, where a question mark would have no
answer.

The HTML sheet loads plotly.js from a CDN, so **viewing** it needs network access
— it is written fine offline but renders blank without one. Static exports
(`--static png|pdf|svg`) are self-contained and are the right choice for
archiving or for an air-gapped machine.

The annotation format is detected from the file's content, not its extension —
`.bed` files converted from `.out` are common enough that extensions cannot be
trusted. Override with `--annot-format` if detection gets it wrong.

### Input routes

Four, and they compose:

| Route | Flags |
|---|---|
| annotation you trust | `--annot FILE --consensus FASTA --family NAME` |
| genome search (v1 path, not yet implemented) | `--blastn` |
| Stockholm seed | `--stk FILE [--family ID]` |
| **both together** | `--annot FILE --stk FILE --family NAME` |

The last is not a contest between the two. The seed says which copies were
*chosen*; the annotation says which *exist*. Given both, the seed supplies the
consensus and the seed-QC panels while panels 1–2 show every annotated genomic
copy, and the sheet states the comparison — `panels 1-2 show 515 annotated
genomic copies; the seed uses 44`. That ratio is the QC signal.

### From a Stockholm seed

A Dfam seed alignment carries the copies, their genomic loci, the alignment and
the consensus in one file, so nothing else is needed:

```bash
teaid --stk families.stk --family rnd-1_family-257 --seed-qc -o sheets/
```

`--seed-qc` appends the seed-QC panels below the grid. Panel 6 follows the shape
of [Dfam's own seed-alignment track](https://dfam.org/family/DF000001423/browser):
a coverage band split into sequences that match the consensus and sequences that
differ, over a pileup of the individual seed sequences drawn as their aligned
runs, so an internal deletion shows as a gap in a lane rather than being smoothed
over. Dfam's 3-sequence floor is a rule across the band and stretches falling
short are shaded. Panel 8 shows the seed's `#=GF TP` as a supplied claim to be
checked, not as a conclusion.

Splitting coverage by agreement matters: a column with 43 aligned sequences of
which 25 disagree is not 43 sequences of support, and depth alone hides that.

`--seed-qc` is off by default, so the default sheet is the same four quadrants
whatever the input was.

Sequence identifiers are Smitten format, in either the 4-part
(`GCA_951799975.1:OX637595.1:15848-16090_+`) or the 2-part
(`OY720097.1:14692470-14693460_+`) shape; both are read, and their 1-based
closed coordinates are converted to the package's 0-based half-open convention
on the way in.

`--pipeline` reverses the input priority to seed-first, for callers whose primary
artefact is a seed.

Run `teaid --help` for the full option list.

### A note on divergence

`perc_div` does not mean the same thing in every file. For homology-based tools
it is divergence from a library consensus (a proxy for the age of an insertion);
for tandem-repeat finders it is array homogeneity; some tools report nothing at
all. TE-Aid tracks which of the three it has, labels the axis accordingly, and
**never plots an absent divergence as zero** — a zero there would read as a
pristine, very recent insertion. Copies with no divergence are left out of panel
1 and counted in the axis label. Use `--divergence-kind` when you know what
produced the file.

### Exit codes

Non-zero codes are distinct so a pipeline can tell why a family failed, and the
ones TE-Aid diagnoses itself also print a stable slug on stderr:

```
teaid: error [no-family]: family 'x' not in families.fa
```

| Code | stderr slug | Meaning |
|---|---|---|
| 0 | — | a sheet was written (possibly with panels dropped — see below) |
| 1 | *none* | uncaught internal error; a traceback |
| 2 | *none* | argument error, straight from `argparse` |
| 3 | `no-input` | an input file is missing or unreadable |
| 4 | `no-family` | the family is not in the consensus library |
| 5 | `no-evidence` | the family has nothing plottable |
| 6 | `bad-seed` | the Stockholm file could not be parsed |

Parse the slug, not the prose — and note that codes 1 and 2 carry **no** slug, so
branch on the code first. **Exit 0 does not mean a complete sheet:** a missing
`blastn`, `getorf` or `bathsearch` drops the affected panel with a
`teaid: warning: …` line and still exits 0. Callers that need a full sheet must
check stderr for warnings.

[`docs/INTEGRATION.md`](docs/INTEGRATION.md) §3 is the full contract.

## Building on TE-Aid

TE-Aid is deliberately general — nothing in it is specific to any pipeline,
consortium or assembly. If you are writing a tool that calls it, or forking it to
add project-specific behaviour, read **[docs/INTEGRATION.md](docs/INTEGRATION.md)**:
input routes, the fail-soft contract, the coordinate and divergence conventions
that bite, using the package as a library, and how to fork so your fork can still
track upstream.

## Development

```bash
pip install -e '.[dev]'
pytest                        # ~250 tests, about 30 seconds
```

The suite never touches your real protein cache: `tests/conftest.py` redirects
`TEAID_CACHE` to a temporary directory unless you set it yourself. To exercise
the real library instead, point it at a prebuilt one:

```bash
TEAID_CACHE=dev-data/protein-cache pytest
```

A quick check that a tree is healthy rather than subtly broken — `teaid --version`
should print `teaid 2.0.0.dev0`, and:

```bash
teaid --annot dev-data/GCA_963082875.1.fa.out \
      --consensus dev-data/GCA_963082875.1-families.fa \
      --family ltr-1_family-65 -o /tmp/check
```

should end with exactly:

```
ltr-1_family-65: 16 copies, 0 full length, consensus 2,092 bp, terminal-repeat candidates: 1 LTR-like, 1 TIR-like
```

If BATH is not on `PATH` you also get `warning: protein row skipped: …` and still
exit 0 — that is the fail-soft contract working, not a broken tree.

Test data used during development comes from the
[GenomeArk systematic repeat annotations](https://genomeark.s3.amazonaws.com/index.html?prefix=downstream_analyses/repeats/systematic_annotations/RepeatModeler-v2.0.8/)
(RepeatMasker `.out`, family FASTA, and Stockholm seeds for 482 assemblies).
Downloads belong in `dev-data/`, which is git-ignored.

## Citing

The v1 method is described in
["A beginner's guide to manual curation of transposable elements"](https://doi.org/10.1186/s13100-021-00259-7),
Goubert, Craig, Bilat, Peona, Vogan & Protasio, *Mobile DNA* (2022).

Support: open an issue, or [email the maintainer](mailto:goubert.clement@gmail.com).
