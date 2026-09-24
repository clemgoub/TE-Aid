# Integrating with TE-Aid v2

For anyone building a tool that **calls** TE-Aid, or **forks** it to add
project-specific behaviour. If you are here to modify TE-Aid for one pipeline,
read [Forking](#forking) first.

TE-Aid renders evidence for one TE family so a human can judge it. It is
deliberately general: nothing in this repository is specific to any pipeline,
consortium or assembly. That is what makes it worth building on.

---

## 1. What TE-Aid is, in one paragraph

One invocation produces one evidence sheet for one family: an interactive HTML
page (and optionally PNG/PDF/SVG) laid out as v1's 2×2 quadrant grid — copies vs
divergence, consensus coverage, a square self dot-plot, and structure — with an
optional seed-QC row appended. **It never asserts a classification.** That is a
hard boundary, not a default: in seed-QC mode the sheet flags disagreement
between an externally supplied expected class and the evidence shown, and that
flag only means anything because the sheet itself has no opinion.

The fourth quadrant splits in two (0.42 / 0.58) as soon as there is anything to
draw in it: **panel 4** is what the *dot-plot* implies — terminal-repeat
candidates only — and **panel 5** is what *sequence annotation* implies — ORFs
with their protein-domain hits. The split is by provenance, and a fork that
mixes the two loses the one thing a curator reads the quadrant for: whether
two independent lines of evidence agree.

---

## 2. Input routes

Four routes, designed to compose:

| Route | Flags | Use when |
|---|---|---|
| **a** | `--annot FILE --consensus FASTA --family NAME` | you have an annotation you trust |
| **b** | `--genome FASTA --consensus FASTA` (v1: `-g`, `-q`; also `--blastn`) | you have only a genome and a consensus |
| **c** | `--stk FILE [--family ID]` | you have a Stockholm seed |
| **d** | `--stk FILE` with `--annot FILE` or `--genome FASTA` | you have a seed and its genome — the richest sheet |

`--annot` accepts RepeatMasker `.out`, RepeatMasker GFF, or BED16, gzipped or
not. **The format is detected from content, not extension**, because `.bed`
files converted from `.out` are common; override with `--annot-format`.

### Route d: seed and annotation together

These are **not** alternatives. The seed says which copies were *chosen*; the
annotation says which *exist*. Given both:

- the **seed stays primary** — it supplies the consensus and the seed-QC panels;
- **panels 1–2 show every annotated genomic copy**, not just the seed's;
- the sheet states the comparison (`panels 1-2 show 515 annotated genomic
  copies; the seed uses 44`).

That ratio is the QC signal: a seed built from 44 of 515 copies may be
representative or may have sampled one corner of the family, and only the
comparison distinguishes them. If the family is absent from the annotation,
TE-Aid warns on stderr and falls back to the seed's own sequences — it does not
fail.

`--pipeline` reverses the standalone input priority (`--annot`, `--genome`,
`--stk`) to seed-first, for callers whose primary artefact is a seed. It does
not change route d, where the seed is primary either way. `--annot` and
`--genome` together, without a seed, use the annotation and say on the sheet
that the genome was ignored.

### Route b: v1's genome search

The consensus is searched against the genome with `blastn` (megablast, E ≤ 10e-8
by default, `-e` to change), exactly as v1 did, and every HSP becomes one copy.
That makes this the leakiest route: megablast misses diverged copies, and an
insertion interrupted by another element arrives as two copies, because a
search has no ID column to join them. The sheet's notes say so. Its divergence
axis reads `100 − blastn identity (%)`, not divergence: an uncorrected local
identity is not RepeatMasker's Kimura-corrected value.

The blastn database is built beside the genome on first use
(`makeblastdb -in G -out G`, as v1), so a v1 database is reused; if that
directory is not writable it goes under the cache. `-r` applies v1's
`reduce.cpp` rule, dropping hits nested inside another hit's locus. Without
`blastn` or `makeblastdb` the route has no copies, so it fails with `no-input`
(exit 3) rather than dropping a panel.

---

## 3. The fail-soft contract

A caller processing many families must be able to keep going when one fails and
record *why*. Failures that TE-Aid itself diagnoses carry a distinct exit code
**and** a stable machine-readable slug on stderr:

```
teaid: error [bad-seed]: families.stk: no record with #=GF ID 'absent'
       409 records, e.g. ltr-1_family-1, ltr-1_family-2, …
```

| Exit | Slug | Meaning |
|---|---|---|
| 0 | — | a sheet was written — but see *exit 0 is not "complete"* below |
| 1 | *none* | uncaught internal error; a Python traceback, not a message |
| 2 | `usage` | argument error: an unknown flag, a missing input, a value out of range, or `--family` omitted for a multi-entry library |
| 3 | `no-input` | an input file is missing or unreadable, or the `--genome` route cannot run BLAST+ |
| 4 | `no-family` | the family is not in the consensus library |
| 5 | `no-evidence` | the family exists but carries nothing plottable |
| 6 | `bad-seed` | the Stockholm file could not be parsed |
| 8 | `batch-incomplete` | batch mode (`--all`, `--families`): at least one family failed; `results.tsv` says which and why |

**Only exit 1 has no slug.** It is an unhandled exception and prints a Python
traceback. Treat any failure whose stderr has no `[slug]` bracket as
fatal-unknown and log the stderr verbatim; do not try to parse it. Argument
errors from `argparse` itself carry the `usage` slug too, and they arrive as
`SystemExit(2)` when TE-Aid is called in-process.

A multi-record seed given without `--family` is not an argument error. It reports
`bad-seed` and exit 6, because the file was read and found ambiguous.

**Exit 0 is not "complete".** The sheet is fail-soft by design: if `blastn`,
`getorf` or `bathsearch` is missing or fails, the affected panel is dropped, a
`teaid: warning: …` line goes to stderr with **no slug**, and the run still
exits 0. Missing `blastn` costs the dot-plot and the terminal-repeat candidates
(panels 3 and 4); missing `getorf` costs the ORF rectangles; missing
`bathsearch`, or an unbuildable protein library, costs the domain arrows. A
caller that needs a *full* sheet must scan stderr for `warning:`, because the
exit code will not tell it.

**Parse the slug, not the prose.** Messages will be improved; slugs and codes
will not change without a version bump. Warnings are non-fatal and never affect
the exit code.

TE-Aid does not validate seeds beyond what it needs in order to draw them. If
you are *producing* seeds, gate them with `stk lint` from
[dfam-curator](https://github.com/Dfam-consortium/dfam-curator); TE-Aid is the
human-facing half of the gate, not a replacement for the machine half.

---

## 4. Data conventions

Get these wrong and the errors are silent and plausible-looking.

**Coordinates.** Internally every coordinate is **0-based half-open**, genomic
and consensus alike. Conversions happen once, at the reader:

| Source | Convention | Conversion |
|---|---|---|
| RepeatMasker `.out`, GFF | 1-based fully closed | `start - 1` |
| Smitten identifiers in Stockholm | 1-based fully closed | `start - 1` |
| BED16 genomic columns | 0-based half-open | none |
| BED16 `repeat_start`/`repeat_end` | 1-based inclusive | `start - 1` |

That last row is not a typo: a BED16 row genuinely mixes both conventions, its
genomic columns 0-based and its consensus columns inherited 1-based from the
`.out` it was converted from.

**Smitten identifiers.** Both shapes are read:

```
GCA_951799975.1:OX637595.1:15848-16090_+     assembly, sequence, span, strand
OY720097.1:14692470-14693460_+               sequence, span, strand
```

Emit the 4-part form if you have an assembly accession. The 2-part form is what
RepeatModeler writes and what all 482 GenomeArk seed sets use, so a reader that
handles only the documented 4-part form reads none of them.

**Divergence is not one quantity.** Every record carries a `divergence_kind` of
`consensus`, `array_homogeneity`, or `none`, because `perc_div` means divergence
from a library consensus for homology tools but *array homogeneity* for
tandem-repeat finders, and nothing at all for pure detectors. **Absent
divergence is never plotted as zero** — a zero there reads as a pristine, very
recent insertion, the opposite of "unknown". BED16 does not record which kind it
holds, so it is inferred from the class label; override with
`--divergence-kind`.

`--divergence-kind` **applies to BED16 only**. `.out` and GFF3 state their own
kind, so the flag is silently discarded for them (`readers/__init__.py`) with no
warning — it is not a general override, and a caller that passes it for a `.out`
gets no error and no effect.

**Spanning is not depth.** Panel 2 counts copies that *span* a position,
internal deletions included. Panel 6 counts sequences that contribute an actual
*base* there. The gap between them is the family's internal deletion structure,
and it is large in practice — across 409 GenomeArk seeds, 362 differ, by up to
37 sequences at one position. **If you gate seeds on "depth ≥ 3", decide which
one you mean.** A span-based gate is optimistic: a family can pass it and still
have fewer than three real bases inside a common deletion.

**The protein library is a cache, and callers must manage it.** Panel 5 searches
a BATH pHMM library that is *built on first use*: ~130 sequential HTTPS requests
to InterPro plus two `bathbuild` runs, roughly five minutes and ~240 MB, cached
under `$TEAID_CACHE` (default `~/.teaid/proteins`) and keyed by a hash of the
curated domain table. Three consequences for anyone calling TE-Aid:

- **Point `TEAID_CACHE` at shared, writable storage** so every invocation reuses
  one library rather than each user building their own.
- **`build()` locks while it builds.** Concurrent first runs queue on an
  exclusive `flock` at `<cache>/<signature>.lock`: the first builds, and the rest
  wait and then reuse what it finished. Warming the cache with one serial run is
  no longer required, but it still avoids a long silent wait in the workers. The
  lock needs a writable cache directory; on a read-only one, prebuild the
  library and pass `--proteins`.
- **`--proteins FILE` searches a prebuilt library and skips the build entirely**,
  which suits an offline or read-only environment. It replaces tiers 1–2 only:
  tier 3 still runs against `--repeatpeps` or the discovered RepeatPeps.lib, as in
  any other run. The sheet's provenance line names the library searched, and a
  note says why whenever tier 3 did not run.

`--rebuild-proteins` forces a rebuild; it is the recovery path if a cache is
ever suspect. `--deep` **replaces** tiers 1–2 rather than adding to them, and
turns off tier 3 as well, so a deep run carries no Pfam accessions and panel 5
loses its order-derived colours.

**Fragments.** Rows that share an ID are one copy: BED16 column 16, the `.out`
`ID` column, or the GFF3 `ID` attribute, all read into `Copy.hit_id`. TE-Aid
groups them by default, keyed on chromosome and ID, and counts copies, not
rows, in panels 1 and 2, the full-length count, the subtitle and the stdout
summary. `--no-group-fragments` counts rows instead. A `.out` written with
RepeatMasker `-no_id` has no ID column and is read with every row its own copy.

What a grouped copy measures:

- **coverage and full length** use the union of the consensus its rows align
  to, never the span from first to last row;
- **divergence** is the rows' mean weighted by consensus length, RepeatMasker's
  own rule when it merges annotations;
- **each join** between rows is classified from the genomic and consensus gaps:
  an *insertion* in the copy when the genome carries more (often a nested
  element), a *deletion* when the consensus carries more, and *flush* within
  50 bp. Panel 1 draws the first two.

A producer that writes one constant into the ID column would collapse a family
into one copy. The sheet's notes line names the largest group, which is how
that shows up.

---

## 5. Using TE-Aid as a library

The CLI is a thin wrapper. The pieces are independently usable:

```python
from teaid import analysis, readers, report, theme
from teaid.readers import stockholm
from teaid.sequences import read_fasta

annotation = readers.read("genome.fa.out")          # format auto-detected
family = annotation.for_family("rnd-1_family-257", exact=False)
consensus = read_fasta("families.fa").get("rnd-1_family-257")

coverage = analysis.coverage(family, len(consensus))
full = analysis.full_length(family, len(consensus), threshold=0.9)
hits = analysis.self_blast(consensus.sequence)
terminal = analysis.terminal_repeats(hits, len(consensus))

seed = stockholm.read("families.stk", "rnd-1_family-257")
depth = seed.depth()                 # base-level, per consensus position
copies = seed.to_annotation()        # same record struct as any annotation

data = report.SheetData(                 # these five are required, the rest default
    family="rnd-1_family-257",
    consensus_length=len(consensus),
    copies=family.copies,
    full_length=full,
    coverage=coverage,
)
figure = report.build_figure(data, theme.LIGHT)
report.write_html(figure, "sheet.html", data, theme.LIGHT)
```

Every reader returns the same `Annotation` of `Copy` records, so a new input
format only needs a reader — nothing downstream changes.

### Outputs

The CLI writes into `-o/--output` (default: `.`, created if absent):

```
<output>/<stem>.teaid.html            always
<output>/<stem>.teaid.{png,pdf,svg}   with --static png|pdf|svg (needs kaleido)
<output>/<stem>.selfhits.tsv          with --tables or --all-tables
<output>/<stem>.orfs.tsv              with --tables or --all-tables
<output>/<stem>.protein_hits.tsv      with --tables or --all-tables
<output>/<stem>.copies.tsv            with --all-tables
```

The tables are v1's `-t` and `-T`. Each has a header row and uses **1-based,
fully closed** coordinates. `copies.tsv` has one row per annotation row, with a
`copy` number shared by the rows grouped into one insertion, so it can be read
at either level. `protein_hits.tsv` is the sheet's hit table, `drawn` column
included, and leads with a `consensus` column so tables from many families
concatenate. Self hits keep blastn's orientation: `s_start > s_end` on a reverse
hit.

`<stem>` is the consensus name up to the first `#`, with every character outside
`[A-Za-z0-9-._]` replaced by `_` — so `rnd-1_family-257#LTR/Gypsy` writes
`rnd-1_family-257.teaid.html`. A seed with no `#=GF ID` is named `seed`.
In a single-family run, **collisions overwrite silently**: two families whose
names differ only in a replaced character land on the same file. Batch mode
checks for this up front and writes the later ones as `<stem>.2`, `<stem>.3`,
with a warning and the stem recorded in `results.tsv`.

### Batch mode

`--all` draws every family: each record of a seed (`--stk`), each consensus
with annotated copies (`--annot`), or each consensus (`--genome`). `--families
FILE` draws those named, one per line, reading only the first tab-separated
column, so `tools/select_families.py` output works as it is; `-` reads stdin.
`-j N` draws N families at once, in separate processes.

Each input is read once: the annotation (a 147k-row `.fa.out` took several
seconds a family), the consensus library, every seed record, the genome's
blastn database and the protein library, which is built before any worker
starts. Each worker receives only its own family's rows.

Besides the sheets, a batch writes:

```
<output>/results.tsv   one row per family, with a header
<output>/index.html    a contact sheet over the batch, grouped by TE class
```

`results.tsv` columns: `family stem exit slug seconds copies annotation_rows
full_length class consensus_length warnings error`. `warnings` joins the
family's warning lines with ` | `; `error` holds its error line. One family's
failure or crash does not stop the batch. The run exits 0 only when every family
drew, and 8 (`batch-incomplete`) otherwise. Families `--all` skipped for having
no annotated copies are not failures: the run says how many on stderr and on the
index. `--no-index` skips the contact sheet.

Stdout keeps the single-family contract: one `wrote <path>` line per file,
including the two above, plus each family's summary line. Stderr carries one
progress line per family (`teaid: [3/50] fam: ok (12 s)`) rather than each
family's own warnings, which are in the manifest.

Each file written prints one `wrote <path>` line on **stdout**. Errors and
warnings go to stderr, so stdout stays parseable.

The HTML sheet loads plotly.js **from a CDN** — it needs network access when
*viewed*, not when written. The static exports are self-contained.

---

## 6. Forking

Fork when you need behaviour that would not make sense to someone curating TEs
who has never heard of your pipeline. Project-specific gating, bespoke panels,
house classification vocabularies, and pipeline-shaped inputs all belong in a
fork rather than upstream.

Design your fork so it can track upstream:

- **Add readers, do not edit them.** A new input format is a module in
  `teaid/readers/` returning an `Annotation`. Registering it takes **four**
  edits, and missing any one fails in a different way: add the module; add it to
  `_READERS`; add it to `FORMATS`, or `--annot-format yours` is rejected by
  `choices=readers.FORMATS`; and add a branch to `detect()` **before** the
  ≥16-column BED16 catch-all, which will otherwise claim any wide TSV first.
  Note that `read()` forwards `divergence_kind` to bed16 only, so a new format
  that needs it must be added to that guard too.
- **Add panels, do not renumber them.** `report.build_figure` assembles the grid
  from blocks. Keep the 2×2 quadrants and the square dot-plot — that layout is
  most of why the tool gets used, and a vertical stack that forces scrolling
  loses it.
- **Keep the classification boundary.** If your fork asserts a class, the
  seed-QC disagreement flag stops meaning anything.
- **Put project logic in the caller** where you can. Much of what a pipeline
  needs — selecting families, ordering work, deciding what passes — sits
  naturally outside TE-Aid entirely.

Upstream welcomes back anything genuinely general: a new annotation dialect, a
better terminal-repeat heuristic, a panel that any curator would want.

---

## 7. Not yet available

| Wanted | Blocked on |
|---|---|
| Panel 5's **nucleotide** row (Dfam slice + `nhmmer`) | work-order step 8 |
| Panel 7, consensus vs contributing library entries | an agreed input hook — see below |
| TSD detection, `#=GF TD` output | the flank-analysis path; work-order step 9 |

Panel 5's **protein** row ships: ORFs from `getorf` with frameshift-aware domain
hits from `bathsearch`. It needs BATH and EMBOSS on `PATH`; without them the row
is dropped with a warning and the run still exits 0.

**Panel 7's input hook is deliberately unspecified.** It needs the source
library entries a rebuilt consensus came from, which only a producer has. When
it is designed it should take something generic — a FASTA of candidate source
entries — rather than any one pipeline's internal structure. If you are building
a seed producer, that is the interface to propose.
