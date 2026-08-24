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

---

## 2. Input routes

Four routes, designed to compose:

| Route | Flags | Use when |
|---|---|---|
| **a** | `--annot FILE --consensus FASTA --family NAME` | you have an annotation you trust |
| **b** | `--blastn` *(not implemented)* | you have only a genome and a consensus |
| **c** | `--stk FILE [--family ID]` | you have a Stockholm seed |
| **d** | `--annot FILE --stk FILE --family NAME` | you have both — the richest sheet |

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

`--pipeline` reverses the standalone input priority (`--annot`, `--blastn`,
`--stk`) to seed-first, for callers whose primary artefact is a seed. It does
not change route d, where the seed is primary either way.

---

## 3. The fail-soft contract

A caller processing many families must be able to keep going when one fails and
record *why*. Every failure therefore carries **both** a distinct exit code and a
stable machine-readable slug on stderr:

```
teaid: error [bad-seed]: families.stk: no record with #=GF ID 'absent'
       409 records, e.g. ltr-1_family-1, ltr-1_family-2, …
```

| Exit | Slug | Meaning |
|---|---|---|
| 0 | — | success |
| 2 | `usage` | bad or missing arguments |
| 3 | `no-input` | an input file is missing or unreadable |
| 4 | `no-family` | the family is not in the consensus library |
| 5 | `no-evidence` | the family exists but carries nothing plottable |
| 6 | `bad-seed` | the Stockholm file could not be parsed |

**Parse the slug, not the prose.** Messages will be improved; slugs and codes
will not change without a version bump. Warnings (`teaid: warning: …`) are
non-fatal and do not affect the exit code.

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

**Spanning is not depth.** Panel 2 counts copies that *span* a position,
internal deletions included. Panel 6 counts sequences that contribute an actual
*base* there. The gap between them is the family's internal deletion structure,
and it is large in practice — across 409 GenomeArk seeds, 362 differ, by up to
37 sequences at one position. **If you gate seeds on "depth ≥ 3", decide which
one you mean.** A span-based gate is optimistic: a family can pass it and still
have fewer than three real bases inside a common deletion.

**Fragments.** BED16 column 16 (`hit_id`, the `.out` `ID` column) groups
fragments of one insertion split by an indel or a nested element.
`Annotation.fragment_groups()` exposes the grouping. Counting fragments as
independent copies inflates copy number — the specific failure that motivated
v2's move away from search-driven discovery.

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

data = report.SheetData(family="…", consensus_length=…, copies=…, …)
figure = report.build_figure(data, theme.LIGHT)
report.write_html(figure, "sheet.html", data, theme.LIGHT)
```

Every reader returns the same `Annotation` of `Copy` records, so a new input
format only needs a reader — nothing downstream changes.

---

## 6. Forking

Fork when you need behaviour that would not make sense to someone curating TEs
who has never heard of your pipeline. Project-specific gating, bespoke panels,
house classification vocabularies, and pipeline-shaped inputs all belong in a
fork rather than upstream.

Design your fork so it can track upstream:

- **Add readers, do not edit them.** A new input format is a module in
  `teaid/readers/` returning an `Annotation`; register it in `readers/__init__`.
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
| Panel 5, homology evidence (protein + nucleotide hits) | BATH; work-order step 6 |
| Panel 7, consensus vs contributing library entries | an agreed input hook — see below |
| `--blastn` legacy route | work-order step 10 |
| TSD detection, `#=GF TD` output | the flank-analysis path; work-order step 9 |

**Panel 7's input hook is deliberately unspecified.** It needs the source
library entries a rebuilt consensus came from, which only a producer has. When
it is designed it should take something generic — a FASTA of candidate source
entries — rather than any one pipeline's internal structure. If you are building
a seed producer, that is the interface to propose.
