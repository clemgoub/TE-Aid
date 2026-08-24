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
3. **Stockholm (Dfam seed) input**, for use with a companion seed-building pipeline.
4. **Stronger homology evidence**: frameshift-aware translated pHMM search (BATH),
   and optional nucleotide hits against a taxonomic slice of Dfam.

## Status

| Step | State |
|---|---|
| Package skeleton, CLI, annotation readers (`.out` / GFF3 / BED16) | ✅ done |
| Panels 1–3 (copies vs divergence · coverage · self dot-plot) | ✅ done |
| Interactive HTML sheet + static export, light and dark | ✅ done |
| Panel 4 (structure: ORFs, TIR/LTR candidates) | ✅ done |
| Panel 5 (homology evidence) | ⬜ needs BATH |
| `--stk` reader and `--seed-qc` panels | ⬜ |
| `--blastn` legacy path | ⬜ |

The sheet keeps **v1's 2×2 quadrant layout** — copies vs divergence and coverage
on top, the self dot-plot (square, 1:1) and structure below — so the whole family
is readable at a glance. Panel 5 will split the bottom-right quadrant with panel
4 once there is homology evidence to draw.

## Install

```bash
git clone -b v2 https://github.com/clemgoub/TE-Aid.git
cd TE-Aid
python3 -m venv .venv && source .venv/bin/activate
pip install -e .              # add '[static]' for PNG/PDF export via kaleido
```

Requires Python ≥ 3.10, plus two external tools on `PATH`:

- `blastn` ([NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/)) for the self dot-plot
- `getorf` ([EMBOSS](https://emboss.sourceforge.net/)) for the ORF track

Either may be absent: the sheet still renders, and the affected quadrant says
what is missing rather than disappearing.

## Usage

```bash
teaid --annot genome.fa.out \
      --consensus families.fa \
      --family rnd-1_family-257 \
      --output sheets/
```

This writes `sheets/rnd-1_family-257.teaid.html`. Add `--static png` (or `pdf`,
`svg`) for a publication figure alongside it, and `--theme dark` for a dark sheet.

In the HTML sheet you can drag to pan, scroll to zoom, and double-click a panel
to autoscale it. **Reset view** returns all four quadrants to the full consensus
span at once — unlike Plotly's built-in reset, which autoscales each panel
independently and so leaves them on different x-ranges.

The annotation format is detected from the file's content, not its extension —
`.bed` files converted from `.out` are common enough that extensions cannot be
trusted. Override with `--annot-format` if detection gets it wrong.

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

Non-zero codes are distinct so a pipeline can tell why a family failed:

| Code | Meaning |
|---|---|
| 0 | success |
| 2 | usage error |
| 3 | an input file is missing or unreadable |
| 4 | the family is not in the consensus library |
| 5 | the family has no annotated copies |

## Development

```bash
pip install -e '.[dev]'
pytest
```

Test data used during development comes from the
[GenomeArk systematic repeat annotations](https://genomeark.s3.amazonaws.com/index.html?prefix=downstream_analyses/repeats/systematic_annotations/RepeatModeler-v2.0.8/)
(RepeatMasker `.out`, family FASTA, and Stockholm seeds for 482 assemblies).
Downloads belong in `dev-data/`, which is git-ignored.

## Citing

The v1 method is described in
["A beginner's guide to manual curation of transposable elements"](https://doi.org/10.1186/s13100-021-00259-7),
Goubert, Craig, Bilat, Peona, Vogan & Protasio, *Mobile DNA* (2022).

Support: open an issue, or [email the maintainer](mailto:goubert.clement@gmail.com).
