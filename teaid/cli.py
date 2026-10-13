"""Command line entry point.

Four input routes, designed to compose: ``--annot`` (an annotation you already
trust), ``--genome`` (v1's route: blastn of the consensus against a genome, also
spelled ``--blastn`` or ``-g``), ``--stk`` (a Stockholm seed), and ``--stk``
**with** ``--annot`` or ``--genome``, which is not a contest between the two but
the richest sheet, showing the family as it exists in the genome against the
copies a seed actually used.

Standalone priority among single routes is annot, genome, stk; ``--pipeline``
reverses that to seed-first for callers whose primary artefact is a seed.

**Fail-soft contract.** A caller processing many families must be able to keep
going when one fails and record why, so every failure is machine-readable: a
distinct exit code, and a stderr line carrying a stable slug::

    teaid: error [no-family]: family 'x' not in families.fa

Parse the slug, not the prose. See docs/INTEGRATION.md for the full contract.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from pathlib import Path

from . import (__version__, analysis, classcheck, genome, homology, orfs, proteins,
               readers, report, tables, theme)
from .readers import stockholm
from .records import Annotation, DivergenceKind
from .sequences import Consensus, read_fasta

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_NO_INPUT = 3  # a required file is missing or unreadable
EXIT_NO_FAMILY = 4  # the requested family is not in the annotation or library
EXIT_NO_EVIDENCE = 5  # the family exists but carries nothing plottable
EXIT_BAD_SEED = 6  # the Stockholm file could not be parsed
EXIT_BATCH_INCOMPLETE = 8  # batch mode: at least one family failed; see results.tsv

# Exit code paired with the slug printed on stderr, so a caller can branch on
# either. Keep the two in step.
_SLUG = {
    EXIT_USAGE: "usage",
    EXIT_NO_INPUT: "no-input",
    EXIT_NO_FAMILY: "no-family",
    EXIT_NO_EVIDENCE: "no-evidence",
    EXIT_BAD_SEED: "bad-seed",
    EXIT_BATCH_INCOMPLETE: "batch-incomplete",
}


def fail(code: int, message: str, hint: str | None = None) -> int:
    """Report a failure in the form the pipeline parses, and return its code."""
    print(f"teaid: error [{_SLUG.get(code, 'error')}]: {message}", file=sys.stderr)
    if hint:
        print(f"       {hint}", file=sys.stderr)
    return code


class _Parser(argparse.ArgumentParser):
    """argparse, with its errors reported in the fail-soft form.

    Stock argparse prints ``teaid: error: ...`` with no slug, so a batch caller
    could not tell a mistyped flag from a crash. With this override every error
    TE-Aid diagnoses carries a slug, and only an uncaught exception (exit 1)
    has none. The exit code stays 2 and argparse still raises SystemExit.
    """

    def error(self, message: str):
        self.print_usage(sys.stderr)
        self.exit(EXIT_USAGE, f"teaid: error [{_SLUG[EXIT_USAGE]}]: {message}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = _Parser(
        prog="teaid",
        description=(
            "Build an evidence sheet per TE family, to support manual curation. "
            "TE-Aid renders evidence; it never asserts a classification."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "example:\n"
            "  teaid --annot genome.fa.out --consensus families.fa "
            "--family rnd-1_family-257 -o sheets/\n"
        ),
    )

    source = parser.add_argument_group("input")
    source.add_argument(
        "--annot",
        metavar="FILE",
        help="annotation of genomic copies: RepeatMasker .out, GFF3, or BED16 (may be .gz)",
    )
    source.add_argument(
        "--annot-format",
        choices=readers.FORMATS,
        help="override annotation format detection",
    )
    source.add_argument(
        "--stk",
        metavar="FILE",
        help=(
            "Stockholm seed alignment (Dfam format, may be .gz). Supplies the "
            "copies, their loci, the consensus and the alignment in one file"
        ),
    )
    source.add_argument(
        "-g",
        "--genome",
        "--blastn",
        dest="genome",
        metavar="FASTA",
        help=(
            "reference genome, for the v1 route: copies are found by blastn of the "
            "consensus against it. Leakier than an annotation, and the sheet says so. "
            "The blastn database is built beside the genome on first use, as in v1"
        ),
    )
    source.add_argument(
        "-c",
        "--consensus",
        "-q",
        "--query",
        dest="consensus",
        metavar="FASTA",
        help=(
            "consensus FASTA; a library of many families or a single sequence "
            "(v1: -q/--query). Required with --annot and --genome; with --stk the "
            "seed supplies the consensus"
        ),
    )
    source.add_argument(
        # -f is the full-length threshold, as in v1; the family name moved to -n.
        "-n",
        "--family",
        metavar="NAME",
        help="family to report on; defaults to the only entry in a single-sequence FASTA",
    )
    source.add_argument(
        "--divergence-kind",
        choices=[k.value for k in DivergenceKind],
        help=(
            "what the annotation's divergence column measures, when you know it. "
            "BED16 does not record this; the default infers it from the class label"
        ),
    )

    source.add_argument(
        "--pipeline",
        action="store_true",
        help="reverse the input priority to seed-first, when a seed is the primary input",
    )
    source.add_argument(
        "--seed-qc",
        action="store_true",
        dest="seed_qc",
        help="add the seed-QC panels: alignment depth and the expected class label",
    )

    batch = parser.add_argument_group("batch (many families in one run)")
    batch.add_argument(
        "--all", action="store_true",
        help=("every family: each record of the seed, each consensus with annotated "
              "copies (--annot), or each consensus (--genome)"),
    )
    batch.add_argument(
        "--families", metavar="FILE",
        help=("the families to draw, one per line ('-' reads stdin). Only the first "
              "tab-separated column is read, so tools/select_families.py output works"),
    )
    batch.add_argument(
        "-j", "--jobs", type=int, default=1, metavar="N",
        help="families drawn at once, in separate processes (default 1)",
    )
    batch.add_argument(
        "--no-index", action="store_true",
        help="do not write index.html, the contact sheet over the batch",
    )

    output = parser.add_argument_group("output")
    output.add_argument("-o", "--output", metavar="DIR", default=".", help="output directory")
    output.add_argument(
        "--theme", choices=sorted(theme.THEMES), default="light", help="colour theme"
    )
    output.add_argument(
        "-a", "--alpha", type=float, default=None, metavar="A",
        help="opacity of partial copies in panel 1, 0-1 (default 0.65; v1's -a)",
    )
    output.add_argument(
        "-F", "--full-length-alpha", type=float, default=None, metavar="A",
        help="opacity of full-length copies in panel 1, 0-1 (default 0.9; v1's -F)",
    )
    output.add_argument(
        "-y", "--auto-y", type=float, default=None, metavar="MAX", dest="divergence_max",
        help="fix panel 1's divergence axis to 0-MAX instead of fitting it (v1's -y)",
    )
    output.add_argument(
        "-D", "--emboss-dotmatcher", action="store_true",
        help=("accepted so v1 scripts run, and ignored: panel 3, the self blastn "
              "dot-plot, replaces the dotmatcher plot"),
    )
    output.add_argument(
        "-t",
        "--tables",
        action="store_true",
        help=("also write the self dot-plot, ORF and protein-hit coordinates as TSV, "
              "1-based, as v1's -t did"),
    )
    output.add_argument(
        "-T",
        "--all-tables",
        action="store_true",
        help="as --tables, plus one row per genomic copy row, as v1's -T did",
    )
    output.add_argument(
        "--static",
        metavar="FMT",
        choices=["png", "pdf", "svg"],
        help="also write a static export (requires kaleido)",
    )

    tuning = parser.add_argument_group("tuning")
    tuning.add_argument(
        "-f",
        "--full-length-threshold",
        type=float,
        default=0.9,
        metavar="FRAC",
        help="fraction of the consensus a copy must cover to count as full length (default 0.9)",
    )
    tuning.add_argument(
        "--no-group-fragments",
        action="store_true",
        help="count every annotation row as a copy, instead of joining the rows "
             "that share an ID (RepeatMasker ID, BED16 column 16) into one insertion",
    )
    tuning.add_argument(
        "--word-size",
        type=int,
        default=7,
        help="blastn word size for the self dot-plot (default 7)",
    )
    tuning.add_argument(
        "--self-evalue",
        type=float,
        default=1e-3,
        help="blastn e-value for the self dot-plot (default 1e-3)",
    )
    tuning.add_argument(
        "--no-dotplot", action="store_true", help="skip the self dot-plot (no blastn call)"
    )
    tuning.add_argument(
        "-m",
        "--min-orf",
        type=int,
        default=400,
        metavar="BP",
        help="minimum ORF length in nucleotides (default 400, as in v1)",
    )
    tuning.add_argument(
        "-e",
        "--e-value",
        dest="genome_evalue",
        type=float,
        default=genome.DEFAULT_EVALUE,
        metavar="E",
        help="genome route: blastn E-value for keeping a hit (default 10e-8, as v1)",
    )
    tuning.add_argument(
        "-r",
        "--remove-redundant",
        action="store_true",
        help=("genome route: drop hits whose locus lies inside another hit's, "
              "v1's reduce.cpp rule"),
    )
    tuning.add_argument(
        "-R",
        "--no-reverse-orfs",
        action="store_true",
        help="only report ORFs on the forward strand (getorf -noreverse)",
    )
    tuning.add_argument(
        "--no-orfs", action="store_true", help="skip the ORF track (no getorf call)"
    )

    homology_group = parser.add_argument_group("protein homology (panel 5)")
    homology_group.add_argument(
        "--no-homology", action="store_true",
        help="skip the protein row (no bathsearch call)",
    )
    homology_group.add_argument(
        "--min-model-coverage", type=float, default=0.0, metavar="FRAC",
        help=("leave out tier 1-2 hits spanning less than FRAC of their model (0-1). "
              "Off by default until calibrated; left-out hits stay in the hit table"),
    )
    homology_group.add_argument(
        "--cpus", type=int, default=None, metavar="N",
        help="threads for bathsearch (default: BATH's own, 2)",
    )
    homology_group.add_argument(
        "--proteins", metavar="FILE",
        help="a prebuilt BATH pHMM library to search instead of building tiers 1-2; "
             "tier 3 (blastp vs RepeatPeps) still runs",
    )
    homology_group.add_argument(
        "--repeatpeps", metavar="FILE",
        help="path to RepeatPeps.lib (ships with RepeatMasker). Needed for the "
             "superfamilies Pfam cannot model: piggyBac, Maverick, Crypton, Penelope",
    )
    homology_group.add_argument(
        "--deep", action="store_true",
        help="search all of RepeatPeps as pHMMs rather than only the Pfam-blind "
             "superfamilies. Builds a ~6.4 GB library on first use; see docs/BRIEF_v2.md §5.1a",
    )
    homology_group.add_argument(
        "--homology-evalue", type=float, default=1e-3, metavar="E",
        help="E-value cutoff for the protein search (default 1e-3)",
    )
    homology_group.add_argument(
        "--rebuild-proteins", action="store_true",
        help="rebuild the cached protein library even if it is present",
    )

    parser.add_argument("--version", action="version", version=f"teaid {__version__}")
    return parser


def _protein_homology(args, consensus, found_orfs, notes, sources,
                      prep: _Prepared) -> tuple[list, list]:
    """Protein hits for panel 5, degrading to an empty row rather than failing.

    A missing BATH install or an unbuildable library costs the sheet one panel;
    it must not cost it the whole run, since the other four panels are
    independent of it. Returns (hits kept for drawing, every hit found), and
    appends to ``notes`` and ``sources`` so the sheet says which library was
    searched and, when tier 3 did not run, why.
    """
    from pathlib import Path as _Path

    if prep.protein_error:
        notes.append(f"protein row skipped: {prep.protein_error}")
        return [], []
    try:
        if prep.protein_library is not None:
            library = prep.protein_library
        elif args.proteins:
            # A prebuilt pHMM library replaces tiers 1-2 only. Tier 3 still
            # follows --repeatpeps and the usual discovery, as in any other run.
            library = proteins.Library(
                hmm=_Path(args.proteins),
                repeatpeps=proteins.find_repeatpeps(args.repeatpeps),
                tiers=(f"prebuilt {_Path(args.proteins).name}",),
            )
        else:
            library = proteins.build(
                deep=args.deep,
                repeatpeps=args.repeatpeps,
                force=args.rebuild_proteins,
                log=lambda m: print(f"teaid: {m}", file=sys.stderr),
            )
    except proteins.LibraryError as exc:
        print(f"teaid: warning: protein row skipped: {exc}", file=sys.stderr)
        notes.append(f"protein row skipped: {exc}")
        return [], []

    try:
        hits = homology.search(
            consensus.sequence,
            library,
            name=consensus.bare_name,
            evalue=args.homology_evalue,
            cpus=args.cpus,
        )
    except (homology.BathNotFound, RuntimeError) as exc:
        print(f"teaid: warning: protein row skipped: {exc}", file=sys.stderr)
        notes.append(f"protein row skipped: {exc}")
        return [], []
    sources.append(f"proteins <b>{', '.join(library.tiers)}</b> (bathsearch)")

    # Tier 3: named-element evidence from RepeatPeps, via the ORF peptides.
    # Answers a different question from tiers 1-2 -- which named element this
    # resembles, rather than which domain it encodes -- and costs 17 MB.
    if args.deep:
        notes.append("tier 3 (blastp vs RepeatPeps) off: --deep already searches "
                     "all of RepeatPeps")
    elif args.no_orfs:
        notes.append("tier 3 (blastp vs RepeatPeps) off: it searches the ORF "
                     "peptides, and --no-orfs turned those off")
    elif library.repeatpeps is None:
        notes.append("tier 3 (blastp vs RepeatPeps) off: RepeatPeps.lib not found")
    else:
        hits += homology.search_repeatpeps(found_orfs, library.repeatpeps,
                                           target=consensus.bare_name)
        sources.append(f"<b>{library.repeatpeps.name}</b> (blastp)")

    covered, short = homology.apply_coverage_floor(hits, args.min_model_coverage)
    if short:
        notes.append(f"{len(short)} protein hit{'s' if len(short) != 1 else ''} spanning "
                     f"under {args.min_model_coverage:.0%} of their model left out "
                     f"(--min-model-coverage); they are in the table")
    kept = homology.best_per_region(covered)
    dropped = len(covered) - len(kept)
    if dropped:
        notes.append(f"{dropped} competing protein hit{'s' if dropped != 1 else ''} "
                     f"collapsed — all {len(hits)} are listed under the sheet")
    return kept, hits


def _grouping_note(rows: Annotation, insertions: list, disabled: bool) -> str | None:
    """Say on the sheet how rows became copies, so the count can be checked.

    Naming the largest group matters: a producer that writes one constant into
    the ID column would otherwise collapse a family to a single 'copy' silently.
    """
    if disabled:
        return f"fragment grouping off: each of {len(rows):,} annotation rows counted as a copy"
    if len(insertions) == len(rows):
        return None
    largest = max(i.n_fragments for i in insertions)
    return (f"{len(rows):,} annotation rows grouped into "
            f"{report._count(len(insertions), 'copy', 'copies')} by their shared ID "
            f"(largest: {largest} rows)")


def _seed_qc_fields(seed, enabled: bool) -> dict:
    """Seed-QC inputs for the sheet, or nothing when the panels are off.

    Flag-gated and off by default (see docs/BRIEF_v2.md §2): the default sheet
    is the same four quadrants whether or not the input happened to be a seed.
    """
    if not enabled or seed is None:
        return {}
    import numpy as np

    depth = np.asarray(seed.depth(), dtype=np.int64)
    return {
        "seed_depth": depth,
        # The same positions counted the way panel 2 counts them, so panel 6 can
        # show both and make the difference between 'spans' and 'has a base'
        # visible rather than leaving it to be inferred across two panels.
        # Seed sequences carry no hit_id, so grouping leaves them one per row.
        "seed_span": analysis.coverage(seed.to_annotation().insertions(), len(depth)),
        "seed_mismatches": np.asarray(seed.mismatches(), dtype=np.int64),
        "seed_blocks": seed.aligned_blocks(),
        "expected_class": seed.expected_class,
        "seed_sequence_count": len(seed),
    }


@dataclass(slots=True)
class _Prepared:
    """Inputs already read, handed to one family's run.

    Empty for a single family, which then reads everything itself. Batch mode
    reads each input once, and gives every family its own slice of it.
    """

    annotation: Annotation | None = None  # this family's rows, or the whole file
    library: object | None = None  # a ConsensusLibrary
    seed: object | None = None
    genome_db: Path | None = None
    protein_library: object | None = None
    protein_error: str | None = None
    stem: str | None = None
    # Filled in by the run, for the batch manifest and index.
    facts: dict = dataclass_field(default_factory=dict)


@dataclass(slots=True)
class _Loaded:
    """One family, resolved from whichever input source was chosen."""

    family: str
    consensus: Consensus
    annotation: Annotation
    notes: list[str] = dataclass_field(default_factory=list)
    sources: list[str] = dataclass_field(default_factory=list)
    seed: object | None = None
    # All annotated genomic copies of this family, when an annotation was
    # supplied alongside a seed. Drives panels 1-2 so they show the family as it
    # exists in the genome, leaving panel 6 to show the seed's own sampling.
    genomic: Annotation | None = None


def _load_genomic_context(args, family: str, prep: _Prepared) -> Annotation | int | None:
    """Annotated genomic copies of ``family``, for seed QC against the genome.

    Returns None — not a failure — when the annotation simply has no rows for
    this family. A seed whose family is absent from the annotation is itself
    worth knowing about, but it is not a reason to refuse the sheet.
    """
    annot_path = Path(args.annot)
    annotation = prep.annotation
    if annotation is None:
        annotation = _read_annotation(args)
        if isinstance(annotation, int):
            return annotation

    subset = annotation.for_family(family, exact=False)
    if not len(subset):
        print(
            f"teaid: warning: no copies of {family!r} in {annot_path}; "
            f"panels 1-2 fall back to the seed's own sequences",
            file=sys.stderr,
        )
        return None
    return subset


def _read_annotation(args) -> Annotation | int:
    annot_path = Path(args.annot)
    if not annot_path.exists():
        return fail(EXIT_NO_INPUT, f"annotation not found: {annot_path}")
    kwargs = {}
    if args.divergence_kind:
        kwargs["divergence_kind"] = DivergenceKind(args.divergence_kind)
    try:
        return readers.read(annot_path, fmt=args.annot_format, **kwargs)
    except readers.FormatDetectionError as exc:
        return fail(EXIT_NO_INPUT, str(exc))


def _read_library(args, route: str):
    if not args.consensus:
        return fail(EXIT_USAGE, f"--consensus is required with {route}")
    consensus_path = Path(args.consensus)
    if not consensus_path.exists():
        return fail(EXIT_NO_INPUT, f"consensus FASTA not found: {consensus_path}")
    library = read_fasta(consensus_path)
    if not len(library):
        return fail(EXIT_NO_INPUT, f"no sequences in {consensus_path}")
    return library


def _resolve_consensus(args, route: str, prep: _Prepared) -> tuple[Path, Consensus] | int:
    """The one consensus a route draws, from ``--consensus`` and ``--family``."""
    library = prep.library if prep.library is not None else _read_library(args, route)
    if isinstance(library, int):
        return library
    consensus_path = Path(args.consensus)
    family = args.family
    if family is None:
        if len(library) != 1:
            return fail(
                EXIT_USAGE,
                f"--family is required: {consensus_path} holds {len(library):,} sequences",
            )
        family = next(iter(library)).bare_name

    consensus = library.get(family)
    if consensus is None:
        near = library.suggest(family)
        return fail(
            EXIT_NO_FAMILY,
            f"family {family!r} not in {consensus_path}",
            f"did you mean: {', '.join(near)}" if near else None,
        )
    return consensus_path, consensus


def _load_annotation(args, prep: _Prepared) -> _Loaded | int:
    """The ``--annot`` path: an annotation plus a separate consensus FASTA."""
    annot_path = Path(args.annot)
    if not args.consensus:
        return fail(EXIT_USAGE, "--consensus is required with --annot")
    if prep.annotation is None and not annot_path.exists():
        return fail(EXIT_NO_INPUT, f"annotation not found: {annot_path}")
    resolved = _resolve_consensus(args, "--annot", prep)
    if isinstance(resolved, int):
        return resolved
    consensus_path, consensus = resolved
    family = consensus.bare_name

    annotation = prep.annotation
    if annotation is None:
        annotation = _read_annotation(args)
        if isinstance(annotation, int):
            return annotation

    subset = annotation.for_family(family, exact=False)
    if not len(subset):
        return fail(
            EXIT_NO_EVIDENCE,
            f"no copies of {family!r} in {annot_path} "
            f"({len(annotation):,} rows, {len(annotation.families()):,} families)",
        )

    notes = []
    if annotation.skipped:
        notes.append(f"{len(annotation.skipped):,} unparsable rows skipped")
    sources = [
        f"annotation <b>{annot_path.name}</b> ({annotation.source_format})",
        f"consensus <b>{consensus_path.name}</b>",
    ]
    return _Loaded(consensus.bare_name, consensus, subset, notes, sources)


def _genome_copies(args, consensus: Consensus, prep: _Prepared) -> Annotation | int:
    """Copies of ``consensus`` found by blastn against ``--genome`` (the v1 route)."""
    genome_path = Path(args.genome)
    if not genome_path.exists():
        return fail(EXIT_NO_INPUT, f"genome not found: {genome_path}")
    try:
        database = prep.genome_db or genome.ensure_database(
            genome_path, log=lambda m: print(f"teaid: {m}", file=sys.stderr))
        copies = genome.search(consensus.sequence, consensus.bare_name, database,
                               evalue=args.genome_evalue, family=consensus.bare_name,
                               class_label=consensus.class_label)
    except (genome.BlastnNotFound, RuntimeError) as exc:
        # Without blastn this route has no copies at all, so it is an input
        # failure rather than a dropped panel.
        return fail(EXIT_NO_INPUT, str(exc))
    if args.remove_redundant:
        copies = genome.remove_redundant(copies)
    return copies


def _load_genome(args, prep: _Prepared) -> _Loaded | int:
    """The ``--genome`` path: v1's route, a consensus searched against a genome."""
    resolved = _resolve_consensus(args, "--genome", prep)
    if isinstance(resolved, int):
        return resolved
    consensus_path, consensus = resolved
    copies = _genome_copies(args, consensus, prep)
    if isinstance(copies, int):
        return copies
    if not len(copies):
        return fail(
            EXIT_NO_EVIDENCE,
            f"blastn found no copies of {consensus.bare_name!r} in {args.genome} "
            f"at E <= {args.genome_evalue:g}",
        )
    notes = [
        "copies found by blastn (megablast) against the genome, one per HSP: this "
        "route misses diverged copies and splits interrupted ones, so an annotation "
        "gives a truer count"
    ]
    if args.remove_redundant:
        notes.append("-r: hits nested inside another hit's locus removed, as in v1")
    sources = [
        f"genome <b>{Path(args.genome).name}</b> (blastn, E ≤ {args.genome_evalue:g})",
        f"consensus <b>{consensus_path.name}</b>",
    ]
    return _Loaded(consensus.bare_name, consensus, copies, notes, sources)


def _load_seed(args, prep: _Prepared) -> _Loaded | int:
    """The ``--stk`` path: one Stockholm record supplies everything.

    Copies, loci, consensus and alignment all come from the seed, so no separate
    FASTA is needed and none is consulted — a consensus passed alongside would
    silently disagree with the alignment the depth panel is drawn from.
    """
    seed_path = Path(args.stk)
    if prep.seed is None and not seed_path.exists():
        return fail(EXIT_NO_INPUT, f"seed alignment not found: {seed_path}")

    try:
        seed = prep.seed if prep.seed is not None else stockholm.read(seed_path, args.family)
    except stockholm.MalformedSeedError as exc:
        hint = None
        if args.family:
            try:
                available = stockholm.identifiers(seed_path)
                hint = f"{len(available):,} records, e.g. {', '.join(available[:4])}"
            except stockholm.MalformedSeedError:
                pass
        return fail(EXIT_BAD_SEED, str(exc), hint)

    sequence = seed.consensus()
    if not sequence:
        return fail(
            EXIT_NO_EVIDENCE,
            f"seed {seed.identifier!r} yields an empty consensus "
            f"(no #=GC RF line and no aligned columns)",
        )

    annotation = seed.to_annotation()
    if not len(annotation):
        return fail(
            EXIT_NO_EVIDENCE, f"seed {seed.identifier!r} holds no aligned sequences"
        )

    notes = []
    if args.consensus:
        notes.append("--consensus ignored: the seed supplies its own consensus")
        print(
            "teaid: warning: --consensus ignored; the seed supplies its own consensus",
            file=sys.stderr,
        )
    declared = seed.declared_count
    if declared is not None and declared != len(seed):
        notes.append(f"#=GF SQ declares {declared} sequences but {len(seed)} are present")

    name = seed.identifier or "seed"
    sources = [f"seed <b>{seed_path.name}</b> (Stockholm, {len(seed)} sequences)"]
    return _Loaded(
        name, Consensus(name, sequence), annotation, notes, sources, seed=seed
    )


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not 0 < args.full_length_threshold <= 1:
        parser.error("--full-length-threshold must be in (0, 1]")
    if not 0 <= args.min_model_coverage <= 1:
        parser.error("--min-model-coverage must be in [0, 1]")
    for flag, value in (("--alpha", args.alpha), ("--full-length-alpha", args.full_length_alpha)):
        if value is not None and not 0 <= value <= 1:
            parser.error(f"{flag} must be in [0, 1]")
    if args.divergence_max is not None and args.divergence_max <= 0:
        parser.error("--auto-y must be positive")

    # Input priority: --annot, then --genome (v1's blastn route), then --stk;
    # --pipeline reverses to seed-first. Whichever source is chosen, everything
    # downstream works from the same record struct.
    order = ("stk", "genome", "annot") if args.pipeline else ("annot", "genome", "stk")
    chosen = next((name for name in order if getattr(args, name)), None)
    if chosen is None:
        parser.error("give an input: --annot, --genome (-g) or --stk")

    # The priority above settles which source *drives* the sheet when only one
    # can. A seed given together with genomic copies is not that contest: the
    # seed is the thing being QC'd, so it supplies the consensus and the seed-QC
    # panels while the annotation or genome supplies context for panels 1-2.
    if args.stk and (args.annot or args.genome):
        chosen = "stk"

    if args.all or args.families:
        return _batch(args, chosen)
    return _draw(args, chosen, _Prepared())


def _draw(args, chosen: str, prep: _Prepared) -> int:
    """Everything after route selection, for one family."""
    loaders = {"stk": _load_seed, "annot": _load_annotation, "genome": _load_genome}
    loaded = loaders[chosen](args, prep)
    if isinstance(loaded, int):
        return loaded
    if chosen == "annot" and args.genome:
        loaded.notes.append("--genome ignored: an annotation takes priority over a search")

    # --stk with genomic copies: the seed says which copies were *chosen* to
    # build the family; the annotation (or a genome search) says which copies
    # *exist*. Showing both is the point of seed QC: a seed built from 12 of 340
    # copies may be perfectly good or may have sampled one corner of the family,
    # and only the comparison distinguishes them.
    if chosen == "stk" and (args.annot or args.genome):
        if args.annot:
            genomic = _load_genomic_context(args, loaded.family, prep)
            label = f"annotation <b>{Path(args.annot).name}</b>"
            kind = "annotated genomic copies"
        else:
            genomic = _genome_copies(args, loaded.consensus, prep)
            label = f"genome <b>{Path(args.genome).name}</b>"
            kind = "genomic copies found by blastn"
            if not isinstance(genomic, int) and not len(genomic):
                genomic = None
        if isinstance(genomic, int):
            return genomic
        if genomic is not None:
            loaded.genomic = genomic
            loaded.sources.append(f"{label} ({genomic.source_format})")
            copies = len(genomic.insertions(group=not args.no_group_fragments))
            loaded.notes.append(
                f"panels 1-2 show {copies:,} {kind}; "
                f"the seed uses {len(loaded.annotation):,}"
            )

    family = loaded.family
    consensus = loaded.consensus
    subset = loaded.annotation
    notes = loaded.notes
    seed = loaded.seed
    consensus_length = len(consensus)
    # Panels 1-2 describe the family in the genome. With an annotation supplied
    # alongside a seed that is every annotated copy; otherwise it is whatever
    # the single source gave us.
    shown = loaded.genomic if loaded.genomic is not None else subset
    # The one place the grouping is decided; everything downstream takes these.
    insertions = shown.insertions(group=not args.no_group_fragments)
    grouping = _grouping_note(shown, insertions, args.no_group_fragments)
    if grouping:
        notes.append(grouping)
    if args.emboss_dotmatcher:
        notes.append("-D ignored: the EMBOSS dotmatcher plot is retired; panel 3 is the "
                     "self dot-plot")
        print("teaid: warning: -D ignored; panel 3, the self dot-plot, replaces dotmatcher",
              file=sys.stderr)
    coverage = analysis.coverage(insertions, consensus_length)
    full = analysis.full_length(insertions, consensus_length, args.full_length_threshold)

    self_hits = []
    terminal: dict[str, list] = {}
    if not args.no_dotplot:
        try:
            self_hits = analysis.self_blast(
                consensus.sequence,
                name=consensus.bare_name,
                word_size=args.word_size,
                evalue=args.self_evalue,
            )
            terminal = analysis.terminal_repeats(self_hits, consensus_length)
        except (analysis.BlastNotFound, RuntimeError) as exc:
            notes.append(f"dot-plot skipped: {exc}")
            print(f"teaid: warning: {exc}", file=sys.stderr)

    found_orfs = []
    if not args.no_orfs:
        try:
            found_orfs = orfs.find_orfs(
                consensus.sequence,
                name=consensus.bare_name,
                min_size=args.min_orf,
                reverse=not args.no_reverse_orfs,
            )
        except (orfs.GetorfNotFound, RuntimeError) as exc:
            notes.append(f"ORF track skipped: {exc}")
            print(f"teaid: warning: {exc}", file=sys.stderr)

    protein_hits: list = []
    all_hits: list = []
    class_check = None
    if not args.no_homology:
        protein_hits, all_hits = _protein_homology(
            args, consensus, found_orfs, notes, loaded.sources, prep
        )
        if protein_hits and seed is not None and args.seed_qc:
            class_check = classcheck.check(seed.expected_class, protein_hits)

    data = report.SheetData(
        family=consensus.bare_name,
        consensus_length=consensus_length,
        insertions=insertions,
        full_length=full,
        coverage=coverage,
        self_hits=self_hits,
        terminal_repeats=terminal,
        orfs=found_orfs,
        orf_min_size=args.min_orf,
        class_label=consensus.class_label,
        full_length_threshold=args.full_length_threshold,
        source_format=shown.source_format,
        sources=loaded.sources,
        notes=notes,
        homology=protein_hits,
        homology_all=all_hits if not args.no_homology else [],
        class_check=class_check,
        partial_opacity=0.65 if args.alpha is None else args.alpha,
        full_opacity=0.9 if args.full_length_alpha is None else args.full_length_alpha,
        divergence_max=args.divergence_max,
        **_seed_qc_fields(seed, args.seed_qc),
    )

    figure = report.build_figure(data, theme.get(args.theme))

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = prep.stem or _safe_stem(consensus.bare_name)
    html_path = out_dir / f"{stem}.teaid.html"
    report.write_html(figure, html_path, data, theme.get(args.theme))
    print(f"wrote {html_path}")

    if args.tables or args.all_tables:
        for path in tables.write(data, out_dir, stem, copies=args.all_tables):
            print(f"wrote {path}")

    if args.static:
        static_path = out_dir / f"{stem}.teaid.{args.static}"
        try:
            report.export_static(figure, static_path)
            print(f"wrote {static_path}")
        except TimeoutError as exc:
            print(f"teaid: warning: {exc}; the HTML sheet was written", file=sys.stderr)
        except RuntimeError as exc:
            hint = ("; install kaleido with 'pip install teaid[static]'"
                    if "kaleido" in str(exc).lower() else "")
            print(f"teaid: warning: static export failed ({exc}){hint}", file=sys.stderr)

    rows = f" ({len(shown):,} annotation rows)" if len(shown) != len(insertions) else ""
    summary = (
        f"{consensus.bare_name}: {report._count(len(insertions), 'copy', 'copies')}{rows}, "
        f"{len(full):,} full length, consensus {consensus_length:,} bp"
    )
    if terminal.get("LTR") or terminal.get("TIR"):
        summary += (
            f", terminal-repeat candidates: {len(terminal.get('LTR', []))} LTR-like, "
            f"{len(terminal.get('TIR', []))} TIR-like"
        )
    if found_orfs:
        longest = max(o.length_aa for o in found_orfs)
        summary += f", {len(found_orfs)} ORFs (longest {longest:,} aa)"
    print(summary)
    prep.facts.update(
        copies=len(insertions), rows=len(shown), full_length=len(full),
        class_label=consensus.class_label or "", length=consensus_length,
    )
    return EXIT_OK


def _family_list(args, names: list[str]) -> list[str] | int:
    """The families a batch draws: all of ``names``, or those named in --families."""
    if args.all:
        return names
    source = sys.stdin if args.families == "-" else None
    if source is None:
        path = Path(args.families)
        if not path.exists():
            return fail(EXIT_NO_INPUT, f"family list not found: {path}")
        text = path.read_text()
    else:
        text = source.read()
    wanted = []
    for line in text.splitlines():
        name = line.split("\t", 1)[0].strip()
        if name and not name.startswith("#"):
            wanted.append(name)
    return wanted


def _unique_stems(names: list[str]) -> dict[str, str]:
    """A file stem per family, suffixed where two names would share one.

    One family per invocation could leave collisions to the caller. Once TE-Aid
    runs the loop it has to handle them, or two sheets silently become one.
    """
    stems: dict[str, str] = {}
    seen: dict[str, int] = {}
    for name in names:
        base = _safe_stem(name.split("#", 1)[0])
        seen[base] = seen.get(base, 0) + 1
        stems[name] = base if seen[base] == 1 else f"{base}.{seen[base]}"
    return stems


def _batch_one(args, chosen: str, family: str, prep: _Prepared) -> dict:
    """Draw one family of a batch, capturing its output. Runs in a worker."""
    import contextlib
    import io
    import time
    import traceback

    args = argparse.Namespace(**{**vars(args), "family": family, "all": False, "families": None})
    out, err = io.StringIO(), io.StringIO()
    start = time.perf_counter()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = _draw(args, chosen, prep)
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
        except Exception:  # one family's crash must not end the batch
            traceback.print_exc()
            code = 1
    return dict(family=family, stem=prep.stem, code=code,
                seconds=time.perf_counter() - start,
                stdout=out.getvalue(), stderr=err.getvalue(), facts=prep.facts)


_MANIFEST = ["family", "stem", "exit", "slug", "seconds", "copies", "annotation_rows",
             "full_length", "class", "consensus_length", "warnings", "error"]


def _batch(args, chosen: str) -> int:
    """Draw many families in one run, reading each input once.

    Writes each sheet as usual, plus ``results.tsv`` (one row per family) and
    ``index.html``. Exits 0 only when every family drew; otherwise exit 8,
    ``batch-incomplete``, with the reasons in the manifest.
    """
    import concurrent.futures

    from . import index as index_page
    from .records import _bare_name

    if args.family:
        return fail(EXIT_USAGE, "--family draws one family; drop it to use --all or --families")
    if args.jobs < 1:
        return fail(EXIT_USAGE, "--jobs must be at least 1")
    log = lambda message: print(f"teaid: {message}", file=sys.stderr)  # noqa: E731
    base = _Prepared()
    skipped_note = None

    if chosen == "stk":
        seed_path = Path(args.stk)
        if not seed_path.exists():
            return fail(EXIT_NO_INPUT, f"seed alignment not found: {seed_path}")
        try:
            seeds = {s.identifier or f"seed-{i}": s
                     for i, s in enumerate(stockholm.read_all(seed_path), start=1)}
        except stockholm.MalformedSeedError as exc:
            return fail(EXIT_BAD_SEED, str(exc))
        names = list(seeds)
    else:
        library = _read_library(args, f"--{chosen}")
        if isinstance(library, int):
            return library
        base.library = library
        names = [entry.bare_name for entry in library]

    by_family: dict[str, list] = {}
    annotation = None
    if args.annot:
        annotation = _read_annotation(args)
        if isinstance(annotation, int):
            return annotation
        for copy in annotation.copies:
            by_family.setdefault(_bare_name(copy.family), []).append(copy)
        if chosen == "annot" and args.all:
            absent = [n for n in names if _bare_name(n) not in by_family]
            names = [n for n in names if _bare_name(n) in by_family]
            if absent:
                skipped_note = (f"{len(absent):,} library families have no annotated copies "
                                f"and were not drawn")
    if args.genome and chosen != "annot":
        genome_path = Path(args.genome)
        if not genome_path.exists():
            return fail(EXIT_NO_INPUT, f"genome not found: {genome_path}")
        try:
            base.genome_db = genome.ensure_database(genome_path, log=log)
        except (genome.BlastnNotFound, RuntimeError) as exc:
            return fail(EXIT_NO_INPUT, str(exc))
    if not args.no_homology and not args.proteins:
        # Built once here, so no worker waits on the lock with nothing to show.
        try:
            base.protein_library = proteins.build(deep=args.deep, repeatpeps=args.repeatpeps,
                                                  force=args.rebuild_proteins, log=log)
        except proteins.LibraryError as exc:
            base.protein_error = str(exc)
            log(f"warning: protein row skipped for every family: {exc}")

    families = _family_list(args, names)
    if isinstance(families, int):
        return families
    if not families:
        return fail(EXIT_NO_FAMILY, "no families to draw")
    stems = _unique_stems(families)
    renamed = [f for f in families if stems[f] != _safe_stem(f.split("#", 1)[0])]
    for family in renamed:
        log(f"warning: {family!r} is written as {stems[family]}.teaid.html; "
            f"its name collides with another family's")

    def prep_for(family: str) -> _Prepared:
        prep = _Prepared(library=base.library, genome_db=base.genome_db,
                         protein_library=base.protein_library,
                         protein_error=base.protein_error, stem=stems[family])
        if chosen == "stk":
            prep.seed = seeds.get(family)
        if annotation is not None:
            # Only this family's rows travel to the worker, not the whole file.
            prep.annotation = type(annotation)(
                copies=by_family.get(_bare_name(family), []),
                source_path=annotation.source_path,
                source_format=annotation.source_format,
                skipped=annotation.skipped,
            )
        return prep

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    total = len(families)

    def report_one(result: dict) -> None:
        results.append(result)
        sys.stdout.write(result["stdout"])
        slug = _SLUG.get(result["code"], "error" if result["code"] else "")
        status = "ok" if result["code"] == EXIT_OK else f"failed [{slug}]"
        log(f"[{len(results)}/{total}] {result['family']}: {status} "
            f"({result['seconds']:.0f} s)")

    if args.jobs == 1:
        for family in families:
            if chosen == "stk" and family not in seeds:
                report_one(_missing(family, stems[family], "not a record in the seed"))
                continue
            report_one(_batch_one(args, chosen, family, prep_for(family)))
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = []
            for family in families:
                if chosen == "stk" and family not in seeds:
                    report_one(_missing(family, stems[family], "not a record in the seed"))
                    continue
                futures.append(pool.submit(_batch_one, args, chosen, family, prep_for(family)))
            for future in concurrent.futures.as_completed(futures):
                report_one(future.result())

    order = {family: i for i, family in enumerate(families)}
    results.sort(key=lambda r: order[r["family"]])
    manifest = out_dir / "results.tsv"
    with manifest.open("w") as handle:
        handle.write("\t".join(_MANIFEST) + "\n")
        for r in results:
            handle.write("\t".join(str(v) for v in _manifest_row(r)) + "\n")
    print(f"wrote {manifest}")
    if not args.no_index:
        page = index_page.write(out_dir, results, title=_batch_title(args, chosen),
                                note=skipped_note)
        print(f"wrote {page}")

    failed = [r for r in results if r["code"] != EXIT_OK]
    summary = f"{total - len(failed):,} of {total:,} families drawn"
    if skipped_note:
        summary += f"; {skipped_note}"
    if failed:
        log(f"{summary}; {len(failed):,} failed, see {manifest}")
        return fail(EXIT_BATCH_INCOMPLETE, f"{len(failed):,} of {total:,} families failed; "
                                           f"see {manifest}")
    log(summary)
    return EXIT_OK


def _missing(family: str, stem: str, why: str) -> dict:
    line = f"teaid: error [{_SLUG[EXIT_NO_FAMILY]}]: {family!r} {why}\n"
    return dict(family=family, stem=stem, code=EXIT_NO_FAMILY, seconds=0.0,
                stdout="", stderr=line, facts={})


def _manifest_row(r: dict) -> list:
    facts = r.get("facts") or {}
    lines = [ln.strip() for ln in r["stderr"].splitlines() if ln.strip()]
    warnings = " | ".join(ln.split("warning:", 1)[1].strip() for ln in lines if "warning:" in ln)
    errors = " | ".join(ln for ln in lines if "error" in ln and "warning:" not in ln)
    slug = "" if r["code"] == EXIT_OK else _SLUG.get(r["code"], "error")
    return [r["family"], r["stem"], r["code"], slug, f"{r['seconds']:.1f}",
            facts.get("copies", ""), facts.get("rows", ""), facts.get("full_length", ""),
            facts.get("class_label", ""), facts.get("length", ""),
            warnings.replace("\t", " "), errors.replace("\t", " ")[:500]]


def _batch_title(args, chosen: str) -> str:
    source = {"stk": args.stk, "annot": args.annot, "genome": args.genome}[chosen]
    return f"TE-Aid sheets: {Path(source).name}"


def main_legacy_alias(argv: list[str] | None = None) -> int:
    """Entry point for the historical ``TE-Aid`` command name."""
    # Say exactly what carries over from v1, rather than promising it all.
    print(
        "note: 'TE-Aid' is now 'teaid', and this alias runs it. Every v1 flag is "
        "accepted with its v1 meaning; -D is ignored, because panel 3 replaces the "
        "dotmatcher plot. The family name, which v1 did not need, is -n.",
        file=sys.stderr,
    )
    return main(argv)


def _safe_stem(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-._" else "_" for ch in name)


if __name__ == "__main__":
    raise SystemExit(main())
