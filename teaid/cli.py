"""Command line entry point.

Four input routes, designed to compose: ``--annot`` (an annotation you already
trust), ``--blastn`` (the v1 rediscovery path, not implemented yet), ``--stk``
(a Stockholm seed), and ``--annot`` **with** ``--stk`` — which is not a contest
between the two but the richest sheet, showing the family as it exists in the
genome against the copies a seed actually used.

Standalone priority among single routes is annot, blastn, stk; ``--pipeline``
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

from . import (__version__, analysis, classcheck, homology, orfs, proteins,
               readers, report, theme)
from .readers import stockholm
from .records import Annotation, DivergenceKind
from .sequences import Consensus, read_fasta

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_NO_INPUT = 3  # a required file is missing or unreadable
EXIT_NO_FAMILY = 4  # the requested family is not in the annotation or library
EXIT_NO_EVIDENCE = 5  # the family exists but carries nothing plottable
EXIT_BAD_SEED = 6  # the Stockholm file could not be parsed

# Exit code paired with the slug printed on stderr, so a caller can branch on
# either. Keep the two in step.
_SLUG = {
    EXIT_USAGE: "usage",
    EXIT_NO_INPUT: "no-input",
    EXIT_NO_FAMILY: "no-family",
    EXIT_NO_EVIDENCE: "no-evidence",
    EXIT_BAD_SEED: "bad-seed",
}


def fail(code: int, message: str, hint: str | None = None) -> int:
    """Report a failure in the form the pipeline parses, and return its code."""
    print(f"teaid: error [{_SLUG.get(code, 'error')}]: {message}", file=sys.stderr)
    if hint:
        print(f"       {hint}", file=sys.stderr)
    return code


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="teaid",
        description=(
            "Build an evidence sheet for one TE family, to support manual curation. "
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
        "-c",
        "--consensus",
        metavar="FASTA",
        help=(
            "consensus FASTA; a library of many families or a single sequence. "
            "Required with --annot; with --stk the seed supplies the consensus"
        ),
    )
    source.add_argument(
        "-f",
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

    output = parser.add_argument_group("output")
    output.add_argument("-o", "--output", metavar="DIR", default=".", help="output directory")
    output.add_argument(
        "--theme", choices=sorted(theme.THEMES), default="light", help="colour theme"
    )
    output.add_argument(
        "--static",
        metavar="FMT",
        choices=["png", "pdf", "svg"],
        help="also write a static export (requires kaleido)",
    )

    tuning = parser.add_argument_group("tuning")
    tuning.add_argument(
        "-t",
        "--full-length-threshold",
        type=float,
        default=0.9,
        metavar="FRAC",
        help="fraction of the consensus a copy must span to count as full length (default 0.9)",
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
        "--proteins", metavar="FILE",
        help="a prebuilt BATH pHMM library to search, instead of building one",
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


def _protein_homology(args, consensus, notes) -> tuple[list, str | None]:
    """Protein hits for panel 5, degrading to an empty row rather than failing.

    A missing BATH install or an unbuildable library costs the sheet one panel;
    it must not cost it the whole run, since the other four panels are
    independent of it.
    """
    from pathlib import Path as _Path

    try:
        if args.proteins:
            library = proteins.Library(hmm=_Path(args.proteins), repeatpeps=None,
                                       tiers=("prebuilt",))
        else:
            library = proteins.build(
                deep=args.deep,
                repeatpeps=args.repeatpeps,
                force=args.rebuild_proteins,
                log=lambda m: print(f"teaid: {m}", file=sys.stderr),
            )
    except proteins.LibraryError as exc:
        print(f"teaid: warning: protein row skipped: {exc}", file=sys.stderr)
        return [], f"protein row skipped: {exc}"

    try:
        hits = homology.search(
            consensus.sequence,
            library,
            name=consensus.bare_name,
            evalue=args.homology_evalue,
        )
    except (homology.BathNotFound, RuntimeError) as exc:
        print(f"teaid: warning: protein row skipped: {exc}", file=sys.stderr)
        return [], f"protein row skipped: {exc}"

    kept = homology.best_per_region(hits)
    dropped = len(hits) - len(kept)
    note = None
    if dropped:
        note = f"{dropped} competing protein hit{'s' if dropped != 1 else ''} collapsed"
    return kept, note


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
        "seed_span": analysis.coverage(seed.to_annotation(), len(depth)),
        "seed_mismatches": np.asarray(seed.mismatches(), dtype=np.int64),
        "seed_blocks": seed.aligned_blocks(),
        "expected_class": seed.expected_class,
        "seed_sequence_count": len(seed),
    }


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


def _load_genomic_context(args, family: str) -> Annotation | int | None:
    """Annotated genomic copies of ``family``, for seed QC against the genome.

    Returns None — not a failure — when the annotation simply has no rows for
    this family. A seed whose family is absent from the annotation is itself
    worth knowing about, but it is not a reason to refuse the sheet.
    """
    annot_path = Path(args.annot)
    if not annot_path.exists():
        return fail(EXIT_NO_INPUT, f"annotation not found: {annot_path}")

    kwargs = {}
    if args.divergence_kind:
        kwargs["divergence_kind"] = DivergenceKind(args.divergence_kind)
    try:
        annotation = readers.read(annot_path, fmt=args.annot_format, **kwargs)
    except readers.FormatDetectionError as exc:
        return fail(EXIT_NO_INPUT, str(exc))

    subset = annotation.for_family(family, exact=False)
    if not len(subset):
        print(
            f"teaid: warning: no copies of {family!r} in {annot_path}; "
            f"panels 1-2 fall back to the seed's own sequences",
            file=sys.stderr,
        )
        return None
    return subset


def _load_annotation(args) -> _Loaded | int:
    """The ``--annot`` path: an annotation plus a separate consensus FASTA."""
    if not args.consensus:
        return fail(EXIT_USAGE, "--consensus is required with --annot")

    annot_path, consensus_path = Path(args.annot), Path(args.consensus)
    for path, label in ((annot_path, "annotation"), (consensus_path, "consensus FASTA")):
        if not path.exists():
            return fail(EXIT_NO_INPUT, f"{label} not found: {path}")

    library = read_fasta(consensus_path)
    if not len(library):
        return fail(EXIT_NO_INPUT, f"no sequences in {consensus_path}")

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

    kwargs = {}
    if args.divergence_kind:
        kwargs["divergence_kind"] = DivergenceKind(args.divergence_kind)
    try:
        annotation = readers.read(annot_path, fmt=args.annot_format, **kwargs)
    except readers.FormatDetectionError as exc:
        return fail(EXIT_NO_INPUT, str(exc))

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


def _load_seed(args) -> _Loaded | int:
    """The ``--stk`` path: one Stockholm record supplies everything.

    Copies, loci, consensus and alignment all come from the seed, so no separate
    FASTA is needed and none is consulted — a consensus passed alongside would
    silently disagree with the alignment the depth panel is drawn from.
    """
    seed_path = Path(args.stk)
    if not seed_path.exists():
        return fail(EXIT_NO_INPUT, f"seed alignment not found: {seed_path}")

    try:
        seed = stockholm.read(seed_path, args.family)
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

    # Input priority: --annot, then --blastn, then --stk; --pipeline reverses to
    # seed-first. Whichever source is chosen, everything downstream works from
    # the same record struct.
    order = ("stk", "annot") if args.pipeline else ("annot", "stk")
    chosen = next((name for name in order if getattr(args, name)), None)
    if chosen is None:
        parser.error("give an input: --annot or --stk (--blastn is not implemented yet)")

    # The priority above settles which source *drives* the sheet when only one
    # can. A seed given together with an annotation is not that contest: the
    # seed is the thing being QC'd, so it supplies the consensus and the seed-QC
    # panels while the annotation supplies genomic context for panels 1-2.
    if args.stk and args.annot:
        chosen = "stk"

    loaded = (_load_seed if chosen == "stk" else _load_annotation)(args)
    if isinstance(loaded, int):
        return loaded

    # --stk and --annot are not alternatives when both are given. The seed says
    # which copies were *chosen* to build the family; the annotation says which
    # copies *exist* in the genome. Showing both is the point of seed QC: a seed
    # built from 12 of 340 copies may be perfectly good or may have sampled one
    # corner of the family, and only the comparison distinguishes them.
    if chosen == "stk" and args.annot:
        genomic = _load_genomic_context(args, loaded.family)
        if isinstance(genomic, int):
            return genomic
        if genomic is not None:
            loaded.genomic = genomic
            loaded.sources.append(
                f"annotation <b>{Path(args.annot).name}</b> "
                f"({genomic.source_format})"
            )
            loaded.notes.append(
                f"panels 1-2 show {len(genomic):,} annotated genomic copies; "
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
    coverage = analysis.coverage(shown, consensus_length)
    full = analysis.full_length(shown, consensus_length, args.full_length_threshold)

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

    protein_hits = []
    class_check = None
    if not args.no_homology:
        protein_hits, note = _protein_homology(args, consensus, notes)
        if note:
            notes.append(note)
        if protein_hits and seed is not None and args.seed_qc:
            class_check = classcheck.check(seed.expected_class, protein_hits)

    data = report.SheetData(
        family=consensus.bare_name,
        consensus_length=consensus_length,
        copies=list(shown),
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
        class_check=class_check,
        **_seed_qc_fields(seed, args.seed_qc),
    )

    figure = report.build_figure(data, theme.get(args.theme))

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = _safe_stem(consensus.bare_name)
    html_path = out_dir / f"{stem}.teaid.html"
    report.write_html(figure, html_path, data, theme.get(args.theme))
    print(f"wrote {html_path}")

    if args.static:
        static_path = out_dir / f"{stem}.teaid.{args.static}"
        try:
            report.strip_help_markers(figure).write_image(static_path, scale=2)
            print(f"wrote {static_path}")
        except Exception as exc:  # kaleido raises a variety of types
            print(
                f"teaid: static export failed ({exc}); install kaleido with "
                f"'pip install teaid[static]'",
                file=sys.stderr,
            )

    summary = (
        f"{consensus.bare_name}: {len(shown):,} copies, {len(full):,} full length, "
        f"consensus {consensus_length:,} bp"
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
    return EXIT_OK


def main_legacy_alias(argv: list[str] | None = None) -> int:
    """Entry point for the historical ``TE-Aid`` command name."""
    print(
        "note: 'TE-Aid' is now 'teaid'. This alias will keep working; "
        "please update your scripts.",
        file=sys.stderr,
    )
    return main(argv)


def _safe_stem(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-._" else "_" for ch in name)


if __name__ == "__main__":
    raise SystemExit(main())
