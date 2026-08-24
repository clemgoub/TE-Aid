"""Command line entry point.

Standalone input priority is ``--annot`` (an annotation you already trust),
then ``--blastn`` (the v1 rediscovery path), then ``--stk``. Only ``--annot`` is
implemented at this stage.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__, analysis, readers, report, theme
from .records import DivergenceKind
from .sequences import read_fasta

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_NO_INPUT = 3  # a required file is missing or unreadable
EXIT_NO_FAMILY = 4  # the requested family is not in the annotation or library
EXIT_NO_EVIDENCE = 5  # the family exists but carries nothing plottable


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
        "-c",
        "--consensus",
        metavar="FASTA",
        required=True,
        help="consensus FASTA; a library of many families or a single sequence",
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

    parser.add_argument("--version", action="version", version=f"teaid {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.annot:
        parser.error("--annot is required (the --blastn and --stk paths are not yet implemented)")

    annot_path, consensus_path = Path(args.annot), Path(args.consensus)
    for path, label in ((annot_path, "annotation"), (consensus_path, "consensus FASTA")):
        if not path.exists():
            print(f"teaid: {label} not found: {path}", file=sys.stderr)
            return EXIT_NO_INPUT

    if not 0 < args.full_length_threshold <= 1:
        parser.error("--full-length-threshold must be in (0, 1]")

    library = read_fasta(consensus_path)
    if not len(library):
        print(f"teaid: no sequences in {consensus_path}", file=sys.stderr)
        return EXIT_NO_INPUT

    family = args.family
    if family is None:
        if len(library) != 1:
            print(
                f"teaid: --family is required: {consensus_path} holds "
                f"{len(library):,} sequences",
                file=sys.stderr,
            )
            return EXIT_USAGE
        family = next(iter(library)).bare_name

    consensus = library.get(family)
    if consensus is None:
        print(f"teaid: family {family!r} not in {consensus_path}", file=sys.stderr)
        near = library.suggest(family)
        if near:
            print(f"       did you mean: {', '.join(near)}", file=sys.stderr)
        return EXIT_NO_FAMILY

    kwargs = {}
    if args.divergence_kind:
        kwargs["divergence_kind"] = DivergenceKind(args.divergence_kind)
    try:
        annotation = readers.read(annot_path, fmt=args.annot_format, **kwargs)
    except readers.FormatDetectionError as exc:
        print(f"teaid: {exc}", file=sys.stderr)
        return EXIT_NO_INPUT

    subset = annotation.for_family(family, exact=False)
    notes: list[str] = []
    if annotation.skipped:
        notes.append(f"{len(annotation.skipped):,} unparsable rows skipped")

    if not len(subset):
        print(
            f"teaid: no copies of {family!r} in {annot_path} "
            f"({len(annotation):,} rows, {len(annotation.families()):,} families)",
            file=sys.stderr,
        )
        return EXIT_NO_EVIDENCE

    consensus_length = len(consensus)
    coverage = analysis.coverage(subset, consensus_length)
    full = analysis.full_length(subset, consensus_length, args.full_length_threshold)

    self_hits = []
    if not args.no_dotplot:
        try:
            self_hits = analysis.self_blast(
                consensus.sequence,
                name=consensus.bare_name,
                word_size=args.word_size,
                evalue=args.self_evalue,
            )
        except (analysis.BlastNotFound, RuntimeError) as exc:
            notes.append(f"dot-plot skipped: {exc}")
            print(f"teaid: warning: {exc}", file=sys.stderr)

    data = report.SheetData(
        family=consensus.bare_name,
        consensus_length=consensus_length,
        copies=list(subset),
        full_length=full,
        coverage=coverage,
        self_hits=self_hits,
        class_label=consensus.class_label,
        full_length_threshold=args.full_length_threshold,
        source_format=annotation.source_format,
        notes=notes,
    )

    figure = report.build_figure(data, theme.get(args.theme))

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = _safe_stem(consensus.bare_name)
    html_path = out_dir / f"{stem}.teaid.html"
    figure.write_html(
        html_path,
        include_plotlyjs="cdn",
        full_html=True,
        config={"scrollZoom": True, "displaylogo": False, "responsive": True},
    )
    print(f"wrote {html_path}")

    if args.static:
        static_path = out_dir / f"{stem}.teaid.{args.static}"
        try:
            figure.write_image(static_path, scale=2)
            print(f"wrote {static_path}")
        except Exception as exc:  # kaleido raises a variety of types
            print(
                f"teaid: static export failed ({exc}); install kaleido with "
                f"'pip install teaid[static]'",
                file=sys.stderr,
            )

    summary = (
        f"{consensus.bare_name}: {len(subset):,} copies, {len(full):,} full length, "
        f"consensus {consensus_length:,} bp"
    )
    if self_hits:
        terminal = analysis.terminal_repeats(self_hits, consensus_length)
        if terminal["LTR"] or terminal["TIR"]:
            summary += (
                f", terminal-repeat candidates: {len(terminal['LTR'])} LTR-like, "
                f"{len(terminal['TIR'])} TIR-like"
            )
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
