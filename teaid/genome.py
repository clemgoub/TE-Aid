"""The v1 input route: copies found by ``blastn`` of the consensus against a genome.

v1 took a consensus and a genome and nothing else, and this module keeps that
working (``-q/--query`` with ``-g/--genome``, or ``--blastn GENOME``). It is the
leakiest route v2 offers, and the sheet says so. Megablast misses diverged
copies. Each HSP is one copy, because a search has no ID column to join the
pieces of an interrupted insertion. An annotation the user already trusts gives
a better sheet from the same genome.

What it copies from v1, deliberately:

- the database is built beside the genome (``makeblastdb -in G -out G``), so an
  existing v1 database is reused as it is;
- ``blastn`` runs with its default task (megablast) at v1's E-value, 10e-8;
- ``-r/--remove-redundant`` is v1's ``reduce.cpp`` rule, ported unchanged.

What it does not copy: v1 plotted ``100 - pident`` as "divergence to consensus".
That is an uncorrected identity from one local alignment, not RepeatMasker's
Kimura-corrected divergence, so it is stored as
:attr:`DivergenceKind.BLASTN_IDENTITY` and panel 1's axis names it for what it is.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import tempfile
from pathlib import Path

from .records import Annotation, Copy, DivergenceKind

# v1's default, written the way v1 wrote it (TE-Aid shell script, EVALUE).
DEFAULT_EVALUE = 10e-8

_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"


class BlastnNotFound(RuntimeError):
    pass


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise BlastnNotFound(
            f"{tool} not found on PATH; the genome route needs NCBI BLAST+ "
            "(https://blast.ncbi.nlm.nih.gov/)"
        )
    return path


def _is_database(prefix: Path) -> bool:
    if shutil.which("blastdbcmd") is None:
        return any(prefix.with_name(prefix.name + ext).exists() for ext in (".nin", ".nal"))
    result = subprocess.run(["blastdbcmd", "-db", str(prefix), "-info"],
                            capture_output=True, text=True, check=False)
    return result.returncode == 0


def ensure_database(genome: Path, *, log=print) -> Path:
    """The blastn database for ``genome``, built on first use.

    Built beside the genome, as v1 did, so a v1 database is reused. When that
    directory is not writable it goes to the cache instead, keyed on the
    genome's path, size and modification time.
    """
    from .proteins import cache_root

    if _is_database(genome):
        return genome
    _require("makeblastdb")
    stat = genome.stat()
    key = hashlib.sha256(f"{genome.resolve()}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    fallback = cache_root() / "genomes" / key.hexdigest()[:16] / genome.name
    if _is_database(fallback):
        return fallback
    reason = "no writable place for it"
    for prefix in (genome, fallback):
        try:
            prefix.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            reason = str(exc)
            continue
        log(f"building a blastn database for {genome.name} at {prefix.parent} (once per genome)")
        # No -parse_seqids, as in v1: it rejects some headers v1 accepted.
        result = subprocess.run(
            ["makeblastdb", "-in", str(genome), "-out", str(prefix), "-dbtype", "nucl"],
            capture_output=True, text=True, check=False,
        )
        if result.returncode == 0 and _is_database(prefix):
            return prefix
        reason = result.stderr.strip()[:300] or f"exit {result.returncode}"
    raise RuntimeError(f"makeblastdb could not build a database for {genome}: {reason}")


def search(sequence: str, name: str, database: Path, *, evalue: float = DEFAULT_EVALUE,
           family: str | None = None, class_label: str | None = None) -> Annotation:
    """Genomic copies of one consensus, one per HSP, as an :class:`Annotation`."""
    blastn = _require("blastn")
    with tempfile.TemporaryDirectory(prefix="teaid-genome-") as tmp:
        query = Path(tmp) / "consensus.fa"
        query.write_text(f">{name}\n{sequence}\n")
        result = subprocess.run(
            [blastn, "-query", str(query), "-db", str(database), "-evalue", str(evalue),
             "-outfmt", _OUTFMT],
            capture_output=True, text=True, check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"blastn failed: {result.stderr.strip()[:300]}")
    annotation = parse(result.stdout, family=family or name, class_label=class_label,
                       consensus_length=len(sequence))
    annotation.source_path = str(database)
    return annotation


def parse(text: str, *, family: str, class_label: str | None,
          consensus_length: int) -> Annotation:
    """``-outfmt 6`` rows as copies. Coordinates become 0-based half-open."""
    annotation = Annotation(source_format="blastn")
    for lineno, line in enumerate(text.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) < 12:
            if line.strip():
                annotation.skipped.append((lineno, f"expected 12 fields, got {len(fields)}"))
            continue
        try:
            pident = float(fields[2])
            q_start, q_end = int(fields[6]), int(fields[7])
            s_start, s_end = int(fields[8]), int(fields[9])
        except ValueError as exc:
            annotation.skipped.append((lineno, str(exc)))
            continue
        strand = "+" if s_start <= s_end else "-"
        annotation.copies.append(Copy(
            chrom=fields[1],
            start=min(s_start, s_end) - 1,
            end=max(s_start, s_end),
            strand=strand,
            family=family,
            class_label=class_label,
            divergence=100.0 - pident,
            divergence_kind=DivergenceKind.BLASTN_IDENTITY,
            consensus_start=min(q_start, q_end) - 1,
            consensus_end=max(q_start, q_end),
            consensus_left=max(0, consensus_length - max(q_start, q_end)),
            source_line=lineno,
        ))
    return annotation


def remove_redundant(annotation: Annotation) -> Annotation:
    """v1's ``-r``: drop a hit whose genomic locus lies inside a kept hit's locus.

    Ported from ``reduce.cpp`` without change. Hits are taken in order of
    sequence, start ascending and end descending. A hit is dropped when a kept
    hit on the same sequence contains its locus, or has the identical locus and
    a query start no larger. Scores play no part.
    """
    order = sorted(annotation.copies, key=lambda c: (c.chrom, c.start, -c.end))
    kept: list[Copy] = []
    for copy in order:
        redundant = False
        for other in kept:
            if other.chrom != copy.chrom or other.start > copy.start or other.end < copy.end:
                continue
            if other.start < copy.start or other.end > copy.end:
                redundant = True
                break
            if (other.consensus_start or 0) <= (copy.consensus_start or 0):
                redundant = True
                break
        if not redundant:
            kept.append(copy)
    return Annotation(copies=kept, source_path=annotation.source_path,
                      source_format=annotation.source_format, skipped=annotation.skipped)
