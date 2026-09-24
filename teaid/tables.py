"""Coordinate tables beside the sheet: v1's ``-t`` and ``-T``.

v1 wrote the self dot-plot, ORF and protein-hit coordinates as text files with
``-t``, and added the genomic hit table with ``-T``, for anyone scripting on top
of a sheet. The HTML sheet already offers the protein hits as a TSV download,
but a pipeline cannot click a button, so these are written as files.

Every table has a header row and uses **1-based, fully closed** coordinates, as
v1's did and as a human reading them expects; internally everything is 0-based
half-open (see ``records``).
"""

from __future__ import annotations

import csv
from pathlib import Path

from .report import SheetData, protein_hit_rows


def _write(path: Path, header: list[str], rows: list[list]) -> Path:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def write(data: SheetData, out_dir: Path, stem: str, *, copies: bool = False) -> list[Path]:
    """Write the tables for one sheet; ``copies`` adds the genomic copy table."""
    written = []
    # Self hits keep blastn's orientation: s_start > s_end on a reverse hit, as
    # in v1's table. The stored pair is (sstart - 1, send) either way.
    self_rows = [
        [h.q_start + 1, h.q_end, h.s_start + 1, h.s_end, "-" if h.is_reverse else "+",
         f"{h.identity:.1f}", f"{h.evalue:.2g}", f"{h.bitscore:.1f}"]
        for h in data.self_hits if not h.is_trivial_diagonal
    ]
    written.append(_write(out_dir / f"{stem}.selfhits.tsv",
                          ["q_start", "q_end", "s_start", "s_end", "strand", "identity_pct",
                           "evalue", "bitscore"], self_rows))
    written.append(_write(out_dir / f"{stem}.orfs.tsv",
                          ["start", "end", "strand", "length_aa"],
                          [[*o.coordinates_1based(), o.strand, o.length_aa] for o in data.orfs]))
    hits = protein_hit_rows(data)
    header = list(hits[0]) if hits else [
        "consensus", "name", "class", "accession", "start", "end", "strand", "evalue",
        "score", "identity_pct", "model_coverage_pct", "frameshifts", "stop_codons",
        "source", "drawn"]
    written.append(_write(out_dir / f"{stem}.protein_hits.tsv", header,
                          [[row[k] for k in header] for row in hits]))
    if copies:
        full = {id(i) for i in data.full_length}
        rows = []
        for n, insertion in enumerate(data.insertions, start=1):
            for fragment in insertion.fragments:
                g0, g1 = fragment.genomic_1based()
                c = fragment.consensus_1based() or ("", "")
                rows.append([
                    n, fragment.chrom, g0, g1, fragment.strand, c[0], c[1],
                    "" if fragment.divergence is None else f"{fragment.divergence:.2f}",
                    fragment.divergence_kind.value, fragment.hit_id or "",
                    insertion.n_fragments,
                    "yes" if id(insertion) in full else "no",
                ])
        written.append(_write(out_dir / f"{stem}.copies.tsv",
                              ["copy", "chrom", "start", "end", "strand", "consensus_start",
                               "consensus_end", "divergence", "divergence_kind", "hit_id",
                               "rows_in_copy", "full_length"], rows))
    return written
