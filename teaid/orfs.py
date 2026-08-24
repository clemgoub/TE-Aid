"""Open reading frames, via EMBOSS ``getorf``.

Kept as an external ``getorf`` call rather than reimplemented, matching v1 (see
``docs/BRIEF_v2.md`` §5.3). EMBOSS stays a dependency; only ``dotmatcher`` was
retired, replaced by the self-blastn dot-plot.

``getorf`` writes each ORF as a FASTA record whose header carries the
nucleotide interval on the input sequence::

    >consensus_1 [412 - 1683] Description
    >consensus_2 [2210 - 1109] (REVERSE SENSE) Description

On reverse-sense ORFs the two numbers descend. Both are 1-based fully closed and
are normalised here to 0-based half-open, ascending, with the direction carried
by ``strand`` instead of by coordinate order.

The ORF track and the protein homology row answer different questions: an ORF
shows an intact open frame, while a translated pHMM search finds homology
whether or not the frame survived. A pseudogenised copy shows protein homology
and no ORF, which is informative rather than contradictory.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

# '[start - end]' optionally followed by '(REVERSE SENSE)'
_HEADER = re.compile(r"\[\s*(\d+)\s*-\s*(\d+)\s*\]\s*(\(REVERSE SENSE\))?")


class GetorfNotFound(RuntimeError):
    pass


@dataclass(slots=True)
class ORF:
    """One open reading frame on the consensus, 0-based half-open."""

    start: int
    end: int
    strand: str  # '+' or '-'
    peptide: str = ""

    @property
    def length_nt(self) -> int:
        return self.end - self.start

    @property
    def length_aa(self) -> int:
        return len(self.peptide) if self.peptide else self.length_nt // 3

    def coordinates_1based(self) -> tuple[int, int]:
        return self.start + 1, self.end


def find_orfs(
    sequence: str,
    name: str = "consensus",
    *,
    min_size: int = 400,
    reverse: bool = True,
) -> list[ORF]:
    """Run ``getorf`` over a consensus and return ORFs sorted by position.

    ``min_size`` is in nucleotides and matches v1's ``--min-orf`` default of 400.
    ``getorf`` is invoked with its default ``-find 0`` (translate the regions
    between stop codons), also as in v1.
    """
    if shutil.which("getorf") is None:
        raise GetorfNotFound(
            "getorf not found on PATH; install EMBOSS to draw the ORF track"
        )

    with tempfile.TemporaryDirectory(prefix="teaid-getorf-") as tmp:
        query = Path(tmp) / "query.fa"
        out = Path(tmp) / "orfs.fa"
        query.write_text(f">{name}\n{sequence}\n")

        command = [
            "getorf",
            "-sequence", str(query),
            "-outseq", str(out),
            "-minsize", str(min_size),
            "-reverse" if reverse else "-noreverse",
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"getorf failed: {result.stderr.strip()}")

        text = out.read_text() if out.exists() else ""

    return _parse(text, len(sequence))


def _parse(fasta_text: str, sequence_length: int) -> list[ORF]:
    orfs: list[ORF] = []
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        if header is None:
            return
        match = _HEADER.search(header)
        if match is None:
            return
        first, second = int(match.group(1)), int(match.group(2))
        reverse_sense = match.group(3) is not None
        start, end = (second, first) if first > second else (first, second)
        orf = ORF(
            start=max(0, start - 1),
            end=min(sequence_length, end),
            strand="-" if reverse_sense or first > second else "+",
            peptide="".join(chunks),
        )
        if orf.end > orf.start:
            orfs.append(orf)

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header, chunks = line, []
        else:
            chunks.append(line)
    flush()

    return sorted(orfs, key=lambda o: (o.start, o.end))
