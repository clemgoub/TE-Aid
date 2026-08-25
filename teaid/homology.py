"""Protein homology evidence for panel 5, via BATH.

``bathsearch --fs`` is a frameshift-aware translated pHMM search. That is the
whole reason it replaces v1's ``getorf`` + ``blastp`` for *homology*: an
ORF-finder-then-align approach cannot see a frameshifted or pseudogenised ORF by
construction, and those are exactly the copies that make TE proteins hard to
annotate. The ORF track stays ``getorf``-based and complementary — ORFs show
open frames, this shows homology regardless of frame integrity.

Borrowed from RepeatClassifier's *evidence search*, explicitly not its
classification decision. TE-Aid renders the hits and leaves the verdict alone.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

# bathsearch --tblout columns, in order. Two of them are why BATH is here at all:
# 'shifts' counts frameshifts inside the alignment and 'stops' counts in-frame
# stop codons, so a pseudogenised domain is reported *and* marked as such.
_TBL_FIELDS = (
    "hit_id", "target", "target_acc", "query", "query_acc", "hmm_len",
    "hmm_from", "hmm_to", "seq_len", "ali_from", "ali_to", "evalue",
    "score", "bias", "pid", "shifts", "stops", "description",
)


class BathNotFound(RuntimeError):
    pass


@dataclass(slots=True)
class ProteinHit:
    """One translated pHMM hit on the consensus."""

    query: str  # the domain or protein the model came from
    query_accession: str
    start: int  # 0-based half-open, consensus coordinates
    end: int
    strand: str
    evalue: float
    score: float
    identity: float
    hmm_from: int
    hmm_to: int
    hmm_len: int
    frameshifts: int
    stop_codons: int

    @property
    def is_disrupted(self) -> bool:
        """Carries a frameshift or an in-frame stop — a pseudogenised domain.

        Worth surfacing rather than hiding: it is evidence the element was once
        coding, which a copy with an intact ORF and a copy with none do not
        distinguish between.
        """
        return self.frameshifts > 0 or self.stop_codons > 0

    @property
    def coverage(self) -> float:
        """Fraction of the model the hit spans, 0-1."""
        return (self.hmm_to - self.hmm_from + 1) / self.hmm_len if self.hmm_len else 0.0

    def coordinates_1based(self) -> tuple[int, int]:
        return self.start + 1, self.end


def search(
    sequence: str,
    library,
    name: str = "consensus",
    *,
    evalue: float = 1e-3,
    frameshift_aware: bool = True,
    cpus: int | None = None,
) -> list[ProteinHit]:
    """Run ``bathsearch`` over one consensus and return its hits.

    ``library`` is a :class:`teaid.proteins.Library`.
    """
    if library is None or library.hmm is None:
        return []
    if shutil.which("bathsearch") is None:
        raise BathNotFound(
            "bathsearch not found on PATH; build BATH from "
            "https://github.com/TravisWheelerLab/BATH to draw the protein row"
        )

    with tempfile.TemporaryDirectory(prefix="teaid-bath-") as tmp:
        query = Path(tmp) / "consensus.fa"
        table = Path(tmp) / "hits.tbl"
        query.write_text(f">{name}\n{sequence}\n")

        command = [
            "bathsearch",
            "--tblout", str(table),
            "-o", "/dev/null",
            "-E", str(evalue),
        ]
        if frameshift_aware:
            command.append("--fs")
        if cpus:
            command += ["--cpu", str(cpus)]
        command += [str(library.hmm), str(query)]

        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"bathsearch failed: {result.stderr.strip()[:400]}")
        text = table.read_text() if table.exists() else ""

    return _parse_tblout(text)


def _parse_tblout(text: str) -> list[ProteinHit]:
    hits: list[ProteinHit] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        # The final column is a free-text description that may contain spaces,
        # so it takes the remainder. Splitting one field short instead makes the
        # last numeric column swallow it, and every row is then dropped by the
        # parse guard below -- silently, since a bad row is meant to be skipped.
        fields = line.split(None, len(_TBL_FIELDS) - 1)
        if len(fields) < len(_TBL_FIELDS) - 1:
            continue
        row = dict(zip(_TBL_FIELDS, fields))
        try:
            ali_from, ali_to = int(row["ali_from"]), int(row["ali_to"])
            # bathsearch reports a reverse-strand hit with descending
            # coordinates; the interval is normalised and the direction kept.
            strand = "+" if ali_to >= ali_from else "-"
            start, end = (ali_from, ali_to) if strand == "+" else (ali_to, ali_from)
            hits.append(
                ProteinHit(
                    query=row["query"],
                    query_accession=row["query_acc"],
                    start=start - 1,
                    end=end,
                    strand=strand,
                    evalue=float(row["evalue"]),
                    score=float(row["score"]),
                    identity=float(row["pid"]),
                    hmm_from=int(row["hmm_from"]),
                    hmm_to=int(row["hmm_to"]),
                    hmm_len=int(row["hmm_len"]),
                    frameshifts=int(row["shifts"]),
                    stop_codons=int(row["stops"]),
                )
            )
        except (ValueError, KeyError):
            continue
    hits.sort(key=lambda h: h.evalue)
    return hits


def best_per_region(hits: list[ProteinHit], overlap: float = 0.5) -> list[ProteinHit]:
    """Collapse competing hits, keeping the best-scoring per region.

    The rule is RepeatClassifier's: where models compete over the same stretch
    of consensus, the strongest wins. Cited as the source in the docs, and
    borrowed for *evidence selection only* — never for the classification
    decision RepeatClassifier goes on to make.
    """
    kept: list[ProteinHit] = []
    for hit in sorted(hits, key=lambda h: (-h.score, h.evalue)):
        span = hit.end - hit.start
        if span <= 0:
            continue
        clash = False
        for other in kept:
            shared = min(hit.end, other.end) - max(hit.start, other.start)
            if shared > 0 and shared / min(span, other.end - other.start) >= overlap:
                clash = True
                break
        if not clash:
            kept.append(hit)
    kept.sort(key=lambda h: h.start)
    return kept
