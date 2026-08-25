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
    # Which tier the model came from, so the sheet can say what kind of evidence
    # a hit is: a conserved domain, or a named element from RepeatPeps.
    source: str = "pfam"

    @property
    def is_disrupted(self) -> bool:
        """Carries a frameshift or an in-frame stop — a pseudogenised domain.

        Worth surfacing rather than hiding: it is evidence the element was once
        coding, which a copy with an intact ORF and a copy with none do not
        distinguish between.
        """
        return self.frameshifts > 0 or self.stop_codons > 0

    @property
    def display_name(self) -> str:
        """A short label for a lane.

        Tier-1 hits are Pfam domains and already read well (``RVT_1``). Tier-2
        hits are RepeatPeps proteins, whose model name is the whole FASTA header
        (``ACROBAT1_tnp#DNA/PiggyBac``); the class after the ``#`` is redundant
        with the panel and pushes the useful half off the axis.
        """
        return self.query.split("#", 1)[0] or self.query

    @property
    def source_class(self) -> str | None:
        """The class label a RepeatPeps-derived model carries, if any."""
        return self.query.split("#", 1)[1] if "#" in self.query else None

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


def search_repeatpeps(
    orfs,
    repeatpeps: Path,
    *,
    evalue: float = 1e-5,
) -> list[ProteinHit]:
    """Tier 3: ``blastp`` of ORF peptides against RepeatPeps, as v1 did.

    Cheap (17 MB, no pHMMs) and it answers a different question from tiers 1-2:
    those say *which domain* a region encodes, this says *which named element*
    it most resembles. Only intact ORFs can hit here, by construction — that
    limitation is exactly why the frameshift-aware tiers exist alongside it.

    Peptide coordinates are converted back to the consensus by v1's rule:
    ``consensus = orf_start + 3 * (peptide_position - 1)``.
    """
    if not orfs or repeatpeps is None or not Path(repeatpeps).exists():
        return []
    if shutil.which("blastp") is None or shutil.which("makeblastdb") is None:
        return []

    with tempfile.TemporaryDirectory(prefix="teaid-blastp-") as tmp:
        tmpdir = Path(tmp)
        query = tmpdir / "orfs.faa"
        with query.open("w") as handle:
            for index, orf in enumerate(orfs):
                if orf.peptide:
                    handle.write(f">orf{index}\n{orf.peptide}\n")
        if query.stat().st_size == 0:
            return []

        database = tmpdir / "peps"
        made = subprocess.run(
            ["makeblastdb", "-in", str(repeatpeps), "-out", str(database), "-dbtype", "prot"],
            capture_output=True, text=True, check=False,
        )
        if made.returncode != 0:
            return []

        result = subprocess.run(
            ["blastp", "-query", str(query), "-db", str(database), "-evalue", str(evalue),
             "-outfmt", "6 qseqid sseqid pident qstart qend evalue bitscore", "-max_target_seqs", "5"],
            capture_output=True, text=True, check=False,
        )
        if result.returncode != 0:
            return []

    # v1's best-hit-per-ORF rule: highest bitscore wins.
    best: dict[str, list[str]] = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 7:
            continue
        current = best.get(fields[0])
        if current is None or float(fields[6]) > float(current[6]):
            best[fields[0]] = fields

    hits: list[ProteinHit] = []
    for key, fields in best.items():
        orf = orfs[int(key[3:])]
        q_start, q_end = int(fields[3]), int(fields[4])
        if orf.strand == "+":
            start = orf.start + 3 * (q_start - 1)
            end = orf.start + 3 * q_end
        else:
            end = orf.end - 3 * (q_start - 1)
            start = orf.end - 3 * q_end
        hits.append(
            ProteinHit(
                query=fields[1],
                query_accession="-",
                start=max(0, min(start, end)),
                end=max(start, end),
                strand=orf.strand,
                evalue=float(fields[5]),
                score=float(fields[6]),
                identity=float(fields[2]),
                hmm_from=q_start,
                hmm_to=q_end,
                hmm_len=max(1, orf.length_aa),
                frameshifts=0,
                stop_codons=0,
                source="repeatpeps-blastp",
            )
        )
    hits.sort(key=lambda h: h.evalue)
    return hits


def best_per_region(
    hits: list[ProteinHit],
    overlap: float = 0.5,
    *,
    per_source: bool = True,
) -> list[ProteinHit]:
    """Collapse competing hits, keeping the best-scoring per region.

    The rule is RepeatClassifier's: where models compete over the same stretch
    of consensus, the strongest wins. Cited as the source in the docs, and
    borrowed for *evidence selection only* — never for the classification
    decision RepeatClassifier goes on to make.

    **Competition is within a source, not across them.** Two reasons, and either
    alone is sufficient:

    - The scores are not comparable. A ``blastp`` bitscore over a 1,000-residue
      ORF runs into the thousands; a profile-HMM bit score for a 200-position
      domain is of order 100. Ranked together, tier 3 wins every contest by
      construction and every interpretable Pfam domain disappears — which is
      exactly what happened before this was split.
    - They answer different questions. Tiers 1-2 say *which domain* a region
      encodes; tier 3 says *which named element* it most resembles. A region
      showing both is more informative than a region showing whichever number
      happened to be larger.
    """
    if per_source:
        groups: dict[str, list[ProteinHit]] = {}
        for hit in hits:
            groups.setdefault(getattr(hit, "source", "pfam"), []).append(hit)
        merged = [
            h
            for group in groups.values()
            for h in best_per_region(group, overlap, per_source=False)
        ]
        merged.sort(key=lambda h: h.start)
        return merged

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
