"""Consensus sequence access.

Family names rarely agree across files: a RepeatModeler library header reads
``rnd-1_family-257#Unknown`` while the matching RepeatMasker ``.out`` rows say
``rnd-1_family-257``. Everything here matches on the bare name (the part before
``#``, case-folded) so the two line up without the caller having to care.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class Consensus:
    name: str  # as written in the FASTA, '#class' suffix included
    sequence: str

    @property
    def bare_name(self) -> str:
        return self.name.split("#", 1)[0]

    @property
    def class_label(self) -> str | None:
        return self.name.split("#", 1)[1] if "#" in self.name else None

    def __len__(self) -> int:
        return len(self.sequence)


class ConsensusLibrary:
    """A consensus FASTA, indexed by bare family name."""

    def __init__(self, entries: list[Consensus], source_path: str | None = None):
        self._entries = entries
        self.source_path = source_path
        self._by_bare: dict[str, Consensus] = {}
        for entry in entries:
            self._by_bare.setdefault(entry.bare_name.casefold(), entry)

    def __len__(self) -> int:
        return len(self._entries)

    def __iter__(self):
        return iter(self._entries)

    def get(self, family: str) -> Consensus | None:
        return self._by_bare.get(family.split("#", 1)[0].strip().casefold())

    def names(self) -> list[str]:
        return [e.name for e in self._entries]

    def suggest(self, family: str, limit: int = 5) -> list[str]:
        """Near-miss names, for a useful error when a family is not found."""
        import difflib

        want = family.split("#", 1)[0].strip().casefold()
        return difflib.get_close_matches(want, list(self._by_bare), n=limit, cutoff=0.6)


def read_fasta(path: str | Path) -> ConsensusLibrary:
    """Read a (possibly multi-entry) FASTA into a name-indexed library.

    Deliberately not Biopython: this is a flat parse with no alphabet handling,
    it is on the startup path for every run, and it keeps the sequence names
    exactly as written rather than truncating at the first space.
    """
    path = Path(path)
    entries: list[Consensus] = []
    name: str | None = None
    chunks: list[str] = []

    with path.open("r", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    entries.append(Consensus(name, "".join(chunks)))
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            elif name is not None:
                chunks.append(line)

    if name is not None:
        entries.append(Consensus(name, "".join(chunks)))

    return ConsensusLibrary(entries, source_path=str(path))
