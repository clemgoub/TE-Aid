"""The protein query library behind panel 5, in three tiers.

The tiering is a **cost** decision taken deliberately (see ``docs/BRIEF_v2.md``
§5.1a); the numbers behind it are worth keeping close to the code:

- One single-sequence pHMM per RepeatPeps protein is **~6.4 GB** — measured, at
  718 MB for the first 2,022 of 18,011 sequences.
- Reducing redundancy first does not help. cd-hit at 90% identity removes 0.5%
  of RepeatPeps; it is already a non-redundant curated library.
- But RepeatPeps cannot be dropped either, because **Pfam has no model at all**
  for piggyBac, Maverick/Polinton and Crypton, and only a generic GIY-YIG for
  Penelope. Those four superfamilies are 827 RepeatPeps proteins.

So:

===== ============================================ ========= ===================
Tier  What                                         Size      When
===== ============================================ ========= ===================
1     curated Pfam TE domains as pHMMs             ~10 MB    always
2     pHMMs for *only* the superfamilies Pfam       ~200 MB   always
      cannot model
3     ``blastp`` against the whole RepeatPeps       17 MB     fallback
      FASTA, as v1 did
deep  the whole of RepeatPeps as pHMMs             ~6.4 GB   opt-in
===== ============================================ ========= ===================

**This is not a sensitivity claim.** Whether tiers 1–2 recover what ``--deep``
would is what the §5.5 benchmark exists to measure; the framework is built first
so the benchmark has something to run against.

Everything is cached under ``~/.teaid/proteins/``, keyed by a signature of the
inputs, so editing ``te_domains.tsv`` rebuilds rather than silently reusing a
stale library.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import shutil
import subprocess
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

DATA = Path(__file__).parent / "data"
TE_DOMAINS = DATA / "te_domains.tsv"

# Superfamilies with no usable Pfam model, so RepeatPeps is the only cover.
# Keys are the class labels RepeatPeps uses in its FASTA headers.
PFAM_BLIND_SUPERFAMILIES = {
    "DNA/PiggyBac": "no Pfam model; InterPro covers piggyBac only via PANTHER and CDD",
    "DNA/Maverick": "no Pfam model, and no transposon-specific InterPro entry",
    "DNA/Cryp": "no Pfam model; only the generic tyrosine-recombinase fold",
    "LINE/Penelope": "only the generic GIY-YIG PF01541; the PLE-specific model is CDD's",
}

_HMM_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{}/?annotation=hmm"

# Where RepeatPeps is looked for when it is not given explicitly. It ships with
# RepeatMasker rather than being separately downloadable — the URL v1 used is
# dead, since RepeatMasker moved to Dfam-consortium and dropped it from git.
_REPEATPEPS_GUESSES = (
    "~/RepeatPeps.lib",
    "~/Downloads/RepeatMasker/Libraries/RepeatPeps.lib",
    "/usr/local/RepeatMasker/Libraries/RepeatPeps.lib",
    "/opt/RepeatMasker/Libraries/RepeatPeps.lib",
    "/opt/homebrew/share/RepeatMasker/Libraries/RepeatPeps.lib",
)


class LibraryError(RuntimeError):
    """The protein library could not be built or found."""


@dataclass(slots=True)
class Library:
    """A built protein query library, ready for ``bathsearch``."""

    hmm: Path | None  # BATH-format pHMMs (tiers 1-2, or the deep library)
    repeatpeps: Path | None  # RepeatPeps FASTA, for the blastp fallback
    tiers: tuple[str, ...] = ()

    @property
    def usable(self) -> bool:
        return self.hmm is not None or self.repeatpeps is not None


def cache_root() -> Path:
    """``~/.teaid/proteins``, overridable for tests and for shared installs."""
    root = os.environ.get("TEAID_CACHE")
    return Path(root).expanduser() if root else Path.home() / ".teaid" / "proteins"


def curated_accessions() -> list[tuple[str, str]]:
    """(accession, name) for every domain in the curated table."""
    if not TE_DOMAINS.exists():
        raise LibraryError(f"curated domain table missing: {TE_DOMAINS}")
    out = []
    for line in TE_DOMAINS.read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        out.append((parts[0], parts[1]))
    return out


def find_repeatpeps(explicit: str | Path | None = None) -> Path | None:
    """Locate RepeatPeps.lib, which ships inside a RepeatMasker installation."""
    if explicit:
        path = Path(explicit).expanduser()
        if not path.exists():
            raise LibraryError(f"RepeatPeps not found: {path}")
        return path
    for guess in _REPEATPEPS_GUESSES:
        path = Path(guess).expanduser()
        if path.exists():
            return path
    return None


def _signature(deep: bool, repeatpeps: Path | None) -> str:
    """Identity of the inputs, so an edited table or library rebuilds."""
    digest = hashlib.sha256()
    digest.update(TE_DOMAINS.read_bytes() if TE_DOMAINS.exists() else b"")
    digest.update(b"deep" if deep else b"default")
    if repeatpeps and repeatpeps.exists():
        stat = repeatpeps.stat()
        digest.update(f"{repeatpeps.name}:{stat.st_size}".encode())
    return digest.hexdigest()[:16]


def _require(tool: str) -> str:
    found = shutil.which(tool)
    if found is None:
        raise LibraryError(
            f"{tool} not found on PATH. BATH provides bathbuild/bathconvert/bathsearch; "
            f"build it from https://github.com/TravisWheelerLab/BATH"
        )
    return found


def fetch_pfam_hmms(accessions: list[str], dest: Path, *, log=print) -> int:
    """Download Pfam's own HMMER models for the curated accessions.

    Fetched rather than vendored: the models are large, they are versioned
    upstream, and a stale copy in the repo would silently drift from the
    accession list beside it.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with dest.open("w") as out:
        for index, acc in enumerate(accessions, start=1):
            for attempt in range(4):
                try:
                    request = urllib.request.Request(
                        _HMM_URL.format(acc), headers={"Accept": "*/*"}
                    )
                    raw = urllib.request.urlopen(request, timeout=60).read()
                    break
                except (urllib.error.URLError, TimeoutError, OSError):
                    if attempt == 3:
                        raw = b""
                    time.sleep(1.5 * (attempt + 1))
            if not raw:
                log(f"  warning: could not fetch {acc}; skipping")
                continue
            text = gzip.decompress(raw).decode() if raw[:2] == b"\x1f\x8b" else raw.decode()
            out.write(text if text.endswith("\n") else text + "\n")
            written += 1
            if index % 25 == 0:
                log(f"  fetched {index}/{len(accessions)} Pfam models")
    return written


def extract_pfam_blind(repeatpeps: Path, dest: Path) -> int:
    """RepeatPeps entries for the superfamilies Pfam cannot model (tier 2)."""
    keep = tuple(PFAM_BLIND_SUPERFAMILIES)
    written = 0
    with repeatpeps.open(errors="replace") as handle, dest.open("w") as out:
        emit = False
        for line in handle:
            if line.startswith(">"):
                label = line[1:].split("#", 1)[1].split()[0] if "#" in line else ""
                emit = label.startswith(keep)
                written += emit
            if emit:
                out.write(line)
    return written


def build(
    *,
    deep: bool = False,
    repeatpeps: str | Path | None = None,
    force: bool = False,
    log=print,
) -> Library:
    """Build (or reuse) the protein library, returning the paths to search.

    Degrades rather than failing: without RepeatPeps, tier 1 alone is still a
    usable library and the sheet says which tiers it had.
    """
    peps = find_repeatpeps(repeatpeps)
    signature = _signature(deep, peps)
    home = cache_root() / signature
    hmm = home / ("deep.bhmm" if deep else "default.bhmm")
    stamp = home / ("deep.done" if deep else "default.done")
    tiers: list[str] = []

    if stamp.exists() and hmm.exists() and not force:
        tiers = stamp.read_text().split(",")
        return Library(hmm=hmm, repeatpeps=peps, tiers=tuple(t for t in tiers if t))

    home.mkdir(parents=True, exist_ok=True)
    _require("bathconvert")
    _require("bathbuild")
    parts: list[Path] = []

    if deep:
        if peps is None:
            raise LibraryError(
                "--deep needs RepeatPeps.lib; pass --repeatpeps or install RepeatMasker"
            )
        log("building the deep library: all of RepeatPeps as pHMMs (~6.4 GB, slow)")
        deep_hmm = home / "repeatpeps.bhmm"
        _run([_require("bathbuild"), str(deep_hmm), str(peps)])
        parts.append(deep_hmm)
        tiers.append("deep:RepeatPeps")
    else:
        # Tier 1 — the curated Pfam domains.
        accessions = [acc for acc, _ in curated_accessions()]
        raw = home / "pfam.hmm"
        if not raw.exists() or force:
            log(f"fetching {len(accessions)} Pfam models from InterPro")
            got = fetch_pfam_hmms(accessions, raw, log=log)
            log(f"  {got}/{len(accessions)} models fetched")
        pfam = home / "pfam.bhmm"
        if not pfam.exists() or force:
            _run([_require("bathconvert"), str(pfam), str(raw)])
        parts.append(pfam)
        tiers.append(f"pfam:{len(accessions)}")

        # Tier 2 — only the superfamilies Pfam cannot model.
        if peps is not None:
            subset = home / "pfam_blind.fa"
            if not subset.exists() or force:
                n = extract_pfam_blind(peps, subset)
                log(f"extracting {n} RepeatPeps proteins Pfam cannot model")
            blind = home / "pfam_blind.bhmm"
            if not blind.exists() or force:
                log("building tier-2 pHMMs (~200 MB)")
                _run([_require("bathbuild"), str(blind), str(subset)])
            parts.append(blind)
            tiers.append("repeatpeps-gaps")
        else:
            log("RepeatPeps not found: piggyBac, Maverick, Crypton and Penelope "
                "have no Pfam model and will not be searched")

    # A single part is moved into place rather than copied. Concatenating would
    # double the library on disk, which is tolerable at 233 MB and impossible
    # for the ~6.4 GB deep build.
    if len(parts) == 1:
        parts[0].replace(hmm)
    else:
        with hmm.open("wb") as out:
            for part in parts:
                out.write(part.read_bytes())
        # The parts are rebuildable from the small cached inputs beside them
        # (pfam.hmm, pfam_blind.fa), so keeping both costs disk for nothing.
        for part in parts:
            part.unlink(missing_ok=True)
    stamp.write_text(",".join(tiers))
    return Library(hmm=hmm, repeatpeps=peps, tiers=tuple(tiers))


def _run(command: list[str]) -> None:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise LibraryError(f"{Path(command[0]).name} failed: {result.stderr.strip()[:400]}")
