"""Reader for RepeatMasker GFF output.

RepeatMasker emits two dialects under the ``.gff``/``.gff3`` name and both are
in circulation, so both are accepted:

GFF3 (``-gff`` with RepeatMasker >= 4.1)::

    OY720097.1  RepeatMasker  dispersed_repeat  4  8148  6759  +  .  \
        Target=(ACCCTG)n 1 8154

GFF2-style (older releases), same columns but space-separated attributes with
the family name quoted and prefixed::

    OY720097.1  RepeatMasker  similarity  4  8148  6759  +  .  \
        Target "Motif:(ACCCTG)n" 1 8154

Coordinates in columns 4/5 are 1-based fully closed. The ``Target`` attribute
carries the consensus interval, also 1-based closed, and is the only place
consensus coordinates appear — a GFF without it yields records with no
consensus positions, which panels 1 and 2 cannot use.

The class/family label is not part of the standard RepeatMasker GFF output. It
is read from a ``Class``/``class`` attribute when a producer supplies one, and
is otherwise None.
"""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path
from urllib.parse import unquote

from ..records import Annotation, Copy, DivergenceKind

_ARRAY_HOMOGENEITY_CLASSES = frozenset({"simple_repeat", "low_complexity", "satellite"})

# Target=NAME START END  |  Target "Motif:NAME" START END
_TARGET = re.compile(
    r'Target[=\s]+"?(?:Motif:)?(?P<name>[^";\s]+)"?\s+(?P<start>\d+)\s+(?P<end>\d+)'
)


def _open(path: Path) -> io.TextIOBase:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def _attributes(blob: str) -> dict[str, str]:
    """Parse GFF3 ``key=value;`` pairs. GFF2 space-syntax is handled separately."""
    attrs: dict[str, str] = {}
    for part in blob.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, _, value = part.partition("=")
        attrs[key.strip()] = unquote(value.strip().strip('"'))
    return attrs


def read(path: str | Path) -> Annotation:
    path = Path(path)
    annotation = Annotation(source_path=str(path), source_format="gff3")

    with _open(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            try:
                annotation.copies.append(_parse_row(fields, lineno))
            except (ValueError, IndexError) as exc:
                annotation.skipped.append((lineno, str(exc)))

    return annotation


def _parse_row(fields: list[str], lineno: int) -> Copy:
    if len(fields) < 9:
        raise ValueError(f"expected 9 tab-separated columns, got {len(fields)}")

    start = int(fields[3])  # 1-based closed
    end = int(fields[4])
    strand = fields[6] if fields[6] in {"+", "-"} else "."
    attrs = _attributes(fields[8])

    target = _TARGET.search(fields[8])
    if target is None:
        raise ValueError("no Target attribute; consensus coordinates unavailable")
    family = target.group("name")
    c_start, c_end = int(target.group("start")), int(target.group("end"))
    if c_start > c_end:
        c_start, c_end = c_end, c_start

    class_label = attrs.get("Class") or attrs.get("class") or attrs.get("Classification")

    # RepeatMasker GFF has no divergence column; column 6 is the SW score, which
    # is emphatically not a divergence. Recording it as one would put an
    # arbitrary alignment score on an axis labelled 'divergence (%)'.
    divergence = None
    kind = DivergenceKind.NONE
    if "Divergence" in attrs or "perc_div" in attrs:
        raw_div = attrs.get("Divergence") or attrs.get("perc_div")
        try:
            divergence = float(raw_div)
            head = (class_label or "").split("/", 1)[0].casefold()
            kind = (
                DivergenceKind.ARRAY_HOMOGENEITY
                if head in _ARRAY_HOMOGENEITY_CLASSES
                else DivergenceKind.CONSENSUS
            )
        except (TypeError, ValueError):
            divergence = None

    return Copy(
        chrom=fields[0],
        start=start - 1,
        end=end,
        strand=strand,
        family=family,
        class_label=class_label,
        divergence=divergence,
        divergence_kind=kind,
        consensus_start=c_start - 1,
        consensus_end=c_end,
        consensus_left=None,
        hit_id=attrs.get("ID"),
        source_line=lineno,
    )
