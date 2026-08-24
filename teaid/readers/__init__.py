"""Annotation readers, and the sniffing that picks between them.

One shared record struct (:mod:`teaid.records`), three annotation dialects.
Callers should use :func:`read` and let the format be detected, or pass an
explicit ``fmt`` when the extension lies.
"""

from __future__ import annotations

import gzip
from pathlib import Path

from ..records import Annotation
from . import bed16, gff3, rmout

FORMATS = ("rmout", "gff3", "bed16")

_READERS = {"rmout": rmout.read, "gff3": gff3.read, "bed16": bed16.read}


class FormatDetectionError(ValueError):
    """Raised when a file matches no known annotation dialect."""


def _first_data_lines(path: Path, limit: int = 40) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    lines: list[str] = []
    with opener(path, "rt", errors="replace") as handle:  # type: ignore[operator]
        for line in handle:
            if line.strip():
                lines.append(line.rstrip("\n"))
            if len(lines) >= limit:
                break
    return lines


def detect(path: str | Path) -> str:
    """Identify the annotation dialect by content, not by extension.

    Content wins because the extensions are unreliable in the wild: the hub
    ships ``.bed`` files converted from ``.out``, and ``.gff``/``.gff3`` are
    used interchangeably for two different attribute syntaxes.
    """
    path = Path(path)
    lines = _first_data_lines(path)
    if not lines:
        raise FormatDetectionError(f"{path} is empty")

    if any(line.lstrip().startswith("##gff-version") for line in lines):
        return "gff3"

    for line in lines:
        if line.startswith("#"):
            continue
        tabbed = line.split("\t")
        # GFF is 9 tab-separated columns whose 4th and 5th are integers and
        # whose 9th carries a Target attribute.
        if len(tabbed) == 9 and tabbed[3].isdigit() and tabbed[4].isdigit():
            if "Target" in tabbed[8]:
                return "gff3"
        # BED16 is >=16 tab-separated columns with integer coordinates.
        if len(tabbed) >= 16 and tabbed[1].isdigit() and tabbed[2].isdigit():
            return "bed16"
        # RepeatMasker .out is whitespace-aligned: numeric SW score, then a
        # float divergence, with '+' or 'C' in the strand column.
        fields = line.split()
        if len(fields) >= 15 and fields[0].lstrip("-").isdigit() and fields[8] in {"+", "C"}:
            return "rmout"

    # A .out with only header lines in the sample still identifies by its header.
    if any(line.lstrip().startswith("SW") for line in lines[:3]):
        return "rmout"

    raise FormatDetectionError(
        f"{path}: not a recognisable RepeatMasker .out, GFF3, or BED16 file"
    )


def read(path: str | Path, *, fmt: str | None = None, **kwargs) -> Annotation:
    """Read an annotation, detecting the dialect unless ``fmt`` is given."""
    resolved = fmt or detect(path)
    if resolved not in _READERS:
        raise ValueError(f"unknown annotation format {resolved!r}; expected one of {FORMATS}")
    reader = _READERS[resolved]
    if resolved != "bed16":
        kwargs.pop("divergence_kind", None)
    return reader(path, **kwargs)


__all__ = ["read", "detect", "FORMATS", "FormatDetectionError", "rmout", "gff3", "bed16"]
