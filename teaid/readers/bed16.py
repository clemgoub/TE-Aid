"""Reader for BED16, the VGP repeat-hub interchange format.

Column definitions are normative in ``VGP_TEbed/docs/INPUT_FORMAT.md``::

     1 chrom                 9 perc_del
     2 chromStart (0-based) 10 perc_ins
     3 chromEnd (half-open) 11 query_left
     4 name                 12 repeat_class_family
     5 score (0-1000)       13 repeat_start
     6 strand               14 repeat_end
     7 SW_score             15 repeat_left
     8 perc_div             16 hit_id

Tab-separated, already 0-based half-open (unlike ``.out``), optional ``#``
header, optionally gzipped. **Position is authoritative** — a header row's
column names are ignored, and real files disagree about them (some write
``RM_ID`` for column 16).

Two column pairs are easy to transpose and mean different things:
``query_left`` (11) is genome-side, ``repeat_left`` (15) is consensus-side.
This reader only takes 15.

Permissiveness is deliberate. The spec marks columns 5, 6, 12 and 16
non-nullable, but production files from every structural caller violate all
four benignly — ``score=0`` as a stand-in, ``strand=.`` for tandem arrays,
``repeat_class_family=NA`` as REPET's abstention token. A validator built from
the table alone rejects roughly half the shipped inputs, so nullability is
enforced only where a missing value would make the record meaningless.
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

from ..records import Annotation, Copy, DivergenceKind

# Documented null tokens (INPUT_FORMAT §1), plus the empty field.
_NULL = frozenset({"NA", "na", ".", "", "nan", "NaN"})

_HEADER_FIRST_FIELDS = frozenset({"chrom", "chr", "chromosome"})

# BED16 carries no statement of which quantity perc_div holds; that lives in the
# hub's per-tool manifest, which a standalone TE-Aid run does not have. These
# class labels are produced by tandem/array callers whose divergence is
# unit-to-unit array homogeneity rather than divergence from a consensus.
_ARRAY_HOMOGENEITY_CLASSES = frozenset(
    {"satellite", "simple_repeat", "low_complexity", "tandem", "trf"}
)


def _open(path: Path) -> io.TextIOBase:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def _text(token: str) -> str | None:
    return None if token in _NULL else token


def _num(token: str, cast) -> float | int | None:
    if token in _NULL:
        return None
    try:
        return cast(token)
    except ValueError:
        # The spec's own coercion rule: an unparseable token becomes null
        # rather than killing the row.
        return None


def _is_header(fields: list[str]) -> bool:
    return fields[0].startswith("#") or fields[0].casefold() in _HEADER_FIRST_FIELDS


def read(path: str | Path, *, divergence_kind: DivergenceKind | None = None) -> Annotation:
    """Parse a BED16 file into the shared record struct.

    ``divergence_kind`` overrides the per-record heuristic. Pass it when you
    know what produced the file — a structural caller's ``perc_div``, if it has
    one at all, is rarely divergence from a consensus.
    """
    path = Path(path)
    annotation = Annotation(source_path=str(path), source_format="bed16")

    with _open(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            fields = line.split("\t")
            if _is_header(fields):
                continue
            try:
                annotation.copies.append(_parse_row(fields, lineno, divergence_kind))
            except (ValueError, IndexError) as exc:
                annotation.skipped.append((lineno, str(exc)))

    return annotation


def _parse_row(
    fields: list[str], lineno: int, forced_kind: DivergenceKind | None
) -> Copy:
    if len(fields) < 16:
        raise ValueError(f"expected >=16 tab-separated fields, got {len(fields)}")

    start = _num(fields[1], int)
    end = _num(fields[2], int)
    if start is None or end is None:
        raise ValueError("chromStart/chromEnd missing or non-numeric")

    strand = fields[5] if fields[5] in {"+", "-"} else "."
    class_label = _text(fields[11])
    divergence = _num(fields[7], float)

    if divergence is None:
        kind = DivergenceKind.NONE
    elif forced_kind is not None:
        kind = forced_kind
    else:
        kind = _infer_divergence_kind(class_label)

    consensus_start = _num(fields[12], int)
    consensus_end = _num(fields[13], int)
    # repeat_start/repeat_end are consensus-side and, per the .out conversion
    # rule the hub applies upstream, already reordered to ascending. Guard
    # anyway; a producer that skipped the reorder should not transpose silently.
    if (
        consensus_start is not None
        and consensus_end is not None
        and consensus_start > consensus_end
    ):
        consensus_start, consensus_end = consensus_end, consensus_start

    return Copy(
        chrom=fields[0],
        start=start,
        end=end,
        strand=strand,
        family=fields[3],
        class_label=class_label,
        divergence=divergence,
        divergence_kind=kind,
        # BED16's repeat_start is 1-based inclusive (inherited unchanged from
        # the .out it was converted from) while the genomic columns are 0-based;
        # the two conventions genuinely do coexist in one row.
        consensus_start=None if consensus_start is None else consensus_start - 1,
        consensus_end=consensus_end,
        consensus_left=_num(fields[14], int),
        hit_id=_text(fields[15]),
        source_line=lineno,
    )


def _infer_divergence_kind(class_label: str | None) -> DivergenceKind:
    if class_label is None:
        return DivergenceKind.CONSENSUS
    head = class_label.split("/", 1)[0].casefold()
    return (
        DivergenceKind.ARRAY_HOMOGENEITY
        if head in _ARRAY_HOMOGENEITY_CLASSES
        else DivergenceKind.CONSENSUS
    )
