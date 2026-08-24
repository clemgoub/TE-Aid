"""Reader for RepeatMasker ``.out`` annotation files.

The format is fixed-width in spirit but ragged in practice, so fields are taken
by splitting on whitespace, not by column offset. Three header lines are
skipped (two of column titles, one blank).

The layout, per row::

    SW  div  del  ins  query  qbegin  qend  (qleft)  strand  repeat  class/family
        rbegin  rend  (rleft)  ID  [*]

Two traps this reader exists to absorb:

1. ``strand`` is ``+`` or ``C`` (complement), never ``-``.
2. On ``C`` rows the three consensus columns are written **in reverse order**,
   i.e. ``(left) end begin`` instead of ``begin end (left)``. Reading them
   positionally without accounting for this silently transposes the consensus
   coordinates of every reverse-strand copy — which is most of them.

Verified against GenomeArk RepeatModeler-v2.0.8 systematic annotations, e.g.::

    14672  5.0 1.8 0.4  OY720097.1  8303 10606 (23705997) C rnd-1_family-257 \
        Unknown  (0) 1755 1  3

where the family consensus is 1755 bp: begin=1, end=1755, left=0.
"""

from __future__ import annotations

import re
from pathlib import Path

from ..records import Annotation, Copy, DivergenceKind

# RepeatMasker delegates these classes to TRF and dust, which report how uniform
# an array is rather than how far it has drifted from a library consensus. Same
# column, different meaning — see DivergenceKind.
_ARRAY_HOMOGENEITY_CLASSES = frozenset({"simple_repeat", "low_complexity"})

_PAREN = re.compile(r"^\((-?\d+)\)$")


def _int_maybe_paren(token: str) -> int:
    """Parse ``123`` or ``(123)`` — RepeatMasker parenthesises 'remaining'."""
    m = _PAREN.match(token)
    return int(m.group(1)) if m else int(token)


def _divergence_kind(class_label: str) -> DivergenceKind:
    return (
        DivergenceKind.ARRAY_HOMOGENEITY
        if class_label.split("/", 1)[0].casefold() in _ARRAY_HOMOGENEITY_CLASSES
        else DivergenceKind.CONSENSUS
    )


def read(path: str | Path) -> Annotation:
    """Parse a RepeatMasker ``.out`` file into the shared record struct."""
    path = Path(path)
    annotation = Annotation(source_path=str(path), source_format="rmout")

    with path.open("r", errors="replace") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split()
            # Header lines start with the column titles; the SW score column is
            # numeric on every data row, which is a cheaper and more robust test
            # than counting how many header lines to skip.
            if not fields[0].lstrip("-").isdigit():
                continue
            try:
                annotation.copies.append(_parse_row(fields, lineno))
            except (ValueError, IndexError) as exc:
                annotation.skipped.append((lineno, str(exc)))

    return annotation


def _parse_row(fields: list[str], lineno: int) -> Copy:
    if len(fields) < 15:
        raise ValueError(f"expected >=15 fields, got {len(fields)}")

    divergence = float(fields[1])
    chrom = fields[4]
    q_begin = int(fields[5])  # 1-based, fully closed
    q_end = int(fields[6])
    raw_strand = fields[8]
    family = fields[9]
    class_label = fields[10]

    if raw_strand == "+":
        strand = "+"
        c_begin = _int_maybe_paren(fields[11])
        c_end = _int_maybe_paren(fields[12])
        c_left = _int_maybe_paren(fields[13])
    elif raw_strand == "C":
        strand = "-"
        # Reversed on complement rows: (left) end begin
        c_left = _int_maybe_paren(fields[11])
        c_end = _int_maybe_paren(fields[12])
        c_begin = _int_maybe_paren(fields[13])
    else:
        raise ValueError(f"unrecognised strand {raw_strand!r} (expected '+' or 'C')")

    # Belt and braces: whatever the row order claimed, the consensus interval is
    # stored ascending.
    if c_begin > c_end:
        c_begin, c_end = c_end, c_begin

    return Copy(
        chrom=chrom,
        start=q_begin - 1,  # to 0-based half-open
        end=q_end,
        strand=strand,
        family=family,
        class_label=class_label,
        divergence=divergence,
        divergence_kind=_divergence_kind(class_label),
        consensus_start=c_begin - 1,
        consensus_end=c_end,
        consensus_left=c_left,
        is_overlapping_lower_score=fields[-1] == "*",
        source_line=lineno,
    )
