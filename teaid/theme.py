"""Colour and chrome tokens for the evidence sheet.

Three categorical slots carry three fixed roles across every panel, so a colour
means the same thing wherever it appears on the sheet:

- **base** — the ordinary genomic evidence: every annotated copy, the coverage
  trace, same-strand (direct) self-matches.
- **highlight** — the notable subset: full-length copies.
- **inverted** — opposite-strand features: TIR-like self-matches.

Both palettes were validated with the data-viz validator over all pairs
(``--pairs all``): light passes lightness band, chroma, CVD separation (worst
deutan ΔE 9.2) and normal-vision floor (worst ΔE 24.0); dark passes all five
checks including 3:1 contrast. On the light surface ``inverted`` sits at 2.74:1,
below the 3:1 bar, so it is never the sole carrier of meaning — every series is
named in a legend and, where it matters, direct-labelled.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Theme:
    name: str
    surface: str
    paper: str
    text_primary: str
    text_secondary: str
    muted: str
    grid: str
    axis: str
    base: str
    highlight: str
    inverted: str
    # ORF strand, carried over from v1 where a forward ORF had a black outline
    # and a reverse one a red outline. 'Black' becomes near-white on the dark
    # surface: the meaning is 'forward', not 'black'.
    orf_forward: str
    orf_reverse: str

    @property
    def is_dark(self) -> bool:
        return self.name == "dark"


LIGHT = Theme(
    name="light",
    surface="#fcfcfb",
    paper="#f9f9f7",
    text_primary="#0b0b0b",
    text_secondary="#52514e",
    muted="#898781",
    grid="#e1e0d9",
    axis="#c3c2b7",
    base="#2a78d6",
    highlight="#eb6834",
    inverted="#1baf7a",
    orf_forward="#1c1c1a",
    orf_reverse="#d03b3b",
)

DARK = Theme(
    name="dark",
    surface="#1a1a19",
    paper="#0d0d0d",
    text_primary="#ffffff",
    text_secondary="#c3c2b7",
    muted="#898781",
    grid="#2c2c2a",
    axis="#383835",
    base="#3987e5",
    highlight="#d95926",
    inverted="#199e70",
    orf_forward="#e8e8e2",
    orf_reverse="#e8756a",
)

THEMES = {"light": LIGHT, "dark": DARK}

FONT_FAMILY = 'system-ui, -apple-system, "Segoe UI", sans-serif'


def get(name: str) -> Theme:
    try:
        return THEMES[name]
    except KeyError:
        raise ValueError(f"unknown theme {name!r}; expected 'light' or 'dark'") from None
