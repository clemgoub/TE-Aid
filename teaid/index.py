"""``index.html``: a contact sheet over the sheets a batch run wrote.

Grouped by TE class, in a grid, because comparing sheets across a library is
what a batch is for, and one full-width sheet per screen defeats it. Class
colours are v1's, so the grouping reads as "LTR" or "LINE" to a curator before
anything else does. A family that failed is a card carrying its slug and error,
so a partial batch reads as partial rather than silently short.

Nothing here is specific to one library: every number on the page is computed
from the batch's own results. ``tools/build_example_index.py`` builds the dev
set's page from the same class table.
"""

from __future__ import annotations

import html
from collections import defaultdict
from pathlib import Path

from . import tetypes

CLASS_ORDER = ["LTR", "LINE", "SINE", "TIR", "RC", "MAV", "PLE", "DIRS", "Satellite",
               "Unknown"]
CLASS_TITLE = {
    "LTR": "LTR retrotransposons",
    "LINE": "LINEs",
    "SINE": "SINEs",
    "TIR": "TIR (DNA) transposons",
    "RC": "Rolling-circle (Helitron)",
    "MAV": "Maverick / Polinton",
    "PLE": "Penelope-like elements",
    "DIRS": "DIRS / Crypton",
    "Satellite": "Satellite",
    "Unknown": "Unclassified and non-TE repeats",
}

_STYLE = """
:root { --ink:#14140f; --ink-2:#55544c; --ink-3:#8b8a80; --bg:#faf9f5; --card:#fff;
  --rule:#e3e1d7; --accent:#7a5c1e; --bad:#b3261e; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --ink:#f2f1ea; --ink-2:#b6b4a8; --ink-3:#86847a; --bg:#14140f; --card:#1c1c16;
  --rule:#33322a; --accent:#d8bd7d; --bad:#f2958c; } }
* { box-sizing: border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
  font:15px/1.55 ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif; }
.wrap { max-width:1220px; margin:0 auto; padding:40px 16px 80px; }
header { border-bottom:2px solid var(--ink); padding-bottom:16px; }
h1 { font-size:26px; line-height:1.2; margin:0 0 4px; }
.sub { color:var(--ink-2); margin:0; }
.stats { display:flex; gap:28px; flex-wrap:wrap; margin:20px 0 32px;
  font-variant-numeric:tabular-nums; }
.stat b { display:block; font-size:20px; }
.stat span { color:var(--ink-3); font-size:12px; text-transform:uppercase; letter-spacing:.08em; }
h2 { font-size:13px; text-transform:uppercase; letter-spacing:.1em; margin:40px 0 4px;
  display:flex; align-items:center; gap:10px; }
h2 .dot { width:11px; height:11px; border-radius:50%; flex:none; }
h2 .n { color:var(--ink-3); font-weight:400; letter-spacing:0; text-transform:none; }
.rule { height:1px; background:var(--rule); margin-bottom:20px; }
.grid { display:grid; gap:18px; grid-template-columns:repeat(auto-fill, minmax(300px, 1fr)); }
.card { background:var(--card); border:1px solid var(--rule); border-radius:10px;
  overflow:hidden; display:flex; flex-direction:column; }
.card figcaption { padding:10px 12px; display:flex; gap:9px; align-items:baseline;
  flex-wrap:wrap; }
.card.failed { border-color:var(--bad); }
.name { font-family:ui-monospace, "SF Mono", Menlo, monospace; font-size:13px;
  font-weight:600; word-break:break-all; }
.badge { font-size:10.5px; padding:2px 7px; border-radius:999px; color:#14140f;
  font-weight:600; }
.facts { flex-basis:100%; color:var(--ink-2); font-size:12px; font-variant-numeric:tabular-nums; }
.error { flex-basis:100%; color:var(--bad); font-family:ui-monospace, monospace;
  font-size:12px; word-break:break-word; }
.shot { display:block; background:#fff; line-height:0; border-top:1px solid var(--rule); }
.card img { display:block; width:100%; height:auto; }
a.open { color:var(--accent); text-decoration:none; font-size:12px;
  border-bottom:1px solid currentColor; margin-left:auto; }
a:focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
.note { color:var(--ink-2); font-size:13px; margin:0 0 8px; }
"""


def _klass(label: str) -> str:
    klass = tetypes.class_of_repeatmasker_label(label or "Unknown")
    return klass if klass in CLASS_TITLE else "Unknown"


def write(out_dir: Path, results: list[dict], *, title: str, note: str | None = None) -> Path:
    """Write ``out_dir/index.html`` from a batch's results; returns its path."""
    groups: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        groups[_klass((r.get("facts") or {}).get("class_label", ""))].append(r)

    drawn = [r for r in results if r["code"] == 0]
    lengths = [r["facts"]["length"] for r in drawn if r.get("facts", {}).get("length")]
    copies = [r["facts"]["copies"] for r in drawn if "copies" in r.get("facts", {})]
    seconds = sum(r["seconds"] for r in results)

    def span(values, fmt):
        return f"{fmt(min(values))} &ndash; {fmt(max(values))}" if values else "&ndash;"

    stats = [
        (f"{len(drawn):,}/{len(results):,}", "drawn"),
        (str(len([k for k in groups if any(r['code'] == 0 for r in groups[k])])), "TE classes"),
        (span(lengths, lambda v: f"{v:,} bp"), "consensus length"),
        (span(copies, lambda v: f"{v:,}"), "copies"),
        (f"{int(seconds) // 60}m {int(seconds) % 60}s", "total time"),
    ]
    parts = [
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{html.escape(title)}</title><style>{_STYLE}</style></head><body>"
        f'<div class="wrap"><header><h1>{html.escape(title)}</h1>'
        f'<p class="sub">One evidence sheet per family, grouped by TE class.</p></header>'
        '<div class="stats">'
        + "".join(f'<div class="stat"><b>{v}</b><span>{k}</span></div>' for v, k in stats)
        + "</div>"
    ]
    if note:
        parts.append(f'<p class="note">{html.escape(note)}</p>')

    for klass in CLASS_ORDER:
        members = sorted(groups.get(klass, []),
                         key=lambda r: (r.get("facts", {}).get("length", 0), r["family"]))
        if not members:
            continue
        colour = tetypes.colour(klass)
        parts.append(
            f'<h2><span class="dot" style="background:{colour}"></span>'
            f'{html.escape(CLASS_TITLE[klass])}<span class="n">{len(members)}</span></h2>'
            f'<div class="rule"></div><div class="grid">'
        )
        for r in members:
            parts.append(_card(out_dir, r, colour))
        parts.append("</div>")
    parts.append("</div></body></html>")
    page = out_dir / "index.html"
    page.write_text("".join(parts), encoding="utf-8")
    return page


def _card(out_dir: Path, r: dict, colour: str) -> str:
    facts = r.get("facts") or {}
    stem = r["stem"]
    failed = r["code"] != 0
    bits = [f'<figure class="card{" failed" if failed else ""}"><figcaption>',
            f'<span class="name">{html.escape(r["family"])}</span>']
    if facts.get("class_label"):
        bits.append(f'<span class="badge" style="background:{colour}">'
                    f'{html.escape(facts["class_label"])}</span>')
    if not failed:
        bits.append(f'<a class="open" href="{html.escape(stem)}.teaid.html">sheet &rarr;</a>')
        length = facts.get("length", 0)
        size = f"{length:,} bp" if length < 1000 else f"{length / 1000:.1f} kb"
        rows = facts.get("rows", 0)
        row_note = f" ({rows:,} rows)" if rows and rows != facts.get("copies") else ""
        bits.append(f'<span class="facts">{size} &middot; {facts.get("copies", 0):,} copies'
                    f'{row_note} &middot; {facts.get("full_length", 0):,} full length'
                    f' &middot; {r["seconds"]:.0f} s</span>')
    else:
        errors = [ln for ln in r["stderr"].splitlines() if "error" in ln]
        message = errors[-1] if errors else f"exit {r['code']}"
        bits.append(f'<span class="error">{html.escape(message[:300])}</span>')
    bits.append("</figcaption>")
    png = out_dir / f"{stem}.teaid.png"
    if not failed and png.exists():
        bits.append(f'<a class="shot" href="{html.escape(png.name)}">'
                    f'<img loading="lazy" src="{html.escape(png.name)}" '
                    f'alt="evidence sheet for {html.escape(r["family"])}"></a>')
    bits.append("</figure>")
    return "".join(bits)
