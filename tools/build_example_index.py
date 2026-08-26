#!/usr/bin/env python3
"""Build a contact sheet over rendered example sheets.

Companion to ``tools/select_families.py``: that picks the families, this makes
the 50 resulting sheets reviewable in one scroll. A grid, not a column --
comparing designs across TE types is the whole point, and one full-width sheet
per screen defeats it.

Grouped by TE class, with v1's own class colours as the section accent, so the
grouping is carried by position *and* by a colour curators already read as
"LTR" or "LINE" rather than by an arbitrary palette.

Expects two TSVs:

- ``--selection`` -- ``tools/select_families.py`` output (name, label, class,
  length, copies, why)
- ``--results``   -- one row per render: name, exit code, seconds, stdout+stderr
  with newlines collapsed to ``|``

Missing renders are drawn as a failure card carrying the captured output, so a
partial run reads as partial rather than silently short.
"""
from __future__ import annotations

import argparse
import html
from collections import defaultdict
from pathlib import Path

from teaid import tetypes

CLASS_ORDER = ["LTR", "LINE", "SINE", "TIR", "RC", "MAV", "Satellite", "Unknown"]
CLASS_TITLE = {
    "LTR": "LTR retrotransposons",
    "LINE": "LINEs",
    "SINE": "SINEs",
    "TIR": "TIR (DNA) transposons",
    "RC": "Rolling-circle (Helitron)",
    "MAV": "Maverick / Polinton",
    "Satellite": "Satellite",
    "Unknown": "Unclassified and non-TE repeats",
}


def build(selection: str, results: str, out_dir: Path) -> None:
    OUT = out_dir
    meta = {}
    for line in Path(selection).read_text().splitlines():
        name, label, klass, length, copies, why = line.split("\t")
        meta[name] = dict(label=label, klass=klass, length=int(length),
                          copies=int(copies), why=why)

    status = {}
    for line in Path(results).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        name, code, secs, out = parts[0], parts[1], parts[2], parts[3]
        status[name] = dict(code=int(code), secs=int(secs), out=out)

    groups = defaultdict(list)
    for name, m in meta.items():
        groups[m["klass"]].append(name)

    ok = sum(1 for s in status.values() if s["code"] == 0)
    total_s = sum(s["secs"] for s in status.values())

    parts = []
    parts.append(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TE-Aid v2 - 50 example sheets</title>
<style>
:root {{
  --ink: #14140f; --ink-2: #55544c; --ink-3: #8b8a80;
  --bg: #faf9f5; --card: #ffffff; --rule: #e3e1d7; --accent: #7a5c1e;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --ink: #f2f1ea; --ink-2: #b6b4a8; --ink-3: #86847a;
    --bg: #14140f; --card: #1c1c16; --rule: #33322a; --accent: #d8bd7d;
  }}
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0; background: var(--bg); color: var(--ink);
  font: 15px/1.55 ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
}}
.wrap {{ max-width: 1220px; margin: 0 auto; padding: 48px 24px 96px; }}
header {{ border-bottom: 2px solid var(--ink); padding-bottom: 20px; margin-bottom: 8px; }}
h1 {{ font-size: 30px; line-height: 1.15; margin: 0 0 6px; letter-spacing: -.02em; }}
.sub {{ color: var(--ink-2); margin: 0; }}
.stats {{ display: flex; gap: 28px; flex-wrap: wrap; margin: 22px 0 40px;
  font-variant-numeric: tabular-nums; }}
.stat b {{ display: block; font-size: 22px; letter-spacing: -.01em; }}
.stat span {{ color: var(--ink-3); font-size: 12px; text-transform: uppercase;
  letter-spacing: .08em; }}
h2 {{ font-size: 13px; text-transform: uppercase; letter-spacing: .1em;
  margin: 44px 0 4px; display: flex; align-items: center; gap: 10px; }}
h2 .dot {{ width: 11px; height: 11px; border-radius: 50%; flex: none; }}
h2 .n {{ color: var(--ink-3); font-weight: 400; letter-spacing: 0; text-transform: none; }}
.rule {{ height: 1px; background: var(--rule); margin-bottom: 22px; }}
.grid {{ display: grid; gap: 20px;
  grid-template-columns: repeat(auto-fill, minmax(330px, 1fr)); }}
.card {{ background: var(--card); border: 1px solid var(--rule); border-radius: 10px;
  overflow: hidden; display: flex; flex-direction: column; }}
.card figcaption {{ padding: 11px 13px; border-bottom: 1px solid var(--rule);
  display: flex; gap: 9px; align-items: baseline; flex-wrap: wrap; }}
.name {{ font-family: ui-monospace, "SF Mono", Menlo, monospace; font-size: 13px;
  font-weight: 600; }}
.badge {{ font-size: 10.5px; padding: 2px 7px; border-radius: 999px; color: #14140f;
  font-weight: 600; }}
.facts {{ flex-basis: 100%; color: var(--ink-2); font-size: 12px;
  font-variant-numeric: tabular-nums; }}
.why {{ flex-basis: 100%; color: var(--ink-3); font-size: 11.5px; font-style: italic;
  line-height: 1.4; }}
/* The sheet is the point, so it gets the room: click to open it full size. */
.shot {{ display: block; background: #fff; line-height: 0; }}
.card img {{ display: block; width: 100%; height: auto; }}
.fail {{ padding: 14px 13px; color: #e06c62; font-family: ui-monospace, monospace;
  font-size: 12px; word-break: break-word; }}
a.open {{ color: var(--accent); text-decoration: none; font-size: 12px;
  border-bottom: 1px solid currentColor; margin-left: auto; }}
a.shot:focus-visible, a.open:focus-visible {{ outline: 2px solid var(--accent);
  outline-offset: 2px; }}
@media (prefers-reduced-motion: no-preference) {{
  .card {{ transition: border-color .15s ease; }}
  .card:hover {{ border-color: var(--ink-3); }}
}}
footer {{ margin-top: 64px; padding-top: 18px; border-top: 1px solid var(--rule);
  color: var(--ink-3); font-size: 13px; }}
</style></head><body><div class="wrap">
<header>
  <h1>Fifty families, one sheet each</h1>
  <p class="sub">A stratified sample of the RepeatModeler2 library for
  <code>GCA_963082875.1</code> &mdash; every class label the library contains,
  across its full range of consensus length and copy number.</p>
</header>
<div class="stats">
  <div class="stat"><b>{ok}/{len(meta)}</b><span>rendered</span></div>
  <div class="stat"><b>8</b><span>TE classes</span></div>
  <div class="stat"><b>267 bp &ndash; 11.6 kb</b><span>consensus length</span></div>
  <div class="stat"><b>12 &ndash; 10,851</b><span>copies</span></div>
  <div class="stat"><b>{total_s // 60}m {total_s % 60}s</b><span>total CPU</span></div>
</div>
""")

    for klass in CLASS_ORDER:
        names = groups.get(klass, [])
        if not names:
            continue
        names.sort(key=lambda n: meta[n]["length"])
        colour = tetypes.colour(klass)
        parts.append(
            f'<h2><span class="dot" style="background:{colour}"></span>'
            f'{html.escape(CLASS_TITLE[klass])}<span class="n">{len(names)}</span></h2>'
            f'<div class="rule"></div><div class="grid">'
        )
        for name in names:
            m = meta[name]
            st = status.get(name)
            png = OUT / f"{name}.teaid.png"
            kb = f"{m['length']:,} bp" if m["length"] < 1000 else f"{m['length'] / 1000:.1f} kb"
            parts.append('<figure class="card"><figcaption>')
            parts.append(f'<span class="name">{html.escape(name)}</span>')
            parts.append(
                f'<span class="badge" style="background:{colour}">'
                f'{html.escape(m["label"])}</span>'
            )
            parts.append(
                f'<a class="open" href="{name}.teaid.html">interactive &rarr;</a>'
            )
            secs = f" &middot; {st['secs']}s" if st else ""
            parts.append(
                f'<span class="facts">{kb} &middot; {m["copies"]:,} copies{secs}</span>'
            )
            parts.append(f'<span class="why">{html.escape(m["why"])}</span>')
            parts.append("</figcaption>")
            if png.exists():
                parts.append(
                    f'<a class="shot" href="{png.name}" title="open full size">'
                    f'<img loading="lazy" src="{png.name}" '
                    f'alt="evidence sheet for {html.escape(name)}"></a>'
                )
            else:
                msg = st["out"].replace("|", " ") if st else "not rendered"
                parts.append(f'<div class="fail">no static render &mdash; {html.escape(msg[:400])}</div>')
            parts.append("</figure>")
        parts.append("</div>")

    parts.append(
        '<footer>Rendered with the annotation route '
        '(<code>--annot</code> + <code>--consensus</code>), BATH and EMBOSS on '
        'PATH, against the prebuilt tier-1+2 protein library. Click a sheet to '
        'open it full size, or <b>interactive &rarr;</b> for the zoomable HTML. '
        'Regenerate the family list with '
        '<code>tools/select_families.py</code>.</footer>'
        "</div></body></html>"
    )
    (OUT / "index.html").write_text("".join(parts), encoding="utf-8")
    print(f"wrote {OUT / 'index.html'} ({ok}/{len(meta)} rendered)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--selection", required=True,
                        help="tools/select_families.py output")
    parser.add_argument("--results", required=True,
                        help="TSV: name, exit code, seconds, captured output")
    parser.add_argument("--out", default="dev-data/examples",
                        help="directory holding the rendered sheets (default: %(default)s)")
    args = parser.parse_args()
    build(args.selection, args.results, Path(args.out))


if __name__ == "__main__":
    main()
