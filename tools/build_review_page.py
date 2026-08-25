#!/usr/bin/env python3
"""Generate the te_domains review page from the table itself.

The page is built from ``teaid/data/te_domains.tsv`` and the exclusion records in
``build_te_domains.py``, so its counts cannot drift from the data being reviewed.
"""

from __future__ import annotations

import collections
import html
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_te_domains import EXCLUDED, EXCLUDED_CATEGORIES  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
TSV = REPO / "teaid" / "data" / "te_domains.tsv"
IPR = "https://www.ebi.ac.uk/interpro/entry/pfam/{}/"

CLANS = [
    ("CL0219", "RNase_H", 121, 75, "host enzymes — Argonaute and Piwi, Cas9, RuvC, the "
     "spliceosome's PRP8, DNA polymerase subunits, exonucleases"),
    ("CL0027", "RdRP", 15, 13, "viral RNA-dependent RNA polymerases. Only PF00078 and "
     "PF07727 are reverse transcriptases"),
    ("CL0169", "Rep", 24, 22, "plasmid and virus replication proteins — geminivirus, "
     "parvovirus, papillomavirus, MobA relaxases"),
    ("CL0523", "GAG-polyprotein", 16, 6, "domesticated host genes built from gag — PEG10, "
     "Arc, PNMA, RTL1, LDOC"),
]

GAPS = [
    ("piggyBac", "DNA/PiggyBac", 152, "no Pfam model. InterPro covers it only through "
     "PANTHER and CDD"),
    ("Maverick / Polinton", "DNA/Maverick", 349, "no Pfam model, and no InterPro entry that "
     "is transposon-specific"),
    ("Crypton", "DNA/Cryp", 141, "no Pfam model; only the generic tyrosine-recombinase fold"),
    ("Penelope", "LINE/Penelope", 185, "only the generic GIY-YIG endonuclease PF01541; the "
     "PLE-specific model is CDD's, not Pfam's"),
]


def rows():
    lines = TSV.read_text(encoding="utf-8").splitlines()[1:]
    return [dict(zip(("acc", "name", "order", "prov", "why"), ln.split("\t"))) for ln in lines]


def e(text: str) -> str:
    return html.escape(str(text))


def build() -> str:
    data = rows()
    groups: dict[tuple[str, str], list[dict]] = collections.OrderedDict()
    for r in data:
        groups.setdefault((r["order"], r["why"]), []).append(r)

    order_counts = collections.Counter(r["order"] for r in data)
    n_orders = len(order_counts)

    group_html = []
    current_order = None
    for (order, why), members in groups.items():
        if order != current_order:
            group_html.append(
                f'<h3 class="order">{e(order)}'
                f'<span class="order-n">{order_counts[order]} '
                f'domain{"s" if order_counts[order] != 1 else ""}</span></h3>'
            )
            current_order = order
        chips = "".join(
            f'<a class="chip" href="{IPR.format(r["acc"])}" target="_blank" rel="noopener">'
            f'<code>{e(r["acc"])}</code><span>{e(r["name"])}</span></a>'
            for r in members
        )
        group_html.append(
            f'<div class="group"><p class="why">{e(why)}</p>'
            f'<div class="chips">{chips}</div></div>'
        )

    clan_rows = "".join(
        f"<tr><td><code>{e(a)}</code></td><td>{e(n)}</td>"
        f'<td class="num">{t}</td><td class="num bad">{x}</td><td>{e(note)}</td></tr>'
        for a, n, t, x, note in CLANS
    )
    gap_rows = "".join(
        f"<tr><td>{e(fam)}</td><td><code>{e(lab)}</code></td>"
        f'<td class="num">{n}</td><td>{e(note)}</td></tr>'
        for fam, lab, n, note in GAPS
    )
    excl_rows = "".join(
        f'<tr><td><a href="{IPR.format(a)}" target="_blank" rel="noopener">'
        f"<code>{e(a)}</code></a></td><td>{e(why)}</td></tr>"
        for a, why in sorted(EXCLUDED.items())
    )
    cat_rows = "".join(
        f"<div class='cat'><h4>{e(name)}</h4><p>{e(why)}</p></div>"
        for name, why in EXCLUDED_CATEGORIES.items()
    )

    return TEMPLATE.format(
        n_domains=len(data), n_orders=n_orders, n_excluded=len(EXCLUDED),
        n_categories=len(EXCLUDED_CATEGORIES),
        groups="\n".join(group_html), clan_rows=clan_rows, gap_rows=gap_rows,
        excl_rows=excl_rows, cat_rows=cat_rows,
    )


TEMPLATE = """<meta charset="utf-8">
<title>TE Domain Set Review</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap">
<style>
  :root {{
    color-scheme: light;
    --ground:  #f7f8fa;
    --surface: #ffffff;
    --sunk:    #eef1f6;
    --ink:     #14181f;
    --ink-2:   #454e5c;
    --muted:   #6a7484;
    --rule:    #dde2ea;
    --accent:  #2a78d6;
    --accent-soft: #e7f0fb;
    --bad:     #c0392f;
    --bad-soft:#fbeceb;
    --shadow:  0 1px 2px rgba(20,24,31,.05), 0 8px 24px -16px rgba(20,24,31,.28);
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      color-scheme: dark;
      --ground:  #0f1216;
      --surface: #171b21;
      --sunk:    #1d222a;
      --ink:     #eef1f5;
      --ink-2:   #b6bfcd;
      --muted:   #8892a2;
      --rule:    #2a313b;
      --accent:  #63a2ee;
      --accent-soft: #162435;
      --bad:     #e8756a;
      --bad-soft:#2e1b1a;
      --shadow:  0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
    }}
  }}
  :root[data-theme="dark"] {{
    color-scheme: dark;
    --ground:  #0f1216;
    --surface: #171b21;
    --sunk:    #1d222a;
    --ink:     #eef1f5;
    --ink-2:   #b6bfcd;
    --muted:   #8892a2;
    --rule:    #2a313b;
    --accent:  #63a2ee;
    --accent-soft: #162435;
    --bad:     #e8756a;
    --bad-soft:#2e1b1a;
    --shadow:  0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
  }}

  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    background: var(--ground);
    color: var(--ink);
    font-family: "IBM Plex Sans", system-ui, -apple-system, sans-serif;
    font-size: 16px;
    line-height: 1.6;
    -webkit-font-smoothing: antialiased;
  }}
  .wrap {{ max-width: 940px; margin: 0 auto; padding: 56px 24px 96px; }}

  h1, h2, h3, h4 {{ font-family: "IBM Plex Serif", Georgia, serif; text-wrap: balance; margin: 0; }}
  h1 {{ font-size: 2.1rem; font-weight: 600; letter-spacing: -.015em; line-height: 1.2; }}
  h2 {{ font-size: 1.35rem; font-weight: 600; letter-spacing: -.01em; }}
  h3 {{ font-size: 1.02rem; font-weight: 600; }}
  h4 {{ font-size: .95rem; font-weight: 600; }}
  p {{ margin: 0; max-width: 68ch; }}
  code {{ font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: .87em; }}
  a {{ color: var(--accent); }}
  a:focus-visible, .chip:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; border-radius: 3px; }}

  .eyebrow {{
    font-family: "IBM Plex Mono", monospace; font-size: .74rem; font-weight: 500;
    letter-spacing: .1em; text-transform: uppercase; color: var(--muted);
  }}
  header {{ display: flex; flex-direction: column; gap: 12px; padding-bottom: 28px; border-bottom: 1px solid var(--rule); }}
  .lede {{ color: var(--ink-2); font-size: 1.06rem; }}

  .stats {{ display: flex; flex-wrap: wrap; gap: 28px; margin-top: 6px; }}
  .stat b {{ display: block; font-family: "IBM Plex Serif", serif; font-size: 1.6rem; font-weight: 600; line-height: 1.1; font-variant-numeric: tabular-nums; }}
  .stat span {{ font-size: .8rem; color: var(--muted); }}

  section {{ display: flex; flex-direction: column; gap: 18px; padding-top: 44px; }}
  .section-head {{ display: flex; flex-direction: column; gap: 7px; }}

  .card {{
    background: var(--surface); border: 1px solid var(--rule); border-radius: 8px;
    padding: 22px 24px; box-shadow: var(--shadow);
    display: flex; flex-direction: column; gap: 14px;
  }}
  .ask {{ border-left: 3px solid var(--accent); }}
  .ask h3 {{ display: flex; align-items: baseline; gap: 10px; }}
  .ask .tag {{
    font-family: "IBM Plex Mono", monospace; font-size: .7rem; letter-spacing: .08em;
    text-transform: uppercase; color: var(--accent); background: var(--accent-soft);
    padding: 3px 8px; border-radius: 4px; white-space: nowrap;
  }}
  .rec {{ background: var(--sunk); border-radius: 6px; padding: 12px 14px; font-size: .93rem; }}
  .rec b {{ color: var(--ink); }}

  table {{ width: 100%; border-collapse: collapse; font-size: .9rem; }}
  .scroll {{ overflow-x: auto; }}
  th {{
    text-align: left; font-family: "IBM Plex Mono", monospace; font-weight: 500;
    font-size: .72rem; letter-spacing: .07em; text-transform: uppercase;
    color: var(--muted); border-bottom: 1px solid var(--rule); padding: 0 12px 8px 0;
  }}
  td {{ padding: 9px 12px 9px 0; border-bottom: 1px solid var(--rule); vertical-align: top; color: var(--ink-2); }}
  td:first-child, th:first-child {{ padding-left: 0; }}
  tr:last-child td {{ border-bottom: none; }}
  .num {{ font-variant-numeric: tabular-nums; text-align: right; white-space: nowrap; color: var(--ink); }}
  .bad {{ color: var(--bad); font-weight: 500; }}

  .order {{ display: flex; align-items: baseline; gap: 12px; margin-top: 26px; padding-bottom: 7px; border-bottom: 1px solid var(--rule); }}
  .order:first-of-type {{ margin-top: 0; }}
  .order-n {{ font-family: "IBM Plex Mono", monospace; font-size: .74rem; color: var(--muted); font-weight: 400; }}
  .group {{ display: flex; flex-direction: column; gap: 10px; padding: 14px 0 4px; }}
  .why {{ color: var(--muted); font-size: .89rem; max-width: 74ch; }}
  .chips {{ display: flex; flex-wrap: wrap; gap: 7px; }}
  .chip {{
    display: inline-flex; align-items: baseline; gap: 7px; text-decoration: none;
    background: var(--surface); border: 1px solid var(--rule); border-radius: 5px;
    padding: 5px 10px; color: var(--ink-2); font-size: .84rem; transition: border-color .12s, color .12s;
  }}
  .chip code {{ color: var(--accent); font-weight: 500; }}
  .chip:hover {{ border-color: var(--accent); color: var(--ink); }}

  .cat {{ display: flex; flex-direction: column; gap: 5px; padding: 13px 0; border-bottom: 1px solid var(--rule); }}
  .cat:last-child {{ border-bottom: none; }}
  .cat p {{ font-size: .89rem; color: var(--muted); }}

  .excl-table td:first-child {{ width: 108px; white-space: nowrap; }}
  .excl-table code {{ color: var(--bad); }}

  footer {{ margin-top: 56px; padding-top: 22px; border-top: 1px solid var(--rule); color: var(--muted); font-size: .85rem; display: flex; flex-direction: column; gap: 8px; }}

  @media (max-width: 620px) {{
    .wrap {{ padding: 36px 18px 64px; }}
    h1 {{ font-size: 1.7rem; }}
    .stats {{ gap: 18px; }}
  }}
  @media (prefers-reduced-motion: reduce) {{ * {{ transition: none !important; }} }}
</style>

<div class="wrap">
  <header>
    <span class="eyebrow">TE-Aid v2 · work order step 6 · draft for review</span>
    <h1>The curated Pfam TE-domain set</h1>
    <p class="lede">A draft of <code>teaid/data/te_domains.tsv</code>, the protein query set behind
      panel 5. Nothing here is wired into the search path yet — §5.4 asks for your eyes first,
      and two findings below need a decision before the code gets written.</p>
    <div class="stats">
      <div class="stat"><b>{n_domains}</b><span>domains included</span></div>
      <div class="stat"><b>{n_orders}</b><span>TE orders covered</span></div>
      <div class="stat"><b>{n_excluded}</b><span>explicit exclusions</span></div>
      <div class="stat"><b>{n_categories}</b><span>excluded categories</span></div>
    </div>
  </header>

  <section>
    <div class="section-head">
      <span class="eyebrow">What I need from you</span>
      <h2>Three things</h2>
    </div>

    <div class="card ask">
      <h3><span class="tag">review</span> The host-enzyme exclusions</h3>
      <p>The list you asked to see is <a href="#exclusions">further down</a>: 29 accessions
        excluded one by one with the reason, plus 3 whole categories. The CL0219 calls are the
        ones you flagged. The judgement I applied throughout: <em>would a hit here, on a TE
        consensus, be evidence of transposable-element origin?</em></p>
      <div class="rec"><b>Most likely to be contentious:</b> I excluded cellular
        <code>PF00075</code> RNase H but included the RT-specific <code>PF17917</code> and
        <code>PF17919</code> instead, so retroelement RNase H is still covered without the host
        enzyme coming with it. I also excluded <code>PF00098</code> zinc knuckle (real in gag,
        but far commoner in host RNA-binding proteins) and <code>PF00136</code> DNA polymerase B
        (Mavericks carry one, so does every replisome).</div>
    </div>

    <div class="card ask">
      <h3><span class="tag">decide</span> RepeatPeps as pHMMs is ~6.4&nbsp;GB</h3>
      <p>§5.1 plans one single-sequence pHMM per RepeatPeps protein. I started that build and
        stopped it: at 2,022 of 18,011 sequences it had written 718&nbsp;MB, which extrapolates
        to <b>~6.4&nbsp;GB</b> — this machine has 7&nbsp;GB free, so it would have filled the disk.</p>
      <p>The obvious escape hatch does not work. Clustering RepeatPeps with cd-hit at 90%
        identity removes <b>0.5%</b> of it (18,011&nbsp;→&nbsp;17,916); at 80% it removes 3%.
        The library is already non-redundant, so the size is intrinsic: 16.2&nbsp;M residues at
        ~396 bytes each.</p>
      <div class="rec"><b>What I'd suggest, but it is your call:</b> build pHMMs only for the
        superfamilies Pfam cannot model at all (next section) — 870 proteins, <b>~200&nbsp;MB</b>
        instead of 6.4&nbsp;GB — and search the rest of RepeatPeps with <code>blastp</code> as v1
        did. That keeps frameshift-aware search where the conserved catalytic domains are, and
        keeps piggyBac and Maverick covered, without shipping a 6&nbsp;GB model. The alternative
        is to accept the full library as an opt-in <code>--deep</code> arm and let the §5.5
        benchmark settle it.</div>
    </div>

    <div class="card ask">
      <h3><span class="tag">decide</span> Whether prokaryotic IS elements belong</h3>
      <p>I left IS transposases, phage integrases and conjugative-transposon proteins out. They
        are genuine transposases, but not of the elements TE-Aid curates. Including them would
        let bacterial contamination in a eukaryotic assembly show up as a protein hit — arguably
        useful, but a different question.</p>
      <div class="rec"><b>Default taken:</b> excluded. Say the word and they go in behind a
        flag rather than in the default set.</div>
    </div>
  </section>

  <section>
    <div class="section-head">
      <span class="eyebrow">Evidence</span>
      <h2>Why this is an accession list and not a clan selection</h2>
      <p>Your §5.4 argument holds, and it holds for four clans rather than one. Counts are from
        InterPro, fetched at build time.</p>
    </div>
    <div class="card scroll">
      <table>
        <thead><tr><th>Clan</th><th>Name</th><th class="num">Members</th><th class="num">Not TE</th><th>What the rest of it is</th></tr></thead>
        <tbody>{clan_rows}</tbody>
      </table>
    </div>
    <p class="why">And in the other direction: the hAT C-terminal dimerisation region
      (<code>PF05699</code>) and the Helitron-associated <code>PF13952</code> sit in no clan at
      all, so a clan-based selection misses them however it is filtered.</p>
  </section>

  <section>
    <div class="section-head">
      <span class="eyebrow">Evidence</span>
      <h2>Where Pfam has no model, and RepeatPeps is the only cover</h2>
      <p>This is the concrete answer to why the query set needs both halves — and it is what
        makes dropping RepeatPeps outright a bad trade.</p>
    </div>
    <div class="card scroll">
      <table>
        <thead><tr><th>Superfamily</th><th>RepeatPeps label</th><th class="num">Proteins</th><th>Pfam coverage</th></tr></thead>
        <tbody>{gap_rows}</tbody>
      </table>
    </div>
  </section>

  <section id="exclusions">
    <div class="section-head">
      <span class="eyebrow">For review</span>
      <h2>What was excluded, and why</h2>
      <p>The judgement call in each case is whether a hit would be evidence of TE origin or
        would manufacture it.</p>
    </div>
    <div class="card scroll">
      <table class="excl-table">
        <thead><tr><th>Accession</th><th>Reason for exclusion</th></tr></thead>
        <tbody>{excl_rows}</tbody>
      </table>
    </div>
    <div class="card">{cat_rows}</div>
  </section>

  <section>
    <div class="section-head">
      <span class="eyebrow">The table</span>
      <h2>{n_domains} domains, grouped by what justifies them</h2>
      <p>Grouped the way it was curated: each block is one rationale, and every accession in it
        was included for that reason. Names come from InterPro at build time, so a typo or a dead
        accession fails the build rather than shipping. Accessions link out to InterPro.</p>
    </div>
    <div class="card">{groups}</div>
  </section>

  <footer>
    <span>Generated from <code>teaid/data/te_domains.tsv</code> by <code>tools/build_review_page.py</code> — counts cannot drift from the data.</span>
    <span>TE-Aid v2 · branch <code>v2</code> · domain names and clan membership from InterPro/Pfam</span>
  </footer>
</div>
"""


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "te_domains_review.html"
    out.write_text(build(), encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
