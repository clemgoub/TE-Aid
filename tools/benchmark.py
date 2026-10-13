#!/usr/bin/env python3
"""The protein benchmark of bench/PROTOCOL.md: run arms, compare them, report, verify.

    benchmark.py run --arm default --set dmel dev-data/dfam/dmel.fa --out bench/results/RUN
    benchmark.py shuffle dev-data/dfam/dmel.fa --n 5 --seed 7 --out bench/shuffled.fa
    benchmark.py compare bench/results/RUN --set dmel --arm default
    benchmark.py report bench/results/RUN
    benchmark.py verify bench/results/RUN [--repeat bench/results/RUN2]

Arms (bench/PROTOCOL.md, "Arms"):

- ``default``: bathsearch --fs over the tier 1-2 library, the sheet's search;
- ``nofs``: the same without --fs;
- ``pfam-orf``: hmmsearch of the 130 tier-1 models on getorf peptides;
- ``tier3``: blastp of getorf peptides against RepeatPeps;
- ``deep``: bathsearch --fs over all of RepeatPeps as single-sequence models.

Every arm searches at E <= 1 and keeps every hit, so one run serves every
threshold. A run writes ``<out>/<set>/<arm>/hits.tsv`` and ``timing.tsv``, and
records its inputs and tools in ``<out>/manifest.json``.

Writing under bench/results/ is refused until bench/PROTOCOL.md is committed,
unmodified and signed: the protocol says no result may exist before the commit
that signs it, and this is where that is enforced. Runs anywhere else (a
scratch directory, the goby examples) are allowed, for testing the tool.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from teaid import classcheck, homology, orfs, proteins  # noqa: E402

PROTOCOL = ROOT / "bench" / "PROTOCOL.md"
RESULTS = ROOT / "bench" / "results"
UNSIGNED_MARK = "draft, not yet signed"
ARMS = ("default", "nofs", "pfam-orf", "tier3", "deep")
SEARCH_E = 1.0
PRIMARY_E = 1e-3
E_GRID = (1e-3, 1e-5, 1e-10)
FLOOR_GRID = tuple(i / 10 for i in range(10))
OVERLAP = 0.5  # homology.best_per_region's rule, reused for matching regions
# The four superfamilies tier 2 exists for, by the subtype that names them.
PFAM_BLIND = {"piggyBac": ("piggybac",), "Maverick": ("maverick", "polinton"),
              "Crypton": ("crypton",), "Penelope": ("penelope",)}

# The draft's numbers. They are parameters until the protocol is signed, and
# report takes the signed ones from --thresholds.
DRAFT_THRESHOLDS = {
    "contradiction_ceiling": 0.05, "agreeing_kept": 0.90,
    "a1_lower_bound": 0.85, "a2_median_span": 0.80,
    "a2_high_span": 0.90, "a2_high_share": 0.60, "a4_cost_ratio": 2.0,
}

HIT_FIELDS = ("family", "label", "model", "model_class", "accession", "source", "start", "end", "strand",
              "evalue", "score", "hmm_from", "hmm_to", "hmm_len", "coverage",
              "frameshifts", "stop_codons")


# --- Inputs ---------------------------------------------------------------

def read_families(path: Path) -> list[tuple[str, str, str]]:
    """(name, label, sequence) per record of a ``>name#label`` FASTA."""
    out, name, label, chunks = [], None, "", []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, label, "".join(chunks)))
            head = line[1:].split()[0]
            name, _, label = head.partition("#")
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        out.append((name, label, "".join(chunks)))
    return out


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    """A shuffle that keeps every dinucleotide count (Altschul and Erickson, 1985).

    Walks a random Eulerian path through the graph whose edges are the
    sequence's dinucleotides. The last edge out of each letter is drawn so that
    together they form a tree into the final letter (here by Wilson's
    algorithm), the other edges are shuffled, and the walk takes them in order.
    """
    seq = seq.upper()
    if len(seq) < 3:
        return seq
    edges: dict[str, list[str]] = {}
    for a, b in zip(seq, seq[1:]):
        edges.setdefault(a, []).append(b)
    last = seq[-1]
    while True:
        # Wilson's algorithm on the multigraph: a random arborescence into `last`.
        last_edge: dict[str, str] = {}
        in_tree = {last}
        for start in edges:
            u = start
            path: dict[str, str] = {}
            while u not in in_tree:
                path[u] = rng.choice(edges[u])
                u = path[u]
            u = start
            while u not in in_tree:
                last_edge[u] = path[u]
                in_tree.add(u)
                u = path[u]
        order: dict[str, list[str]] = {}
        for a, outs in edges.items():
            rest = list(outs)
            if a in last_edge:
                rest.remove(last_edge[a])
            rng.shuffle(rest)
            order[a] = rest + ([last_edge[a]] if a in last_edge else [])
        walk, u = [seq[0]], seq[0]
        cursor = {a: 0 for a in order}
        while cursor.get(u, 0) < len(order.get(u, [])):
            v = order[u][cursor[u]]
            cursor[u] += 1
            walk.append(v)
            u = v
        if len(walk) == len(seq):
            return "".join(walk)


# --- Arms -----------------------------------------------------------------

def _repeatpeps() -> Path | None:
    return proteins.find_repeatpeps(ROOT / "dev-data" / "RepeatPeps.lib")


def _pfam_orf(seq: str, name: str, library, cpus: int, min_orf: int) -> list[homology.ProteinHit]:
    """hmmsearch of the tier-1 models on getorf peptides, in consensus coordinates."""
    pfam = Path(library.hmm).parent / "pfam.hmm"
    found = orfs.find_orfs(seq, name, min_size=min_orf)
    if not found:
        return []
    with tempfile.TemporaryDirectory(prefix="teaid-bench-") as tmp:
        peptides = Path(tmp) / "orfs.fa"
        peptides.write_text("".join(f">orf{i}\n{o.peptide}\n" for i, o in enumerate(found)))
        table = Path(tmp) / "dom.tbl"
        subprocess.run(["hmmsearch", "--domtblout", str(table), "-E", str(SEARCH_E),
                        "--domE", str(SEARCH_E), "--cpu", str(cpus), "-o", "/dev/null",
                        str(pfam), str(peptides)], check=True, capture_output=True)
        rows = [l.split() for l in table.read_text().splitlines() if l and not l.startswith("#")]
    hits = []
    for f in rows:
        orf = found[int(f[0][3:])]
        ali_from, ali_to = int(f[17]), int(f[18])
        if orf.strand == "+":
            start, end = orf.start + (ali_from - 1) * 3, orf.start + ali_to * 3
        else:
            start, end = orf.end - ali_to * 3, orf.end - (ali_from - 1) * 3
        hits.append(homology.ProteinHit(
            query=f[3], query_accession=f[4], start=start, end=end, strand=orf.strand,
            evalue=float(f[12]), score=float(f[13]), identity=0.0, hmm_from=int(f[15]),
            hmm_to=int(f[16]), hmm_len=int(f[5]), frameshifts=0, stop_codons=0,
            source="pfam-orf", target=name))
    return hits


def run_arm(arm: str, families, cpus: int, min_orf: int):
    """Yield (name, label, hits, seconds) for every family, one arm."""
    library = None
    if arm in ("default", "nofs", "pfam-orf"):
        library = proteins.build(repeatpeps=_repeatpeps(), log=lambda m: print(m, file=sys.stderr))
    elif arm == "deep":
        library = proteins.build(deep=True, repeatpeps=_repeatpeps(),
                                 log=lambda m: print(m, file=sys.stderr))
    peps = _repeatpeps()
    for name, label, seq in families:
        t0 = time.perf_counter()
        if arm in ("default", "deep", "nofs"):
            hits = homology.search(seq, library, name, evalue=SEARCH_E,
                                   frameshift_aware=arm != "nofs", cpus=cpus)
        elif arm == "pfam-orf":
            hits = _pfam_orf(seq, name, library, cpus, min_orf)
        elif arm == "tier3":
            hits = homology.search_repeatpeps(orfs.find_orfs(seq, name, min_size=min_orf),
                                              peps, evalue=SEARCH_E, target=name) if peps else []
        else:
            raise ValueError(f"unknown arm {arm!r}")
        yield name, label, hits, time.perf_counter() - t0


def tool_versions() -> dict:
    def first(cmd, line=0):
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            return (out.stdout or out.stderr).splitlines()[line].strip()
        except Exception as exc:  # recorded, not fatal: the manifest says what was missing
            return f"unavailable: {exc}"
    return {"bathsearch": first(["bathsearch", "-h"], 1), "hmmsearch": first(["hmmsearch", "-h"], 1),
            "blastp": first(["blastp", "-version"]), "getorf": first(["getorf", "-version"])}


def git(*args) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True).stdout.strip()


def protocol_state() -> dict:
    """Whether bench/PROTOCOL.md is committed, unmodified and signed, and when."""
    committed = git("log", "-1", "--format=%H %ct", "--", str(PROTOCOL.relative_to(ROOT)))
    dirty = bool(git("status", "--porcelain", "--", str(PROTOCOL.relative_to(ROOT))))
    signed = PROTOCOL.exists() and UNSIGNED_MARK not in PROTOCOL.read_text()
    commit, _, when = committed.partition(" ")
    return {"commit": commit or None, "committed_at": int(when) if when else None,
            "modified": dirty, "signed": signed}


def guard(out: Path) -> None:
    """Refuse to write a result under bench/results/ before the protocol is signed."""
    try:
        out.resolve().relative_to(RESULTS.resolve())
    except ValueError:
        return
    state = protocol_state()
    problems = [msg for bad, msg in (
        (not state["commit"], "bench/PROTOCOL.md is not committed"),
        (state["modified"], "bench/PROTOCOL.md has uncommitted changes"),
        (not state["signed"], f"bench/PROTOCOL.md still says '{UNSIGNED_MARK}'"),
    ) if bad]
    if problems:
        raise SystemExit("refusing to write under bench/results/: " + "; ".join(problems)
                         + ". Use a directory outside bench/results/ to test the tool.")


def cmd_run(args) -> None:
    out = Path(args.out)
    guard(out)
    set_name, fasta = args.set
    fasta = Path(fasta)
    families = read_families(fasta)
    if args.limit:
        families = families[: args.limit]
    arm_dir = out / set_name / args.arm
    arm_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "teaid_commit": git("rev-parse", "HEAD"), "protocol": protocol_state(),
        "tools": tool_versions(), "search_evalue": SEARCH_E, "sets": {}, "arms": {}}
    manifest["sets"][set_name] = {"fasta": str(fasta), "sha256": sha256(fasta),
                                  "families": len(families), "limit": args.limit}
    with open(arm_dir / "hits.tsv", "w") as hits_out, open(arm_dir / "timing.tsv", "w") as t_out:
        hits_out.write("\t".join(HIT_FIELDS) + "\n")
        t_out.write("family\tseconds\thits\n")
        for name, label, hits, seconds in run_arm(args.arm, families, args.cpus, args.min_orf):
            for h in hits:
                row = dict(family=name, label=label, model=h.display_name,
                           model_class=h.source_class or "", accession=h.query_accession,
                           source=h.source, start=h.start, end=h.end, strand=h.strand,
                           evalue=f"{h.evalue:.3g}", score=f"{h.score:.1f}", hmm_from=h.hmm_from,
                           hmm_to=h.hmm_to, hmm_len=h.hmm_len, coverage=f"{h.coverage:.4f}",
                           frameshifts=h.frameshifts, stop_codons=h.stop_codons)
                hits_out.write("\t".join(str(row[k]) for k in HIT_FIELDS) + "\n")
            t_out.write(f"{name}\t{seconds:.3f}\t{len(hits)}\n")
            print(f"{set_name}/{args.arm}: {name} {len(hits)} hits {seconds:.1f} s", file=sys.stderr)
    manifest["arms"].setdefault(set_name, {})[args.arm] = {
        "cpus": args.cpus, "min_orf": args.min_orf,
        "finished": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


def cmd_shuffle(args) -> None:
    rng = random.Random(args.seed)
    families = read_families(Path(args.fasta))
    picked = rng.sample(families, min(args.n, len(families)))
    with open(args.out, "w") as out:
        for name, label, seq in picked:
            out.write(f">{name}_shuf{args.seed}#{label}\n{dinucleotide_shuffle(seq, rng)}\n")
    print(f"wrote {len(picked)} shuffled sequences (seed {args.seed}) to {args.out}", file=sys.stderr)


# --- Analysis -------------------------------------------------------------

@dataclass(frozen=True)
class Hit:
    """One row of hits.tsv, with the fields homology's filters and classcheck read."""

    family: str
    label: str
    model: str
    model_class: str | None  # a RepeatPeps protein's class label; None for Pfam domains
    accession: str
    source: str
    start: int
    end: int
    evalue: float
    score: float
    coverage: float

    @property
    def query_accession(self) -> str:
        return self.accession

    @property
    def source_class(self) -> str | None:
        return self.model_class


def read_hits(path: Path) -> dict[str, list[Hit]]:
    """Hits per family, in file order, plus every family that had none (timing.tsv)."""
    by_family: dict[str, list[Hit]] = {}
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    for line in lines[1:]:
        r = dict(zip(head, line.split("\t")))
        by_family.setdefault(r["family"], []).append(Hit(
            r["family"], r["label"], r["model"], r.get("model_class") or None, r["accession"],
            r["source"], int(r["start"]), int(r["end"]), float(r["evalue"]), float(r["score"]),
            float(r["coverage"])))
    timing = path.parent / "timing.tsv"
    if timing.exists():
        for line in timing.read_text().splitlines()[1:]:
            by_family.setdefault(line.split("\t")[0], [])
    return by_family


def drawn(hits: list[Hit], evalue: float, floor: float) -> list[Hit]:
    """What the sheet would draw for one family: E filter, coverage floor, best per region."""
    kept = [h for h in hits if h.evalue <= evalue]
    kept, _ = homology.apply_coverage_floor(kept, floor)
    return homology.best_per_region(kept, OVERLAP)


def shared(a_start, a_end, b_start, b_end) -> float:
    shorter = min(a_end - a_start, b_end - b_start)
    return max(0, min(a_end, b_end) - max(a_start, b_start)) / shorter if shorter > 0 else 0.0


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def compare(ref: dict[str, list[Hit]], arm: dict[str, list[Hit]], evalue: float = PRIMARY_E) -> list[dict]:
    """One row per reference region: found, span and pieces, as bench/PROTOCOL.md defines them."""
    rows = []
    for family, ref_hits in ref.items():
        regions = drawn(ref_hits, evalue, 0.0)
        candidates = drawn(arm.get(family, []), evalue, 0.0)
        for region in sorted(regions, key=lambda h: h.start):
            matches = [h for h in candidates
                       if shared(region.start, region.end, h.start, h.end) >= OVERLAP]
            covered, cursor = 0, region.start
            for h in sorted(matches, key=lambda h: h.start):
                lo, hi = max(h.start, cursor), min(h.end, region.end)
                if hi > lo:
                    covered += hi - lo
                    cursor = hi
            rows.append({"family": family, "label": region.label, "start": region.start,
                         "end": region.end, "model": region.model, "found": bool(matches),
                         "span": covered / (region.end - region.start), "pieces": len(matches)})
    return rows


def judged(hits: list[Hit], orders) -> tuple[int, int]:
    """(agreeing, contradicting) hits, judged against their family's label."""
    agree = contra = 0
    for h in hits:
        verdict = classcheck.judge(h.label, h, orders).verdict
        agree += verdict == "agree"
        contra += verdict == "contradict"
    return agree, contra


def specificity(hits: dict[str, list[Hit]], th: dict) -> dict:
    """The Question 1 grid for one set, and the floor its decision rule picks per E-value."""
    orders = classcheck.domain_orders()
    grid, chosen = [], {}
    for evalue in E_GRID:
        base = None
        for floor in FLOOR_GRID:
            a = c = n = 0
            for fam_hits in hits.values():
                d = drawn(fam_hits, evalue, floor)
                n += len(d)
                fa, fc = judged(d, orders)
                a, c = a + fa, c + fc
            base = a if base is None else base
            row = {"evalue": evalue, "floor": floor, "drawn": n, "agree": a, "contradict": c,
                   "rate": c / (a + c) if a + c else 0.0, "kept": a / base if base else 1.0}
            grid.append(row)
            if evalue not in chosen and row["rate"] <= th["contradiction_ceiling"] + 1e-12 \
                    and row["kept"] >= th["agreeing_kept"] - 1e-12:
                chosen[evalue] = floor
    return {"grid": grid, "chosen": {str(e): chosen.get(e) for e in E_GRID}}


def pfam_blind_group(label: str) -> str | None:
    sub = label.partition("/")[2].casefold() or label.casefold()
    for group, prefixes in PFAM_BLIND.items():
        if sub.startswith(prefixes):
            return group
    return None


def criteria(regions: list[dict], cost: dict, th: dict) -> dict:
    n, k = len(regions), sum(r["found"] for r in regions)
    lo, hi = wilson(k, n)
    spans = [r["span"] for r in regions if r["found"]]
    median = statistics.median(spans) if spans else 0.0
    high = sum(s >= th["a2_high_span"] for s in spans) / len(spans) if spans else 0.0
    blind = [r for r in regions if pfam_blind_group(r["label"])]
    missed = [r for r in blind if not r["found"]]
    ratio = cost.get("ratio")
    return {
        "A1": {"regions": n, "found": k, "wilson": [lo, hi], "pass": lo >= th["a1_lower_bound"]},
        "A2": {"median_span": median, "share_high": high,
               "pass": median >= th["a2_median_span"] and high >= th["a2_high_share"]},
        "A3": {"regions": len(blind), "missed": len(missed),
               "by_group": {g: sum(pfam_blind_group(r["label"]) == g for r in blind) for g in PFAM_BLIND},
               "pass": bool(blind) and not missed},
        "A4": {**cost, "pass": ratio is not None and ratio >= th["a4_cost_ratio"]},
    }


def median_seconds(arm_dir: Path) -> float | None:
    t = arm_dir / "timing.tsv"
    if not t.exists():
        return None
    secs = [float(l.split("\t")[1]) for l in t.read_text().splitlines()[1:]]
    return statistics.median(secs) if secs else None


def load_thresholds(path: str | None) -> dict:
    th = dict(DRAFT_THRESHOLDS)
    if path:
        th.update(json.loads(Path(path).read_text()))
    return th


def cmd_compare(args) -> None:
    run = Path(args.run)
    ref = read_hits(run / args.set / args.ref / "hits.tsv")
    arm = read_hits(run / args.set / args.arm / "hits.tsv")
    rows = compare(ref, arm, args.evalue)
    out = run / args.set / f"regions.{args.ref}-vs-{args.arm}.tsv"
    with open(out, "w") as f:
        f.write("family\tlabel\tstart\tend\tmodel\tfound\tspan\tpieces\n")
        for r in rows:
            f.write(f"{r['family']}\t{r['label']}\t{r['start']}\t{r['end']}\t{r['model']}\t"
                    f"{int(r['found'])}\t{r['span']:.4f}\t{r['pieces']}\n")
    print(f"wrote {len(rows)} regions to {out}", file=sys.stderr)


def cmd_report(args) -> None:
    run = Path(args.run)
    th = load_thresholds(args.thresholds)
    manifest = json.loads((run / "manifest.json").read_text())
    report = {"thresholds": th, "sets": {}}
    lines = [f"# Benchmark verdict: {run.name}", "",
             f"TE-Aid {manifest.get('teaid_commit', '?')[:10]}, protocol "
             f"{(manifest.get('protocol') or {}).get('commit') or 'uncommitted'}, "
             f"signed: {(manifest.get('protocol') or {}).get('signed')}.", "",
             "Thresholds: " + ", ".join(f"{k} {v}" for k, v in th.items()), ""]
    for set_name in manifest["sets"]:
        set_dir = run / set_name
        entry = {}
        if set_name != "shuffled" and (set_dir / "default" / "hits.tsv").exists():
            entry["specificity"] = specificity(read_hits(set_dir / "default" / "hits.tsv"), th)
        if (set_dir / "deep" / "hits.tsv").exists() and (set_dir / "default" / "hits.tsv").exists():
            regions = compare(read_hits(set_dir / "deep" / "hits.tsv"),
                              read_hits(set_dir / "default" / "hits.tsv"))
            deep_s, default_s = median_seconds(set_dir / "deep"), median_seconds(set_dir / "default")
            cost = {"deep_median_s": deep_s, "default_median_s": default_s,
                    "ratio": deep_s / default_s if deep_s and default_s else None}
            entry["criteria"] = criteria(regions, cost, th)
        report["sets"][set_name] = entry
        lines += [f"## {set_name}", ""]
        if "specificity" in entry:
            lines += ["Question 1, floor chosen per E-value: " + ", ".join(
                f"E <= {e}: {f if f is not None else 'none'}" for e, f in entry["specificity"]["chosen"].items()), ""]
        if "criteria" in entry:
            for name, c in entry["criteria"].items():
                lines.append(f"- {name}: {'pass' if c['pass'] else 'FAIL'} "
                             + json.dumps({k: v for k, v in c.items() if k != 'pass'}))
            lines.append("")
    (run / "verdict.json").write_text(json.dumps(report, indent=2) + "\n")
    (run / "verdict.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def cmd_verify(args) -> int:
    run = Path(args.run)
    manifest = json.loads((run / "manifest.json").read_text())
    checks: list[tuple[str, bool, str]] = []
    state = protocol_state()
    created = datetime.datetime.fromisoformat(manifest["created"]).timestamp()
    signed_first = bool(state["signed"] and state["committed_at"] and state["committed_at"] <= created)
    checks.append(("protocol signed and committed before the run", signed_first,
                   f"protocol commit {state['commit']}, run created {manifest['created']}"))
    for set_name, arms in manifest.get("arms", {}).items():
        for arm in arms:
            hits = read_hits(run / set_name / arm / "hits.tsv")
            rows = compare(hits, hits)
            ok = all(r["found"] and abs(r["span"] - 1) < 1e-9 for r in rows)
            checks.append((f"{set_name}/{arm} matched against itself", ok,
                           f"{len(rows)} regions, recall {sum(r['found'] for r in rows)}/{len(rows)}"))
    shuffled = manifest.get("arms", {}).get("shuffled", {})
    for arm in shuffled:
        hits = read_hits(run / "shuffled" / arm / "hits.tsv")
        n = sum(1 for hs in hits.values() for h in hs if h.evalue <= PRIMARY_E)
        checks.append((f"shuffled/{arm}: no hits at E <= {PRIMARY_E}", n == 0, f"{n} hits"))
    if not shuffled:
        checks.append(("negative control run", False, "no 'shuffled' set in this run"))
    if args.repeat:
        other = Path(args.repeat)
        for set_name, arms in manifest.get("arms", {}).items():
            for arm in arms:
                a, b = run / set_name / arm / "hits.tsv", other / set_name / arm / "hits.tsv"
                same = b.exists() and a.read_text() == b.read_text()
                checks.append((f"{set_name}/{arm} identical on a second run", same, str(b)))
    for name, ok, detail in checks:
        print(f"{'pass' if ok else 'FAIL'}  {name}  ({detail})")
    return 0 if all(ok for _, ok, _ in checks) else 1


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="run one arm over one family set")
    p.add_argument("--arm", choices=ARMS, required=True)
    p.add_argument("--set", nargs=2, metavar=("NAME", "FASTA"), required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cpus", type=int, default=2)
    p.add_argument("--min-orf", type=int, default=400, help="getorf minimum, nt (the CLI's default)")
    p.add_argument("--limit", type=int, help="first N families only, for testing")
    p = sub.add_parser("shuffle", help="dinucleotide-shuffled consensus sequences, the negative control")
    p.add_argument("fasta")
    p.add_argument("--n", type=int, default=5)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("compare", help="match reference regions against an arm")
    p.add_argument("run")
    p.add_argument("--set", required=True)
    p.add_argument("--ref", default="deep")
    p.add_argument("--arm", default="default")
    p.add_argument("--evalue", type=float, default=PRIMARY_E)
    p = sub.add_parser("report", help="Question 1 grid and criteria A1-A4")
    p.add_argument("run")
    p.add_argument("--thresholds", help="JSON with the signed numbers; the draft's otherwise")
    p = sub.add_parser("verify", help="the protocol's checks before any number counts")
    p.add_argument("run")
    p.add_argument("--repeat", help="a second run of the same arms, for determinism")
    args = parser.parse_args(argv)
    if args.cmd == "run":
        cmd_run(args)
    elif args.cmd == "shuffle":
        cmd_shuffle(args)
    elif args.cmd == "compare":
        cmd_compare(args)
    elif args.cmd == "report":
        cmd_report(args)
    elif args.cmd == "verify":
        return cmd_verify(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
