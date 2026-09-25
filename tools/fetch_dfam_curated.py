#!/usr/bin/env python3
"""Fetch Dfam's curated families for one clade, through the Dfam REST API.

We calibrate the protein thresholds (docs/BRIEF_v2.md §5.1a) against families
whose classification a curator set, and the GenomeArk goby library has none.
Dfam's curated families are classified, and the API returns them with their
consensus sequences in a few megabytes, with no FamDB partition to download.

Writes three files beside ``--out``:

- ``<out>.fa`` — consensus sequences, headers ``>name#type/subtype accession`` so
  that ``teaid`` reads the class label as it does for any RepeatModeler library;
- ``<out>.tsv`` — accession, name, length, type, subtype, the full Dfam
  classification path, and the number of curated coding sequences;
- ``<out>.provenance.json`` — the Dfam release, the query, the date, the family
  count and the FASTA's sha256, so a calibration run can cite what it used.

Usage::

    tools/fetch_dfam_curated.py --clade 7227 --out dev-data/dfam/dmel
    tools/fetch_dfam_curated.py --clade 9606 --out dev-data/dfam/human

``--clade`` takes an NCBI Taxonomy ID. Ancestral families are included
(``clade_relatives=ancestors``), because a genome's repeats include families
defined higher in its lineage. Uncurated ("DR") families are not.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://dfam.org/api"
PAGE = 100


def get(path: str, **params) -> dict:
    url = f"{API}/{path}" + (f"?{urllib.parse.urlencode(params)}" if params else "")
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return json.load(response)
        except OSError as exc:
            if attempt == 4:
                raise SystemExit(f"fetch failed after 5 tries: {url}: {exc}")
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--clade", required=True, help="NCBI Taxonomy ID, e.g. 7227")
    parser.add_argument("--out", required=True, help="output prefix, e.g. dev-data/dfam/dmel")
    args = parser.parse_args()

    query = dict(clade=args.clade, clade_relatives="ancestors", format="full")
    version = get("version")
    records, start, total = [], 0, None
    while total is None or start < total:
        page = get("families", **query, limit=PAGE, start=start)
        total = page["total_count"]
        records += page["results"]
        start += PAGE
        print(f"  {len(records):,}/{total:,}", file=sys.stderr)
        time.sleep(0.5)  # be polite to a public API
    if len(records) != total:
        raise SystemExit(f"expected {total} families, received {len(records)}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fasta_lines, rows = [], []
    for r in records:
        label = "/".join(x for x in (r.get("repeat_type_name"), r.get("repeat_subtype_name")) if x)
        fasta_lines += [f">{r['name']}#{label or 'Unknown'} {r['accession']}", r["consensus_sequence"]]
        rows.append([r["accession"], r["name"], r["length"], r.get("repeat_type_name") or "",
                     r.get("repeat_subtype_name") or "", r["classification"],
                     len(r.get("coding_seqs") or [])])
    fasta = "\n".join(fasta_lines) + "\n"
    (out.with_suffix(".fa")).write_text(fasta)
    with out.with_suffix(".tsv").open("w") as handle:
        handle.write("accession\tname\tlength\ttype\tsubtype\tclassification\tcoding_seqs\n")
        for row in rows:
            handle.write("\t".join(str(v) for v in row) + "\n")
    provenance = dict(
        dfam_version=version.get("dfam_version"), dfam_release_date=version.get("release_date"),
        api=API, query=query, fetched=datetime.date.today().isoformat(), families=len(records),
        fasta_sha256=hashlib.sha256(fasta.encode()).hexdigest(),
    )
    out.with_suffix(".provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote {len(records):,} families to {out.with_suffix('.fa')}", file=sys.stderr)


if __name__ == "__main__":
    main()
