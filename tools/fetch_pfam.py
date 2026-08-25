"""Fetch Pfam clan membership and name searches from InterPro into a cache.

Usage: fetch_pfam.py [cache_dir]   (default: current directory)

Feeds tools/build_te_domains.py. Cached so re-runs are cheap and so the
curated table can always be rebuilt from the same evidence.
"""
import json, urllib.request, urllib.parse, time, os
CACHE = (__import__("sys").argv[1] if len(__import__("sys").argv)>1 else ".") + "/pfam_cache.json"
def get(url):
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            body = urllib.request.urlopen(req, timeout=90).read()
            if not body.strip():
                return {"results": [], "next": None}   # 204: no matches
            return json.loads(body)
        except Exception as e:
            last = e; time.sleep(2 * (attempt + 1))
    raise last
def page_all(url, cap=600):
    out = []
    while url and len(out) < cap:
        d = get(url); out += d["results"]; url = d.get("next"); time.sleep(0.2)
    return out
def slim(r):
    m = r["metadata"]
    n = m["name"]
    if isinstance(n, dict):
        return {"acc": m["accession"], "short": n.get("short") or "", "name": n.get("name") or ""}
    return {"acc": m["accession"], "short": "", "name": n}

data = json.load(open(CACHE)) if os.path.exists(CACHE) else {"clans": {}, "search": {}}
for clan in ["CL0219","CL0027","CL0523","CL0169","CL0329"]:
    if clan in data["clans"]: continue
    data["clans"][clan] = [slim(r) for r in page_all(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/set/pfam/{clan}/?page_size=200")]
    print(f"{clan}: {len(data['clans'][clan])}", flush=True); json.dump(data, open(CACHE,"w"))

TERMS = ["transposase","transposon","retrotransposon","integrase","reverse transcriptase",
         "retroviral","gag protein","helitron","tyrosine recombinase","DDE",
         "piggyBac","mariner","hAT element","Mutator","Harbinger","CACTA","Crypton","Maverick",
         "Polinton","Penelope","ribonuclease H","aspartyl protease","zinc knuckle","rolling circle"]
for t in TERMS:
    if t in data["search"]: continue
    data["search"][t] = [slim(r) for r in page_all(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/?search={urllib.parse.quote(t)}&page_size=200")]
    print(f"search {t!r}: {len(data['search'][t])}", flush=True); json.dump(data, open(CACHE,"w"))
print("DONE")
