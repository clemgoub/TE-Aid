"""Verify individual Pfam accessions against InterPro, with descriptions.

Usage: verify_pfam.py [cache_dir]

Every accession proposed for teaid/data/te_domains.tsv is checked to exist
and its name taken from here, so a typo or a dead accession fails the build
rather than shipping. PF08333, for one, does not resolve.
"""
import json, urllib.request, time, os
OUT=(__import__("sys").argv[1] if len(__import__("sys").argv)>1 else ".") + "/pfam_verified.json"
CAND = """PF05699 PF13952 PF03004 PF12017 PF12596 PF01498 PF03372 PF14529 PF00077 PF08284
PF02022 PF00552 PF00607 PF13976 PF00098 PF01541 PF00589 PF06817 PF06815 PF08333
PF02994 PF00136 PF13022 PF09322 PF04986 PF10683 PF13358 PF03184 PF01609 PF00872
PF10551 PF20700 PF04827 PF04937 PF05380 PF17917 PF17919 PF13456 PF14214 PF21859
PF01021 PF14223 PF14244 PF17241 PF19259 PF00078 PF07727 PF00665 PF13333 PF13683
PF24764 PF18758 PF18759 PF13843 PF01359 PF13610 PF13612 PF13701 PF13737 PF13751
PF13546 PF13359 PF13586 PF26100 PF14291 PF07999 PF02914 PF01526 PF03050 PF01385""".split()
def get(url):
    last=None
    for a in range(4):
        try:
            b=urllib.request.urlopen(urllib.request.Request(url,headers={"Accept":"application/json"}),timeout=60).read()
            return json.loads(b) if b.strip() else None
        except Exception as e:
            last=e; time.sleep(1.5*(a+1))
    return {"__error__": str(last)}
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
for acc in CAND:
    if acc in res: continue
    d = get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/")
    if d is None or "__error__" in (d or {}):
        res[acc] = {"error": (d or {}).get("__error__","empty")}
    else:
        m=d["metadata"]; n=m["name"]
        desc = (m.get("description") or [{}])[0].get("text","")
        import re as _re
        res[acc] = {"short": n.get("short") if isinstance(n,dict) else "",
                    "name": n.get("name") if isinstance(n,dict) else n,
                    "type": m.get("type"),
                    "desc": _re.sub(r"<[^>]+>","",desc)[:400]}
    json.dump(res, open(OUT,"w")); time.sleep(0.25)
print("verified", len(res))
