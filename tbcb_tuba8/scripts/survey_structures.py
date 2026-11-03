#!/usr/bin/env python3
"""Check existence/metadata of candidate PDB entries and search RCSB for the TBCB/TBCE/tubulin system."""
import json, urllib.request, urllib.parse, sys

def jget(url, data=None):
    req = urllib.request.Request(url, data=(json.dumps(data).encode() if data else None),
                                 headers={"Content-Type":"application/json","User-Agent":"tbcb-research/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r: return json.loads(r.read().decode())

CANDIDATES = ["9M1M","9EEB","9EDR","1JFF","5IJ0","3RYF","6S8K","1WHG","2YQ6","4ICV","5N5N","7SJ7"]
print("=== candidate entries ===")
for pid in CANDIDATES:
    try:
        e = jget(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
        t = e.get("struct",{}).get("title","")[:72]
        meth = ",".join(x.get("method","") for x in e.get("exptl",[]))
        res = e.get("rcsb_entry_info",{}).get("resolution_combined",["-"])
        dep = e.get("rcsb_accession_info",{}).get("initial_release_date","")[:10]
        print(f"  {pid}  {meth:22s} res={str(res[0]):>6s}  {dep}  {t}")
    except Exception as ex:
        print(f"  {pid}  NOT FOUND / error ({type(ex).__name__})")

print("\n=== full-text search: tubulin cofactor complexes ===")
for term in ["tubulin cofactor B","tubulin cofactor E","TBCD TBCE ARL2","tubulin folding cofactor"]:
    q = {"query":{"type":"terminal","service":"full_text","parameters":{"value":term}},
         "return_type":"entry","request_options":{"paginate":{"start":0,"rows":12},
         "results_verbosity":"compact"}}
    try:
        r = jget("https://search.rcsb.org/rcsbsearch/v2/query", q)
        ids = r.get("result_set",[])
        print(f"  '{term}': n={r.get('total_count')} -> {ids[:12]}")
    except Exception as ex:
        print(f"  '{term}': error {type(ex).__name__} {ex}")
