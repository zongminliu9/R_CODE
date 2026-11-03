#!/usr/bin/env python3
import json, urllib.request
def jget(url, data=None):
    req = urllib.request.Request(url, data=(json.dumps(data).encode() if data else None),
        headers={"Content-Type":"application/json","User-Agent":"tbcb-research/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r: return json.loads(r.read().decode())

print("=== newly-found entries ===")
for pid in ["9ZNW","9ZNX","9ZNY","9ZNZ","9M1J","9M1K","9M1L","9ZNV","1T0Y","2KJ6","1V6E","2KJR","4B6M","4B6W","1QSD","2YUH","3MXZ","1H7C","5AJ8","2L3L","4CQI"]:
    try:
        e = jget(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
        t = e.get("struct",{}).get("title","")[:80]
        meth = ",".join(x.get("method","")[:9] for x in e.get("exptl",[]))
        res = e.get("rcsb_entry_info",{}).get("resolution_combined",["-"])
        npol = e.get("rcsb_entry_info",{}).get("polymer_entity_count",0)
        print(f"  {pid} {meth:10s} res={str(res[0]):>6s} nent={npol}  {t}")
    except Exception as ex: print(f"  {pid} err {type(ex).__name__}")

print("\n=== search: CAP-Gly + EEY/tyrosine peptide complexes ===")
for term in ["CAP-Gly EEY","CAP-Gly domain peptide complex","p150glued CAP-Gly tubulin peptide","CLIP-170 CAP-Gly"]:
    q={"query":{"type":"terminal","service":"full_text","parameters":{"value":term}},
       "return_type":"entry","request_options":{"paginate":{"start":0,"rows":14},"results_verbosity":"compact"}}
    try:
        r=jget("https://search.rcsb.org/rcsbsearch/v2/query",q)
        print(f"  '{term}': n={r.get('total_count')} -> {r.get('result_set',[])[:14]}")
    except Exception as ex: print(f"  '{term}': err {ex}")
