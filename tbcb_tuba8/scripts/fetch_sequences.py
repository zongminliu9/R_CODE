#!/usr/bin/env python3
"""Resolve + download reviewed UniProt entries for the TBCB/tubulin system (mouse primary, human secondary).
Pure stdlib: no deps. Writes FASTA + JSON metadata into data/seq and data/meta."""
import json, os, sys, time, urllib.parse, urllib.request

OUT_SEQ, OUT_META = "data/seq", "data/meta"
os.makedirs(OUT_SEQ, exist_ok=True); os.makedirs(OUT_META, exist_ok=True)

GENES = ["Tuba8","Tuba1a","Tuba1b","Tuba1c","Tuba3a","Tuba4a","Tubb2a","Tubb2b",
         "Tubb3","Tubb4a","Tubb4b","Tubb5","Tbcb","Tbce","Tbca","Tbcc","Tbcd",
         "Arl2","Cct1","Tcp1","Cct2","Cct5","Pfdn1","Pfdn2","Vash1","Ttl"]
ORGS = {"mouse": 10090, "human": 9606}

def get(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "tbcb-tuba8-research/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode()
        except Exception as e:
            if i == tries-1: raise
            time.sleep(2**i)

def search(gene, taxid):
    q = f'(gene_exact:{gene}) AND (organism_id:{taxid}) AND (reviewed:true)'
    fields = "accession,id,gene_primary,gene_names,protein_name,length,sequence,organism_name,cc_subcellular_location,ft_mod_res"
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q)
           + "&fields=" + fields + "&format=json&size=10")
    return json.loads(get(url)).get("results", [])

index = {}
for org, taxid in ORGS.items():
    for gene in GENES:
        try:
            res = search(gene, taxid)
        except Exception as e:
            print(f"  !! {org}/{gene}: {type(e).__name__}", file=sys.stderr); continue
        if not res:
            print(f"  -- {org}/{gene}: no reviewed entry"); continue
        # prefer entry whose primary gene matches exactly (case-insensitive)
        res.sort(key=lambda e: (e.get("genes",[{}])[0].get("geneName",{}).get("value","").lower() != gene.lower(),
                                -e.get("sequence",{}).get("length",0)))
        e = res[0]
        acc = e["primaryAccession"]
        seq = e["sequence"]["value"]
        prim = e.get("genes",[{}])[0].get("geneName",{}).get("value", gene)
        name = e.get("proteinDescription",{}).get("recommendedName",{}).get("fullName",{}).get("value","")
        key = f"{org}_{prim}"
        index[key] = {"accession": acc, "entry": e.get("uniProtkbId"), "gene": prim,
                      "organism": org, "protein_name": name, "length": len(seq), "sequence": seq}
        with open(f"{OUT_SEQ}/{key}.fasta","w") as fh:
            fh.write(f">{key}|{acc}|{prim}|{org} {name}\n")
            for i in range(0,len(seq),60): fh.write(seq[i:i+60]+"\n")
        print(f"  ok {key:18s} {acc:8s} len={len(seq):4d}  {name[:58]}")
        with open(f"{OUT_META}/{key}.json","w") as fh: json.dump(e, fh, indent=1)

with open(f"{OUT_META}/sequence_index.json","w") as fh: json.dump(index, fh, indent=1)
print(f"\n{len(index)} entries -> {OUT_META}/sequence_index.json")
