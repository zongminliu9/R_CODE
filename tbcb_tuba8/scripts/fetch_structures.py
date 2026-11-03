#!/usr/bin/env python3
"""Download PDB/mmCIF for the working structure set + AlphaFold models. Reports chain/ligand composition."""
import json, os, urllib.request, gzip, io, sys
os.makedirs("data/pdb", exist_ok=True); os.makedirs("data/afdb", exist_ok=True)

PDB = {
 # cofactor-tubulin machinery (late maturation)
 "9M1M":"TBC-DEC-Arl2-ab-tubulin GDP 2.21A","9M1L":"TBC-DE-Arl2-ab-tubulin GTP 2.55A",
 "9M1J":"TBC-DE-Arl2-b-tubulin 2.14A","9EEB":"TBCD/E/G + tubulin dimer 3.7A",
 "9EDR":"TBCD/E/G/C + tubulin, TBCC N-term 3.8A",
 # TBCB domains (mouse!) and homologs
 "1WHG":"mouse TBCB CAP-Gly (NMR)","1V6E":"mouse TBCB UBL N-term (NMR)",
 "4B6M":"T.brucei TBCB CAP-Gly 1.59A","4B6W":"T.brucei TBCB architecture 2.35A",
 "1T0Y":"TBCB UBL (NMR)","4ICV":"human TBCE UBL 1.45A",
 # CAP-Gly : EEY/F peptide complexes  (binding-mode templates)
 "2HL3":"CAP-Gly + EEY peptide","2HQH":"CLIP-170 CAP-Gly complex","2HKN":"CLIP-170 CAP-Gly complex",
 "9YNE":"p150glued CAP-Gly + tubulin peptide","9YNH":"p150glued CAP-Gly + tubulin peptide",
 "8PR0":"p150glued CAP-Gly + tubulin peptide","3RDV":"CAP-Gly complex",
 # tubulin references
 "1JFF":"ab-tubulin 3.5A (canonical numbering)","3RYF":"GTP-tubulin:RB3-SLD 2.52A",
 "5IJ0":"human TUBB3 microtubule","7SJ7":"recombinant human MT 3.8A",
}
AFDB = {"Q9JJZ2":"mouse TUBA8","Q9D1E6":"mouse TBCB","Q8CIV8":"mouse TBCE",
        "P68369":"mouse TUBA1A","P68368":"mouse TUBA4A","P05213":"mouse TUBA1B",
        "Q99426":"human TBCB","Q9NY65":"human TUBA8"}

def dl(url, path, gz=False):
    req = urllib.request.Request(url, headers={"User-Agent":"tbcb-research/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r: raw = r.read()
    if gz: raw = gzip.decompress(raw)
    open(path,"wb").write(raw); return len(raw)

print("=== PDB (mmCIF) ===")
ok=[]
for pid, desc in PDB.items():
    p = f"data/pdb/{pid}.cif"
    if os.path.exists(p) and os.path.getsize(p)>1000: print(f"  cached {pid}"); ok.append(pid); continue
    try:
        n = dl(f"https://files.rcsb.org/download/{pid}.cif.gz", p, gz=True)
        print(f"  ok {pid}  {n/1e6:6.2f} MB  {desc}"); ok.append(pid)
    except Exception as e:
        print(f"  !! {pid}  {type(e).__name__} {e}")

print("\n=== AlphaFold DB ===")
for acc, desc in AFDB.items():
    p = f"data/afdb/{acc}.pdb"
    if os.path.exists(p) and os.path.getsize(p)>1000: print(f"  cached {acc}"); continue
    try:
        meta = json.loads(urllib.request.urlopen(urllib.request.Request(
            f"https://alphafold.ebi.ac.uk/api/prediction/{acc}",
            headers={"User-Agent":"tbcb-research/1.0"}), timeout=60).read().decode())
        n = dl(meta[0]["pdbUrl"], p)
        print(f"  ok {acc}  {n/1e6:5.2f} MB  v{meta[0].get('latestVersion')}  {desc}")
    except Exception as e:
        print(f"  !! {acc}  {type(e).__name__} {e}")

json.dump({"pdb":ok,"afdb":list(AFDB)}, open("data/meta/structure_index.json","w"), indent=1)
