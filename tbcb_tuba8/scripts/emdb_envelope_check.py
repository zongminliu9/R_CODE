#!/usr/bin/env python3
"""Coarse but independent test of the AF3 TBCB placement against EMD-2447.

EMD-2447 is the only experimental structure of a TBCB-containing complex (alpha-tubulin + TBCE +
TBCB). At 21 A nominal resolution on a 32^3 grid with 4.66 A voxels it cannot say which residues
touch; it constrains overall dimensions, how many lobes the assembly has, and whether a candidate
sub-complex can sit inside it with contiguous room left for the third component.

That last point is the actual test here. If the AF3 TBCB+alpha sub-complex is too large or the
wrong shape, the remaining envelope volume will not be TBCE-sized or will be fragmented, and the
AF3 architecture is contradicted. If it fits with a single TBCE-sized lobe left over, the AF3
architecture is merely NOT contradicted -- at this resolution that is the strongest statement
available, and it is not the same as support.
"""
import os, json, struct, warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy import ndimage
import biotite.structure as struc, biotite.structure.io.pdbx as pdbx

MAP = "data/emdb/emd_2447.map"
with open(MAP,"rb") as f: hdr = f.read(1024)
nx,ny,nz,mode = struct.unpack("<4i", hdr[0:16])
mx,my,mz = struct.unpack("<3i", hdr[28:40])
cella = struct.unpack("<3f", hdr[40:52])
vox = cella[0]/mx
d = np.fromfile(MAP, dtype=np.float32, offset=1024, count=nx*ny*nz).reshape(nz,ny,nx)

MASS = dict(alpha=50.1, TBCE=59.3, TBCB=27.4)
TOT = sum(MASS.values())
target_vol = TOT*1000*1.21                      # 1.21 A^3 per Da
# threshold that encloses the expected mass
ths = np.linspace(d.max()*0.2, d.max()*0.95, 400)
vols = np.array([ (d>t).sum()*vox**3 for t in ths ])
t_mass = ths[np.argmin(np.abs(vols-target_vol))]
mask = d > t_mass
print(f"map: {nx}^3 grid, {vox:.2f} A voxels, nominal 21 A")
print(f"mass-consistent threshold {t_mass:.4f} -> {mask.sum()} voxels, "
      f"{mask.sum()*vox**3:.0f} A^3 (target {target_vol:.0f} for {TOT:.1f} kDa)")

lab, nlab = ndimage.label(mask)
sizes = ndimage.sum(mask, lab, range(1, nlab+1))
print(f"connected components at that threshold: {nlab}  "
      f"sizes(A^3) {[int(s*vox**3) for s in sorted(sizes, reverse=True)[:5]]}")
big = (lab == (1+int(np.argmax(sizes))))

def principal_dims(binary_mask, voxel):
    pts = np.argwhere(binary_mask).astype(float)*voxel
    pts -= pts.mean(0)
    w = np.linalg.eigvalsh(np.cov(pts.T))[::-1]
    # extent along each principal axis (full width, 2*sqrt(3*lambda) approximates a uniform body)
    return np.sqrt(w)*2*np.sqrt(3), float(np.sqrt(w.sum()))
dims, rg = principal_dims(big, vox)
print(f"envelope principal extents {dims[0]:.0f} x {dims[1]:.0f} x {dims[2]:.0f} A, "
      f"Rg {rg:.1f} A")

def coords_of(cif, chains=None):
    arr = pdbx.get_structure(pdbx.CIFFile.read(cif), model=1, use_author_fields=True)
    arr = arr[struc.filter_amino_acids(arr)]
    if chains is not None: arr = arr[np.isin(arr.chain_id, chains)]
    return arr

def atom_dims(arr):
    p = arr.coord.astype(float) - arr.coord.mean(0)
    w = np.linalg.eigvalsh(np.cov(p.T))[::-1]
    return np.sqrt(w)*2*np.sqrt(3), float(np.sqrt(w.sum()))

AF3 = "/home/user/uploads_fold/e50dbe95-fold_tbcb_tuba8/fold_tbcb_tuba8_model_0.cif"
cplx = coords_of(AF3)
chs = {c: int((cplx.chain_id==c).sum()) for c in np.unique(cplx.chain_id)}
chT = min(chs, key=chs.get); chA = max(chs, key=chs.get)
print(f"\nAF3 complex: chain {chT} = TBCB, chain {chA} = alpha-tubulin")
for nm, sel in (("TBCB+alpha (AF3)", None), ("alpha alone", [chA]), ("TBCB alone", [chT])):
    a = cplx if sel is None else cplx[np.isin(cplx.chain_id, sel)]
    dm, r = atom_dims(a)
    print(f"  {nm:20s} extents {dm[0]:6.0f} x {dm[1]:5.0f} x {dm[2]:5.0f} A   Rg {r:5.1f} A")

# volume bookkeeping: what is left for TBCE after TBCB+alpha?
v_cplx = (MASS["alpha"]+MASS["TBCB"])*1000*1.21
v_left = mask.sum()*vox**3 - v_cplx
print(f"\nvolume bookkeeping at the mass-consistent threshold:")
print(f"  envelope          {mask.sum()*vox**3:8.0f} A^3")
print(f"  TBCB+alpha needs  {v_cplx:8.0f} A^3")
print(f"  left for TBCE     {v_left:8.0f} A^3   (TBCE expects "
      f"{MASS['TBCE']*1000*1.21:.0f} A^3, ratio {v_left/(MASS['TBCE']*1000*1.21):.2f})")

print(f"\nwhat this test can and cannot conclude:")
print(f"  voxel {vox:.2f} A and 21 A resolution -> the envelope has "
      f"{mask.sum()} voxels for {TOT:.0f} kDa, i.e. about "
      f"{mask.sum()/ (TOT*1000/110):.1f} voxels per residue.")
print(f"  A residue-level interface assignment needs <5 A resolution. This map is 4x coarser,")
print(f"  so it CANNOT test which alpha-tubulin surface TBCB occupies. It can only test gross")
print(f"  architecture, and the volume bookkeeping above is the whole of that test.")
json.dump(dict(voxel_A=float(vox), threshold=float(t_mass),
               envelope_volume=float(mask.sum()*vox**3), n_components=int(nlab),
               envelope_extents=[float(x) for x in dims], envelope_Rg=rg,
               volume_left_for_TBCE=float(v_left),
               ratio_to_expected=float(v_left/(MASS['TBCE']*1000*1.21))),
          open("results/af3/emdb2447_envelope.json","w"), indent=1)
