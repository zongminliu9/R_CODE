"""Assemble the force-field files for the TBCB/alpha-tubulin intermediate-state campaign.

Protein  : Amber ff14SB, ISOLDE's ParmEd conversion (ff/amberff14SB.xml, unprefixed Amber types so the
           GTP template can share classes with it).  Checked against OpenMM's own protein.ff14SB.xml by
           a00b_validate_ff.py before anything is run.
GTP      : Meagher, Redman & Carlson, J Comput Chem 24:1016 (2003) polyphosphate set, as converted by
           ISOLDE (bryce_set.xml, residue BRYCE_GTP; source phos.lib + frcmod.phos).  Only the template and
           the phosphate-specific terms are taken from it.  Every other bonded/LJ term comes from ff14SB's
           parm10, which reproduces what tleap gives with leaprc.protein.ff14SB + GTP.prep + frcmod.phos
           (bryce_set.xml also carries FAD/heme frcmod terms that would otherwise override parm10, e.g. the
           CB-NA-C-O improper and the CB-CB bond; those are deliberately NOT copied).
Water    : TIP3P, Joung-Cheatham Na+/Cl- (copied from OpenMM amber14/tip3p.xml).
Mg2+     : Li, Roberts, Chakravorty & Merz, JCTC 9:2733 (2013) 12-6 IOD set for TIP3P
           (Rmin/2 1.395 A, eps 0.01491 kcal/mol), charge +2.  ISOLDE's tip3p_IOD_multivalent.xml lists the
           same LJ but with charge 1.0, which is wrong, so it is written here explicitly.
Outputs  : ff/gtp_meagher.xml, ff/water_ions_tip3p_iod.xml, ff/gtp_bonds.xml, ff/gtp_hydrogens.xml
"""
import os, sys, math
import xml.etree.ElementTree as ET

FF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ff')
A = ET.parse(os.path.join(FF, 'amberff14SB.xml')).getroot()
B = ET.parse(os.path.join(FF, 'bryce_set.xml')).getroot()
T = ET.parse(os.path.join(FF, 'openmm_amber14_tip3p.xml')).getroot()

# ---------------------------------------------------------------- GTP template
src = [r for r in B.find('Residues') if r.get('name') == 'BRYCE_GTP']
assert len(src) == 1
src = src[0]
atoms = [(a.get('name'), a.get('type'), float(a.get('charge'))) for a in src.findall('Atom')]
bonds = [(b.get('atomName1'), b.get('atomName2')) for b in src.findall('Bond')]
qtot = sum(q for _, _, q in atoms)
assert abs(qtot + 4.0) < 1e-6, qtot
assert len(atoms) == 44 and len(bonds) == 46

ff = ET.Element('ForceField')
info = ET.SubElement(ff, 'Info')
ET.SubElement(info, 'Source').text = ('GTP(4-): Meagher, Redman & Carlson J Comput Chem 24:1016 (2003); template and '
                                       'phosphate-specific terms from ISOLDE bryce_set.xml (BRYCE_GTP), remaining terms ff14SB/parm10')
at = ET.SubElement(ff, 'AtomTypes')
# the only new atom type is the terminal phosphate oxygen O3; everything else is an existing ff14SB type
o3 = [x for x in B.find('AtomTypes') if x.get('name') == 'O3']
assert len(o3) == 1
ET.SubElement(at, 'Type', name='O3', **{'class': 'O3'}, element='O', mass='16.0')
known = {x.get('name') for x in A.find('AtomTypes')}
for _, t, _ in atoms:
    assert t in known or t == 'O3', f'type {t} missing from ff14SB'
rs = ET.SubElement(ff, 'Residues')
r = ET.SubElement(rs, 'Residue', name='GTP')
for n, t, q in atoms:
    ET.SubElement(r, 'Atom', name=n, type=t, charge=repr(q))
for a1, a2 in bonds:
    ET.SubElement(r, 'Bond', atomName1=a1, atomName2=a2)
# NOTE: BRYCE_GTP carries two ExternalBond records (O1G, C2'); a free GTP has none, so they are dropped.

def pick(tag, child, want):
    out = []
    for x in B.find(tag):
        if x.tag != child:
            continue
        t = tuple(x.get(f'type{i}') or x.get(f'class{i}') or '' for i in range(1, len(want[0]) + 1))
        if t in want or t[::-1] in want:
            out.append(x)
    return out

hb = ET.SubElement(ff, 'HarmonicBondForce')
for x in pick('HarmonicBondForce', 'Bond', [('O3', 'P')]):
    hb.append(x)
ha = ET.SubElement(ff, 'HarmonicAngleForce')
for x in pick('HarmonicAngleForce', 'Angle', [('O3', 'P', 'O3'), ('OS', 'P', 'O3')]):
    ha.append(x)
pt = ET.SubElement(ff, 'PeriodicTorsionForce')
PHOS = [('O2', 'P', 'OS', 'CT'), ('O2', 'P', 'OS', 'P'), ('O3', 'P', 'OS', 'P'), ('P', 'OS', 'P', 'OS'), ('H1', 'CT', 'OS', 'P')]
got = pick('PeriodicTorsionForce', 'Proper', PHOS)
assert len(got) == len(PHOS), [g.attrib for g in got]
for x in got:
    pt.append(x)
nb = ET.SubElement(ff, 'NonbondedForce', coulomb14scale=A.find('NonbondedForce').get('coulomb14scale'),
                   lj14scale=A.find('NonbondedForce').get('lj14scale'))
ET.SubElement(nb, 'UseAttributeFromResidue', name='charge')
o3nb = [x for x in B.find('NonbondedForce') if x.get('type') == 'O3']
assert len(o3nb) == 1
ET.SubElement(nb, 'Atom', type='O3', sigma=o3nb[0].get('sigma'), epsilon=o3nb[0].get('epsilon'))
ET.indent(ff)
ET.ElementTree(ff).write(os.path.join(FF, 'gtp_meagher.xml'))

# bond and hydrogen definitions, so PDB/PDBx GTP residues get bonds and Modeller.addHydrogens can build H
bd = ET.Element('Residues'); rr = ET.SubElement(bd, 'Residue', name='GTP')
for a1, a2 in bonds:
    ET.SubElement(rr, 'Bond', **{'from': a1, 'to': a2})
ET.indent(bd); ET.ElementTree(bd).write(os.path.join(FF, 'gtp_bonds.xml'))
el = {n: t for n, t, _ in atoms}
hd = ET.Element('Residues'); rh = ET.SubElement(hd, 'Residue', name='GTP')
for a1, a2 in bonds:
    for h, p in ((a1, a2), (a2, a1)):
        if el[h] in ('H1', 'H2', 'H5', 'H', 'HO'):
            ET.SubElement(rh, 'H', name=h, parent=p)
nH = len(rh.findall('H'))
assert nH == 12, nH
ET.indent(hd); ET.ElementTree(hd).write(os.path.join(FF, 'gtp_hydrogens.xml'))

# ---------------------------------------------------------------- water + ions
w = ET.Element('ForceField')
ET.SubElement(ET.SubElement(w, 'Info'), 'Source').text = ('TIP3P and Joung-Cheatham Na+/Cl- from OpenMM amber14/tip3p.xml; '
                                                          'Mg2+ Li-Merz 2013 12-6 IOD (TIP3P), charge +2')
keep_types = {'tip3p-O', 'tip3p-H', 'tip3p_standard-Na+', 'tip3p_standard-Cl-'}
wt = ET.SubElement(w, 'AtomTypes')
for x in T.find('AtomTypes'):
    if x.get('name') in keep_types:
        wt.append(x)
ET.SubElement(wt, 'Type', name='mg_iod-Mg2+', **{'class': 'mg_iod-Mg2+'}, element='Mg', mass='24.305')
wr = ET.SubElement(w, 'Residues')
for x in T.find('Residues'):
    if x.get('name') in ('HOH', 'NA', 'CL'):
        wr.append(x)
mg = ET.SubElement(wr, 'Residue', name='MG'); ET.SubElement(mg, 'Atom', name='MG', type='mg_iod-Mg2+', charge='2.0')
for tag in ('HarmonicBondForce', 'HarmonicAngleForce'):
    if T.find(tag) is not None:
        w.append(T.find(tag))
wnb = ET.SubElement(w, 'NonbondedForce', coulomb14scale=T.find('NonbondedForce').get('coulomb14scale'),
                    lj14scale=T.find('NonbondedForce').get('lj14scale'))
for x in T.find('NonbondedForce'):
    if x.tag == 'UseAttributeFromResidue' or x.get('type') in keep_types:
        wnb.append(x)
RMIN2_A, EPS_KCAL = 1.395, 0.01491
sigma_nm = 2 * RMIN2_A / 2 ** (1 / 6) / 10.0
ET.SubElement(wnb, 'Atom', type='mg_iod-Mg2+', sigma=repr(sigma_nm), epsilon=repr(EPS_KCAL * 4.184))
if T.find('Script') is not None:
    pass
ET.indent(w); ET.ElementTree(w).write(os.path.join(FF, 'water_ions_tip3p_iod.xml'))

print('GTP: 44 atoms, net charge', round(qtot, 4), '| phosphate-specific terms: bonds', len(hb), 'angles', len(ha),
      'propers', len(pt), '| H defs', nH)
print('Mg2+ IOD sigma', round(sigma_nm, 6), 'nm eps', round(EPS_KCAL * 4.184, 6), 'kJ/mol')
