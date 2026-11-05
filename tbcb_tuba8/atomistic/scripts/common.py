"""Shared pieces: force field, system construction, staged equilibration, robust checkpoint I/O, CV helpers."""
import os, sys, json, gzip, time, shutil, math
import numpy as np
import openmm as mm
from openmm import app, unit

HERE = os.path.dirname(os.path.abspath(__file__))
FF = os.path.join(HERE, '..', 'ff')
NM = unit.nanometer
KJ = unit.kilojoule_per_mole
SOLVENT = {'HOH', 'NA', 'CL', 'WAT'}
_bonds_loaded = [False]


def forcefield():
    return app.ForceField(os.path.join(FF, 'amberff14SB.xml'), os.path.join(FF, 'gtp_meagher.xml'),
                          os.path.join(FF, 'water_ions_tip3p_iod.xml'))


def load_cif(path):
    if not _bonds_loaded[0]:
        app.Topology.loadBondDefinitions(os.path.join(FF, 'gtp_bonds.xml')); _bonds_loaded[0] = True
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        px = app.PDBxFile(fh)
    return px.topology, px.positions


def write_cif(path, top, pos):
    tmp = path + '.tmp'
    with gzip.open(tmp, 'wt') as fh:
        app.PDBxFile.writeFile(top, pos, fh, keepIds=True)
    os.replace(tmp, path)


def make_system(top, hmr=True, constraints=app.HBonds, ff=None):
    ff = ff or forcefield()
    kw = dict(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * NM, constraints=constraints, rigidWater=True,
              ewaldErrorTolerance=5e-4)
    if hmr:
        kw['hydrogenMass'] = 3.0 * unit.amu
    return ff.createSystem(top, **kw)


def get_platform(name=None):
    for n in ([name] if name else []) + ['CUDA', 'OpenCL', 'CPU']:
        try:
            p = mm.Platform.getPlatformByName(n)
            return p, ({'Precision': 'mixed'} if n in ('CUDA', 'OpenCL') else {})
        except Exception:
            continue
    raise RuntimeError('no platform')


def solute_indices(top):
    return [a.index for a in top.atoms() if a.residue.name not in SOLVENT]


def backbone_indices(top):
    return [a.index for a in top.atoms() if a.name in ('N', 'CA', 'C', 'O') and a.residue.name not in SOLVENT | {'GTP', 'MG'}]


def heavy_solute_indices(top):
    return [a.index for a in top.atoms() if a.residue.name not in SOLVENT and a.element is not None and a.element.symbol != 'H']


def mg_shell_pairs(top, pos, cutoff=0.26):
    P = np.array(pos.value_in_unit(NM)) if hasattr(pos, 'value_in_unit') else np.asarray(pos)
    atoms = list(top.atoms())
    pairs = []
    for m in [a.index for a in atoms if a.residue.name == 'MG']:
        d = np.linalg.norm(P - P[m], axis=1)
        for i in np.where((d > 0) & (d < cutoff))[0]:
            if atoms[i].element.symbol == 'O':
                pairs.append((int(m), int(i)))
    return pairs


def add_restraints(system, top, pos, idx, k_name='k_pos', mg_pairs=None):
    """Position restraints (periodic distance) on idx with global strength k_name, and flat-bottom Mg-O."""
    f = mm.CustomExternalForce(f'{k_name}*periodicdistance(x, y, z, x0, y0, z0)^2')
    f.addGlobalParameter(k_name, 0.0)
    for p in ('x0', 'y0', 'z0'):
        f.addPerParticleParameter(p)
    P = pos.value_in_unit(NM) if hasattr(pos, 'value_in_unit') else pos
    for i in idx:
        f.addParticle(int(i), [float(x) for x in P[i]])
    system.addForce(f)
    if mg_pairs:
        g = mm.CustomBondForce('k_mg*step(r-0.23)*(r-0.23)^2')
        g.addGlobalParameter('k_mg', 0.0)
        g.setUsesPeriodicBoundaryConditions(True)
        for a, b in mg_pairs:
            g.addBond(a, b, [])
        system.addForce(g)
    return f


def density_g_cm3(system, box):
    mass = sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles()))
    a, b, c = [np.array(v.value_in_unit(NM)) for v in box]
    vol = abs(np.dot(a, np.cross(b, c))) * 1e-21  # cm^3
    return mass * 1.66053906660e-24 / vol


def equilibrate(top, pos, T=310.0, seed=1, platform=None, log=print, npt=True, total_free_ns=1.0, scale=1.0,
                system_fn=None):
    """Staged equilibration.  Returns (positions, velocities, box, record).  The system used here carries the
    restraint forces and a barostat; callers build their own production system and load the returned state."""
    system = system_fn(top) if system_fn else make_system(top, hmr=True)
    bb = backbone_indices(top); heavy = heavy_solute_indices(top)
    pairs = mg_shell_pairs(top, pos)
    add_restraints(system, top, pos, heavy, 'k_heavy', pairs)
    add_restraints(system, top, pos, bb, 'k_bb')
    if npt:
        system.addForce(mm.MonteCarloBarostat(1 * unit.bar, T * unit.kelvin, 25))
    integ = mm.LangevinMiddleIntegrator(50 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integ.setRandomNumberSeed(seed)
    plat, props = get_platform(platform)
    ctx = mm.Context(system, integ, plat, props)
    ctx.setPositions(pos); ctx.setPeriodicBoxVectors(*top.getPeriodicBoxVectors())
    has_mg = bool(pairs)
    def setp(name, v):
        if name != 'k_mg' or has_mg:
            ctx.setParameter(name, v)
    setp('k_heavy', 1000.0); setp('k_mg', 50000.0); setp('k_bb', 0.0)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 2000)
    rec = {'seed': seed, 'T': T, 'stages': []}
    P0 = np.array(ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(NM))
    def stage(name, ps, dt, Tset, kh, kb, kmg):
        integ.setStepSize(dt * unit.picoseconds); integ.setTemperature(Tset * unit.kelvin)
        setp('k_heavy', kh); setp('k_bb', kb); setp('k_mg', kmg)
        integ.step(max(1, int(round(ps * scale / dt))))
        st = ctx.getState(getPositions=True, getEnergy=True)
        P = np.array(st.getPositions(asNumpy=True).value_in_unit(NM))
        A_, B_ = P[bb] - P[bb].mean(0), P0[bb] - P0[bb].mean(0)
        U_, S_, Vt_ = np.linalg.svd(A_.T @ B_); d_ = np.sign(np.linalg.det(Vt_.T @ U_.T))
        R_ = Vt_.T @ np.diag([1, 1, d_]) @ U_.T
        rmsd = float(np.sqrt(np.mean(np.sum(((R_ @ A_.T).T - B_) ** 2, 1)))) * 10
        box = st.getPeriodicBoxVectors()
        r = dict(stage=name, ps=ps * scale, T=Tset, k_heavy=kh, k_bb=kb, k_mg=kmg, E_kcal=st.getPotentialEnergy().value_in_unit(KJ) / 4.184,
                 backbone_rmsd_A=round(rmsd, 3), density=round(density_g_cm3(system, box), 4),
                 mg_dist_A=[round(float(np.linalg.norm(P[a] - P[b])) * 10, 3) for a, b in pairs])
        rec['stages'].append(r); log(json.dumps(r))
    for i, Tset in enumerate(np.linspace(50, T, 6)):
        stage(f'heat{i}', 20, 0.002, float(Tset), 1000.0, 0.0, 50000.0)
    stage('npt_heavy1000', 400, 0.004, T, 1000.0, 0.0, 50000.0)
    stage('npt_bb500', 400, 0.004, T, 0.0, 500.0, 50000.0)
    stage('npt_bb100', 400, 0.004, T, 0.0, 100.0, 5000.0)
    stage('npt_bb10', 400, 0.004, T, 0.0, 10.0, 0.0)
    stage('npt_free', total_free_ns * 1000, 0.004, T, 0.0, 0.0, 0.0)
    st = ctx.getState(getPositions=True, getVelocities=True)
    return (st.getPositions(asNumpy=True), st.getVelocities(asNumpy=True), st.getPeriodicBoxVectors(), rec)


def push(src, dst, tries=5):
    for t in range(tries):
        try:
            tmp = dst + '.tmp'
            shutil.copy2(src, tmp)
            assert os.path.getsize(tmp) == os.path.getsize(src)
            os.replace(tmp, dst)
            return True
        except Exception as ex:
            print(f'push {src} -> {dst} failed: {ex}', flush=True); time.sleep(3)
    return False


def save_state_npz(path, pos, vel, box, extra=None):
    tmp = path[:-4] + '.tmp.npz'
    d = dict(pos=np.asarray(pos.value_in_unit(NM) if hasattr(pos, 'value_in_unit') else pos, dtype=np.float64),
             vel=np.asarray(vel.value_in_unit(NM / unit.picosecond) if hasattr(vel, 'value_in_unit') else vel, dtype=np.float64),
             box=np.array([np.array(v.value_in_unit(NM)) if hasattr(v, 'value_in_unit') else np.asarray(v) for v in box]))
    if extra:
        d.update({k: np.asarray(v) for k, v in extra.items()})
    np.savez(tmp, **d)
    z = np.load(tmp); assert z['pos'].shape == d['pos'].shape
    os.replace(tmp, path)


def load_state_npz(path):
    z = np.load(path, allow_pickle=True)
    box = [mm.Vec3(*v) * NM for v in z['box']]
    return z['pos'], z['vel'], box, {k: z[k] for k in z.files if k not in ('pos', 'vel', 'box')}


def dihedral(p0, p1, p2, p3):
    b0 = p0 - p1; b1 = p2 - p1; b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n; w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w); y = np.dot(np.cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))
