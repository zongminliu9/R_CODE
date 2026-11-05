"""Hamiltonian replica exchange on one GPU context, shared by REST2, the alchemical point-mutation legs and the
umbrella-sampling ladders.  All thermodynamic states share one System and differ only in global context parameters,
so every state's energy is evaluated on each replica every iteration (full u_kn for MBAR at no bookkeeping cost).

Robustness features that this project has needed before:
 * the checkpoint holds every replica's positions/velocities in float64, the state permutation, the RNG state and the
   whole energy record; it is written to node-local disk and only then copied to the shared directory, read back and
   verified (count of iterations) before the copy replaces the previous one;
 * the solvated topology is pinned on disk by the caller, so a resumed run always has the same atom count;
 * a wall-clock budget stops the loop cleanly and the next job in the dependency chain resumes;
 * nothing is reported as finished unless the record has the requested number of iterations.
"""
import os, sys, time, json, shutil, math
import numpy as np
import openmm as mm
from openmm import app, unit

KB = 0.0083144626  # kJ/mol/K


class HREX:
    def __init__(self, system, positions, box, states, temperature_K=310.0, timestep_ps=0.004, friction=1.0,
                 steps_per_iter=250, platform='CUDA', shared_dir='.', local_dir='/tmp', tag='hrex',
                 cv_fn=None, save_states=(0,), save_idx=None, save_topology=None, seed=12345, precision='mixed',
                 swap_scheme='neighbors'):
        self.states = states
        self.K = len(states)
        self.T = temperature_K
        self.kT = KB * temperature_K
        self.n = steps_per_iter
        self.shared = shared_dir; self.local = local_dir; self.tag = tag
        os.makedirs(shared_dir, exist_ok=True); os.makedirs(local_dir, exist_ok=True)
        self.cv_fn = cv_fn
        self.save_states = list(save_states)
        self.save_idx = np.asarray(save_idx) if save_idx is not None else None
        self.save_top = save_topology
        self.swap_scheme = swap_scheme
        self.integ = mm.LangevinMiddleIntegrator(temperature_K * unit.kelvin, friction / unit.picosecond,
                                                 timestep_ps * unit.picosecond)
        self.integ.setRandomNumberSeed(seed % 2147483647)
        plat, props = None, {}
        for name in ([platform] if platform else []) + ['CUDA', 'OpenCL', 'CPU']:
            try:
                plat = mm.Platform.getPlatformByName(name)
                props = {'Precision': precision} if name in ('CUDA', 'OpenCL') else {}
                self.ctx = mm.Context(system, self.integ, plat, props)
                break
            except Exception as ex:
                print(f'platform {name} unavailable: {ex}', flush=True)
                plat = None
        assert plat is not None
        self.platform = plat.getName()
        if self.platform in ('CPU', 'Reference') and (os.environ.get('CUDA_VISIBLE_DEVICES') or os.environ.get('SLURM_JOB_GPUS')):
            raise RuntimeError('GPU job fell back to the CPU platform; refusing to run (a GPU context is probably still held)')
        def _vec(v):
            v = v.value_in_unit(unit.nanometer) if hasattr(v, 'value_in_unit') else v
            return mm.Vec3(*[float(x) for x in v])
        self.box = [_vec(v) for v in box] * unit.nanometer
        self.ctx.setPeriodicBoxVectors(*self.box)
        self.params = sorted(states[0].keys())
        for s in states:
            assert sorted(s.keys()) == self.params
        self.rng = np.random.default_rng(seed)
        self.ck_shared = os.path.join(shared_dir, f'{tag}_ck.npz')
        self.ck_local = os.path.join(local_dir, f'{tag}_ck.npz')
        self.rec_u, self.rec_s, self.rec_cv = [], [], []
        self.n_acc = np.zeros(self.K - 1); self.n_try = np.zeros(self.K - 1)
        self.iteration = 0
        self.pos = None
        self.positions0 = positions

    # ------------------------------------------------------------------ state handling
    def set_state(self, k):
        for p in self.params:
            self.ctx.setParameter(p, self.states[k][p])

    def energies(self):
        out = np.empty(self.K)
        for k in range(self.K):
            self.set_state(k)
            out[k] = self.ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / self.kT
        return out

    # ------------------------------------------------------------------ checkpoint
    def _load(self, path):
        z = np.load(path, allow_pickle=True)
        self.pos = z['pos']; self.vel = z['vel']; self.s = z['s']; self.iteration = int(z['iteration'])
        self.rec_u = list(z['rec_u']); self.rec_s = list(z['rec_s'])
        self.rec_cv = list(z['rec_cv']) if 'rec_cv' in z.files else []
        self.n_acc = z['n_acc']; self.n_try = z['n_try']
        self.rng.bit_generator.state = json.loads(str(z['rng']))
        return True

    def resume_or_init(self, init_velocities=True, sequential=False, warm_steps=2500, min_iter=1000):
        """sequential=True: build the starting replicas by walking along the state ladder - replica k is minimised
        and run for warm_steps in state k starting from replica k-1's configuration - so that no state starts from
        a configuration that is badly strained in its Hamiltonian (e.g. a large side chain switched on in a pocket
        shaped for a small one)."""
        cands = [p for p in (self.ck_local, self.ck_shared) if os.path.exists(p)]
        best, it = None, -1
        for p in cands:
            try:
                z = np.load(p, allow_pickle=True); i = int(z['iteration'])
                if i > it:
                    best, it = p, i
            except Exception as ex:
                print(f'unreadable checkpoint {p}: {ex}', flush=True)
        if best is not None:
            self._load(best)
            assert self.pos.shape[1] == self.ctx.getSystem().getNumParticles(), \
                f'checkpoint has {self.pos.shape[1]} atoms, system {self.ctx.getSystem().getNumParticles()}'
            print(f'resumed from {best} at iteration {self.iteration}', flush=True)
            return
        P = np.array(self.positions0.value_in_unit(unit.nanometer) if hasattr(self.positions0, 'value_in_unit')
                     else self.positions0)
        self.pos = P.copy() if P.ndim == 3 else np.repeat(P[None], self.K, axis=0)   # per-replica starts allowed
        self.vel = np.zeros_like(self.pos)
        self.s = np.arange(self.K)
        if sequential:
            X = self.pos[0].copy()
            for k in range(self.K):
                self.set_state(k); self.ctx.setPositions(X)
                mm.LocalEnergyMinimizer.minimize(self.ctx, 10.0, min_iter)
                self.ctx.setVelocitiesToTemperature(self.T * unit.kelvin, int(self.rng.integers(1, 2 ** 30)))
                self.integ.step(warm_steps)
                st = self.ctx.getState(getPositions=True, getVelocities=True)
                X = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                if not np.all(np.isfinite(X)):
                    raise RuntimeError(f'non-finite positions while warming state {k}')
                self.pos[k] = X; self.vel[k] = st.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
            print(f'sequential warm start through {self.K} states done', flush=True)
            return
        if init_velocities:
            for r in range(self.K):
                self.set_state(self.s[r]); self.ctx.setPositions(self.pos[r])
                self.ctx.setVelocitiesToTemperature(self.T * unit.kelvin, int(self.rng.integers(1, 2 ** 30)))
                self.vel[r] = self.ctx.getState(getVelocities=True).getVelocities(asNumpy=True).value_in_unit(
                    unit.nanometer / unit.picosecond)

    def checkpoint(self):
        tmp = self.ck_local[:-4] + '.tmp.npz'
        np.savez(tmp, pos=self.pos, vel=self.vel, s=self.s, iteration=self.iteration,
                 rec_u=np.array(self.rec_u), rec_s=np.array(self.rec_s), rec_cv=np.array(self.rec_cv, dtype=object),
                 n_acc=self.n_acc, n_try=self.n_try, rng=json.dumps(self.rng.bit_generator.state))
        os.replace(tmp, self.ck_local)
        for attempt in range(5):
            try:
                stmp = self.ck_shared[:-4] + '.tmp.npz'
                shutil.copy2(self.ck_local, stmp)
                z = np.load(stmp, allow_pickle=True)
                assert int(z['iteration']) == self.iteration and z['pos'].shape == self.pos.shape
                os.replace(stmp, self.ck_shared)
                return True
            except Exception as ex:
                print(f'push failed ({ex}); retry {attempt + 1}', flush=True)
                time.sleep(5)
        return False

    # ------------------------------------------------------------------ main loop
    def run(self, n_iter_total, budget_min, ck_every_min=10.0, dcd_every=1, log_every=50):
        t0 = time.time(); tck = t0
        dcds = {}
        if self.save_idx is not None and self.save_top is not None:
            for k in self.save_states:
                fn = os.path.join(self.local, f'{self.tag}_s{k}_from{self.iteration:06d}.dcd')
                dcds[k] = (fn, app.DCDFile(open(fn, 'wb'), self.save_top, 0.001 * self.n * 4))
        while self.iteration < n_iter_total:
            u = np.empty((self.K, self.K)); cvs = []
            for r in range(self.K):
                self.set_state(self.s[r])
                self.ctx.setPositions(self.pos[r]); self.ctx.setVelocities(self.vel[r])
                self.integ.step(self.n)
                st = self.ctx.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
                self.pos[r] = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                self.vel[r] = st.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
                if not np.all(np.isfinite(self.pos[r])):
                    raise RuntimeError(f'non-finite positions in replica {r} (state {self.s[r]}) at iteration {self.iteration}')
                u[r] = self.energies()
                if self.cv_fn is not None:
                    cvs.append(self.cv_fn(self.pos[r]))
            # swaps on the state permutation
            inv = np.empty(self.K, int); inv[self.s] = np.arange(self.K)   # inv[state] = replica
            if self.swap_scheme == 'neighbors':
                start = self.iteration % 2
                pairs = [(k, k + 1) for k in range(start, self.K - 1, 2)]
            else:
                pairs = [(k, k + 1) for k in range(self.K - 1)]
                self.rng.shuffle(pairs)
            for a, b in pairs:
                i, j = inv[a], inv[b]
                dlog = -(u[i, b] + u[j, a] - u[i, a] - u[j, b])
                self.n_try[a] += 1
                if dlog >= 0 or self.rng.random() < math.exp(dlog):
                    self.s[i], self.s[j] = b, a; inv[a], inv[b] = j, i
                    self.n_acc[a] += 1
            self.rec_u.append(u.astype(np.float64)); self.rec_s.append(self.s.copy())
            if self.cv_fn is not None:
                self.rec_cv.append(cvs)
            self.iteration += 1
            for k, (fn, d) in dcds.items():
                if self.iteration % dcd_every == 0:
                    r = inv[k] if False else int(np.where(self.s == k)[0][0])
                    d.writeModel(self.pos[r][self.save_idx] * unit.nanometer, periodicBoxVectors=self.box)
            if self.iteration % log_every == 0:
                acc = self.n_acc / np.maximum(self.n_try, 1)
                el = (time.time() - t0) / 60
                print(f'it {self.iteration}/{n_iter_total}  {el:.1f} min  acc min {acc.min():.2f} '
                      f'median {np.median(acc):.2f}  platform {self.platform}', flush=True)
            now = time.time()
            if (now - tck) / 60 > ck_every_min:
                self.checkpoint(); self._push_dcds(dcds); tck = now
            if (now - t0) / 60 > budget_min:
                break
        self.checkpoint(); self._push_dcds(dcds, close=True)
        return self.iteration >= n_iter_total

    def _push_dcds(self, dcds, close=False):
        for k, (fn, d) in list(dcds.items()):
            try:
                d._file.flush() if hasattr(d, '_file') else None
            except Exception:
                pass
            for attempt in range(3):
                try:
                    shutil.copy2(fn, os.path.join(self.shared, os.path.basename(fn)))
                    break
                except Exception as ex:
                    print(f'dcd push failed {ex}', flush=True); time.sleep(3)

    def acceptance(self):
        return (self.n_acc / np.maximum(self.n_try, 1)).tolist()
