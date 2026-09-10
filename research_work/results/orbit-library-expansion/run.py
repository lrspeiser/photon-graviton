"""Expand missing populated training states, then verify numerical path coverage."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import time
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache'
OUT = CACHE / 'orbit-library-expansion'
OUT.mkdir(parents=True, exist_ok=True)
NAMES = ['x_kpc', 'y_kpc', 'z_kpc', 'vx_kms', 'vy_kms', 'vz_kms']
TIMES = np.linspace(0, .25, 501)

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

cov = load('expansion_coverage', HERE.parent / 'orbit-population-coverage/run.py')

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def save(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')

def sample():
    d = pd.read_parquet(CACHE / 'stellar-catalogs/stellar-errors-conditional_minus_0.017.parquet',
                        filters=[('holdout_role', '==', 'training')])
    assert len(d) == 77927 and (d.holdout_role == 'training').all()
    side = pd.read_parquet(CACHE / 'training-orbit-launches/seed-eligibility.parquet')
    d = d.merge(side, on='source_id', validate='one_to_one').sort_values('source_id').reset_index(drop=True)
    assert len(d) == 27884
    angle = np.radians(json.loads((HERE.parent/'stellar-orbit-support/protocol.json').read_text())['bar_angle_degrees'])
    x, v = cov.cylindrical_to_cartesian(d, angle)
    for j, name in enumerate(NAMES): d[name] = np.c_[x, v][:, j]
    return d

def prepare():
    d = sample()
    d = d[(d.seed_eligibility == 'passes_position_screen') & ~d.selected_as_new_seed].reset_index(drop=True)
    assert len(d) == 27642
    audit = pd.read_parquet(CACHE / 'orbit-population-coverage/training-coverage.parquet')
    d = d.merge(audit[['source_id', 'min_velocity_gap_within_0.5_kpc']], on='source_id', validate='one_to_one')
    assert d.source_id.is_monotonic_increasing
    uncovered = d['min_velocity_gap_within_0.5_kpc'].to_numpy() > 50
    states = d[NAMES].to_numpy()
    scaled = states / np.array([.5]*3 + [50.]*3)
    tree = cKDTree(scaled)
    neighbors = []
    for i in range(len(d)):
        ids = np.array(tree.query_ball_point(scaled[i], np.sqrt(2)+1e-12), dtype=int)
        delta = states[ids] - states[i]
        good = (np.linalg.norm(delta[:, :3], axis=1) <= .5) & (np.linalg.norm(delta[:, 3:], axis=1) <= 50)
        neighbors.append(ids[good])
    chosen, records = [], []
    for ir, (lo, hi) in enumerate([(.5, 3.5), (3.5, 5), (5, 9)]):
        for iz, (zl, zh) in enumerate([(0, .2), (.2, .5), (.5, 1.5)]):
            mask = (d.mean_R_kpc >= lo) & ((d.mean_R_kpc <= hi) if hi == 9 else (d.mean_R_kpc < hi))
            mask &= (abs(d.mean_z_kpc) >= zl) & ((abs(d.mean_z_kpc) <= zh) if zh == 1.5 else (abs(d.mean_z_kpc) < zh))
            candidates = np.flatnonzero(mask)
            remaining = uncovered & mask.to_numpy()
            for step in range(4):
                gains = np.array([np.count_nonzero(remaining[neighbors[i]]) if i not in chosen else -1 for i in candidates])
                best = int(np.argmax(gains))
                if gains[best] <= 0: break
                i = int(candidates[best]); chosen.append(i)
                records.append(dict(seed_index=len(chosen)-1, source_id=str(d.iloc[i].source_id),
                    R_stratum=ir, height_stratum=iz, newly_covered_initial_targets=int(gains[best])))
                remaining[neighbors[i]] = False
    seeds = d.iloc[chosen][['source_id','holdout_role','window_residual_arcsec'] + NAMES].copy()
    provenance = pd.read_parquet(CACHE/'stellar-orbit-support/stellar-input-provenance.parquet',
                                columns=['source_id','parallax_interval_disjoint_5'])
    seeds = seeds.merge(provenance, on='source_id', validate='one_to_one')
    seeds['R_stratum'] = [r['R_stratum'] for r in records]
    seeds['height_stratum'] = [r['height_stratum'] for r in records]
    assert seeds.source_id.is_unique and (seeds.holdout_role == 'training').all()
    assert seeds.source_id.astype(str).tolist() == [r['source_id'] for r in records]
    seeds.to_parquet(OUT/'launches.parquet', index=False)
    save('selection.json', dict(scope='Training-only populated missing-state selection, before integration.',
         new_seeds=len(seeds), records=records, distance_disagreement_flags=int(seeds.parallax_interval_disjoint_5.sum()),
         old_seed_count=72, eligible_nonseed_targets=len(d),
         initial_uncovered_targets=int(uncovered.sum()),
         launches_sha256=digest(OUT/'launches.parquet'), protocol_sha256=digest(HERE/'protocol.md'),
         input_hashes={str(p.relative_to(ROOT)):digest(p) for p in [cov.FILES['moments'], cov.FILES['eligibility'],
                     CACHE/'orbit-population-coverage/training-coverage.parquet']}, holdouts_opened=False))
    print('Selected', len(seeds), 'new launches; distance flags', int(seeds.parallax_interval_disjoint_5.sum()), flush=True)

def integrate():
    driver = load('expansion_orbit_driver', HERE.parent/'full-bar-orbits/run.py')
    field = driver.Field('full')
    seeds = pd.read_parquet(OUT/'launches.parquet')
    assert digest(OUT/'launches.parquet') == json.loads((HERE/'selection.json').read_text())['launches_sha256']
    rows = []
    hashes = {str(p.relative_to(ROOT)):digest(p) for p in field.paths}
    for j, y0 in enumerate(seeds[NAMES].to_numpy()):
        started = time.monotonic(); attempts = []; previous = None; retained = None; error = None
        def rhs(t, y):
            x, p = y[None, :3], y[None, 3:]
            a = field.evaluate(x)[1]
            return np.c_[p-driver.OMEGA*driver.cross(x), a-driver.OMEGA*driver.cross(p)].ravel()
        for tol in [2e-9, 2e-11, 2e-13]:
            try:
                sol = solve_ivp(rhs, [0, .25], y0, t_eval=TIMES, method='DOP853', rtol=tol, atol=tol*.01)
                if not sol.success: raise RuntimeError(sol.message)
            except (RuntimeError, ValueError) as exc:
                error = str(exc); break
            y = sol.y.T
            item = dict(rtol=tol, nfev=sol.nfev)
            if previous is not None:
                dx = float(np.max(np.linalg.norm(y[:, :3]-previous[:, :3], axis=1)))
                dv = float(np.max(np.linalg.norm(y[:, 3:]-previous[:, 3:], axis=1)))
                item.update(position_difference_kpc=dx, velocity_difference_kms=dv, passed=dx<1e-4 and dv<.01)
            attempts.append(item); retained = y
            if item.get('passed', False): break
            previous = y
        row = dict(seed_index=j, source_id=str(seeds.iloc[j].source_id), attempts=attempts, error=error,
                   distance_disagreement_flag=bool(seeds.iloc[j].parallax_interval_disjoint_5), numerical_pass=False)
        if error is None and retained is not None:
            pot = np.concatenate([field.evaluate(q[:, :3])[0] for q in np.array_split(retained, 32)])
            x, p = retained[:, :3], retained[:, 3:]
            J = .5*np.sum(p*p, axis=1)+pot-driver.OMEGA*(x[:, 0]*p[:, 1]-x[:, 1]*p[:, 0])
            drift = float(np.max(abs(J-J[0]))/220**2)
            row.update(Jacobi_drift_over_220_squared=drift, numerical_pass=bool(attempts[-1].get('passed', False) and drift<1e-5))
            np.savez_compressed(OUT/f'orbit-{j}.npz', times=TIMES, trajectory=retained, source_id=seeds.iloc[j].source_id)
        row['seconds'] = time.monotonic()-started
        rows.append(row)
        save('integration.json', dict(completed=len(rows)==len(seeds), rows=rows, input_field_hashes=hashes,
             driver_sha256=digest(HERE.parent/'full-bar-orbits/run.py'), code_sha256=digest(Path(__file__)),
             launches_sha256=digest(OUT/'launches.parquet'), holdouts_opened=False))
        print(f'{j+1}/{len(seeds)} numerical={row["numerical_pass"]} {row["seconds"]:.1f}s {error or ""}', flush=True)

if __name__ == '__main__':
    {'prepare':prepare, 'integrate':integrate}[sys.argv[1]]()
