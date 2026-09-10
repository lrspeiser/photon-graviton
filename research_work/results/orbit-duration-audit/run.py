"""Extend six declared library probes; compare finite-time occupation, not data fits."""
from pathlib import Path
import hashlib
import importlib.util
import json
import time
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache/orbit-duration-audit'
CACHE.mkdir(parents=True, exist_ok=True)
DRIVER = HERE.parent / 'full-bar-orbits/run.py'
spec = importlib.util.spec_from_file_location('duration_driver', DRIVER)
orbit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(orbit)
SOURCE = ROOT / 'research_work/data-cache/training-orbit-launches/full-field-orbits.npz'
SEEDS = HERE.parent / 'training-orbit-launches/launches.json'
# First selected seed of the low/high-height strata at each of three radii.
# Keep the two known distance flags; selection is not based on orbit residuals.
INDICES = [0, 16, 24, 40, 48, 64]
TIMES = np.linspace(0, 1, 2001)
R_EDGES = [0, .5, 3.5, 5, 9, 20, 40, 80, np.inf]
Z_EDGES = [-np.inf, -1.5, -.5, -.2, 0, .2, .5, 1.5, np.inf]
PHI_EDGES = np.linspace(-np.pi, np.pi, 9)

def digest(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def occupation(y, angular=True):
    x = y[:, :3]
    points = np.c_[np.hypot(x[:, 0], x[:, 1]), x[:, 2], np.arctan2(x[:, 1], x[:, 0])]
    bins = [R_EDGES, Z_EDGES, PHI_EDGES] if angular else [R_EDGES, Z_EDGES]
    h = np.histogramdd(points if angular else points[:, :2], bins=bins)[0]
    assert h.sum() == len(y)
    return h.ravel() / len(y)

def tv(a, b, angular=True):
    return float(.5 * np.abs(occupation(a, angular) - occupation(b, angular)).sum())

def save(value, name):
    (HERE / name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')

def main():
    field = orbit.Field('full')
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in [SOURCE, SEEDS, DRIVER, Path(__file__), *field.paths]}
    with np.load(SOURCE) as archive:
        initial = archive['trajectories'][0]
    seeds = json.loads(SEEDS.read_text())
    rows = []
    for index in INDICES:
        start = time.monotonic()
        attempts, previous, retained, error = [], None, None, None
        def rhs(t, y):
            x, p = y[None, :3], y[None, 3:]
            a = field.evaluate(x)[1]
            return np.c_[p-orbit.OMEGA*orbit.cross(x), a-orbit.OMEGA*orbit.cross(p)].ravel()
        for tol in [2e-9, 2e-11, 2e-13]:
            try:
                sol = solve_ivp(rhs, [0, 1], initial[index], method='DOP853', t_eval=TIMES,
                                rtol=tol, atol=tol*.01)
                if not sol.success:
                    raise RuntimeError(sol.message)
            except (ValueError, RuntimeError) as exc:
                error = str(exc)
                break
            retained = sol.y.T
            attempt = dict(rtol=tol, nfev=sol.nfev)
            if previous is not None:
                dx = np.linalg.norm(retained[:, :3]-previous[:, :3], axis=1).max()
                dv = np.linalg.norm(retained[:, 3:]-previous[:, 3:], axis=1).max()
                attempt.update(position_difference_kpc=float(dx), velocity_difference_kms=float(dv),
                               passed=bool(dx < 1e-4 and dv < .01))
            attempts.append(attempt)
            print(index, attempt, flush=True)
            if attempt.get('passed', False):
                break
            previous = retained
        row = dict(seed_index=index, source_id=seeds[index]['source_id'],
                   distance_disagreement_flag=seeds[index]['parallax_interval_disjoint_5'],
                   attempts=attempts, error=error, numerical_pass=False)
        if retained is not None and error is None:
            np.savez_compressed(CACHE / f'orbit-{index}.npz', times=TIMES, trajectory=retained)
            pot = np.concatenate([field.evaluate(part[:, :3])[0] for part in np.array_split(retained, 128)])
            x, p = retained[:, :3], retained[:, 3:]
            J = .5*(p*p).sum(axis=1)+pot-orbit.OMEGA*(x[:, 0]*p[:, 1]-x[:, 1]*p[:, 0])
            drift = float(np.abs(J-J[0]).max()/220**2)
            row.update(Jacobi_drift_over_220_squared=drift,
                       numerical_pass=bool(attempts[-1].get('passed', False) and drift < 1e-5),
                       maximum_radius_kpc=float(np.linalg.norm(x, axis=1).max()),
                       maximum_absolute_height_kpc=float(abs(x[:, 2]).max()))
            # Equal-duration adjacent windows; exclude each window's initial sample.
            row['occupation_windows'] = []
            for stop in [500, 1000, 2000]:
                half = stop//2
                a, b = retained[1:half+1], retained[half+1:stop+1]
                row['occupation_windows'].append(dict(duration_kpc_per_kms=float(TIMES[stop]),
                    spatial_TV=tv(a, b), R_z_TV=tv(a, b, False),
                    half_sampling_spatial_TV=tv(a[1::2], b[1::2])))
            row['first_quarter_vs_full_TV'] = tv(retained[1:501], retained[1:])
            row['last_two_quarters_TV'] = tv(retained[1001:1501], retained[1501:2001])
            row['trajectory_sha256'] = digest(CACHE / f'orbit-{index}.npz')
        row['seconds'] = time.monotonic()-start
        rows.append(row)
        save(dict(rows=rows, completed=False), 'progress.json')
    for p in [SOURCE, SEEDS, DRIVER, Path(__file__), *field.paths]:
        assert digest(p) == hashes[str(p.relative_to(ROOT))]
    save(dict(scope='Six duration probes, not stationary-population validation or a likelihood.',
              indices=INDICES, selection='First seed in low/high-height stratum at each of three radii.',
              duration_kpc_per_kms=1., rows=rows, source_hashes=hashes,
              spatial_bins=dict(R_kpc=[str(v) for v in R_EDGES], z_kpc=[str(v) for v in Z_EDGES], phi_rad=PHI_EDGES.tolist()),
              total_variation_definition='0.5 sum_bins abs(p_early - p_late); descriptive bin-dependent distance, not a p-value.',
              holdouts_opened=False, population_fit=False, source_physics_derived=False,
              completed=True), 'results.json')
    print('Completed', sum(r['numerical_pass'] for r in rows), 'numerical passes', flush=True)

if __name__ == '__main__':
    main()
