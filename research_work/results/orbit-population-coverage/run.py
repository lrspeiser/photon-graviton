"""Training-only finite-library coverage, not a stellar likelihood or gravity score."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache'
OUT = CACHE / 'orbit-population-coverage'
OUT.mkdir(parents=True, exist_ok=True)
FILES = {
    'moments': CACHE / 'stellar-catalogs/stellar-errors-conditional_minus_0.017.parquet',
    'eligibility': CACHE / 'training-orbit-launches/seed-eligibility.parquet',
    'paths': CACHE / 'training-orbit-launches/full-field-orbits.npz',
    'seeds': CACHE / 'training-orbit-launches/launches.parquet',
    'field_check': HERE.parent / 'training-orbit-launches/field-check.json',
    'protocol': HERE.parent / 'stellar-orbit-support/protocol.json',
}
RADII = [.25, .5, 1.0]
SPEEDS = [25., 50., 100.]

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def cylindrical_to_cartesian(d, angle):
    phi = d.mean_phi_rad.to_numpy() - angle
    c, s = np.cos(phi), np.sin(phi)
    R = d.mean_R_kpc.to_numpy()
    vr, vp = d.mean_vR_kms.to_numpy(), d.mean_vphi_kms.to_numpy()
    return (np.c_[R*c, R*s, d.mean_z_kpc],
            np.c_[vr*c-vp*s, vr*s+vp*c, d.mean_vz_kms])

def coverage(x, v, states):
    tree = cKDTree(states[:, :3])
    gaps = np.full((len(x), len(RADII)), np.inf)
    counts = np.zeros_like(gaps, dtype=int)
    for i, (pos, vel) in enumerate(zip(x, v)):
        indices = tree.query_ball_point(pos, RADII[-1])
        if not indices:
            continue
        candidates = states[indices]
        dx = np.linalg.norm(candidates[:, :3] - pos, axis=1)
        dv = np.linalg.norm(candidates[:, 3:] - vel, axis=1)
        for j, radius in enumerate(RADII):
            inside = dx <= radius
            counts[i, j] = inside.sum()
            if inside.any():
                gaps[i, j] = dv[inside].min()
    return gaps, counts

def summary(mask, gaps, counts):
    rows = []
    for j, radius in enumerate(RADII):
        spatial = counts[mask, j] > 0
        for speed in SPEEDS:
            rows.append(dict(radius_kpc=radius, velocity_radius_kms=speed,
                             stars=int(mask.sum()),
                             spatially_covered=int(spatial.sum()),
                             position_and_velocity_covered=int((gaps[mask, j] <= speed).sum())))
    return rows

def main():
    hashes = {str(path.relative_to(ROOT)): digest(path) for path in FILES.values()}
    archived = json.loads(FILES['field_check'].read_text())
    path_key = str(FILES['paths'].relative_to(ROOT))
    assert hashes[path_key] == archived['input_hashes'][path_key]
    side = pd.read_parquet(FILES['eligibility'])
    d = pd.read_parquet(FILES['moments'], filters=[('holdout_role', '==', 'training')])
    assert len(d) == 77927 and (d.holdout_role == 'training').all()
    d = d.merge(side, on='source_id', validate='one_to_one').sort_values('source_id').reset_index(drop=True)
    assert len(d) == 27884
    # Exclude launch stars from this descriptive comparison, not just their t=0 samples.
    d = d[~d.selected_as_new_seed].reset_index(drop=True)
    assert len(d) == 27812
    angle = np.radians(json.loads(FILES['protocol'].read_text())['bar_angle_degrees'])
    x, v = cylindrical_to_cartesian(d, angle)
    with np.load(FILES['paths']) as archive:
        states = archive['trajectories']
        assert states.shape == (501, 72, 6) and np.isfinite(states).all()
        seed_ids = archive['source_id']
    seeds = pd.read_parquet(FILES['seeds'])
    assert seeds.source_id.tolist() == seed_ids.tolist()
    assert (seeds.holdout_role == 'training').all()
    np.testing.assert_array_equal(states[0], seeds[['x_kpc', 'y_kpc', 'z_kpc', 'vx_kms', 'vy_kms', 'vz_kms']].to_numpy())
    assert not d.source_id.isin(seed_ids).any()
    gaps, counts = coverage(x, v, states[1:].reshape(-1, 6))
    # Exact brute-force comparison on fixed indices verifies the tree search.
    flat = states[1:].reshape(-1, 6)
    for i in np.linspace(0, len(d)-1, 12, dtype=int):
        dx = np.linalg.norm(flat[:, :3] - x[i], axis=1)
        dv = np.linalg.norm(flat[:, 3:] - v[i], axis=1)
        for j, radius in enumerate(RADII):
            inside = dx <= radius
            expected = dv[inside].min() if inside.any() else np.inf
            assert counts[i, j] == inside.sum() and gaps[i, j] == expected
    assert np.all(np.diff(counts, axis=1) >= 0)
    assert np.all(gaps[:, 1:] <= gaps[:, :-1])
    good = (d.seed_eligibility == 'passes_position_screen').to_numpy()
    cells = []
    for lo, hi in [(.5, 3.5), (3.5, 5), (5, 9)]:
        for zl, zh in [(0, .2), (.2, .5), (.5, 1.5)]:
            mask = good & (d.mean_R_kpc >= lo) & ((d.mean_R_kpc <= hi) if hi == 9 else (d.mean_R_kpc < hi))
            mask &= (abs(d.mean_z_kpc) >= zl) & ((abs(d.mean_z_kpc) <= zh) if zh == 1.5 else (abs(d.mean_z_kpc) < zh))
            cells.append(dict(R_range=[lo, hi], absolute_z_range=[zl, zh],
                              stars=int(mask.sum()), spatially_covered=int((counts[mask, 1] > 0).sum()),
                              position_and_velocity_covered=int((gaps[mask, 1] <= 50).sum())))
    assert sum(row['stars'] for row in cells) == good.sum()
    # Sensitivity to finite saved-time spacing: use every second saved sample.
    coarse_gaps, coarse_counts = coverage(x, v, states[2::2].reshape(-1, 6))
    assert np.all(coarse_counts <= counts)
    assert np.all(coarse_gaps >= gaps)
    audit = d[['source_id', 'seed_eligibility']].copy()
    for j, radius in enumerate(RADII):
        audit[f'min_velocity_gap_within_{radius}_kpc'] = gaps[:, j]
        audit[f'samples_within_{radius}_kpc'] = counts[:, j]
    audit.to_parquet(OUT / 'training-coverage.parquet', index=False)
    result = dict(scope='Finite orbit library coverage of training posterior means; not a likelihood or gravity rejection.',
                  training_stars_excluding_seeds=len(d), position_screen_passes=int(good.sum()),
                  saved_path_states=36000, launch_time_samples_excluded=True,
                  launch_stars_excluded=True, all_training=summary(np.ones(len(d), bool), gaps, counts),
                  position_screened_training=summary(good, gaps, counts),
                  half_time_sampling=summary(good, coarse_gaps, coarse_counts), cells=cells,
                  brute_force_comparisons=12, holdouts_opened=False,
                  measurement_uncertainties_integrated=False, population_weights_fitted=False,
                  cutoffs_are_diagnostic_not_error_bars=True, source_hashes=hashes,
                  code_sha256=digest(Path(__file__)), audit_sha256=digest(OUT / 'training-coverage.parquet'))
    for path in FILES.values():
        assert digest(path) == hashes[str(path.relative_to(ROOT))]
    (HERE / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(cells=cells, screened=result['position_screened_training']), indent=2))

if __name__ == '__main__':
    main()
