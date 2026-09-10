"""Field-resolution and occupation audit for completed duration probes."""
import json
from pathlib import Path
import numpy as np
import run

HERE = Path(__file__).resolve().parent
d = json.loads((HERE / 'results.json').read_text())
assert d['completed']
for path, expected in d['source_hashes'].items():
    assert run.digest(run.ROOT/path) == expected
full = run.orbit.full
axispath = full.CACHE/'axisymmetric-reference-80-L128-N512.npz'
coarsepath, finepath = full.CACHE/'fine.npz', full.CACHE/'finer.npz'
axis, coarse, fine = [full.load_field(p) for p in [axispath, coarsepath, finepath]]
hashes = {str(p.relative_to(run.ROOT)): run.digest(p) for p in [Path(__file__), HERE/'results.json', axispath, coarsepath, finepath]}
rows = []
for row in d['rows']:
    if row['error'] is not None:
        rows.append(dict(seed_index=row['seed_index'], assessed=False, reason=row['error']))
        continue
    path = run.CACHE/f"orbit-{row['seed_index']}.npz"
    assert run.digest(path) == row['trajectory_sha256']
    hashes[str(path.relative_to(run.ROOT))] = run.digest(path)
    with np.load(path) as f:
        y = f['trajectory']
        np.testing.assert_array_equal(f['times'], run.TIMES)
    errors = []
    for part in np.array_split(y[:, :3], 128):
        a, b, c = axis.evaluate(part)[1], coarse.evaluate(part)[1], fine.evaluate(part)[1]
        errors.extend((np.linalg.norm(c-b, axis=1)/np.linalg.norm(a+c, axis=1)).tolist())
    angle = np.unwrap(np.arctan2(y[:, 1], y[:, 0]))
    radius = np.hypot(y[:, 0], y[:, 1])
    rows.append(dict(seed_index=row['seed_index'], assessed=True, positions=len(y),
                     max_extra_force_refinement=float(max(errors)), field_pass=bool(max(errors)<.01),
                     net_bar_frame_turns=float((angle[-1]-angle[0])/(2*np.pi)),
                     bar_frame_angle_range_turns=float(np.ptp(angle)/(2*np.pi)),
                     cylindrical_radius_range_kpc=[float(radius.min()), float(radius.max())],
                     numerical_pass=row['numerical_pass']))
    assert run.tv(y[1:501], y[1:501]) == 0
    # Independently count histogram cells using digitized indices.
    x = y[1:501, :3]
    indices = [np.digitize(a, edges)-1 for a, edges in zip(
        [np.hypot(x[:, 0], x[:, 1]), x[:, 2], np.arctan2(x[:, 1], x[:, 0])],
        [run.R_EDGES, run.Z_EDGES, run.PHI_EDGES])]
    linear = np.ravel_multi_index(indices, (8, 8, 8))
    check = np.bincount(linear, minlength=512)/500
    np.testing.assert_array_equal(check, run.occupation(y[1:501]))
    print('field', row['seed_index'], max(errors), flush=True)
with np.load(run.SOURCE) as f:
    original = f['trajectories']
baseline = []
for index in range(72):
    y = original[:, index]
    baseline.append(dict(seed_index=index, spatial_TV=run.tv(y[1:251], y[251:501]),
                         R_z_TV=run.tv(y[1:251], y[251:501], False)))
for path, expected in hashes.items():
    assert run.digest(run.ROOT/path) == expected
run.save(dict(rows=rows, baseline_72=baseline, input_hashes=hashes,
              numerical_passes=sum(r.get('numerical_pass', False) for r in rows),
              field_passes=sum(r.get('field_pass', False) for r in rows),
              occupation_is_bin_dependent=True, population_stationarity_proved=False,
              holdouts_opened=False), 'assessment.json')
