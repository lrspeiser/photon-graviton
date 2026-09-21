"""SM-2: exposed-data, fixed-law AQUAL benchmark, not a new theory.

Run four disjoint shards with --shard 0..3 --shards 4, then --aggregate DIR.
Every attempted galaxy is retained; unresolved sources are never silently dropped.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path
import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'path-memory'))
import pm2a as P
import fields as F
sys.path.insert(0, str(HERE))
import model as M

NAMES = ('newtonian', 'bekenstein_root', 'SM1_screened', 'simple_mu')
A0 = M.A0_SI * M.KPC_M / 1e6
EQUATIONS = dict(zip(NAMES, (F.NEWTONIAN, F.COMPLETION_I, M.equation_for_repository(), F.SIMPLE)))


def write_new(path, value):
    with Path(path).open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance():
    inputs = [ROOT / 'temporal_candidate_audit/data/Rotmod_LTG.zip',
              ROOT / 'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt',
              HERE.parent / 'companion-extensions/mond-inventory-results.json']
    sources = [Path(mod.__file__) for mod in (P, P.BAR, P.I, F, F.Q, M)] + [Path(__file__)]
    return dict(executed_git_sha=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                declared_baseline='abd026f324c6bf90f71a68f0b0d6a4f8b889c6c3',
                python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
                source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources})


def relative(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), .01 * max(float(np.max(np.abs(b))), 1e-300))))


def solve_galaxy(d):
    started = time.monotonic()
    src = P.build_source(d)
    sigma = P.sigma_of(src)
    comp = src['comp']
    rows = d['rotmod'][d['rotmod'][:, 0] > 0]
    radii, observed = rows[:, 0], rows[:, 1]
    total, height = src['m_total_primary'], .1 * src['rd']
    rin = .002 * min(src['rd'], float(radii.min()))
    rm = np.sqrt(F.G * total / A0)
    rout = max(80 * src['rd'], 20 * src['h_gas'], 10 * float(radii.max()), 80 * rm)
    nr_boundary = int(np.ceil(384 * np.log(2 * rout / rin) / np.log(rout / rin)))
    inner_disk = P.thick_disk_inner_mass(sigma, rin, height, nodes=128)
    records = []
    for label, nr, nth, outer in [('coarse', 192, 32, rout), ('fine', 384, 64, rout),
                                  ('boundary', nr_boundary, 64, 2 * rout)]:
        grid = F.Q.Grid(rin, outer, nr, nth)
        cell, _discarded_cylinder_inner = F.Q.exponential_disk_cell_masses(grid, sigma, height, nq_r=8, nq_u=12)
        inner = inner_disk
        if comp['m_bulge'] is not None:
            bulge, bulge_inner = F.Q.spherical_cell_masses(grid, comp['m_bulge'])
            cell = cell + bulge
            inner += bulge_inner
        if not (grid.rf[1] < radii.min() and radii.max() < grid.rf[-2]):
            raise ValueError('Requested force readout would extrapolate')
        for name, eq in EQUATIONS.items():
            if label == 'boundary' and name != 'SM1_screened':
                continue
            solver = F.Solver(grid, cell.copy(), inner, eq)
            scale = 0. if name == 'newtonian' else (8.563e-11 * M.KPC_M / 1e6 if name == 'simple_mu' else A0)
            solver.solve(scale)
            acceleration = np.asarray(solver.radial_gradient_at(radii))
            valid = np.isfinite(acceleration) & (acceleration > 0)
            velocity = np.sqrt(radii * np.maximum(acceleration, 0))
            records.append(dict(grid=label, equation=name, nr=nr, nth=nth, rin=rin, rout=outer,
                converged=bool(solver.converged), residual=float(solver.residual_off_pin),
                iteration_change=float(solver.iteration_change), iterations=len(solver.history),
                mass_relative_error=float(solver.m_tot / total - 1),
                radial_acceleration=acceleration.tolist(),
                velocity_kms=[float(v) if ok else None for v, ok in zip(velocity, valid)],
                nonpositive_or_invalid_rows=int((~valid).sum()),
                RMSE_kms=float(np.sqrt(np.mean((velocity - observed) ** 2))) if valid.all() else None))
    assert src['hash'] == P.source_fingerprint(src), 'Source mutated during calculation'
    def row(grid, name):
        return next(r for r in records if r['grid'] == grid and r['equation'] == name)
    refinements = {n: relative(row('coarse', n)['radial_acceleration'], row('fine', n)['radial_acceleration']) for n in NAMES}
    boundary = relative(row('boundary', 'SM1_screened')['radial_acceleration'], row('fine', 'SM1_screened')['radial_acceleration'])
    valid_circular = all(r['nonpositive_or_invalid_rows'] == 0 for r in records)
    mass_pass = max(abs(r['mass_relative_error']) for r in records) < 1e-3
    numerical = all(r['converged'] for r in records) and mass_pass and max(refinements.values()) < .02 and boundary < .005
    raw_v2 = rows[:, 3] * np.abs(rows[:, 3]) + .5 * rows[:, 4] ** 2 + .7 * rows[:, 5] ** 2
    tabulated = np.sqrt(np.maximum(raw_v2, 0))
    gn = np.asarray(row('fine', 'newtonian')['radial_acceleration'])
    alg = A0 * M.h(np.maximum(gn, 0) / A0)
    difference = relative(row('fine', 'SM1_screened')['radial_acceleration'], alg)
    return dict(galaxy=d['name'], split=d['split'], points=len(radii), radii_kpc=radii.tolist(),
        observed_kms=observed.tolist(), observational_error_kms=rows[:, 2].tolist(),
        source_fingerprint=src['hash'], total_source_Msun=total, height_kpc=height,
        stellar_tail_fraction=src['extrapolated_total_fraction'],
        tabulated_baryon_velocity_kms=tabulated.tolist(),
        tabulated_baryon_RMSE_kms=float(np.sqrt(np.mean((tabulated - observed) ** 2))),
        reconstruction_velocity_RMS_difference_kms=float(np.sqrt(np.mean((np.sqrt(np.maximum(radii * gn, 0)) - tabulated) ** 2))),
        tabulated_clipped_rows=int((raw_v2 < 0).sum()), solves=records, force_refinement=refinements,
        screened_boundary_change=boundary, screened_PDE_vs_same_density_algebraic=difference,
        numerical_certified=bool(numerical), all_models_circular=valid_circular,
        matched_certified=bool(numerical and valid_circular), mass_pass=bool(mass_pass),
        fine_RMSE_kms={n: row('fine', n)['RMSE_kms'] for n in NAMES}, seconds=time.monotonic() - started)


def run_shard(out, shard, nshards):
    if out.exists():
        raise FileExistsError('Use a fresh shard output directory')
    out.mkdir(parents=True)
    all_galaxies = sorted(P.I.sparc_galaxies(), key=lambda d: d['name'])
    assert len(all_galaxies) == 149 and len({d['name'] for d in all_galaxies}) == 149
    selected = [d for i, d in enumerate(all_galaxies) if i % nshards == shard]
    records = []
    for d in selected:
        try:
            record = solve_galaxy(d)
        except Exception:
            record = dict(galaxy=d['name'], split=d['split'], error=traceback.format_exc(), matched_certified=False)
        write_new(out / (d['name'] + '.json'), record)
        records.append(record)
        print('SM2_GALAXY ' + json.dumps({k: record[k] for k in ('galaxy', 'matched_certified', 'fine_RMSE_kms', 'error') if k in record}), flush=True)
    write_new(out / 'shard.json', dict(shard=shard, shards=nshards, expected_galaxies=[d['name'] for d in selected],
                                      provenance=provenance(), rows=records))


def bootstrap_delta(rows, left, right):
    delta = np.array([r['fine_RMSE_kms'][left] - r['fine_RMSE_kms'][right] for r in rows])
    if not len(delta):
        return None
    rng = np.random.default_rng(210926)
    samples = delta[rng.integers(0, len(delta), size=(10000, len(delta)))].mean(1)
    return dict(mean_kms=float(delta.mean()), bootstrap_95_percent_interval_kms=np.quantile(samples, [.025, .975]).tolist(),
                units='galaxy RMSE difference, not significance including observational systematics')


def summarize(rows):
    result = {}
    for split in ('train', 'validation', 'test', 'all'):
        selected = [r for r in rows if split == 'all' or r['split'] == split]
        models = {}
        for name in NAMES:
            complete = [r for r in selected if r['fine_RMSE_kms'][name] is not None]
            values = np.array([r['fine_RMSE_kms'][name] for r in complete])
            weights = np.array([r['points'] for r in complete])
            models[name] = dict(complete_galaxies=len(complete),
                mean_galaxy_RMSE_kms=float(values.mean()) if len(values) else None,
                median_galaxy_RMSE_kms=float(np.median(values)) if len(values) else None,
                pooled_RMSE_kms=float(np.sqrt(np.average(values ** 2, weights=weights))) if len(values) else None)
        paired = [r for r in selected if all(r['fine_RMSE_kms'][n] is not None for n in NAMES)]
        result[split] = dict(galaxies=len(selected), points=sum(r['points'] for r in selected), models=models,
            paired_galaxies=len(paired), screened_minus_root=bootstrap_delta(paired, 'SM1_screened', 'bekenstein_root'),
            screened_minus_simple=bootstrap_delta(paired, 'SM1_screened', 'simple_mu'))
    return result


def aggregate(source, out):
    files = sorted(source.rglob('shard.json'))
    payloads = [json.loads(p.read_text()) for p in files]
    if len(payloads) != 4 or {p['shard'] for p in payloads} != set(range(4)):
        raise ValueError('Exactly four distinct completed shards are required')
    commits = {p['provenance']['executed_git_sha'] for p in payloads}
    if len(commits) != 1:
        raise ValueError('Mixed code versions are not a single experiment')
    rows = sorted([r for p in payloads for r in p['rows']], key=lambda r: r['galaxy'])
    expected = {g for p in payloads for g in p['expected_galaxies']}
    if len(rows) != 149 or len(expected) != 149 or {r['galaxy'] for r in rows} != expected:
        raise ValueError('Missing or duplicate galaxies')
    attempted = [r for r in rows if 'error' not in r]
    certified = [r for r in attempted if r['matched_certified']]
    if out.exists():
        raise FileExistsError('Use a fresh aggregate output directory')
    out.mkdir(parents=True)
    summary = dict(experiment='SM-2', executed_git_sha=next(iter(commits)),
        scope='same-density AQUAL benchmark on exposed data; no microscopic novelty or complete gravity claim',
        attempted_galaxies=149, completed_galaxies=len(attempted), matched_certified_galaxies=len(certified),
        unresolved_galaxies=[r['galaxy'] for r in rows if not r['matched_certified']],
        implementation_errors=[r for r in rows if 'error' in r],
        all_source_scores_provisional=summarize(attempted), common_certified_subset_scores=summarize(certified),
        tabulated_vs_reconstructed_newtonian=dict(
            mean_tabulated_RMSE_kms=float(np.mean([r['tabulated_baryon_RMSE_kms'] for r in attempted])) if attempted else None,
            mean_reconstruction_velocity_RMS_difference_kms=float(np.mean([r['reconstruction_velocity_RMS_difference_kms'] for r in attempted])) if attempted else None),
        all_numerical_certified=len(certified) == 149)
    write_new(out / 'summary.json', summary)
    write_new(out / 'all-galaxies.json', dict(provenance=payloads[0]['provenance'], rows=rows))
    print('SM2_SUMMARY=' + json.dumps(summary), flush=True)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--shard', type=int, default=0)
    parser.add_argument('--shards', type=int, default=4)
    parser.add_argument('--aggregate', type=Path)
    args = parser.parse_args()
    if args.aggregate:
        aggregate(args.aggregate, args.output_dir)
    else:
        if args.shards != 4 or not 0 <= args.shard < 4:
            raise ValueError('The declared campaign consists of four shards')
        run_shard(args.output_dir, args.shard, args.shards)


if __name__ == '__main__':
    main()
