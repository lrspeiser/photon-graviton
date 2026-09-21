"""Run SM-1 without overwriting earlier research or fitting any parameters.

python research_work/results/screened-memory/run.py --output-dir /tmp/sm1
Add --sparc --disks in the full repository. These flags require archived inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import sys
import zipfile
from pathlib import Path
import numpy as np
import model as M

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def algebra_checks():
    x = np.geomspace(1e-14, 1e14, 3001)
    inverse_error = float(np.max(np.abs(M.inverse(M.h(x)) / x - 1)))
    z = np.geomspace(1e-6, 1e2, 81)
    eps = 1e-5
    derivative = (M.W(z * (1 + eps)) - M.W(z * (1 - eps))) / (2 * eps * z)
    energy_error = float(np.max(np.abs(derivative / M.inverse(z) - 1)))
    nodes = np.linspace(-2, 2, 49)
    dx = nodes[1] - nodes[0]
    rng = np.random.default_rng(9271)
    U = rng.normal(0, .04, len(nodes))
    U[[0, -1]] = 0
    q, masses = np.array([-.37, .55]), np.array([.2, .3])

    def energy(u, positions):
        return dx * np.sum(M.W(np.abs(np.diff(u) / dx))) + sum(
            masses[i] * (M.particle_weights(nodes, positions[i])[0] @ u) for i in range(2))

    gradient = np.diff(U) / dx
    force = np.diff(M.mu(np.abs(gradient)) * gradient)
    for i in range(2):
        force -= masses[i] * M.particle_weights(nodes, q[i])[0][1:-1]
    eps = 1e-6
    numeric = []
    for j in range(1, len(nodes) - 1):
        up, um = U.copy(), U.copy()
        up[j] += eps
        um[j] -= eps
        numeric.append((energy(up, q) - energy(um, q)) / (2 * eps))
    field_error = float(np.max(np.abs(np.array(numeric) + force)) / np.max(np.abs(force)))
    pforce, pnum = [], []
    for i in range(2):
        pforce.append(-masses[i] * (M.particle_weights(nodes, q[i])[1] @ U))
        qp, qm = q.copy(), q.copy()
        qp[i] += eps
        qm[i] -= eps
        pnum.append((energy(U, qp) - energy(U, qm)) / (2 * eps))
    particle_error = float(np.max(np.abs(np.array(pnum) + pforce)) / np.max(np.abs(pforce)))
    rows = dict(inverse_relative_error=inverse_error, energy_derivative_relative_error=energy_error,
                field_hamiltonian_gradient_error=field_error, particle_hamiltonian_gradient_error=particle_error,
                analytic_hprime_lower_bound=M.DERIVATIVE_LOWER_BOUND,
                sampled_hprime_minimum=float(M.hprime(x).min()),
                passed=bool(inverse_error < 1e-11 and max(energy_error, field_error, particle_error) < 1e-6))
    return rows


def dynamics():
    runs = [M.toy_run(dt=dt, gamma=gamma) for gamma in (0., 1.) for dt in (.004, .002, .001)]
    groups = [runs[:3], runs[3:]]
    passed = all(g[2]['normalized_energy_error'] < 1e-4
                 and g[2]['normalized_energy_error'] < g[1]['normalized_energy_error'] < g[0]['normalized_energy_error']
                 and all(r['heat_nondecreasing'] for r in g) for g in groups)
    # Explicitly exploratory extension after the declared t=2 experiment showed
    # no resolved body displacement. It is not recast as the preregistered run.
    extended = [M.toy_run(dt=.002, duration=8., gamma=0.),
                M.toy_run(dt=.001, duration=8., gamma=0.),
                M.toy_run(dt=.001, duration=8., gamma=0., n=321),
                M.toy_run(dt=.001, duration=8., gamma=1.)]
    return dict(primary_runs=runs, primary_passed=bool(passed), exploratory_duration_8_runs=extended,
                extension_scope='Post-primary exploration in the same finite Dirichlet box; not a front-speed, stable-orbit, or astrophysical measurement')


def solar():
    rows = []
    for name, au in [('Mercury', .387098), ('Earth', 1.), ('Saturn', 9.5826), ('Neptune', 30.07)]:
        gn = 1.32712440018e20 / (au * 149597870700.) ** 2
        x = gn / M.A0_SI
        rows.append(dict(body=name, radius_AU=au, gN_m_s2=gn,
                         unscreened_excess_m_s2=float(M.A0_SI * np.sqrt(x)),
                         screened_excess_m_s2=float(M.A0_SI * M.excess(x)),
                         fractional_excess=float(M.excess(x) / x)))
    return dict(scope='isolated spherical Sun illustration only; external-field and ephemeris constraints NOT tested', rows=rows)


def sparc():
    guide = ROOT / 'research_work/results/companion-extensions/mond-inventory-results.json'
    raw = ROOT / 'temporal_candidate_audit/data/Rotmod_LTG.zip'
    inventory = json.loads(guide.read_text())['rows']
    records = []
    with zipfile.ZipFile(raw) as archive:
        for row in inventory:
            name = row['galaxy']
            data = np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name + '_rotmod.dat'))))
            data = data[data[:, 0] > 0]
            r, observed = data[:, 0], data[:, 1]
            raw_v2 = data[:, 3] * np.abs(data[:, 3]) + .5 * data[:, 4] ** 2 + .7 * data[:, 5] ** 2
            v2 = np.maximum(raw_v2, 0)
            gn = v2 / r * 1e6 / M.KPC_M
            a = M.A0_SI
            acceleration = dict(baryons=gn, unscreened_root=gn + np.sqrt(a * gn),
                                screened_root=a * M.h(gn / a),
                                simple_MOND=.5 * (gn + np.sqrt(gn * gn + 4 * 8.563e-11 * gn)))
            predictions = {key: np.sqrt(value * r * M.KPC_M / 1e6) for key, value in acceleration.items()}
            records.append(dict(galaxy=name, split=row['split'], points=len(r), clipped_baryon_rows=int((raw_v2 < 0).sum()),
                                rmse_kms={key: float(np.sqrt(np.mean((value - observed) ** 2))) for key, value in predictions.items()}))
    assert len(records) == 149
    summary = {}
    for split in sorted(set(r['split'] for r in records)):
        selected = [r for r in records if r['split'] == split]
        summary[split] = dict(galaxies=len(selected), points=sum(r['points'] for r in selected), models={})
        for name in records[0]['rmse_kms']:
            values = np.array([r['rmse_kms'][name] for r in selected])
            sizes = np.array([r['points'] for r in selected])
            summary[split]['models'][name] = dict(mean_galaxy_RMSE_kms=float(values.mean()),
                median_galaxy_RMSE_kms=float(np.median(values)), pooled_RMSE_kms=float(np.sqrt(np.average(values ** 2, weights=sizes))))
    return dict(scope='exposed-data algebraic screening, NOT the nonspherical AQUAL prediction; no parameters fitted',
                a0_SI=M.A0_SI, summary=summary, per_galaxy=records,
                input_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (guide, raw)})


def disks():
    sys.path.insert(0, str(HERE.parent / 'path-memory'))
    import fields as F
    a = M.A0_SI * M.KPC_M / 1e6
    equations = dict(newtonian=F.NEWTONIAN, unscreened=F.COMPLETION_I, screened=M.equation_for_repository())
    radii = np.array([.3, .5, 1., 2., 4., 8., 12.])
    gx, gw = np.polynomial.legendre.leggauss(128)
    theta, tw = (gx + 1) * np.pi / 4, gw * np.pi / 4
    output = []
    for compactness in (.1, 10.):
        total = compactness * a / F.G
        sigma = lambda r: total / (2 * np.pi) * np.exp(-r)
        rin, height = .002, .15
        inner = float(2 * np.pi * rin ** 2 * np.sum(tw * np.sin(theta) * np.cos(theta) *
                sigma(rin * np.sin(theta)) * (-np.expm1(-rin * np.cos(theta) / height))))
        tests = []
        for label, nr, nth, rout in [('coarse', 128, 24, 80.), ('fine', 256, 48, 80.),
                                     ('boundary', 273, 48, 160.)]:
            grid = F.Q.Grid(rin, rout, nr, nth)
            cell, _wrong_inner = F.Q.exponential_disk_cell_masses(grid, sigma, height, nq_r=8, nq_u=12)
            for name, equation in equations.items():
                if label == 'boundary' and name != 'screened':
                    continue
                solver = F.Solver(grid, cell.copy(), inner, equation)
                solver.solve(a if name != 'newtonian' else 0.)
                tests.append(dict(grid=label, equation=name, nr=nr, nth=nth, rout=rout,
                    converged=bool(solver.converged), residual=float(solver.residual_off_pin),
                    iterations=len(solver.history), mass_relative_error=float(solver.m_tot / total - 1),
                    radial_acceleration=solver.radial_gradient_at(radii).tolist()))
        def values(grid, equation):
            return np.array(next(t for t in tests if t['grid'] == grid and t['equation'] == equation)['radial_acceleration'])
        refinement = float(np.max(np.abs(values('coarse', 'screened') / values('fine', 'screened') - 1)))
        boundary = float(np.max(np.abs(values('boundary', 'screened') / values('fine', 'screened') - 1)))
        departure = values('fine', 'screened') / (a * M.h(values('fine', 'newtonian') / a)) - 1
        output.append(dict(compactness=compactness, radii=radii.tolist(), solves=tests,
            screened_refinement_relative=refinement, screened_boundary_relative=boundary,
            screened_PDE_vs_algebraic_relative=departure.tolist(),
            passed=bool(all(t['converged'] for t in tests) and refinement < .02 and boundary < .005)))
    return dict(scope='synthetic thick exponential disks, not SPARC source reconstructions or lensing', cases=output)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--sparc', action='store_true')
    parser.add_argument('--disks', action='store_true')
    args = parser.parse_args()
    if (args.output_dir / 'results.json').exists():
        raise FileExistsError('Use a fresh output directory; previous evidence is never overwritten')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = dict(experiment='SM-1', baseline='28f154242aa6e6ae010a97af8124d01cc63fed3b',
                   a0_SI=M.A0_SI, algebra=algebra_checks(), dynamics=dynamics(), solar=solar())
    if args.sparc:
        results['sparc'] = sparc()
    if args.disks:
        results['disks'] = disks()
    results['source_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in (HERE / 'model.py', HERE / 'run.py')}
    (args.output_dir / 'results.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    summary = {k: v for k, v in results.items() if k not in ('sparc', 'disks')}
    if 'sparc' in results:
        summary['sparc'] = {k: v for k, v in results['sparc'].items() if k != 'per_galaxy'}
    if 'disks' in results:
        summary['disks'] = results['disks']
    print('SM1_SUMMARY=' + json.dumps(summary, allow_nan=False), flush=True)
    numerical_pass = results['algebra']['passed'] and results['dynamics']['primary_passed']
    if 'disks' in results:
        numerical_pass &= all(c['passed'] for c in results['disks']['cases'])
    if not numerical_pass:
        raise SystemExit('SM-1 numerical gate failed; evidence saved unchanged')


if __name__ == '__main__':
    main()
