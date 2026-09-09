"""Closed Gaussian cloud/photon Hamiltonian: excitation is not occupation growth."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = Path(os.environ.get('PHOTON_GRAVITON_RESULTS', ROOT / 'research_work/generated')) / HERE.name


def optical(x, radius, eta):
    excess = eta * np.exp(-(x / radius)**2) / radius**3
    return 1 + excess, -2*x*excess/radius**2, excess*(-3/radius + 2*x*x/radius**3)


def run(case, solver, cfg):
    count = case['pairs']
    inertia = cfg['inertia']
    eta = cfg['optical_coupling']
    x0 = -cfg['source_distance'] - cfg['pair_spacing']*np.arange(count)
    p0 = case['energy_per_photon'] * optical(x0, cfg['initial_radius'], eta)[0]
    # One representative ray of each reflection-symmetric pair, weighted by two.
    initial = np.r_[cfg['initial_radius'], 0., x0, p0, 0.]

    def rhs(t, y):
        radius, momentum = y[:2]
        assert radius > 0
        x = y[2:2+count]
        p = y[2+count:2+2*count]
        n, nx, nr = optical(x, radius, eta)
        speed = momentum/inertia
        force = 2*np.sum(p*nr/n**2)
        return np.r_[speed, 1/radius**3 - 1/radius**2 + force,
                     1/n, p*nx/n**2, force*speed]

    times = np.linspace(0, cfg['end_time'], 2001)
    sol = solve_ivp(rhs, (0, cfg['end_time']), initial, method='DOP853', t_eval=times, **solver)
    assert sol.success, sol.message
    radius, momentum = sol.y[:2]
    x = sol.y[2:2+count]
    p = sol.y[2+count:2+2*count]
    n = optical(x, radius, eta)[0]
    photon = 2*np.sum(p/n, axis=0)
    cloud = momentum**2/(2*inertia) + 0.5/radius**2 - 1/radius
    # Constant seed rest energy (2/3)*inertia is reported separately; subtracting
    # it from the numerical ledger improves sensitivity to the tiny energy transfer.
    total = photon + cloud
    energy_error = float(np.max(abs(total-total[0])))
    work_error = float(np.max(abs(cloud-cloud[0]-sol.y[-1])))
    assert energy_error < cfg['checks']['energy_absolute_error'], energy_error
    assert work_error < cfg['checks']['work_absolute_error'], work_error
    assert np.all(p > 0) and np.all(n > 0)
    assert count == 0 or float(x[:, -1].min()) > 8*float(radius[-1])
    if count == 0:
        assert np.max(abs(radius-1)) == 0
    return {'name': case['name'], 'solver': solver,
            'initial_photon_energy': float(photon[0]), 'final_photon_energy': float(photon[-1]),
            'cloud_energy_gain': float(cloud[-1]-cloud[0]),
            'integrated_photon_work_on_cloud': float(sol.y[-1, -1]),
            'initial_cloud_binding_energy': float(cloud[0]), 'final_cloud_binding_plus_motion': float(cloud[-1]),
            'seed_rest_energy': 2*inertia/3,
            'maximum_energy_absolute_error': energy_error, 'maximum_work_absolute_error': work_error,
            'minimum_radius': float(radius.min()), 'maximum_radius': float(radius.max()),
            'final_radius': float(radius[-1]), 'final_radial_momentum': float(momentum[-1]),
            'occupation_change': 0, 'occupation_change_basis': 'Exact U(1) symmetry; fixed N variational ansatz, not independent numerical evidence.',
            'times': times[::20].tolist(), 'radius_history': radius[::20].tolist(),
            'cloud_energy_history': cloud[::20].tolist()}


def main():
    cfg = json.loads((HERE/'protocol.json').read_text())
    results = []
    for case in cfg['cases']:
        pair = [run(case, solver, cfg) for solver in cfg['solvers']]
        delta = max(abs(pair[0][key]-pair[1][key]) for key in
                    ['cloud_energy_gain', 'final_radius', 'final_radial_momentum'])
        assert delta < cfg['checks']['refinement_absolute_error'], delta
        results.append({'case': case, 'runs': pair, 'maximum_refinement_absolute_difference': delta})
        print(case['name'], 'cloud gain', pair[-1]['cloud_energy_gain'], 'R range',
              pair[-1]['minimum_radius'], pair[-1]['maximum_radius'], 'refinement', delta)
    OUT.mkdir(parents=True, exist_ok=True)
    result = {'scope': cfg['scope'], 'checks_pass': True, 'cases': results,
              'source_hashes': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [HERE/'check.py', HERE/'protocol.json']}}
    (OUT/'bound-cloud-exchange-results.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    main()
