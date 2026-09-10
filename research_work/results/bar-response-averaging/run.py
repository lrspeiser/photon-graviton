"""Does nonlinear response commute with azimuth averaging of the real bar?"""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORBIT = HERE.parent / 'rotating-bar-orbits'
sys.path.insert(0, str(ORBIT))
spec = importlib.util.spec_from_file_location('bar_orbit_operator', ORBIT / 'run.py')
operator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(operator)
from field import G

PARAMFILE = HERE.parent / 'joint-galaxy-audit/results.json'
PARAM = json.loads(PARAMFILE.read_text())['sparc']['parameters']
A, POWER = PARAM['A'], PARAM['p']
ASTAR = PARAM['a_star_m_s2'] * 3.085677581491367e19 / 1e6


def response(acc, power=POWER):
    g = np.linalg.norm(acc, axis=-1)
    assert np.min(g) > 0
    return A * (g / ASTAR)[..., None] ** (power - 1) * acc


def sample(field, R, z, nphi):
    phi = np.arange(nphi) * 2 * np.pi / nphi
    xyz = np.c_[R * np.cos(phi), R * np.sin(phi), np.full(nphi, z)]
    _, acc = field.evaluate(xyz)
    cyl = np.c_[acc[:, 0] * np.cos(phi) + acc[:, 1] * np.sin(phi),
                -acc[:, 0] * np.sin(phi) + acc[:, 1] * np.cos(phi), acc[:, 2]]
    mean = cyl.mean(axis=0)
    full = response(cyl)
    mean_response = full.mean(axis=0)
    response_mean = response(mean)
    norm = np.linalg.norm(response_mean)
    commutator = np.linalg.norm(mean_response - response_mean) / norm
    rms = np.sqrt(np.mean(np.sum((full - mean_response) ** 2, axis=1))) / norm
    peak = np.max(np.linalg.norm(full - response_mean, axis=1)) / norm
    # Vector components must be averaged in the rotating cylindrical basis.
    # A Cartesian vector average around a ring would erase radial gravity.
    linear_control = np.linalg.norm(response(cyl, power=1).mean(axis=0) - response(mean, power=1)) / norm
    assert linear_control < 1e-12
    assert abs(mean_response[1]) / norm < 1e-10
    return dict(R_kpc=R, z_kpc=z, nphi=nphi,
                mean_baryonic_acceleration_cyl=mean.tolist(),
                response_of_mean_cyl=response_mean.tolist(),
                mean_of_response_cyl=mean_response.tolist(),
                fractional_mean_response_difference=float(commutator),
                fractional_azimuthal_response_rms=float(rms),
                fractional_azimuthal_response_peak=float(peak),
                linear_response_commutation_control=float(linear_control))


def main():
    disk = operator.DiskGrid(256, 129, operator.disks(24))
    fields = {o: operator.Field(disk, order=o) for o in [40, 64]}
    allrows = []
    refinement = []
    field_checks = []
    for R in [.5, 1., 2., 3., 5., 8., 12., 20.]:
        for z in [0., .1, .8, 1.5]:
            pair = [sample(fields[64], R, z, n) for n in [128, 256]]
            allrows.append(pair[-1])
            keys = ['fractional_mean_response_difference', 'fractional_azimuthal_response_rms']
            diff = max(abs(pair[1][k] - pair[0][k]) for k in keys)
            assert diff < 1e-7
            refinement.append(diff)
            low = sample(fields[40], R, z, 256)
            field_checks.append(dict(R_kpc=R, z_kpc=z,
                mean_difference_change=abs(low[keys[0]] - pair[-1][keys[0]]),
                azimuth_rms_change=abs(low[keys[1]] - pair[-1][keys[1]])))
        print('Completed R=', R, flush=True)
    # Independent exact m=0 bar control, keeping all radial/angular l modes.
    axis = operator.Field(disk, order=64)
    bar = axis.bar
    mask = bar.m == 0
    axis.bar = operator.FastMultipole(bar.r, bar.coeff[:, mask], bar.l[mask], bar.m[mask])
    checks = []
    for R, z in [(1., .1), (3., .8), (8., 0.)]:
        r = sample(axis, R, z, 128)
        full = next(row for row in allrows if row['R_kpc'] == R and row['z_kpc'] == z)
        mean_error = np.linalg.norm(np.array(full['mean_baryonic_acceleration_cyl']) - r['mean_baryonic_acceleration_cyl']) / np.linalg.norm(r['mean_baryonic_acceleration_cyl'])
        assert mean_error < 1e-10
        assert r['fractional_mean_response_difference'] < 1e-10
        assert r['fractional_azimuthal_response_rms'] < 1e-10
        checks.append(dict(R_kpc=R, z_kpc=z, exact_m0_average_error=float(mean_error),
                           response_commutation_error=r['fractional_mean_response_difference']))
    # The azimuthal derivative averages to zero. Differentiate mean cylindrical
    # B to measure the m=0 Poisson source, rather than equating B with gravity.
    def source(R, z, h, nphi):
        around = [sample(fields[64], rr, zz, nphi) for rr, zz in
                  [(R+h,z), (R-h,z), (R,z+h), (R,z-h)]]
        sources_out = []
        for key in ['mean_of_response_cyl', 'response_of_mean_cyl']:
            div = ((R+h)*around[0][key][0] - (R-h)*around[1][key][0]) / (2*h*R)
            div += (around[2][key][2] - around[3][key][2]) / (2*h)
            # g_extra=-grad(Phi_extra), so rho_extra=-div(B)/(4*pi*G).
            sources_out.append(-div / (4*np.pi*G))
        return np.asarray(sources_out)
    density_rows = []
    for R in [.5, 1., 3., 8., 20.]:
        for z in [0., .1, .8, 1.5]:
            coarse = source(R,z,.002,128)
            fine = source(R,z,.001,128)
            angular = source(R,z,.001,256)
            scale = max(np.max(abs(angular)), 1.)
            step_error = float(np.max(abs(fine-coarse))/scale)
            angular_error = float(np.max(abs(angular-fine))/scale)
            assert step_error < .003 and angular_error < 1e-4
            density_rows.append(dict(R_kpc=R,z_kpc=z,
                mean_full_response_effective_density_Msun_per_kpc3=float(angular[0]),
                response_after_averaging_effective_density_Msun_per_kpc3=float(angular[1]),
                fractional_source_difference=float(abs(angular[0]-angular[1])/scale),
                radial_vertical_step_change_fraction=step_error,
                azimuth_refinement_change_fraction=angular_error))
    sources = [Path(__file__), PARAMFILE, ORBIT / 'run.py', ORBIT / 'fast_multipole.py',
               HERE.parent / 'bar-field-foundation/field.py', disk.path]
    sources += [ROOT / f'research_work/data-cache/bar-field/bar-L{o}.npz' for o in [40, 64]]
    sources += [ROOT / 'research_work/data-cache/bar-field/nuclei-L16.npz']
    result = dict(scope='Source-vector averaging audit; not a Poisson solve, acceleration prediction or stellar fit.',
                  response_parameters=PARAM, units='kpc, km/s; acceleration (km/s)^2/kpc',
                  response='B=A*(|g_b|/a_star)^(p-1)*g_b; div(g_extra)=div(B)',
                  baseline='Existing conservative order-24 disk spline, nuclei L16, bar L64 and central mass; no halo.',
                  angular_average_resolution=[128,256], max_azimuth_refinement_change=max(refinement),
                  rows=allrows, bar_order_sensitivity=field_checks, axisymmetric_controls=checks,
                  azimuth_mean_poisson_source_checks=density_rows,
                  input_hashes={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  holdouts_opened=False, gravity_fit_performed=False,
                  poisson_source_divergence_computed=True, full_poisson_potential_solved=False,
                  strong_case_established=False)
    (HERE / 'results.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    for label, subset in [('bulge', [r for r in allrows if r['R_kpc'] <= 3]),
                          ('disk', [r for r in allrows if r['R_kpc'] >= 5])]:
        print(label, 'max mean difference', max(r['fractional_mean_response_difference'] for r in subset),
              'max azimuth rms', max(r['fractional_azimuthal_response_rms'] for r in subset))
    print('Angular convergence', max(refinement))
    print('Worst source difference', max(density_rows, key=lambda r:r['fractional_source_difference']))


if __name__ == '__main__':
    main()
