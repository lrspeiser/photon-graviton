"""Leading ordinary spin-two soft radiation conditional on elastic photon scattering.

No cosmological matter density, new population, or hard scattering rate is assumed.
"""
from pathlib import Path
import json
import os
import numpy as np
from scipy.constants import G, hbar, c, electron_volt
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(os.environ.get('PHOTON_GRAVITON_RESULTS', ROOT/'research_work/generated'))/'soft-graviton'
ETA = np.array([-1., -1., 1., 1.])


def hard_momenta(theta, mass):
    """E=1, exact elastic gamma+rest target -> gamma+target, signature +---."""
    direction = np.array([np.sin(theta), 0., np.cos(theta)])
    final_energy = 1/(1+(1-np.cos(theta))/mass)
    k = np.array([1., 0., 0., 1.])
    p = np.array([mass, 0., 0., 0.])
    kp = final_energy*np.r_[1., direction]
    return np.array([k, p, kp, p+k-kp])


def tt_kernel(momenta, directions):
    """Sum of squared two unit-normalized physical tensor polarizations.

    Project each spatial momentum first to avoid subtraction of large longitudinal
    components near a photon direction. Exactly collinear nodes are avoided.
    """
    tensor = np.zeros((len(directions), 3, 3))
    for sign, p in zip(ETA, momenta):
        cosine = directions @ p[1:]
        transverse = p[1:] - cosine[:, None]*directions
        denominator = p[0]-cosine
        tensor += sign*np.einsum('ni,nj->nij', transverse, transverse)/denominator[:, None, None]
    return np.einsum('nij,nij->n', tensor, tensor)-.5*np.trace(tensor, axis1=1, axis2=2)**2


def angular_integral(theta, order, mass=None):
    mu, weights = np.polynomial.legendre.leggauss(order)
    phi = (np.arange(2*order)+.5)*np.pi/order
    mu_grid, phi_grid = np.meshgrid(mu, phi, indexing='ij')
    radius = np.sqrt(1-mu_grid**2)
    directions = np.stack([radius*np.cos(phi_grid), radius*np.sin(phi_grid), mu_grid], axis=-1).reshape(-1, 3)
    if mass is None:
        # Strict M/E -> infinity TT limit: target terms vanish, E'=E.
        momenta = np.array([[1.,0.,0.,1.], [1.,0.,0.,0.],
                            [1.,np.sin(theta),0.,np.cos(theta)], [1.,0.,0.,0.]])
    else:
        momenta = hard_momenta(theta, mass)
    kernel = tt_kernel(momenta, directions).reshape(order, 2*order)
    assert np.min(kernel)>-1e-12
    return float(np.sum(weights[:, None]*kernel)*np.pi/order)


def main():
    rng = np.random.default_rng(8102026)
    ward_errors, shell_errors, polarization_errors = [], [], []
    for _ in range(100):
        momenta = hard_momenta(rng.uniform(.05,3.09), rng.uniform(10,1000))
        direction = rng.normal(size=3); direction /= np.linalg.norm(direction)
        qlow = np.r_[1., -direction]
        denominators = momenta @ qlow
        tensor = np.einsum('a,ai,aj,a->ij', ETA, momenta, momenta, 1/denominators)
        ward_errors.append(float(np.max(abs(qlow @ tensor))))
        mass = momenta[1,0]
        shell_errors.append(float(abs(momenta[3,0]**2-np.dot(momenta[3,1:],momenta[3,1:])-mass**2)/mass**2))
        a = np.cross(direction, [0.,0.,1.]); a /= np.linalg.norm(a)
        b = np.cross(direction,a)
        plus = (np.outer(a,a)-np.outer(b,b))/np.sqrt(2)
        cross = (np.outer(a,b)+np.outer(b,a))/np.sqrt(2)
        direct = np.sum(tensor[1:,1:]*plus)**2+np.sum(tensor[1:,1:]*cross)**2
        projected = tt_kernel(momenta, direction[None,:])[0]
        polarization_errors.append(float(abs(direct-projected)))
    assert max(ward_errors)<1e-10 and max(shell_errors)<1e-12
    assert max(polarization_errors)<1e-10
    angles = []
    for degrees in [0., 10., 30., 90., 180.]:
        theta = np.deg2rad(degrees)
        coarse = angular_integral(theta, 128)
        fine = angular_integral(theta, 256)
        assert abs(fine-coarse)<.005*max(1.,fine)
        assert fine <= 32*np.pi*(1+1e-12)
        angles.append({'photon_angle_degrees': degrees, 'I_heavy': fine,
                       'quadrature_absolute_change': abs(fine-coarse)})
    assert angles[0]['I_heavy']==0.
    # Independent analytic backscatter: tensor kernel=2*cos(graviton angle)^2.
    assert abs(angles[-1]['I_heavy']-8*np.pi/3)<1e-10
    finite_mass = []
    heavy = angular_integral(np.pi/2, 256)
    for mass in [1e2, 1e3, 1e4]:
        result = angular_integral(np.pi/2, 256, mass)
        finite_mass.append({'M_over_E': mass, 'I': result, 'distance_from_heavy': abs(result-heavy)})
    assert all(b['distance_from_heavy']<a['distance_from_heavy'] for a,b in zip(finite_mass[:-1],finite_mass[1:]))
    planck_eV = np.sqrt(hbar*c**5/G)/electron_volt
    bound_C = (16/np.pi)/planck_eV**2  # E=1 eV, strict heavy limit
    delta = 1e-3  # declared illustrative soft ceiling, not a measured threshold
    alpha = .0002488993286382367  # reused empirical alpha / Mpc, not Hubble expansion
    moments = []
    for minimum in [1e-6,1e-12,1e-30]:
        loglo, loghi = np.log(minimum), np.log(delta)
        computed = [quad(lambda t: np.exp(j*t), loglo, loghi, epsabs=1e-15)[0] for j in [0,1,2]]
        exact = [np.log(delta/minimum), delta-minimum, (delta**2-minimum**2)/2]
        assert np.allclose(computed, exact, rtol=1e-10, atol=1e-15)
        moments.append({'fractional_soft_floor': minimum, 'P_over_C': computed[0],
                        'energy_fraction_over_C': computed[1], 'second_moment_over_C': computed[2]})
    results = {
        'scope': 'Leading soft ordinary-graviton radiation per elastic hard event; not a full cross section, redshift theory, or observational validation.',
        'conventions': {'signature': '+---', 'kappa_squared': '32 pi G', 'polarization_norm': 1,
                        'units_in_derivation': 'hbar=c=1', 'Planck_energy_eV': float(planck_eV)},
        'checks': {'random_kinematics_cases': 100, 'max_Ward_residual': max(ward_errors),
                   'max_relative_target_shell_error': max(shell_errors),
                   'max_polarization_sum_error': max(polarization_errors),
                   'forward_leading_soft_term_zero': True, 'backscatter_I_equals_8pi_over_3': True},
        'angular_results': angles, 'finite_mass_convergence': finite_mass,
        'infrared_moments': moments,
        'illustration': {'photon_energy_eV': 1., 'soft_ceiling_fraction': delta,
                          'upper_C_heavy': bound_C, 'upper_soft_energy_fraction_per_hard_event': bound_C*delta,
                          'empirical_alpha_per_Mpc': alpha,
                          'necessary_hard_events_per_Mpc_from_soft_budget_bound': alpha/(bound_C*delta),
                          'warning': 'Necessary soft-energy throughput diagnostic only. Such a huge scattering rate invalidates transparent independent propagation; no gas density or successful cumulative-redshift law is inferred.'},
        'unresolved': ['Finite-energy two-photon one-graviton amplitude and phase space',
                       'Physical target and hard differential scattering rate',
                       'Inclusive real/virtual infrared treatment',
                       'Photon spectrum, direction and duration after propagation',
                       'Capture, stable storage and joint motion/lensing response']}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'soft-graviton-results.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(results['illustration'], indent=2))


if __name__ == '__main__':
    main()
