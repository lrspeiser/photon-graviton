"""Numerical gates for CL-F1. No archived result is overwritten."""
import json
import numpy as np
from scipy.special import erf

import steady_field as SF
from cluster_lensing import BaryonShells, Geometry, WrittenClusterLens, observables, solve_lens


def relative(a, b):
    return float(np.max(np.abs(np.asarray(a) - b)) / max(float(np.max(np.abs(b))), 1e-300))


def gaussian_mass(r, mass, s):
    x = np.asarray(r) / s
    return mass * (erf(x / np.sqrt(2)) - np.sqrt(2 / np.pi) * x * np.exp(-x*x / 2))


def gaussian_lens(n, npj):
    r = np.geomspace(.02, 800., n)
    return WrittenClusterLens(BaryonShells(r, gaussian_mass(r, 2e14, 100.)), [80., 200.], [2e-8, 1e-8], npj)


def gaussian_exact(R):
    M, s = 2e14, 100.
    mb = M * (-np.expm1(-R**2 / (2*s*s)))
    sb = M / (2*np.pi*s*s) * np.exp(-R**2 / (2*s*s))
    mm, sm = np.zeros_like(R), np.zeros_like(R)
    for w, a in zip([80., 200.], [2e-8, 1e-8]):
        q = s*s + w*w
        F = M*w*w/q * np.exp(-R*R / (2*q))
        mm += a * np.sqrt(2*np.pi)*w*R*R*F / (2*SF.G*q)
        sm += a * np.sqrt(2*np.pi)*w*F / (4*np.pi*SF.G) * (2/q - R*R/q**2)
    return mb, sb, mm, sm


class PointLens:
    """Independent analytic source for the root solver, not an active cluster model."""
    def __init__(self, mass):
        self.mass = mass

    def observables(self, r, geom):
        r = np.asarray(r, float)
        p = dict(radius_kpc=r, alpha_rad=4*SF.G*self.mass/(SF.C_KMS**2*r),
                 sigma_Msun_kpc2=np.zeros_like(r), delta_sigma_Msun_kpc2=self.mass/(np.pi*r*r))
        return observables(p, geom)


def checks():
    gates = {}
    R = np.array([.01, 1., 10., 50., 100., 200.])
    M, a = 1e14, 100.
    sphere = BaryonShells([a], [M])
    mass, sigma = sphere.project(R)
    exact_m = M * (1 - np.maximum(1 - R*R/(a*a), 0.)**1.5)
    exact_s = 3*M/(2*np.pi*a*a) * np.sqrt(np.maximum(1 - R*R/(a*a), 0.))
    zero = WrittenClusterLens(sphere).profile(R)
    err = max(relative(mass, exact_m), relative(sigma, exact_s),
              relative(zero['alpha_rad'], 4*SF.G*mass/(SF.C_KMS**2*R)))
    gates['uniform_sphere_and_zero_response'] = dict(error=err, passed=err < 1e-10,
        doubled_mass_control_rejected=relative(2*zero['alpha_rad'], zero['alpha_rad']) > .99)

    R = np.geomspace(1., 2500., 100)
    mb, sb, mm, sm = gaussian_exact(R)
    errors = []
    for n, npj in [(1200, 2400), (2400, 4800)]:
        p = gaussian_lens(n, npj).profile(R)
        errors.append(max(relative(p['baryon_projected_mass_Msun'], mb),
                          relative(p['field_equivalent_projected_mass_Msun'], mm),
                          relative(p['field_equivalent_surface_density_Msun_kpc2'], sm),
                          relative(p['alpha_rad'], SF.deflection_from_projected_mass(mb+mm, R))))
    # Independently computed erroneous projection without sqrt(2*pi)*w.
    w, q = 80., 100.**2 + 80.**2
    F = 2e14*w*w/q*np.exp(-R*R/(2*q))
    bad = R*R*F/(2*SF.G*q)
    good = np.sqrt(2*np.pi)*w*bad
    gates['gaussian_convolution'] = dict(errors=errors,
        passed=bool(errors[-1] < 2e-4 and errors[-1] < errors[0]),
        missing_projection_factor_control_rejected=relative(bad, good) > .9,
        signed_field_density_preserved=bool(np.min(p['field_equivalent_surface_density_Msun_kpc2']) < 0))

    # Direct 3D Gaussian convolution differentiated analytically, then integrated
    # along the ray by SF's independent force route. The 2D path uses shell input.
    def g(r):
        out = SF.G*gaussian_mass(r, 2e14, 100.)/r**2
        for w, amp in zip([80., 200.], [2e-8, 1e-8]):
            q = 100.**2 + w*w
            out += amp*2e14*w**3/q**1.5 * np.exp(-r*r/(2*q))*r/q
        return out
    radii = np.array([5., 80., 300., 900.])
    direct = np.array([SF.deflection_from_g(g, b) for b in radii])
    proj = gaussian_lens(2400, 4800).profile(radii)['alpha_rad']
    err = relative(proj, direct)
    gates['independent_3d_deflection'] = dict(error=err, passed=err < 2e-4)

    geometry = Geometry(2e5, 1e6, 8e5, 'fictional static Euclidean test')
    lens = PointLens(1e14)
    theta_e = np.sqrt(4*SF.G*lens.mass/SF.C_KMS**2*geometry.ratio/geometry.lens_kpc)*SF.ARCSEC
    analytic_images = np.array([(5 - np.sqrt(25+4*theta_e**2))/2,
                               (5 + np.sqrt(25+4*theta_e**2))/2])
    solved = solve_lens(lens, geometry, np.geomspace(.001, 1e5, 500))
    root_err = relative([x['theta_arcsec'] for x in solved['images']], analytic_images)
    ring_err = abs(solved['critical_curves']['lambda_t'][0]['theta_arcsec']/theta_e - 1)
    residual = max(abs(x['residual_arcsec']) for x in solved['images'])
    aligned = solve_lens(lens, geometry, np.geomspace(.001, 1e5, 500), beta_arcsec=0.)
    gates['images_and_critical_curves'] = dict(image_error=root_err, ring_error=ring_err,
        max_residual_arcsec=residual, passed=bool(root_err < 1e-8 and ring_err < 1e-8
        and residual < 1e-9 and len(aligned['images']) == 0 and aligned['aligned_source_rings_only']))

    invalids = [lambda: Geometry(1, 0, 1, 'bad'), lambda: Geometry(2, 1, 1, 'bad'),
                lambda: BaryonShells([1, 2], [3, 2]), lambda: BaryonShells([1, 1], [1, 2]),
                lambda: WrittenClusterLens(sphere, [0], [1]),
                lambda: WrittenClusterLens(sphere, [1], [-1]),
                lambda: WrittenClusterLens(sphere, [np.nan], [1])]
    rejected = 0
    for trial in invalids:
        try:
            trial()
        except ValueError:
            rejected += 1
    gates['invalid_inputs'] = dict(rejected=rejected, expected=len(invalids), passed=rejected == len(invalids))
    passed = all(g['passed'] for g in gates.values()) and all(
        v for g in gates.values() for k, v in g.items() if k.endswith('_control_rejected') or k == 'signed_field_density_preserved')
    return dict(gates=gates, numerical_verification_passed=bool(passed))


if __name__ == '__main__':
    result = checks()
    print(json.dumps(result, indent=2, allow_nan=False))
    raise SystemExit(0 if result['numerical_verification_passed'] else 1)
