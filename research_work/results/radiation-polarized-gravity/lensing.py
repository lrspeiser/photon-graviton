"""Spherical lensing for RPG-1 (T3).

Einstein radii of the six SLACS lenses from population stellar masses. Mass follows the archived
Sersic light profiles, deprojected as spheres, in PF-1's static Euclidean geometry. The declared
response bends light with the full field (Phi = Psi). In spherical symmetry AQUAL is exact in the
algebraic form g = nu(g_N/a) g_N. The Newtonian-baryon bracket uses g_N.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import brentq
from scipy.special import gamma

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent if (HERE.parent/'capture-to-orbit').exists() else \
    Path('C:/Users/henry/Documents/Codex/photon-graviton/research_work/results')
sys.path.insert(0, str(HERE))
import aqual as Q  # noqa: E402

ARCSEC = 206264.806
ALPHA0 = 0.0002488993286382367                       # per Mpc
LENSES = ('J0037-0942', 'J1112+0826', 'J1204+0358', 'J1402+6321', 'J1621+3931', 'J1630+4520')
BRANCH = 'energy_loss_and_event_stretch'             # PF-1 flux law: D_L = D(1+z)
INPUTS = dict(observations=RESULTS/'lensing-data-readiness/lens-observations-and-image-models.json',
              geometry=RESULTS/'lensing-data-readiness/conditional-geometry.json',
              light=RESULTS/'slacs-light-profile-audit/results.json',
              masses=RESULTS/'lens-photometric-audit/normalization-sensitivity.json')


def light_profile(components, Dl):
    """Enclosed-light fraction of spherically deprojected Sersic components on a radial grid (kpc).
    nu(r) = -(1/pi) int_0^inf I'(r cosh u) du; totals from the analytic Sersic integral."""
    comps = [dict(R=c['R_arcsec']*Dl/ARCSEC, n=c['n'], b=c['bn'], A=c['amp_at_R']) for c in components]
    r = np.geomspace(1e-4*min(c['R'] for c in comps), 1e3*max(c['R'] for c in comps), 6000)
    u, w = np.polynomial.legendre.leggauss(256)
    density, totals = np.zeros_like(r), []
    for c in comps:
        top = np.arccosh(np.maximum(c['R']*(1 + 100/c['b'])**c['n']/r, 1.))
        uu, ww = .5*(u + 1)*top[:, None], .5*w*top[:, None]
        R = r[:, None]*np.cosh(uu)
        x = (R/c['R'])**(1/c['n'])
        dI = -c['A']*(c['b']/c['n'])*x/R*np.exp(-c['b']*(x - 1))
        density -= (dI*ww).sum(1)/np.pi
        totals.append(2*np.pi*c['A']*c['R']**2*c['n']*np.exp(c['b'])*gamma(2*c['n'])/c['b']**(2*c['n']))
    enclosed = cumulative_trapezoid(4*np.pi*r*r*density, r, initial=0.)/sum(totals)
    return r, enclosed, [t/sum(totals) for t in totals]


def einstein_radius(mass, r, frac, Dl, ratio, a=None):
    """theta_E (arcsec), b_E (kpc) and g_N(b_E) solving b/Dl = (Dls/Ds) alpha(b), with
    alpha = (4/c^2) int_0^{pi/2} g(b sec t) b sec t dt. A subcritical lens has no ring: zeros."""
    lr = np.log(r)

    def g(x):
        gN = Q.G*mass*np.interp(np.log(x), lr, frac, left=0., right=frac[-1])/(x*x)
        return gN if a is None else Q.nu(gN/a)*gN

    def bend(b):
        return 4/Q.C_KMS**2*quad(lambda t: g(b/np.cos(t))*b/np.cos(t), 0, np.pi/2, limit=200, epsabs=0, epsrel=1e-8)[0]

    def f(b):
        return ratio*bend(b) - b/Dl

    if f(1e-3) <= 0:
        return 0., 0., 0.
    b = brentq(f, 1e-3, 300., xtol=1e-12)
    return b/Dl*ARCSEC, b, float(Q.G*mass*np.interp(np.log(b), lr, frac)/b**2)


def mass_factor(mass, r, frac, Dl, ratio, a, target):
    """Multiplier of the stellar mass that gives the observed Einstein radius (theta rises monotonically with
    mass and is zero below criticality)."""
    def theta(f):
        return einstein_radius(f*mass, r, frac, Dl, ratio, a)[0]
    lo, hi = 1., 1.
    while theta(hi) < target and hi < 1e3:
        hi *= 2
    while theta(lo) > target and lo > 1e-3:
        lo /= 2
    return brentq(lambda f: theta(f) - target, lo, hi, xtol=1e-10)


def load():
    obs = {r['Name']: r for r in json.loads(INPUTS['observations'].read_text())}
    geo = {r['Name']: r for r in json.loads(INPUTS['geometry'].read_text())}
    light = {r['Name']: r for r in json.loads(INPUTS['light'].read_text())['rows']}
    masses = {}
    for r in json.loads(INPUTS['masses'].read_text()):
        if r['propagation_branch'] == BRANCH:
            masses.setdefault((r['Name'], r['imf']), r['conditional_log10_stellar_mass'])
    return obs, geo, light, masses


def run(a_star):
    """All six lenses: Newtonian bracket, RPG-1 with constant a*, RPG-1 with a*(z_l) = (1+z_l) a*."""
    obs, geo, light, masses = load()
    rows = []
    for name in LENSES:
        o, g = obs[name], geo[name]
        Dl, ratio = g['conditional_Dl_Mpc']*1000, g['conditional_Dls_over_Ds']
        r, frac, parts = light_profile(light[name]['components'], Dl)
        row = dict(name=name, z_lens=o['zFG'], z_source=o['zBG'], theta_observed_arcsec=o['bSIE'], Dl_kpc=Dl,
                   Dls_over_Ds=ratio, Dl_check_ln1pz_over_alpha=float(np.log1p(o['zFG'])/ALPHA0*1000/Dl),
                   numerical_light_enclosed=float(frac[-1]),
                   light_fraction_check=[float(p/c['total_light_fraction'])
                                         for p, c in zip(parts, light[name]['components'])])
        for imf in ('Chabrier', 'Salpeter'):
            M = 10**masses[(name, imf)]
            cases = {'newtonian_baryons': None, 'rpg1_constant_a': a_star, 'rpg1_linked_a': a_star*(1 + o['zFG'])}
            res = {}
            for label, a in cases.items():
                theta, b, gN = einstein_radius(M, r, frac, Dl, ratio, a)
                res[label] = dict(theta_arcsec=theta, ratio_to_observed=theta/o['bSIE'], b_E_kpc=b,
                                  gN_at_bE_over_a_star=gN/a_star, subcritical=bool(theta == 0))
            row[imf] = dict(log10_stellar_mass=float(np.log10(M)), cases=res,
                            rpg1_mass_factor_for_observed_theta=float(mass_factor(M, r, frac, Dl, ratio, a_star,
                                                                                  o['bSIE'])),
                            newtonian_mass_factor_for_observed_theta=float(mass_factor(M, r, frac, Dl, ratio, None,
                                                                                       o['bSIE'])))
        rows.append(row)
    return rows
