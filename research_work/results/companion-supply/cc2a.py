"""CC-2 stage 1: regime and supply audit of CF-1's companion bath (see protocol.md).

    python cc2a.py [--canonical] [--output-dir DIR]

Uses CF-1's archived potentials, seed profiles, kernel tables and required exposures A* = (sigma/m) rho_inf, and
the same interaction (elastic, equal-mass, isotropic in the centre-of-mass frame, constant sigma/m) applied to
every pair of companions. The reference time (10 Gyr) is an assumption; the cosmic mean density is a comparison
unit only. Fails fast and logs its progress.
"""
import json
import math
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
sys.path.insert(0, str(RESULTS/'gravitational-focusing'))
import cf1  # noqa: E402  (also puts the shared result folders on the path)
import focus as F  # noqa: E402
import evidence_io  # noqa: E402

G, PER_GYR, T_REF = F.G, cf1.PER_GYR, cf1.T_REF
SIGMA_M = (.1, 1., 10., 100., 1000.)           # cm^2/g, trial scales
DISPERSIONS = (100., 200., 300., 500., 1000.)  # km/s, one-dimensional, declared grid
WAVE_FRACTION = .1                              # particle treatment: wavelength <= 0.1 of the half-mass radius
T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:7.1f}s] {msg}', flush=True)


def _unit(rng, n):
    d = rng.normal(size=(n, 3))
    return d/np.linalg.norm(d, axis=1)[:, None]


def seedless_analytic(u, v_esc, n=4000):
    """Two unbound companions, both at speed w = sqrt(u^2 + v_esc^2), independent random directions.
    With phi the angle between them: relative speed g = 2 w sin(phi/2); centre-of-mass speed V = w cos(phi/2);
    V g = w^2 sin(phi). One ends bound with probability P(phi) = max(0, 1 - u^2/(w^2 sin phi)); two never can.
    Returns <g P> and <g> over isotropic phi (weight sin(phi)/2)."""
    w2 = u*u + v_esc*v_esc
    x, wt = np.polynomial.legendre.leggauss(n)
    phi = .5*math.pi*(x + 1)
    dens = .5*np.sin(phi)*.5*math.pi*wt
    g = 2*math.sqrt(w2)*np.sin(phi/2)
    P = np.maximum(0., 1 - u*u/(w2*np.maximum(np.sin(phi), 1e-300)))
    return float(np.sum(dens*g*P)), float(np.sum(dens*g))


def seedless_mc(u, v_esc, n=400000, seed=7):
    """V1 and V3: the same by sampling collisions."""
    rng = np.random.default_rng(seed)
    w = math.sqrt(u*u + v_esc*v_esc)
    v1, v2 = w*_unit(rng, n), w*_unit(rng, n)
    vcm, g = .5*(v1 + v2), np.linalg.norm(v1 - v2, axis=1)
    e = _unit(rng, n)
    v1p, v2p = vcm + .5*g[:, None]*e, vcm - .5*g[:, None]*e
    k1, k2 = np.einsum('ij,ij->i', v1p, v1p), np.einsum('ij,ij->i', v2p, v2p)
    bound = (k1 < v_esc**2) | (k2 < v_esc**2)
    both = (k1 < v_esc**2) & (k2 < v_esc**2)
    dE = float(np.max(np.abs(k1 + k2 - 2*w*w))/(w*w))
    dP = float(np.max(np.abs(v1p + v2p - v1 - v2))/w)
    return float(np.mean(g*bound)), float(np.mean(g)), int(both.sum()), dE, dP


def bound_bound(v_esc, s, n=200000, seed=11):
    """Two seed companions drawn from an isotropic Maxwellian (1-D dispersion s) truncated at v_esc:
    returns <g P(one escapes)> and <g>."""
    rng = np.random.default_rng(seed)

    def draw(m):
        v = rng.normal(scale=s, size=(m, 3))
        bad = np.einsum('ij,ij->i', v, v) >= v_esc*v_esc
        while bad.any():
            v[bad] = rng.normal(scale=s, size=(int(bad.sum()), 3))
            bad = np.einsum('ij,ij->i', v, v) >= v_esc*v_esc
        return v
    v1, v2 = draw(n), draw(n)
    vcm, g = .5*(v1 + v2), np.linalg.norm(v1 - v2, axis=1)
    e = _unit(rng, n)
    v1p, v2p = vcm + .5*g[:, None]*e, vcm - .5*g[:, None]*e
    esc = (np.einsum('ij,ij->i', v1p, v1p) > v_esc**2) | (np.einsum('ij,ij->i', v2p, v2p) > v_esc**2)
    return float(np.mean(g*esc)), float(np.mean(g))


def maxwell_kernel(u_tab, K_tab, disp, n=4000):
    """Net retention for a Maxwellian bath of 1-D dispersion disp: K_M = int f(u) K(u) du over the tabulated range
    (log-interpolated), with the probability mass outside the table reported."""
    lu = np.log(u_tab)
    u = np.geomspace(u_tab[0], u_tab[-1], n)
    f = math.sqrt(2/math.pi)*u*u/disp**3*np.exp(-u*u/(2*disp*disp))
    K = np.interp(np.log(u), lu, K_tab)
    inside = float(np.trapezoid(f, u))
    return float(np.trapezoid(f*K, u)), 1 - inside


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    arch = json.loads((RESULTS/'gravitational-focusing/cf1-results.json').read_text(encoding='utf-8'))
    out = dict(scope='CC-2 stage 1: regime and supply audit of CF-1 at its required supply; no source, no reservoir, '
                     'no fit. Reference time 10 Gyr is an assumption; the cosmic mean density is a comparison unit.',
               validation={}, systems={})
    # V1 and V3: seedless capture, analytic against Monte Carlo, at five depths
    v1 = {}
    worst_e = worst_p = 0.
    for ratio in (.5, 1., 2., 4., 8.):
        a, g_a = seedless_analytic(300., 300.*ratio)
        m, g_m, both, dE, dP = seedless_mc(300., 300.*ratio)
        v1[f'v_esc/u={ratio:g}'] = dict(analytic=a/g_a, monte_carlo=m/g_m, relative=abs((m/g_m)/(a/g_a) - 1), two_bound=both)
        worst_e, worst_p = max(worst_e, dE), max(worst_p, dP)
    out['validation']['V1'] = dict(rows=v1, passed=bool(max(r['relative'] for r in v1.values()) < 1e-2
                                                         and all(r['two_bound'] == 0 for r in v1.values())))
    out['validation']['V3'] = dict(energy=worst_e, momentum=worst_p, passed=bool(worst_e < 1e-12 and worst_p < 1e-12))
    log(f"V1 {out['validation']['V1']['passed']}  V3 {out['validation']['V3']}")
    worst_v2 = 0.
    for sysd in cf1.systems():
        name = sysd['name']
        if sysd['kind'] == 'lens host' and 'J1630' not in name:
            continue                                   # one lens host is enough here; the rest are in CF-1
        a = arch['systems'][name]
        pot = F.Potential(sysd['r'], sysd['M'], sysd['R_b'], name)
        rho_b, s = cf1.seed_profile(pot)
        Mb, rh, Rb = a['mass_within_R_b_Msun'], a['r_half_kpc'], pot.R_b
        tab_u, tab_K = np.array(a['kernel_table']['u_kms']), np.array(a['kernel_table']['K_kms'])
        row = dict(kind=sysd['kind'], R_b_kpc=Rb, r_half_kpc=rh, trials={})
        for key, t in a['trials'].items():
            u, A = float(key), t['exposure_for_retained_equal_to_baryons']
            tr = dict(K_kms=t['K_kms'])
            # A3 (always), A1/A2 only where CF-1 found a required exposure
            tr['travel_in_T_kpc'] = u*T_REF*PER_GYR              # 1 km/s = 1.0227 kpc/Gyr
            tr['travel_in_T_over_R_b'] = tr['travel_in_T_kpc']/Rb
            tr['crossing_R_b_Gyr'] = Rb/u/PER_GYR
            tr['crossing_r_half_Gyr'] = rh/u/PER_GYR
            lam = {f'{m:g} eV': 2*math.pi*cf1.HBAR*cf1.C_MS**2/(m*cf1.EV_J)/(cf1.KPC_M*1e3)/u for m in (1.34e-24, 1e-22)}
            tr['wavelength_kpc'] = lam
            tr['min_mass_eV_for_particle_regime'] = 2*math.pi*cf1.HBAR*cf1.C_MS**2/cf1.EV_J/(cf1.KPC_M*1e3)/u/(WAVE_FRACTION*rh)
            if A:
                tr['A_star_per_kpc'] = A
                tr['scattering_length_kpc'] = 1/A
                tr['mean_free_path_single_speed_kpc'] = .75/A
                tr['mean_free_path_maxwellian_kpc'] = 1/(math.sqrt(2)*A)
                tr['optical_depth_R_b'] = Rb*A
                tr['optical_depth_r_half'] = rh*A
                tr['collision_time_Gyr'] = 1/(A*u*PER_GYR)
                tr['wavelength_over_mean_free_path'] = {k: v*A/.75 for k, v in lam.items()}
                grav = {}
                for sm in SIGMA_M:
                    rho = A/(sm*cf1.CM2_PER_G)
                    grav[f'{sm:g} cm2/g'] = dict(rho_inf_Msun_kpc3=rho, rho_inf_over_mean=rho/cf1.RHO_MEAN,
                                                 bath_mass_in_R_b_over_baryons=4/3*math.pi*Rb**3*rho/Mb,
                                                 jeans_length_over_R_b=u*math.sqrt(math.pi/(G*rho))/Rb)
                tr['bath_gravity'] = grav
                # A6: seedless capture and bound-bound evaporation against seed-driven capture at A*
                radii = np.geomspace(max(pot.r[1], 1e-2*rh), min(20*rh, .9*Rb), 24)
                lr = np.log(radii)
                ve = np.array([math.sqrt(pot.v_esc2(r)) for r in radii])
                zeta2 = 1 + ve**2/u**2
                gP = np.array([seedless_analytic(u, v)[0] for v in ve])
                shell = 4*np.pi*radii**3                       # dV per d ln r
                I_seedless = float(np.trapezoid(.5*zeta2*gP*shell, lr))          # int 1/2 zeta^2 <gP> dV, kpc^3 km/s
                sr = np.interp(lr, pot.lr, s)
                rb = np.interp(lr, pot.lr, rho_b)
                gE = np.array([bound_bound(v, max(x, 1e-6*v), n=60000)[0] for v, x in zip(ve, sr)])
                I_evap = float(np.trapezoid(.5*(cf1.F_SEED*rb)**2*gE*shell, lr))  # Msun^2/kpc^3 km/s
                seed_mass = cf1.F_SEED*Mb
                ch = {}
                for sm in SIGMA_M:
                    smk = sm*cf1.CM2_PER_G
                    rho = A/smk
                    capture = A*t['K_kms']*seed_mass                              # (sigma/m) rho K M_seed, Msun km/s/kpc
                    ch[f'{sm:g} cm2/g'] = dict(seedless_over_seed_capture=smk*rho*rho*I_seedless/capture,
                                               evaporation_over_seed_capture=smk*I_evap/capture)
                tr['channels'] = ch
            row['trials'][key] = tr
        # A5: mixtures and Maxwellian baths
        K300, K3000 = a['trials']['300']['K_kms'], a['trials']['3000']['K_kms']
        row['mixture_fast_fraction_reversing_retention'] = (K300/(K300 - K3000) if K300 > 0 > K3000 else None)
        row['mixture_K_at_1pct_fast'] = .99*K300 + .01*K3000
        row['maxwellian'] = {}
        for d in DISPERSIONS:
            km, outside = maxwell_kernel(tab_u, tab_K, d)
            row['maxwellian'][f'{d:g}'] = dict(K_kms=km, probability_outside_table=outside)
        # V2: a narrowed Maxwellian-like bath around a tabulated speed reproduces the single-speed kernel
        uk = tab_u[len(tab_u)//2]
        uu = np.geomspace(uk*.999, uk*1.001, 201)
        ww = np.exp(-((uu - uk)/(uk*2e-4))**2)
        narrow = float(np.trapezoid(ww*np.interp(np.log(uu), np.log(tab_u), tab_K), uu)/np.trapezoid(ww, uu))
        worst_v2 = max(worst_v2, abs(narrow/np.interp(math.log(uk), np.log(tab_u), tab_K) - 1))
        out['systems'][name] = row
        t300 = row['trials']['300']
        log(f"{name}: tau(R_b)={t300.get('optical_depth_R_b')}  uT/R_b={t300['travel_in_T_over_R_b']:.3g}  "
            f"fast fraction reversing={row['mixture_fast_fraction_reversing_retention']}")
        if time.time() - T0 > 30*60:
            raise RuntimeError('wall-clock budget exhausted')
    out['validation']['V2'] = dict(max_relative=worst_v2, passed=bool(worst_v2 < 1e-3))
    out['all_validation_passed'] = all(v['passed'] for v in out['validation'].values())
    out['runtime_seconds'] = time.time() - T0
    text = json.dumps(out, indent=1, default=float)
    log(f"validation passed={out['all_validation_passed']}")
    status = evidence_io.finish(args, 'companion-supply-cc2a', text, HERE/'cc2a-results.json', ignore={'/runtime_seconds'})
    return status if out['all_validation_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
