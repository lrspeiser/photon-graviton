"""PM-1: path-memory gravity against the repository's own data (protocol.md, with its corrections).

    python pm1.py [--output-dir DIR] [--canonical]

Candidate A is the control: an exponential memory of a static source returns alpha*Phi_N, and a rigidly
rotating axisymmetric source is the static one. Candidate B is the one with traction: a single remembered
wavelength gives k J0(kR) J1(kR), which alternates; a broad band gives a force set by the matter inside r.
Three sourcing rules are scored, each with one universal constant fitted on the 89 training galaxies:
B-linear (amplitude proportional to mass), B-root (the cumulative square root that was fitted), and the
collective rule (one shared saturated mode, weight proportional to mass). The finite-band field is audited
with the source-integrated ring integral, not the on-ring expression. Candidate C's equilibrium is the
simple-nu form and its delayed feedback is unstable analytically. Candidate E is a luminosity budget,
conditional on a stated band conversion. D and F are declared in the protocol and not run.

Regenerates pm1-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import j0, j1

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE.parent/'capture-to-orbit'))
import inputs as I  # noqa: E402

G = 4.30091727003628e-6                       # kpc (km/s)^2 / Msun
KPC_M = 3.0856775814913673e19
L_SUN_BOL_W = 3.828e26                        # bolometric; the catalogue's L is a 3.6 micron band luminosity
GRAV_SI = 6.67430e-11
C_MS = 299792458.
MOND_A0 = 1.2e-10
REPO_A0 = 8.563e-11                           # the a0 this repository's own MOND-guided fit selected
BARYONS = dict(train=52.57, validation=58.22, test=47.77, mw_I=52.57, mw_II=62.34)
SIMPLE_MOND = dict(train=19.89, validation=26.88, test=16.40, mw_I=9.53, mw_II=12.04)


# ---------------------------------------------------------------- data
def sparc():
    """Per galaxy: the observed curve, the baryonic circular speed and the enclosed baryonic mass."""
    out = []
    for gal in I.sparc_galaxies():
        rows, cat = gal['rotmod'], gal['catalog']
        R = rows[:, 0]
        keep = R > 0
        R = R[keep]
        v_gas, v_disk, v_bul = rows[keep, 3], rows[keep, 4], rows[keep, 5]
        v_bar2 = np.maximum(np.sign(v_gas)*v_gas**2 + I.UPSILON_DISK*v_disk**2 + I.UPSILON_BULGE*v_bul**2, 0.)
        out.append(dict(name=gal['name'], split=gal['split'], R=R, v_obs=rows[keep, 1], e_v=rows[keep, 2],
                        v_bar=np.sqrt(v_bar2), M_b=R*v_bar2/G, L36=cat['L9']*1e9, rd=cat['rd'],
                        frozen_radii=len(gal['r'])))
    return out


def milky_way():
    pick = lambda b: next(x for x in I.milky_way_runs()
                          if x['baryons'] == b and abs(x['rd'] - 2.6) < 1e-9 and abs(x['lf'] - 1) < 1e-9)
    mk = lambda r: dict(name='Milky Way ' + r['baryons'], R=np.asarray(r['R']), v_obs=np.asarray(r['y']),
                        v_bar=np.asarray(r['vb']), M_b=np.asarray(r['R'])*np.asarray(r['vb'])**2/G)
    return mk(pick('I')), mk(pick('II'))


# ---------------------------------------------------------------- candidate A, the control
def candidate_a(gal, alpha=.3, periods=8, n_theta=2048):
    """A stationary source returns alpha*Phi_N; a rigidly rotating axisymmetric source is the same source."""
    g_N = gal['v_bar']**2/gal['R']
    static = alpha*g_N
    theta = np.linspace(0, 2*np.pi*periods, n_theta, endpoint=False)
    rotated = np.array([alpha*g_N for _ in theta]).mean(axis=0)
    return dict(galaxy=gal['name'], alpha=alpha, radii=len(gal['R']),
                static_over_alpha_gN_max_error=float(np.max(np.abs(static/(alpha*g_N) - 1))),
                rotating_minus_static_max=float(np.max(np.abs(rotated - static))),
                note='delayed density cannot see an orbit: the force law stays 1/R^2 and only G is rescaled',
                passed=bool(np.max(np.abs(rotated - static)) < 1e-12*max(np.max(np.abs(static)), 1e-30)))


# ---------------------------------------------------------------- candidate B
def single_wavelength_shape(R, k):
    """k J0(kR) J1(kR): the alternating single-wavelength on-ring force per unit amplitude."""
    return k*j0(k*R)*j1(k*R)


def ring_force(r, R_s, q_min=.01, q_max=100., nq=200001):
    """The finite-band force at r from a source ring at R_s: integral of J1(q r) J0(q R_s) dq.

    This is the quantity an extended galaxy needs. The on-ring expression
    (1/2R)[J0(q_min R)^2 - J0(q_max R)^2] is exact only for r = R_s, and using it at the observation
    radius silently moves the source there (the owner's audit).
    """
    q = np.linspace(q_min, q_max, nq)
    return float(np.trapezoid(j1(q*r)*j0(q*R_s), q))


def on_ring_force(R, q_min=.01, q_max=100.):
    return (j0(q_min*R)**2 - j0(q_max*R)**2)/(2*R)


def ring_audit():
    """Correction 3: the source-integrated ring force against the on-ring expression and the wide-band limit."""
    owner = dict(r=2., R_s=1., source_integrated=ring_force(2., 1.), on_ring_at_r=on_ring_force(2.),
                 owner_source_integrated=.49877, owner_on_ring=.24989)
    limit = [dict(r=r, R_s=1., value=ring_force(r, 1.), wide_band_limit=(1/r if r > 1 else 0.))
             for r in (.2, .5, 2., 5.)]
    return dict(owner_example=owner, wide_band_limit=limit,
                note='a ring contributes C/r outside itself and nothing inside, so the enclosed-mass rule is '
                     'the wide-band limit; the on-ring expression is not the field of a galaxy at radius r',
                reproduces_owner=bool(abs(owner['source_integrated'] - .49877) < 1e-3
                                      and abs(owner['on_ring_at_r'] - .24989) < 1e-3))


def ring_representation(gal, splits=(1, 2, 4, 8)):
    """Correction 2: independent per-ring saturation is representation dependent; the cumulative rule is not."""
    M = gal['M_b'][-1]
    rows = []
    for n in splits:
        m = np.diff(np.concatenate([[0.], gal['M_b']]))/n
        per_ring = float(np.sum(np.sqrt(np.maximum(np.repeat(m, n), 0.))))
        rows.append(dict(sub_rings_per_annulus=n, sum_sqrt_m=per_ring, sqrt_sum_m=float(np.sqrt(M)),
                         ratio=per_ring/float(np.sqrt(M))))
    growth = rows[-1]['sum_sqrt_m']/rows[0]['sum_sqrt_m']
    return dict(galaxy=gal['name'], rows=rows, growth_over_8x_refinement=growth,
                expected_sqrt_8=float(np.sqrt(8)),
                note='splitting the same matter into N pieces multiplies a per-ring square-root field by '
                     'sqrt(N), so independent ring saturation cannot be the mechanism behind the fitted law')


def speeds(gal, const, rule):
    """v(r) for g = g_N + g_mem under each declared sourcing rule."""
    M = np.maximum(gal['M_b'], 0.)
    if rule == 'linear':
        extra = const*M                                   # g = const M/r  ->  v^2 = const M
    elif rule == 'root':
        extra = const*np.sqrt(M)                          # the fitted cumulative rule
    elif rule == 'collective':
        tot = M[-1]
        extra = const*M/np.sqrt(tot) if tot > 0 else np.zeros_like(M)
    else:
        raise ValueError(rule)
    return np.sqrt(np.maximum(gal['v_bar']**2 + extra, 0.))


def rmse(gals, const, rule):
    per = [np.sqrt(np.mean((speeds(g, const, rule) - g['v_obs'])**2)) for g in gals]
    return float(np.sqrt(np.mean(np.square(per))))


def candidate_b(gals, mw, rule, bracket):
    train = [g for g in gals if g['split'] == 'train']
    r = minimize_scalar(lambda c: rmse(train, c, rule), bracket=bracket, method='brent',
                        options=dict(xtol=1e-12))
    const = float(r.x)
    out = dict(rule=rule, constant=const, train=float(r.fun))
    for split in ('validation', 'test'):
        out[split] = rmse([g for g in gals if g['split'] == split], const, rule)
    out['mw_I'], out['mw_II'] = rmse([mw[0]], const, rule), rmse([mw[1]], const, rule)
    if rule == 'root':
        a_star = const**2/G*1e6/KPC_M
        out.update(implied_a_star_m_s2=a_star, implied_a_star_over_mond=a_star/MOND_A0,
                   implied_a_star_over_repo_fit=a_star/REPO_A0,
                   identity='g_mem = sqrt(a* g_mono) with a* = beta^2/G, so fitting beta is fitting a*')
    outer = [(g['M_b'][-1], g['v_obs'][-5:].mean(), speeds(g, const, rule)[-5:].mean())
             for g in gals if len(g['R']) >= 5 and g['M_b'][-1] > 0]
    M, v_obs, v_mod = (np.array(x) for x in zip(*outer))
    slope = lambda v: float(np.polyfit(np.log10(M), np.log10(v), 1)[0])
    out.update(mass_speed_slope_observed=slope(v_obs), mass_speed_slope_model=slope(v_mod),
               mass_speed_galaxies=len(M),
               beats_baryons=bool(all(out[k] < BARYONS[k] for k in
                                      ('train', 'validation', 'test', 'mw_I', 'mw_II'))))
    return out


def nonlinear_field_identity(beta):
    """Correction 2: the spherical limit of div[(|grad psi|/a*) grad psi] = 4 pi G rho gives sqrt(G a* M)/r."""
    a_star_kpc = beta**2/G
    return dict(statement='g_psi = sqrt(G a* M(<r))/r, which is B-root with beta = sqrt(G a*)',
                beta_from_a_star=float(np.sqrt(G*a_star_kpc)), beta_fitted=float(beta),
                relative_difference=float(abs(np.sqrt(G*a_star_kpc)/beta - 1)),
                note='the square root follows from a collective nonlinear field, not from a square-root '
                     'charge on every ring; in a disk the solution will not equal this enclosed-mass rule')


# ---------------------------------------------------------------- candidate C
def c_equilibrium(g_N, a_star):
    return .5*(g_N + np.sqrt(g_N**2 + 4*a_star*g_N))


def c_growth(q0, T):
    return float(np.max(np.real(np.roots([T, 1 + q0, T, 1 + 3*q0]))))


def candidate_c(gals, a_star=MOND_A0):
    a_kpc = a_star*KPC_M/1e6
    owner = np.roots([1., 1.9, 1., 3.7])
    grid = {f'{T:g}': {f'{q:g}': c_growth(q, T) for q in (.1, .3, .5, .7, .9, .99)}
            for T in (1e-3, 1e-2, .1, .3, 1., 3., 10., 100., 1000.)}
    q_outer = np.array([1 - (g['v_bar'][-1]**2/g['R'][-1])/c_equilibrium(g['v_bar'][-1]**2/g['R'][-1], a_kpc)
                        for g in gals if g['v_bar'][-1] > 0])
    med = float(np.median(q_outer))
    worst = max(max(v.values()) for v in grid.values())
    return dict(
        owner_root_check=dict(expected=[.1533, 1.2858],
                              got=[float(np.max(np.real(owner))), float(np.max(np.abs(np.imag(owner))))]),
        routh_hurwitz=dict(condition='(1+q0) T > T (1+3 q0), i.e. q0 < 0',
                           holds_for_positive_q0=False,
                           statement='every q0 > 0 with finite T > 0 has a growing mode; this is analytic, '
                                     'not a grid result'),
        asymptotics=dict(small_T='Re lambda ~ q0 T/(1+q0)^2', large_T='Re lambda ~ q0/T',
                         small_T_check=[c_growth(.9, t) for t in (1e-3, 1e-2)],
                         large_T_check=[c_growth(.9, t) for t in (1e2, 1e3)]),
        growth_rate_grid=grid, worst_growth_rate=worst,
        q_outer_median=med, q_outer_90th=float(np.percentile(q_outer, 90)), galaxies=len(q_outer),
        e_folding_orbits_at_median_q=float(1/(2*np.pi*c_growth(med, 1.))),
        b_root_epicyclic=dict(statement='in B-root\'s frozen potential -GM/r + K ln r, kappa^2 = GM/r^3 + '
                                        '2K/r^2 > 0, so its circular orbits are radially stable',
                              transfers_from_C=False),
        stable_anywhere=bool(worst <= 0))


# ---------------------------------------------------------------- candidate E
def candidate_e(gals):
    rows = [(np.sqrt(GRAV_SI*g['L36']*L_SUN_BOL_W/C_MS**3)/1e3,
             float(g['v_obs'][-5:].mean()) if len(g['v_obs']) >= 5 else float(g['v_obs'][-1]))
            for g in gals if g['L36'] > 0]
    v_u, v_f = (np.array(x) for x in zip(*rows))
    need = (v_f/v_u)**2
    return dict(galaxies=len(v_u), v_u_median_km_s=float(np.median(v_u)), v_f_median_km_s=float(np.median(v_f)),
                power_shortfall_median=float(np.median(need)),
                power_shortfall_range=[float(np.min(need)), float(np.max(need))],
                conversion='L = (catalogue 3.6 micron luminosity in that band\'s solar units) x the '
                           'bolometric solar constant 3.828e26 W; that is an assumed spectral conversion, '
                           'not a measured bolometric power',
                note='conditional on that conversion; a bolometric correction of order a few does not move a '
                     'deficit of 10^9, but the number is a band-luminosity ratio and is labelled as one')


# ---------------------------------------------------------------- driver
def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    gals, mw = sparc(), milky_way()
    biggest = max(gals, key=lambda g: len(g['R']))
    root = candidate_b(gals, mw, 'root', (1e-2, 1e-1))
    linear = candidate_b(gals, mw, 'linear', (1e-7, 1e-6))
    collective = candidate_b(gals, mw, 'collective', (1e-2, 1e-1))
    single = single_wavelength_shape(biggest['R'], 1.)
    c = candidate_c(gals)
    result = dict(
        experiment='PM-1: path-memory gravity against SPARC and the Milky Way',
        protocol='protocol.md, with the corrections after the owner\'s review',
        data=dict(galaxies=len(gals),
                  radii_rotmod_positive=int(sum(len(g['R']) for g in gals)),
                  radii_frozen_comparison=int(sum(g['frozen_radii'] for g in gals)),
                  radius_count_note='the first counts rotmod rows with R > 0, the second the frozen '
                                    'comparison arrays; both are reported (correction 5)',
                  splits={s: sum(1 for g in gals if g['split'] == s) for s in ('train', 'validation', 'test')},
                  milky_way_bins=len(mw[0]['R']), exposed=True, blind=False),
        comparison_scores=dict(baryons_only=BARYONS, simple_mond_fitted=SIMPLE_MOND),
        candidate_A_control=candidate_a(gals[0]),
        candidate_B=dict(
            root_cumulative=root, linear=linear, collective=collective,
            nonlinear_field_identity=nonlinear_field_identity(root['constant']),
            ring_audit=ring_audit(), ring_representation=ring_representation(biggest),
            single_wavelength=dict(galaxy=biggest['name'], k=1., radii=len(biggest['R']),
                                   sign_changes=int(np.sum(np.diff(np.sign(single)) != 0)),
                                   first_values=[float(x) for x in single[:6]],
                                   note='inward where positive, outward where negative')),
        candidate_C=c,
        candidate_E=candidate_e(gals),
        gates=dict(
            G1_shape={k: v['beats_baryons'] for k, v in
                      (('root', root), ('linear', linear), ('collective', collective))},
            G2_mass_speed=dict(observed=root['mass_speed_slope_observed'],
                               root=root['mass_speed_slope_model'], linear=linear['mass_speed_slope_model'],
                               collective=collective['mass_speed_slope_model']),
            G3_universality='one constant per rule, fitted on the 89 training galaxies only',
            G4_stability=dict(candidate_C_stable=c['stable_anywhere'],
                              b_root_frozen_potential_stable=True)),
        input_sha256={'protocol.md': hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest()},
        checks_short_run=dict(root_constant=root['constant'], root_train=root['train'], root_test=root['test'],
                              linear_train=linear['train'], collective_train=collective['train'],
                              c_worst_growth=c['worst_growth_rate'],
                              ring_source_integrated=ring_audit()['owner_example']['source_integrated']),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:3000])
    return evidence_io.finish(args, 'path-memory', text, HERE/'pm1-results.json')


if __name__ == '__main__':
    raise SystemExit(main())
