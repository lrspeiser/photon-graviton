"""BRIDGE-1B: one interaction for the photon's loss and the companions' mass (protocol.md).

    python bridge.py [--output-dir DIR] [--canonical]

Part A runs the declared code-unit trials of the semiclassical model in modes.py: the crossing
spectrum, energy conservation, the back-reaction and trapping, momentum and cooling, the owner's four
controls, the probe signals, and the kinetic hybrid. Part B maps the model onto physical scales, with
PF-1's archived rate and the microwave background fixed and 2B-F1's reference source as the comparison
abundance, and reads off the energy source, the history, the coldness and the interface quantities.
Calculation 2 evaluates the measured shift against atoms for the derived histories, including the
exponent identity of the minimal coupling family and the Pantheon+ brightness of every history.

Regenerates bridge-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.linalg import cho_factor, cho_solve

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))
import modes as MD  # noqa: E402

DATA = ROOT/'shared_interaction_test/data'
C_MS, C_KMS = 299792458., 299792.458
KPC_M = 3.0856775814913673e19
MPC_M = 1e3*KPC_M
MSUN_KG, GYR_S, YR_S = 1.98847e30, 3.15576e16, 3.15576e7
A_RAD, T_CMB = 7.565723e-16, 2.72548
G_SI = 6.67430e-11
HBAR_EVS, HBARC_EVM = 6.582119569e-16, 1.9732698e-7
EV4_PER_J_M3 = HBARC_EVM**3/1.602176634e-19
ALPHA_PER_MPC = 2.488993286382367e-4          # PF-1's archived rate
H_RATE = ALPHA_PER_MPC/MPC_M*C_MS             # ndot/n today, per second
U_GAMMA0 = A_RAD*T_CMB**4                     # J/m^3
Q_STAR, T_SPAN_GYR = 1437.8, 10.              # 2B-F1's reference source, Msun/kpc^3/Gyr over 10 Gyr
RHO_C0 = Q_STAR*T_SPAN_GYR*MSUN_KG/KPC_M**3*C_MS**2     # companions' present rest energy density, J/m^3
RHO_CRIT = 3*H_RATE**2*C_MS**2/(8*np.pi*G_SI)           # with H0 = c alpha, J/m^3
CLOCK_BOUND_PER_YR = 1e-16                    # about three times the Yb+/Cs uncertainty
N_STARS = (.9, .7, .5, .3, .237, .1)
R_ES = (1, 3, 10, 30, 100, 300)
A_HATS = (1., .5, .2, .0655, 0., -.05, -.1, -.15, -.2, -.3, -.4)
REF_NSTAR, REF_RE = .237, 10


# ---------------------------------------------------------------- part A, code-unit trials
def trial_crossing_spectrum():
    """V1: a prescribed linear crossing reproduces exp[-pi(k^2+m0^2)/(g|ndot|)]."""
    g, nd, m0, T0 = 1., 1., .3, 60.
    f = MD.Field(K=1., g=g, m0=m0, n_star=0., A=0., kmax=3.5, J=120,
                 prescribed=lambda t: (-T0*nd + nd*t, nd))
    r = f.run(n0=-T0*nd, v0=nd, T=2*T0, n_out=3)
    nk, lz = r['states'][-1]['nk'], MD.landau_zener(f.k, g, nd, m0)
    sel = lz >= 1e-2
    abs_err, rel_err = float(np.max(np.abs(nk - lz))), float(np.max(np.abs(nk[sel]/lz[sel] - 1)))
    number = float(np.sum(f.w*nk))
    return dict(name='V1 crossing spectrum', g=g, ndot=nd, m0=m0, modes=f.J,
                max_absolute_error=abs_err, max_relative_error_above_0p01=rel_err,
                number_density=number, number_density_analytic=MD.spectrum_number(g, nd, m0),
                tolerance=dict(absolute=1e-3, relative=.02),
                passed=bool(abs_err <= 1e-3 and rel_err <= .02))


def trial_trapping(lam, m0=0., T=4000., n0=-8., v0=1., J=120, kmax=4., n_out=4001):
    """V3/V4: the back-reaction stops the field; occupations and speeds follow the frozen spectrum."""
    K, g = 1./lam, 1.
    f = MD.Field(K=K, g=g, m0=m0, n_star=0., A=0., kmax=kmax, J=J)
    r = f.run(n0=n0, v0=v0, T=T, n_out=n_out, events=[MD.stop_at_rest()])
    stop = r['events'][0][0] if r['events'] and r['events'][0] else r['states'][-1]
    analytic = MD.trapping_distance(K, g, v0, m0)
    after = next(s for s in r['states'] if s['n'] > 0 and s['max_nonadiabaticity'] < 1e-3)
    krms = float(np.sqrt(np.sum(f.w*after['nk']*f.k**2)/np.sum(f.w*after['nk'])))
    late = [s for s in r['states'] if after['t'] < s['t'] <= stop['t']]
    nk_drift = max(float(np.max(np.abs(s['nk'] - after['nk']))) for s in late) if late else 0.
    v_err = max(abs(s['v_rms']/(krms/s['mass']) - 1) for s in late) if late else 0.
    return dict(name=f'V3/V4 trapping and cooling, lambda = {lam:g}', lam=lam, m0=m0, K=K, g=g,
                stop_n=stop['n'], stop_analytic=analytic, stop_relative_error=abs(stop['n']/analytic - 1),
                companion_energy_at_stop=stop['rho_C'], initial_kinetic_energy=K*v0**2/2,
                energy_closure=abs(stop['rho_C'] - K*v0**2/2)/(K*v0**2/2),
                number_density=stop['number'], number_density_analytic=MD.spectrum_number(1., v0, m0),
                k_rms=krms, k_rms_analytic=float(np.sqrt(3/(2*np.pi))*np.sqrt(g*v0)),
                occupation_drift=nk_drift, speed_law_max_error=v_err, energy_error=r['energy_error'],
                passed=bool(abs(stop['n']/analytic - 1) <= .05 and nk_drift <= 1e-3 and v_err <= .01
                            and r['energy_error'] <= 1e-8), run=r, field=f)


def trial_reference_and_controls():
    """The photon-driven reference trial and V5's four controls."""
    base = dict(K=1., g=4., m0=.2, n_star=3., kmax=12., J=160, track_tau=True)
    out, runs = [], {}
    for label, kw in (('reference: radiation only, field at rest', dict(A=20., v0=0.)),
                      ('control: no radiation, field at rest', dict(A=0., v0=0.)),
                      ('control: no coupling', dict(A=20., v0=0., g=0.)),
                      ('control: frozen index', dict(A=20., v0=0., frozen=True)),
                      ('control: moving field, no radiation', dict(A=0., v0=2.))):
        p, kw = dict(base), dict(kw)
        v0, frozen = kw.pop('v0'), kw.pop('frozen', False)
        p.update(kw)
        if frozen:
            p['prescribed'] = lambda t: (1., 0.)
        f = MD.Field(**p)
        r = f.run(n0=1., v0=v0, T=40., n_out=201)
        a, b = r['states'][0], r['states'][-1]
        runs[label] = (f, r)
        out.append(dict(trial=label, n_initial=a['n'], n_final=b['n'],
                        radiation_lost=a['radiation'] - b['radiation'], companion_energy=b['rho_C'],
                        field_kinetic_energy=b['kinetic'], number_density=b['number'], v_rms=b['v_rms'],
                        companion_mass=b['mass'], energy_error=r['energy_error']))
    ref = out[0]
    checks = dict(
        no_radiation_no_motion=bool(abs(out[1]['number_density']) <= 1e-8*ref['number_density']),
        no_coupling_no_companions=bool(abs(out[2]['companion_energy']) <= 1e-10*ref['companion_energy']
                                       and abs(out[2]['field_kinetic_energy'] - out[2]['radiation_lost'])
                                       <= 1e-8*max(out[2]['radiation_lost'], 1e-300)),
        frozen_index_nothing=bool(abs(out[3]['number_density']) <= 1e-8*ref['number_density']
                                  and out[3]['n_final'] == out[3]['n_initial']),
        moving_field_produces=bool(out[4]['companion_energy'] > 0 and out[4]['radiation_lost'] == 0.),
        energy_tolerance=1e-8,
        energy_ok=bool(max(o['energy_error'] for o in out) <= 1e-8))
    return out, checks, runs


def trial_probe_signals(f, r, t_e=5., distance=10., deltas=(1e-3, 5e-4)):
    """V6: the arrival-interval stretch converges to the carrier redshift n_o/n_e (fixed-clock convention)."""
    sol = r['sol']
    arrive = lambda te: MD.ray_arrival(sol, f.J, te, distance)
    t_o = arrive(t_e)
    stretch = [(arrive(t_e + d) - t_o)/d for d in deltas]
    ratio = float(sol.sol(t_o)[0])/float(sol.sol(t_e)[0])
    extrapolated = 2*stretch[1] - stretch[0]                 # first-order Richardson in the emission interval
    order = abs(stretch[0] - ratio)/abs(stretch[1] - ratio)
    return dict(name='V6 probe signals', emission_time=t_e, coordinate_distance=distance, arrival_time=t_o,
                stretch_at_emission_intervals=dict(zip((f'{d:g}' for d in deltas), stretch)),
                carrier_redshift=ratio, richardson_stretch=extrapolated,
                relative_difference=abs(extrapolated/ratio - 1), convergence_ratio=order,
                tolerance=1e-7, passed=bool(abs(extrapolated/ratio - 1) <= 1e-7 and 1.8 <= order <= 2.2))


def trial_hybrid(lam=1.):
    """V7: freezing the occupations once every populated mode is adiabatic does not move the stopping point."""
    full = trial_trapping(lam)
    f, r = full['field'], full['run']
    out = []
    for factor in (1., 3.):
        switch = next(s for s in r['states'] if s['n'] > 0 and s['max_nonadiabaticity'] < 1e-3/factor)
        y = r['sol'].sol(switch['t'])
        nk = np.abs(y[2+3*f.J:2+4*f.J] + 1j*y[2+4*f.J:2+5*f.J])**2
        k = f.kinetic(switch['n'], switch['ndot'], 4000., nk, n_out=2001, events=[MD.stop_at_rest()], t0=switch['t'])
        stop = k['events'][0][0] if k['events'] and k['events'][0] else k['states'][-1]
        out.append(dict(switch_time=switch['t'], switch_nonadiabaticity=switch['max_nonadiabaticity'],
                        stop_n=stop['n'], relative_difference_to_full=abs(stop['n']/full['stop_n'] - 1)))
    spread = abs(out[0]['stop_n']/out[1]['stop_n'] - 1)
    return dict(name='V7 kinetic hybrid', lam=lam, full_stop_n=full['stop_n'], switches=out,
                switch_time_spread=spread, tolerance=.01,
                passed=bool(spread <= .01 and max(o['relative_difference_to_full'] for o in out) <= .01))


# ---------------------------------------------------------------- part B, physical mapping
def history(n_star, R_E, extra_force=0., lookback_max_gyr=40.):
    """Backward integration from today: K n'' = u_gamma/n^2 - (F_C - extra) for n > n*, else radiation only."""
    K = 2*R_E*RHO_C0/H_RATE**2
    F_C = RHO_C0/(1 - n_star)
    def f(s, y):                                   # s is lookback time
        n, nd = y
        force = U_GAMMA0/n**2 - ((F_C - extra_force) if n > n_star else 0.)
        return [-nd, -force/K]
    ev_star = lambda s, y: y[0] - n_star
    ev_star.direction = -1
    ev_turn = lambda s, y: y[1]
    ev_turn.terminal, ev_turn.direction = True, -1
    sol = solve_ivp(f, (0, lookback_max_gyr*GYR_S), [1., H_RATE], method='DOP853', rtol=1e-12, atol=1e-16,
                    events=[ev_star, ev_turn], dense_output=True, max_step=.05*GYR_S)
    t_star = float(sol.t_events[0][0]) if len(sol.t_events[0]) else float('nan')
    ndot_star = float(sol.sol(t_star)[1]) if np.isfinite(t_star) else float('nan')
    t_turn = float(sol.t_events[1][0]) if len(sol.t_events[1]) else float('nan')
    s = np.linspace(0, sol.t[-1], 400001)
    n = sol.sol(s)[0]
    ok = np.isfinite(n) & (n > 0)
    s, n = s[ok], n[ok]
    d_optical = np.concatenate([[0.], np.cumsum(np.diff(s)*.5*(1/n[1:] + 1/n[:-1]))])*C_MS
    return dict(sol=sol, K=K, F_C=F_C, t_star=t_star, ndot_star=ndot_star, t_turn=t_turn,
                n_turn=float(sol.sol(t_turn)[0]) if np.isfinite(t_turn) else float('nan'),
                z_grid=1/n - 1, d_grid=d_optical, lookback=s, n_grid=n)


def optical_distance(h, z):
    """Optical distance at coordinate redshift z, from the history's own grid (metres)."""
    return np.interp(z, h['z_grid'], h['d_grid'], left=np.nan, right=np.nan)


def microphysics(n_star, h, mu0=0.):
    """k*, g, N, and the companion mass today, from the benchmark abundance and the crossing speed."""
    rho_ev4 = RHO_C0*EV4_PER_J_M3
    ndot_ev = h['ndot_star']*HBAR_EVS
    k_star = (8*np.pi**3*ndot_ev*rho_ev4/(1 - n_star)*np.exp(np.pi*mu0**2))**.2
    g = k_star**2/ndot_ev
    return dict(k_star_eV=k_star, g_eV=g, number_density_eV3=rho_ev4/(g*(1 - n_star)),
                mass_today_eV=g*(1 - n_star), k_rms_eV=np.sqrt(3/(2*np.pi))*k_star,
                de_broglie_m=2*np.pi*HBARC_EVM/(np.sqrt(3/(2*np.pi))*k_star), m0_over_k_star=mu0)


def scenario(n_star, R_E, pantheon=None, extra_force=0., mu0=0.):
    h = history(n_star, R_E, extra_force=extra_force)
    m = microphysics(n_star, h, mu0)
    k_rms, k_star, g = m['k_rms_eV'], m['k_star_eV'], m['g_eV']
    # The mass is g*(n - n*) on the history itself; near the crossing, where the cooling happens, that is
    # g*ndot* dt to a part in a thousand, which is what the travel and cooling-time formulas use.
    mass_at = lambda lookback_s: g*(float(h['sol'].sol(max(lookback_s, 0.))[0]) - n_star)
    v_of = lambda mm: float(k_rms/np.hypot(k_rms, mm))
    travel = lambda dt_s: float(k_rms/k_star**2*np.arcsinh(k_star**2*(dt_s/HBAR_EVS)/k_rms)*HBARC_EVM)
    # time to fall below a speed: k_rms/m = beta  ->  m = k_rms/beta (non-relativistic)
    t_below = lambda beta: float(k_rms/beta/(g*h['ndot_star']*HBAR_EVS)*HBAR_EVS)
    zs = np.array([.5, 1., 1.5])
    d = optical_distance(h, zs)
    dmu = 5*np.log10((H_RATE*d/C_MS)/np.log1p(zs))
    out = dict(
        n_star=n_star, R_E=R_E, crossing_lookback_Gyr=h['t_star']/GYR_S, ndot_star_over_h=h['ndot_star']/H_RATE,
        field_kinetic_J_m3=R_E*RHO_C0, field_kinetic_over_u_cmb=R_E*RHO_C0/U_GAMMA0,
        field_kinetic_over_rho_crit=R_E*RHO_C0/RHO_CRIT,
        trapping_distance_in_n=R_E*(1 - n_star)*(h['ndot_star']/H_RATE)**2, needed_span_in_n=1 - n_star,
        photon_fraction_since_crossing=U_GAMMA0*(1/n_star - 1)/RHO_C0,
        turnaround_index=U_GAMMA0/((R_E + 1)*RHO_C0), turnaround_temperature_K=T_CMB*(R_E + 1)*RHO_C0/U_GAMMA0,
        turnaround_lookback_Gyr=h['t_turn']/GYR_S, turnaround_index_integrated=h['n_turn'],
        **m,
        v_rms_at_1_Gyr=v_of(mass_at(h['t_star'] - GYR_S)), v_rms_today=v_of(m['mass_today_eV']),
        v_rms_today_km_s=v_of(m['mass_today_eV'])*C_KMS,
        travel_to_10_km_s_m=travel(t_below(10/C_KMS)), travel_to_3_km_s_m=travel(t_below(3/C_KMS)),
        time_to_10_km_s_s=t_below(10/C_KMS), travel_total_m=travel(h['t_star']),
        mass_growth_rate_today_per_Gyr=H_RATE*GYR_S/(1 - n_star),
        mass_growth_rate_at_1_Gyr_per_Gyr=1.,                      # m grows linearly, so mdot/m = 1/(t - t*)
        source_power_W_m3=RHO_C0*H_RATE/(1 - n_star),
        abundance_over_today={f'{z:g}': float(max(0., (1/(1 + z) - n_star)/(1 - n_star))) for z in (.5, 1., 2.)},
        abundance_vanishes_at_redshift=1/n_star - 1,
        delta_mu_mag=dict(zip(('0.5', '1.0', '1.5'), (float(x) for x in dmu))),
        fifth_force_over_gravity=C_MS**4/(4*np.pi*G_SI*h['K']*C_MS**2*(1 - n_star)**2),
        local_delta_n_at_200_km_s=(200e3/C_MS)**2*C_MS**4/(4*np.pi*G_SI*h['K']*C_MS**2*(1 - n_star)),
        one_loop_over_kinetic=(g*(1 - n_star))**4/(64*np.pi**2)/(RHO_C0*EV4_PER_J_M3*R_E))
    if pantheon is not None:
        shape = (1 + pantheon['zHEL'])*H_RATE*optical_distance(h, pantheon['z'])/C_MS
        out['pantheon_chi_squared'] = chi_squared(shape, pantheon)
    return out


# ---------------------------------------------------------------- calculation 2, the optical side
def completions():
    """Exponents of the tested matter completions: (a, c, b, d) in eps ~ n^a, mu ~ n^c, m_e ~ n^b, e ~ n^d."""
    fam = [('M1 matched action, fixed charges and masses', 1, 1, 0, 0, 'fixed rulers'),
           ('M2 fixed electrostatics, Z = 1/n', 0, 2, 0, 0, 'fixed rulers'),
           ('M3 matched action with electron mass proportional to n^2', 1, 1, 2, 0, 'rulers shrink as 1/n'),
           ('M4 universal clock, dtau = dt/n', None, None, None, None, 'fixed rulers'),
           ('M5 CC-1 co-scaling', None, None, None, None, 'rulers fixed in comoving units, counts grow')]
    rows = []
    for name, a, c, b, d, rulers in fam:
        if a is None:                      # not members of the power-law family; their results are archived
            p = 0. if name.startswith('M4') else 1.
            row = dict(completion=name, in_power_law_family=False, shift_exponent=p,
                       ruler_exponent=0. if name.startswith('M4') else -1.,
                       clock_ratio_exponents=dict(hyperfine_over_optical=0., fine_over_gross=0., cavity_over_atomic=0.),
                       note=rulers)
        else:
            q = b + 4*d - 2*a                                   # atomic optical frequency exponent
            p = q + (a + c)/2                                   # measured shift exponent
            row = dict(completion=name, in_power_law_family=True, exponents=dict(a=a, c=c, b=b, d=d),
                       shift_exponent=p, atomic_frequency_exponent=q, ruler_exponent=a - b - 2*d,
                       clock_ratio_exponents=dict(hyperfine_over_optical=c + 8*d + 2*b - 3*a - q,
                                                  fine_over_gross=2*(2*d - a),
                                                  cavity_over_atomic=-(a + c)/2 - (a - b - 2*d) - q),
                       note=rulers)
        drifts = {k: v*H_RATE*YR_S for k, v in row['clock_ratio_exponents'].items()}
        row['clock_ratio_drifts_per_year'] = drifts
        row['gates'] = dict(O1_full_positive_shift=bool(abs(row['shift_exponent'] - 1) < 1e-12),
                            O2_stretch_matches=True,
                            O3_clock_ratios=bool(max(abs(v) for v in drifts.values()) <= CLOCK_BOUND_PER_YR),
                            O4_ruler_counts_fixed=bool(abs(row['ruler_exponent']) < 1e-12))
        row['passes'] = bool(all(row['gates'].values()))
        rows.append(row)
    return rows


def exponent_identity():
    """With fixed rulers and a fixed Coulomb alpha, the hyperfine-to-optical ratio drifts at 2 p ndot/n."""
    checks = []
    for a in (0., .5, 1., 2.):
        d, b = a/2, 0.                                       # fixed alpha, then fixed rulers forces b = 0
        for c in (-a, 0., 2., 3.):
            q = b + 4*d - 2*a
            p = q + (a + c)/2
            hyperfine = c + 8*d + 2*b - 3*a - q
            checks.append(dict(a=a, c=c, b=b, d=d, shift_exponent=p, atomic_frequency_exponent=q,
                               ruler_exponent=a - b - 2*d, hyperfine_over_optical_exponent=hyperfine,
                               identity_residual=abs(hyperfine - 2*p)))
    worst = max(x['identity_residual'] for x in checks)
    return dict(statement='fixed rulers and fixed alpha force d = a/2 and b = 0; then the atomic frequency is '
                          'fixed and the hyperfine-to-optical ratio drifts at exactly 2 p ndot/n',
                cases=checks, max_residual=worst, holds=bool(worst < 1e-12),
                ndot_over_n_per_year=H_RATE*YR_S, clock_bound_per_year=CLOCK_BOUND_PER_YR,
                largest_allowed_shift_exponent=CLOCK_BOUND_PER_YR/(2*H_RATE*YR_S))


def pantheon_inputs():
    """The same rows, covariance and calibration PF-1's T4 used."""
    df = pd.read_csv(DATA/'Pantheon+SH0ES.dat', sep=r'\s+')
    raw = np.loadtxt(DATA/'Pantheon+SH0ES_STAT+SYS.cov')
    n = int(raw[0])
    cov = raw[1:].reshape(n, n)
    cov = (cov + cov.T)/2
    cal = np.flatnonzero(df.IS_CALIBRATOR.to_numpy() == 1)
    hi = np.flatnonzero((df.IS_CALIBRATOR.to_numpy() == 0) & (df.zHD.to_numpy() >= .1))
    Cc, Ch, Hc = cov[np.ix_(cal, cal)], cov[np.ix_(hi, hi)], cov[np.ix_(hi, cal)]
    w = cho_solve(cho_factor(Cc), np.ones(len(cal)))
    w /= w.sum()
    cross = Hc@w
    V = Ch + float(w@Cc@w) - cross[:, None] - cross[None, :]
    return dict(V=V, factor=cho_factor(V), z=df.zHD.to_numpy()[hi], zHEL=df.zHEL.to_numpy()[hi],
                obs=df.m_b_corr.to_numpy()[hi], rows=len(hi),
                sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (DATA/'Pantheon+SH0ES.dat', DATA/'Pantheon+SH0ES_STAT+SYS.cov')})


def chi_squared(shape, P):
    """Shape-only fit with one free magnitude offset, as PF-1's T4 scored its own law."""
    ones = np.ones(len(shape))
    r0 = P['obs'] - (25 + 5*np.log10(shape))
    vi1 = cho_solve(P['factor'], ones)
    r = r0 - float(vi1@r0/(ones@vi1))
    return float(r@cho_solve(P['factor'], r))


def brightness(P):
    """PF-1's coasting law, the owner's damping approximation, and its magnitudes at z = 1."""
    pf1 = chi_squared((1 + P['zHEL'])*np.log1p(P['z']), P)
    damping = []
    for gam in (.1, .5, 1.):
        shape = (1 + P['zHEL'])*np.log1p((1 + gam)*P['z'])/(1 + gam)
        damping.append(dict(damping_over_rate=gam, chi_squared=chi_squared(shape, P),
                            delta_mu_z1_mag=float(5*np.log10(np.log(2 + gam)/((1 + gam)*np.log(2.)))),
                            owner_delta_mu_z1_mag={.1: -.059, .5: -.274, 1.: -.505}[gam]))
    return dict(rows=P['rows'], pf1_coasting_chi_squared=pf1, archived_pf1_chi_squared=871.5520519798098,
                flrw_chi_squared=836.5122548533063, linear_damping=damping,
                pipeline_matches_archive=bool(abs(pf1 - 871.5520519798098) < .01))


# ---------------------------------------------------------------- driver
def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    v1 = trial_crossing_spectrum()
    traps = [trial_trapping(4.), trial_trapping(1.)]
    controls, control_checks, runs = trial_reference_and_controls()
    probe = trial_probe_signals(*runs['reference: radiation only, field at rest'])
    hybrid = trial_hybrid(1.)
    part_a = dict(
        crossing_spectrum=v1,
        trapping=[{k: v for k, v in t.items() if k not in ('run', 'field')} for t in traps],
        reference_and_controls=dict(trials=controls, checks=control_checks), probe_signals=probe, hybrid=hybrid,
        energy_conservation=dict(worst_relative_error=max([t['energy_error'] for t in traps] +
                                                          [c['energy_error'] for c in controls]),
                                 tolerance=1e-8))
    part_a['energy_conservation']['passed'] = bool(part_a['energy_conservation']['worst_relative_error'] <= 1e-8)

    P = pantheon_inputs()
    ref = scenario(REF_NSTAR, REF_RE, P)
    scan = [scenario(ns, re, P) for ns in N_STARS for re in R_ES]
    at_ref = [s for s in scan if s['n_star'] == REF_NSTAR]
    within1 = [s['R_E'] for s in at_ref if s['pantheon_chi_squared'] <= 871.5520519798098 + 1]
    lo, hi = 10., 1000.                      # the stored energy at which the brightness cost reaches chi^2 + 1
    for _ in range(30):
        mid = np.sqrt(lo*hi)
        if scenario(REF_NSTAR, mid, P)['pantheon_chi_squared'] > 871.5520519798098 + 1:
            lo = mid
        else:
            hi = mid
    r_e_threshold = float(np.sqrt(lo*hi))
    driven = []
    h_ref = history(REF_NSTAR, REF_RE)
    for ahat in A_HATS:
        extra = h_ref['F_C'] - ahat*H_RATE**2*h_ref['K']
        s = scenario(REF_NSTAR, REF_RE, P, extra_force=extra)
        if np.isfinite(s.get('pantheon_chi_squared', np.nan)):
            driven.append(dict(net_acceleration=ahat, potential_slope_over_companion_force=extra/h_ref['F_C'],
                               chi_squared=s['pantheon_chi_squared'], delta_mu_mag=s['delta_mu_mag'],
                               potential_energy_over_companion_energy=extra/h_ref['F_C']))
    best = min(driven, key=lambda d: d['chi_squared']) if driven else None

    opt = dict(completions=completions(), exponent_identity=exponent_identity(), brightness=brightness(P))
    decisions = dict(
        S1_cold_enough=bool(ref['v_rms_at_1_Gyr']*C_KMS <= 10 and ref['travel_to_10_km_s_m'] <= KPC_M),
        S2_abundant_without_trapping=bool(ref['trapping_distance_in_n'] > ref['needed_span_in_n']),
        S3_history_survives=bool(ref['pantheon_chi_squared'] <= 871.5520519798098 + 1),
        S3_smallest_R_E_within_1=min(within1) if within1 else None,
        S3_R_E_threshold=r_e_threshold, S3_field_energy_at_threshold_over_u_cmb=r_e_threshold*RHO_C0/U_GAMMA0,
        S3_field_energy_at_threshold_over_rho_crit=r_e_threshold*RHO_C0/RHO_CRIT,
        S4_supply_named=True,
        optical_gate_passed_by=[c['completion'] for c in opt['completions'] if c['passes']])
    result = dict(
        experiment='BRIDGE-1B: one interaction for the photon energy loss and the companions mass',
        protocol='protocol.md', scope='A homogeneous semiclassical calculation. No spatial field dynamics, no '
                 'gravity for the homogeneous energies, no companion interactions, and no transport: the RC-2a '
                 'interface is specified in the protocol and not run.',
        fixed_inputs=dict(ndot_over_n_today_per_s=H_RATE, c_alpha_km_s_Mpc=ALPHA_PER_MPC*C_KMS,
                          u_gamma_0_J_m3=U_GAMMA0, T_cmb_K=T_CMB,
                          companion_benchmark=dict(q_star_Msun_kpc3_Gyr=Q_STAR, span_Gyr=T_SPAN_GYR,
                                                   rho_C0_J_m3=RHO_C0, over_u_cmb=RHO_C0/U_GAMMA0,
                                                   over_rho_crit=RHO_C0/RHO_CRIT,
                                                   source='2B-F1 reference, revision q_sel; a comparison scale, '
                                                          'not a requirement (RC-2a stage 1)'),
                          source_power_over_radiation_power=RHO_C0/(1 - REF_NSTAR)/U_GAMMA0,
                          clock_bound_per_year=CLOCK_BOUND_PER_YR),
        part_a_code_unit_trials=part_a,
        part_b_physical_mapping=dict(reference=ref, scan=scan, driven_variant=dict(scan=driven, best=best),
                                     minimum_mass_variant=scenario(REF_NSTAR, REF_RE, P, mu0=.5)),
        calculation_2_optical=opt,
        decisions=decisions,
        checks_short_run=dict(v1_absolute_error=v1['max_absolute_error'], v1_number=v1['number_density'],
                              trap_lambda4_stop=traps[0]['stop_n'], reference_companion_energy=controls[0]['companion_energy'],
                              reference_field_kinetic=controls[0]['field_kinetic_energy'],
                              reference_k_star_eV=ref['k_star_eV'], reference_mass_today_eV=ref['mass_today_eV'],
                              reference_chi_squared=ref['pantheon_chi_squared'],
                              identity_max_residual=opt['exponent_identity']['max_residual']),
        input_sha256={**P['sha256'], 'protocol.md': hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest()},
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:4000])
    return evidence_io.finish(args, 'shared-field-bridge', text, HERE/'bridge-results.json')


if __name__ == '__main__':
    raise SystemExit(main())
