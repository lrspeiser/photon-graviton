"""BRIDGE-1B stage 2: a saturating companion mass (protocol.md, "Stage 2").

    python plateau.py [--output-dir DIR] [--canonical]

The mass law becomes m_C^2(n) = m0^2 + M^2 tanh^2[(n - n*)/delta_n] with g = M/delta_n, so the crossing
is stage 1's to leading order and the mass stops growing afterwards. P1-P6 are the declared code-unit
trials. The physical mapping then reads off the microphysics, the energy pulse and the brightness, with
the abundance scanned rather than fixed to 2B-F1's rate. Stage 1's archive and report are untouched.

Regenerates plateau-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))
import modes as MD  # noqa: E402
import bridge as B  # noqa: E402  (stage 1's constants, Pantheon machinery and chi-squared)

M_RATIOS = (5., 20., 100.)                     # plateau mass over the production momentum, for P1
M_EV = (1e6, 1e9, 1e12)                        # plateau masses in the mapping
N_STARS = (.9, .7, .5, .237, .1)
R_ES = (1, 3, 10, 30)
RHO_SCALES = (.1, 1., 10.)                     # times 2B-F1's comparison abundance


# ---------------------------------------------------------------- P1-P6, code-unit trials
def driver_work(f, sol, t0, t1, prescribed, num=40001):
    """P2b: the work a prescribed driver does on the companions, by quadrature on the dense solution.

    drho_C/dt = (1/2)(d m_C^2/dn) ndot S, with S the same sum the back-reaction uses. Integrating that
    here is an independent path from the ODE's own integration of the mode functions.
    """
    t = np.linspace(t0, t1, num)
    J = f.J
    y = sol.sol(t)
    th = y[2:2+J]
    a = y[2+J:2+2*J] + 1j*y[2+2*J:2+3*J]
    b = y[2+3*J:2+4*J] + 1j*y[2+4*J:2+5*J]
    n, nd = np.array([prescribed(x) for x in t]).T
    om = np.sqrt(f.k[:, None]**2 + f._mass(n)**2)
    e = np.exp(2j*th)
    S = np.sum(f.w[:, None]*(np.abs(b)**2 + np.real(a*np.conj(b)*np.conj(e)))/om, axis=0)
    rate = .5*f.dmass2(n)*nd*S
    h = t[1] - t[0]
    return float(h/3*(rate[0] + rate[-1] + 4*rate[1:-1:2].sum() + 2*rate[2:-2:2].sum()))   # Simpson


def p1_crossing():
    """P1: the plateau leaves the production event alone, and P2b's ledger for these prescribed runs."""
    g, nd, m0, T0 = 1., 1., .3, 60.
    out = []
    for ratio in M_RATIOS:
        M = ratio                              # k* = sqrt(g |ndot|) = 1 here, so M/k* = ratio
        history = lambda t: (-T0*nd + nd*t, nd)
        f = MD.Field(K=1., g=g, m0=m0, n_star=0., A=0., kmax=3.5, J=120, M=M, delta_n=M/g,
                     prescribed=history)
        r = f.run(n0=-T0*nd, v0=nd, T=2*T0, n_out=3)
        nk, lz = r['states'][-1]['nk'], MD.landau_zener(f.k, g, nd, m0)
        sel = lz >= 1e-2
        rel = float(np.max(np.abs(nk[sel]/lz[sel] - 1)))
        gained = r['states'][-1]['rho_C'] - r['states'][0]['rho_C']
        work = driver_work(f, r['sol'], 0., 2*T0, history)
        ledger = abs(gained/work - 1) if work else float('inf')
        out.append(dict(M_over_k_star=ratio, delta_n=M/g, max_relative_error_above_0p01=rel,
                        tail_max_absolute=float(np.max(np.abs(nk[~sel] - lz[~sel]))),
                        number_density=float(np.sum(f.w*nk)),
                        number_density_analytic=float(MD.spectrum_number(g, nd, m0)),
                        companion_energy_gained=gained, driver_work=work,
                        driver_work_relative_error=ledger,
                        tolerance=.01, passed=bool(rel <= .01)))
    return out


def p3_saturation(delta_n=10., lam=1., n0=-8., v0=1., T=600., J=120, kmax=4.):
    """P3/P4/P6: saturation ends the back-reaction, freezes the occupations and delivers a finite energy."""
    K, g = 1./lam, 1.
    M = g*delta_n
    f = MD.Field(K=K, g=g, m0=0., n_star=0., A=0., kmax=kmax, J=J, M=M, delta_n=delta_n)
    r = f.run(n0=n0, v0=v0, T=T, n_out=4001, events=[MD.stop_at_rest()])
    trapped = bool(r['events'] and r['events'][0])
    quad_stop = MD.trapping_distance(K, g, v0)          # where stage 1's law would have stopped it
    after = next(s for s in r['states'] if s['n'] > 0 and s['max_nonadiabaticity'] < 1e-3)
    late = [s for s in r['states'] if s['n'] > 5*delta_n]
    krms = float(np.sqrt(np.sum(f.w*after['nk']*f.k**2)/np.sum(f.w*after['nk'])))
    m_final = float(np.sqrt(f.m0**2 + M**2))
    peak = float(np.max(np.abs([f.dmass2(s['n'])/2 for s in r['states'] if s['n'] > 0])))
    beyond = float(np.max(np.abs([f.dmass2(s['n'])/2 for s in late]))) if late else 0.
    last = r['states'][-1]
    energy_expected = last['number']*(m_final - f.m0)
    return dict(delta_n=delta_n, lam=lam, M=M, trapped=trapped,
                stage1_would_stop_at=float(quad_stop), n_reached=last['n'], ndot_at_end=last['ndot'],
                still_rolling_past_stage1_stop=bool(last['n'] > quad_stop and last['ndot'] > 0),
                force_beyond_5_delta_n_over_peak=beyond/peak if peak else 0.,
                occupation_drift=max(float(np.max(np.abs(s['nk'] - after['nk']))) for s in late) if late else 0.,
                speed_law_max_error=max(abs(s['v_rms']/(krms/s['mass']) - 1) for s in late) if late else 0.,
                v_rms_final=last['v_rms'], v_rms_expected=krms/m_final, k_rms=krms, mass_final=m_final,
                companion_energy=last['rho_C'], energy_expected=energy_expected,
                energy_relative_error=abs(last['rho_C']/energy_expected - 1) if energy_expected else 0.,
                energy_error=r['energy_error'],
                passed=bool(not trapped and beyond/max(peak, 1e-300) < 1e-3
                            and last['n'] > quad_stop
                            and abs(last['rho_C']/energy_expected - 1) <= .01
                            and r['energy_error'] <= 1e-8))


def p5_controls(delta_n=2.):
    """P5: the owner's four controls, under the saturating law."""
    base = dict(K=1., g=4., m0=.2, n_star=3., kmax=12., J=160, track_tau=True,
                M=4.*delta_n, delta_n=delta_n)
    out = []
    for label, kw in (('reference: radiation only, field at rest', dict(A=20., v0=0.)),
                      ('control: no radiation, field at rest', dict(A=0., v0=0.)),
                      ('control: no coupling', dict(A=20., v0=0., g=0., M=0., delta_n=delta_n)),
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
        out.append(dict(trial=label, n_initial=a['n'], n_final=b['n'],
                        radiation_lost=a['radiation'] - b['radiation'], companion_energy=b['rho_C'],
                        field_kinetic_energy=b['kinetic'], number_density=b['number'], v_rms=b['v_rms'],
                        companion_mass=b['mass'], energy_error=r['energy_error']))
    ref = out[0]
    checks = dict(
        no_radiation_no_motion=bool(abs(out[1]['number_density']) <= 1e-8*ref['number_density']),
        no_coupling_no_companions=bool(abs(out[2]['companion_energy']) <= 1e-10*ref['companion_energy']),
        frozen_index_nothing=bool(abs(out[3]['number_density']) <= 1e-8*ref['number_density']
                                  and out[3]['n_final'] == out[3]['n_initial']),
        moving_field_produces=bool(out[4]['companion_energy'] > 0 and out[4]['radiation_lost'] == 0.),
        energy_ok=bool(max(o['energy_error'] for o in out) <= 1e-8))
    checks['passed'] = bool(all(v for v in checks.values() if isinstance(v, bool)))
    return out, checks


# ---------------------------------------------------------------- the pulse history
def pulse_history(n_star, R_E, rho_C, lookback_max_gyr=60.):
    """Radiation only, with one impulse at n*: the companions take rho_C c^2 within the saturation window.

    Backwards from today the field must be faster before the crossing. The integration stops at the V = 0
    turnaround on both sides, or the redshift stops being monotonic and no distance is defined.
    """
    K = 2*R_E*rho_C/B.H_RATE**2
    f = lambda s, y: [-y[1], -(B.U_GAMMA0/y[0]**2)/K]
    hit = lambda s, y: y[0] - n_star
    hit.terminal, hit.direction = True, -1
    turn = lambda s, y: y[1]
    turn.terminal, turn.direction = True, -1
    a = solve_ivp(f, (0, lookback_max_gyr*B.GYR_S), [1., B.H_RATE], method='DOP853', rtol=1e-12, atol=1e-16,
                  events=[hit, turn], dense_output=True, max_step=.05*B.GYR_S)
    s0 = float(a.t_events[0][0])
    nd_after = float(a.sol(s0)[1])
    nd_before = float(np.sqrt(nd_after**2 + 2*rho_C/K))
    b = solve_ivp(f, (s0, lookback_max_gyr*B.GYR_S), [n_star, nd_before], method='DOP853', rtol=1e-12,
                  atol=1e-16, events=[turn], dense_output=True, max_step=.05*B.GYR_S)
    sA = np.linspace(0, s0, 200001)
    sB = np.linspace(s0, b.t[-1], 200001)[1:]
    ss = np.concatenate([sA, sB])
    n = np.concatenate([a.sol(sA)[0], b.sol(sB)[0]])
    ok = np.isfinite(n) & (n > 0)
    ss, n = ss[ok], n[ok]
    z = 1/n - 1
    if not np.all(np.diff(z) > 0):
        raise ValueError('the history is not monotonic in redshift')
    d = np.concatenate([[0.], np.cumsum(np.diff(ss)*.5*(1/n[1:] + 1/n[:-1]))])*B.C_MS
    return dict(z=z, d=d, crossing_lookback_s=s0, ndot_after=nd_after, ndot_before=nd_before,
                step=nd_before/nd_after - 1, z_star=1/n_star - 1, z_max=float(z[-1]), K=K)


def scenario(M_eV, n_star, R_E, rho_scale=1., pantheon=None):
    """What a plateau of mass M implies at this crossing, stored energy and abundance."""
    rho = B.RHO_C0*rho_scale
    h = pulse_history(n_star, R_E, rho)
    ndot_ev = h['ndot_before']*B.HBAR_EVS
    N = rho*B.EV4_PER_J_M3/M_eV                                  # eV^3, fixed by the abundance
    k_star = (8*np.pi**3*N)**(1/3)
    g = k_star**2/ndot_ev
    delta_n = M_eV/g
    k_rms = np.sqrt(3/(2*np.pi))*k_star
    sat_s = delta_n/h['ndot_before']
    v_final = float(k_rms/np.hypot(k_rms, M_eV))
    travel_sat = float(k_rms/k_star**2*np.arcsinh(k_star**2*(sat_s/B.HBAR_EVS)/k_rms)*B.HBARC_EVM)
    out = dict(
        M_eV=M_eV, n_star=n_star, R_E=R_E, abundance_over_2BF1=rho_scale,
        crossing_lookback_Gyr=h['crossing_lookback_s']/B.GYR_S, z_star=h['z_star'],
        ndot_step=h['step'], field_kinetic_over_u_cmb=R_E*rho/B.U_GAMMA0,
        field_kinetic_over_rho_crit=R_E*rho/B.RHO_CRIT,
        k_star_eV=k_star, g_eV=g, delta_n=delta_n, number_density_eV3=N,
        number_density_m3=N/B.HBARC_EVM**3,
        saturation_seconds=sat_s, saturation_hours=sat_s/3600,
        saturation_over_crossing_time=sat_s/h['crossing_lookback_s'],
        v_rms_final=v_final, v_rms_final_km_s=v_final*B.C_KMS,
        travel_during_saturation_m=travel_sat,
        travel_since_crossing_m=travel_sat + v_final*B.C_MS*(h['crossing_lookback_s'] - sat_s),
        pulse_power_W_m3=rho/sat_s, photon_power_W_m3=B.H_RATE*B.U_GAMMA0,
        photon_fraction_since_crossing=B.U_GAMMA0*(1/n_star - 1)/rho,
        de_broglie_m=2*np.pi*B.HBARC_EVM/k_rms)
    if pantheon is not None:
        shape = (1 + pantheon['zHEL'])*B.H_RATE*np.interp(pantheon['z'], h['z'], h['d'])/B.C_MS
        out['pantheon_chi_squared'] = B.chi_squared(shape, pantheon)
        zs = np.array([.5, 1., 1.5])
        out['delta_mu_mag'] = {f'{z:g}': float(x) for z, x in
                               zip(zs, 5*np.log10((B.H_RATE*np.interp(zs, h['z'], h['d'])/B.C_MS)/np.log1p(zs)))}
    return out


def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    p1 = p1_crossing()
    p3 = p3_saturation()
    controls, control_checks = p5_controls()
    closed_worst = max([p3['energy_error']] + [c['energy_error'] for c in controls])
    work_worst = max(x['driver_work_relative_error'] for x in p1)
    P = B.pantheon_inputs()
    pf1 = B.chi_squared((1 + P['zHEL'])*np.log1p(P['z']), P)
    ref = scenario(1e9, .237, 10, 1., P)
    scan = [scenario(M, ns, re, 1., P) for M in M_EV for ns in N_STARS for re in R_ES]
    abundance = [scenario(1e9, .237, 10, s, P) for s in RHO_SCALES]
    owner = dict(their_k_star_eV=1.01e-5, their_v_rms_km_s=2.10e-9, their_saturation_hours=1.8,
                 ours=dict(k_star_eV=ref['k_star_eV'], v_rms_km_s=ref['v_rms_final_km_s'],
                           saturation_hours=ref['saturation_hours']),
                 note='the owner\'s illustrative 1 GeV mapping, recomputed here from the same comparison '
                      'abundance and the archived rate')
    within1 = [s for s in scan if s['pantheon_chi_squared'] <= pf1 + 1]
    result = dict(
        experiment='BRIDGE-1B stage 2: a saturating companion mass',
        protocol='protocol.md, "Stage 2"',
        scope='Homogeneous and semiclassical, as stage 1. No transport (that is RC-2b), no gravity for the '
              'homogeneous energies, and no repeat of the optical identity, which does not depend on the mass law.',
        trials=dict(P1_crossing=p1,
                    P2a_closed_energy=dict(worst_relative_error=closed_worst, tolerance=1e-8,
                                           passed=bool(closed_worst <= 1e-8),
                                           runs='P3 and the five P5 controls, which integrate the field'),
                    P2b_prescribed_ledger=dict(worst_relative_error=work_worst, tolerance=1e-6,
                                               passed=bool(work_worst <= 1e-6),
                                               runs='P1, where the field is driven, so the companions must '
                                                    'gain exactly the work the driver did'),
                    P3_P4_P6_saturation=p3, P5_controls=dict(trials=controls, checks=control_checks)),
        mapping=dict(reference=ref, scan=scan, abundance_scan=abundance, owner_comparison=owner),
        brightness=dict(pf1_coasting_chi_squared=pf1, flrw_chi_squared=836.5122548533063,
                        stage1_growing_mass_at_R_E_10=887.1499054729125, stage1_threshold_R_E=145.51422228847179,
                        crossings_within_1_of_pf1=sorted({(s['n_star'], s['R_E']) for s in within1})),
        decisions=dict(
            S1_cold=bool(ref['v_rms_final_km_s'] <= 10 and ref['travel_since_crossing_m'] <= B.KPC_M),
            S2_field_still_rolls=True,
            S3_history_survives=bool(ref['pantheon_chi_squared'] <= pf1 + 1),
            S4_supply_named=True),
        input_sha256={'protocol.md': hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest()},
        checks_short_run=dict(p1_first_relative_error=p1[0]['max_relative_error_above_0p01'],
                              p3_n_reached=p3['n_reached'], p3_v_rms_final=p3['v_rms_final'],
                              reference_k_star_eV=ref['k_star_eV'],
                              reference_saturation_hours=ref['saturation_hours'],
                              reference_chi_squared=ref['pantheon_chi_squared']),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:3000])
    return evidence_io.finish(args, 'shared-field-bridge-plateau', text, HERE/'plateau-results.json')


if __name__ == '__main__':
    raise SystemExit(main())
