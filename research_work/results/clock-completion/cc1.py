"""CC-1: the co-scaling completion as a consistency calculation (see protocol.md).

This is a separately labeled branch. It does not claim that the nonexpanding premise is met, that the past was
nonsingular, or that other completions are excluded. Matter and light couple to g_m = -dt^2 + n(t)^2 dx^2
(c = 1); the checks N1-N6 use the tolerances declared in protocol.md, and O1-O8 are derived here.

    python cc1.py [--canonical] [--output-dir DIR]

About 3 minutes. The run fails fast on a non-finite state and logs its progress.
"""
import json
import math
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
import atoms as A  # noqa: E402
import ledger as L  # noqa: E402

T_CLASSICAL, T_QUANTUM = 2e4, 4e4
T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:7.1f}s] {msg}', flush=True)


def pad_frequency(times, lo, hi):
    """Angular frequency from the zero crossings inside [lo, hi] (a straight-line fit of crossing index)."""
    t = times[(times > lo) & (times < hi)]
    return math.pi*np.polyfit(np.arange(len(t)), t, 1)[0]**-1


def free_mode(prof, law, k=1., pad=200.):
    """A spatially homogeneous mode across the ramp, integrated in the frame's own time.
    law 'field':       d/dt(n qdot) = -(k^2/n) q            (PF-1's wave law in t)
    law 'conformal':   q'' = -k^2 q                          (the same law in eta)
    law 'tensor':      d/dt(n^3 hdot) = -n k^2 h             (standard tensor waves in g_m, in t)
    law 'P2':          PF-1's postulate P2 for gravitational waves, identical to 'field'.
    Returns the frequency and amplitude in the pads before and after the ramp."""
    if law == 'conformal':
        rhs, nodes, n = (lambda s, y: [y[1], -k*k*y[0]]), [-pad, 0., prof.eta_T, prof.eta_T + pad], (lambda s: 1.)
        p = 1
    else:
        p = 3 if law == 'tensor' else 1
        nodes, n = [-pad, 0., prof.T, prof.T + 2*pad], prof.n
        rhs = lambda t, y: [y[1]/prof.n(t)**p, -prof.n(t)**(p - 2)*k*k*y[0]]

    def zero(t, y):
        return y[0]
    (tz,), _ = A.integrate(rhs, [1., 0.], nodes, [zero])
    w_pre, w_post = pad_frequency(tz, nodes[0], nodes[1]), pad_frequency(tz, nodes[2], nodes[3])

    def amplitude(t_end):
        from scipy.integrate import solve_ivp
        y = [1., 0.]
        for a, b in zip(nodes[:-1], nodes[1:]):
            if a >= t_end:
                break
            y = solve_ivp(rhs, (a, min(b, t_end)), y, method='DOP853', rtol=A.RTOL, atol=A.ATOL).y[:, -1]
        nn = n(t_end)
        return math.hypot(y[0], y[1]/(nn**p*(k/nn if law != 'conformal' else k)))
    return dict(omega_pre=float(w_pre), omega_post=float(w_post),
                amplitude_pre=float(amplitude(nodes[1])), amplitude_post=float(amplitude(nodes[3])))


def measured_redshift(n1, n2, n3):
    """N4. 1+z = (propagation factor in the frame's time) x (clock factor in the same time)."""
    prof = A.Profile(T_CLASSICAL)
    fld, conf = free_mode(prof, 'field'), free_mode(prof, 'conformal')
    prop_field = fld['omega_pre']/fld['omega_post']
    prop_conf = conf['omega_pre']/conf['omega_post']
    rows = {'field frame, classical atom (N1)': prop_field*n1['clock_factor'],
            'field frame, quantum clock (N3)': prop_field*n3['clock_factor'],
            'conformal frame, classical atom (N2)': prop_conf*n2['clock_factor']}
    target = 2.   # n_o/n_e across the ramp
    dev = {key: float(v/target - 1) for key, v in rows.items()}
    return dict(n_o_over_n_e=target, propagation_factor=dict(field=prop_field, conformal=prop_conf),
                clock_factor=dict(N1=n1['clock_factor'], N2=n2['clock_factor'], N3=n3['clock_factor']),
                one_plus_z=rows, deviation=dev, passed=bool(max(abs(v) for v in dev.values()) < 1e-6))


def cavity(ell=1., epochs=21):
    """N6. Rigid ends at proper separation ell sit at comoving 0 and ell/n(t); light moves with deta = dx.
    The ramp is slow enough that ndot ell/n <= 1e-7. For contrast, ends fixed in comoving coordinates."""
    from scipy.integrate import quad
    prof = A.Profile(.5*math.pi*1e7*ell)
    # eta differences are integrated over the short interval itself; differencing eta at t ~ 1e7 would
    # cancel away the precision a 1e-6 test needs
    deta = lambda a, d: quad(lambda s: 1/prof.n(s), a, a + d, epsabs=0, epsrel=1e-13)[0]
    rows, worst = [], 0.
    for t0 in np.linspace(0, prof.T, epochs):
        d1 = brentq(lambda d: deta(t0, d) - ell/prof.n(t0 + d), 0, 4*ell, xtol=1e-14, rtol=1e-15)
        d2 = brentq(lambda d: deta(t0 + d1, d) - ell/prof.n(t0 + d1), 0, 4*ell, xtol=1e-14, rtol=1e-15)
        fixed = brentq(lambda d: deta(t0, d) - 2*ell, 0, 8*ell, xtol=1e-14, rtol=1e-15)
        rows.append(dict(n=prof.n(t0), round_trip_over_2ell=(d1 + d2)/(2*ell), comoving_ends_round_trip_over_2ell=fixed/(2*ell)))
        worst = max(worst, abs((d1 + d2)/(2*ell) - 1))
    return dict(max_ndot_ell_over_n=.5*prof.w*ell, epochs=rows, max_deviation=float(worst), passed=bool(worst < 1e-6))


def pulse_stretch(D=2e4):
    """O3. Two emissions at t_e and t_e + dt arrive where eta(t_o) = eta(t_e) + D; the observer's atoms tick in t."""
    prof = A.Profile(T_CLASSICAL)
    t_o = lambda te: prof.t_of_eta(prof.eta(te) + D)
    te, h = -1e3, 1e-2
    return dict(emitted_before_ramp_at=te, arrival=t_o(te), stretch=(t_o(te + h) - t_o(te - h))/(2*h),
                n_o_over_n_e=prof.n(t_o(te))/prof.n(te))


def distances(z=(.1, .5, 1., 2.)):
    """O6 for a linear n (the V = 0 limit, coasting), in units of c/H0 with H0 = ndot/n today."""
    rows = []
    for zz in z:
        DM = math.log1p(zz)
        rows.append(dict(z=zz, D_M=DM, D_L=(1 + zz)*DM, D_A=DM/(1 + zz), duality=((1 + zz)*DM)/((1 + zz)**2*DM/(1 + zz))))
    return rows


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    prof_c = A.Profile(T_CLASSICAL)
    log('N1 classical atom, field frame')
    n1 = A.classical_field_frame(prof_c)
    log('N2 classical atom, conformal frame')
    n2 = A.classical_conformal_frame(prof_c)
    a, b = n1.pop('crossing_times'), n2.pop('crossing_times_as_t')
    n2['max_crossing_time_difference_from_N1'] = float(np.max(np.abs(a - b))) if len(a) == len(b) else None
    log(f"N1 passed={n1['passed']} dev={n1['max_frequency_deviation']:.2e}/{n1['max_proper_size_deviation']:.2e}; "
        f"N2 passed={n2['passed']} dev={n2['max_deviation_of_omega_eta_over_n']:.2e}")
    log('N3 quantum clock, field frame')
    n3 = A.quantum_clock(A.Profile(T_QUANTUM), log=log)
    log(f"N3 passed={n3['passed']} dev={n3['max_deviation']}")
    n4 = measured_redshift(n1, n2, n3)
    log(f"N4 passed={n4['passed']} dev={n4['deviation']}")
    log('N5 energy ledger')
    n5 = L.ledger()
    log(f"N5 passed={n5['passed']} checks={n5['checks']}")
    n6 = cavity()
    log(f"N6 passed={n6['passed']} dev={n6['max_deviation']:.2e}")
    gw = {law: free_mode(prof_c, law) for law in ('field', 'tensor', 'P2')}
    for g in gw.values():
        g['amplitude_ratio_post_to_pre'] = g['amplitude_post']/g['amplitude_pre']
    derived = dict(
        O1_measured_redshift='1+z = n_o/n_e in both frames (N4); the field frame puts it in propagation, the conformal frame in the clocks',
        O2_clock_ratio=dict(w20_over_w10_pre=n3['clock_ratio_pre'], max_deviation=n3['max_deviation']['ratio']),
        O3_pulse_stretch=pulse_stretch(),
        O4_local_light_speed='n dx/dt = 1 on null curves of g_m at every epoch; see N6 for cavities',
        O5_gravitational_waves=dict(
            modes=gw, arrival_time_difference='zero: both laws share the null cones of g_m',
            siren_distance_over_D_L=dict(tensor='1', P2='1/(1+z)'),
            P2_over_tensor_amplitude=gw['P2']['amplitude_ratio_post_to_pre']/gw['tensor']['amplitude_ratio_post_to_pre']),
        O6_distances_coasting=distances(),
        O8_history_V0=L.history())
    passed = all(x['passed'] for x in (n1, n2, n3, n4, n5, n6))
    result = dict(branch='co-scaling (separately labeled); a consistency calculation, not a claim that the nonexpanding premise is met',
                  N1=n1, N2=n2, N3=n3, N4=n4, N5=n5, N6=n6, derived=derived, all_passed=passed,
                  runtime_seconds=time.time() - T0)
    text = json.dumps(result, indent=1, default=float)
    log(f'all_passed={passed}')
    status = evidence_io.finish(args, 'clock-completion-cc1', text, HERE/'cc1-results.json', ignore={'/runtime_seconds'})
    return status if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
