"""Round 18, step 1: the companion as a flowing medium, from local equations (one dimension, exact).

Round 17 made the companion's wave one-way by deleting the inward couplings by hand. Here the wave is given its own
local equations, in two versions, and everything the review asked for is computed from them: which way disturbances
travel, the three speeds (crests, energy, relative to the stream), the energy and momentum the pieces exchange with
the medium, and the force on a piece.

A. Waves carried by the stream. The companion streams at U and its waves move at c relative to it:
       (d_t + U d_x)^2 psi - c^2 d_x^2 psi = sum_j Q_j(t) delta(x - x_j)          (the 1D analogue-acoustics equation)
   The pieces enter through L_int = sum_j Q_j psi(x_j): piece j feeds the medium the power Q_j d_t psi(x_j) and feels
   the force Q_j d_x psi(x_j). With Q_j = Re(q_j e^{-i w t}) the time-harmonic Green's function is exact:
       U > c:  g(x) = (i / 2wc) (e^{i k_f x} - e^{i k_s x}) for x > 0, 0 for x < 0,   k_f = w/(U+c), k_s = w/(U-c)
       U < c:  g(x) = (i / 2wc) e^{i k_+ x} (x > 0),  (i / 2wc) e^{i k_- x} (x < 0),   k_+- = w/(U +- c)
   The conserved quantities of a steady medium are the pseudo-energy H = pi^2/2 - U pi psi_x + c^2 psi_x^2/2 (pi =
   psi_t + U psi_x) and the pseudo-momentum; per mode the pseudo-energy flux is +-(w^2/c)|A|^2 (minus for the slow
   mode, whose frequency in the stream's frame is negative) and the pseudo-momentum flux is k/w times it.

B. Waves that travel through the matter's frame at c, with the stream absorbing the part that moves against it:
       (d_t - c d_x) u_+ = S - c k_L(x) u_+,   (d_t + c d_x) u_- = S - c k_R(x) u_-   (u_+ moves left, u_- right)
   with psi_t = (u_+ + u_-)/2 and psi_x = (u_+ - u_-)/(2c) (the undamped equations are psi_tt - c^2 psi_xx = S); an
   outflow from x = 0 damps left-movers at x > 0 and right-movers at x < 0, amplitude e^{-kappa} per unit length
   travelled. What is absorbed, energy and momentum, goes to the stream (booked).

Checks: both Green's functions are verified against a direct time-domain simulation of the local equations with
point-like pieces, forces and powers measured as time averages of Q d_x psi and Q d_t psi.
"""
import argparse, json, time
from pathlib import Path
import numpy as np

W = 1.0          # the pieces' frequency
C = 1.0          # the wave speed (relative to the stream in A, relative to matter in B)


# --------------------------------------------------------------------------------------- A: waves carried by the stream
def g_carried(x, U, c=C, w=W):
    x = np.asarray(x, float)
    pre = 1j / (2 * w * c)
    if U > c:
        kf, ks = w / (U + c), w / (U - c)
        return np.where(x > 0, pre * (np.exp(1j * kf * x) - np.exp(1j * ks * x)), 0.0 + 0j)
    kp, km = w / (U + c), w / (U - c)
    return np.where(x >= 0, pre * np.exp(1j * kp * x), pre * np.exp(1j * km * x))


def dg_carried(x, U, c=C, w=W):
    x = np.asarray(x, float)
    pre = 1j / (2 * w * c)
    if U > c:
        kf, ks = w / (U + c), w / (U - c)
        return np.where(x > 0, pre * (1j * kf * np.exp(1j * kf * x) - 1j * ks * np.exp(1j * ks * x)), 0.0 + 0j)
    kp, km = w / (U + c), w / (U - c)
    return np.where(x >= 0, pre * 1j * kp * np.exp(1j * kp * x), pre * 1j * km * np.exp(1j * km * x))


def self_terms_carried(U, c=C, w=W):
    """Per |q|^2: the pseudo-energy a lone piece feeds the medium, the force on it (+ = downstream), and the power it
    radiates in the stream's own frame (pseudo-energy + U x drag)."""
    if U > c:
        P = 0.0                                                     # g(0) = 0: the two modes' pseudo-energies cancel
        F = 1.0 / (4 * (U ** 2 - c ** 2))                           # average of the one-sided slopes 1/(U^2 - c^2) and 0
    else:
        P = 1.0 / (4 * c)
        F = U / (4 * c * (c ** 2 - U ** 2))
    return P, F, P + U * F


def receiver_carried(U, lead=np.pi / 2, c=C, w=W, dmax=None, n=4001):
    """A receiver at distance d downstream of a source (|q_S| = 1), locked with the given lead to the source's wave
    at its position, |q_R| = 1. Returns the coherent force and the pseudo-energy it feeds, per |q_R||phi|, averaged over
    d across many beat lengths, and their ratio (the 'speed' in P = F v)."""
    if U > c:
        beat = 2 * np.pi / (w / (U - c) - w / (U + c))
    else:
        beat = 2 * np.pi / (w / (U + c))
    dmax = dmax or 20 * beat
    d = np.linspace(0.05 * beat, dmax, n)
    phi = g_carried(d, U, c, w); dphi = dg_carried(d, U, c, w)
    qR = np.exp(-1j * lead) * phi / np.abs(phi)                     # leads the local wave by `lead`
    F = 0.5 * np.real(np.conj(qR) * dphi)                           # coherent force (- = toward the source)
    P = 0.5 * w * np.imag(np.conj(qR) * phi)                        # pseudo-energy fed into the medium
    amp = np.abs(phi)
    # with a quarter-cycle lead F = -(k_f + k_s)|phi|/4 exactly (supersonic) and P = (w/2)|phi|: both follow the height
    return dict(F_per_amp=float(np.mean(F) / np.mean(amp)), P_per_amp=float(np.mean(P) / np.mean(amp)),
                F_over_amp_spread=float(np.std(F / np.maximum(amp, 1e-12)) / abs(np.mean(F / np.maximum(amp, 1e-12)))),
                v_eff=float(np.mean(P) / abs(np.mean(F))), F_mean=float(np.mean(F)), P_mean=float(np.mean(P)),
                amp_mean=float(np.mean(amp)))


def mode_fluxes_carried(U, c=C, w=W):
    """For a lone piece (|q| = 1): pseudo-energy and pseudo-momentum fluxes carried away downstream and upstream."""
    pre = 1 / (2 * w * c)
    if U > c:
        kf, ks = w / (U + c), w / (U - c)
        Ef, Es = (w ** 2 / c) * pre ** 2, -(w ** 2 / c) * pre ** 2   # per-mode pseudo-energy flux, time averaged x 2
        return dict(E_down=0.5 * (Ef + Es), E_up=0.0, P_down=0.5 * (kf / w * Ef + ks / w * Es), P_up=0.0)
    kp, km = w / (U + c), w / (U - c)
    E = (w ** 2 / c) * pre ** 2
    return dict(E_down=0.5 * E, E_up=0.5 * E, P_down=0.5 * kp / w * E, P_up=0.5 * (-km / w) * E)


# ------------------------------------------------------------------------- B: waves the stream absorbs when they go upstream
def g_absorbed(xr, xs, kappa, c=C, w=W):
    """Green's function of the two-way wave with the counter-streaming part absorbed: from a piece at xs to a point xr,
    for an outflow from x = 0. Emission is split equally between the two directions (no recoil); the amplitude of the
    part moving against the stream falls by e^{-kappa x distance travelled inward}."""
    xr = np.asarray(xr, float)
    k = w / c
    d = np.abs(xr - xs)
    right = xr >= xs
    # distance travelled inward: moving right from xs < 0 toward 0, or moving left from xs > 0 toward 0
    inward_R = np.clip(np.minimum(xr, 0.0) - np.minimum(xs, 0.0), 0, None) * (xs < 0)
    inward_L = np.clip(np.maximum(xs, 0.0) - np.maximum(xr, 0.0), 0, None) * (xs > 0)
    att = np.where(right, np.exp(-kappa * inward_R), np.exp(-kappa * inward_L))
    return (1j / (2 * w * c)) * np.exp(1j * k * d) * att


# ------------------------------------------------------------------------------- time-domain checks of the local equations
def simulate_carried(U, pieces, c=C, w=W, L=(-30.0, 60.0), h=0.01, T=160.0, avg_from=100.0, width=0.03):
    """psi on a grid, (d_t + U d_x)^2 psi - c^2 psi_xx = sum Q_j delta_h(x - x_j), as psi_t = pi - U psi_x,
    pi_t = -U pi_x + c^2 psi_xx + S, second-order centred differences, RK4, sponge layers at both ends. Pieces are
    prescribed oscillators Q_j = Re(q_j e^{-iwt}) (ramped on over t < 20). Returns time-averaged force and power."""
    x = np.arange(L[0], L[1], h); n = len(x)
    dt = 0.25 * h / (abs(U) + c)
    sponge = np.zeros(n); edge = 8.0
    sponge += np.where(x < L[0] + edge, ((L[0] + edge - x) / edge) ** 2, 0) * 4.0
    sponge += np.where(x > L[1] - edge, ((x - L[1] + edge) / edge) ** 2, 0) * 4.0
    kern = [np.exp(-0.5 * ((x - p['x']) / width) ** 2) / (np.sqrt(2 * np.pi) * width) for p in pieces]

    def dx(f):
        g = np.empty_like(f); g[1:-1] = (f[2:] - f[:-2]) / (2 * h); g[0] = g[-1] = 0.0; return g

    def dxx(f):
        g = np.empty_like(f); g[1:-1] = (f[2:] - 2 * f[1:-1] + f[:-2]) / h ** 2; g[0] = g[-1] = 0.0; return g

    def rhs(t, psi, pi):
        ramp = min(1.0, t / 20.0)
        S = sum(ramp * np.real(p['q'] * np.exp(-1j * w * t)) * kk for p, kk in zip(pieces, kern))
        return pi - U * dx(psi) - sponge * psi, -U * dx(pi) + c ** 2 * dxx(psi) + S - sponge * pi

    psi = np.zeros(n); pi = np.zeros(n); t = 0.0
    F = np.zeros(len(pieces)); P = np.zeros(len(pieces)); cnt = 0
    while t < T:
        k1 = rhs(t, psi, pi); k2 = rhs(t + dt / 2, psi + dt / 2 * k1[0], pi + dt / 2 * k1[1])
        k3 = rhs(t + dt / 2, psi + dt / 2 * k2[0], pi + dt / 2 * k2[1]); k4 = rhs(t + dt, psi + dt * k3[0], pi + dt * k3[1])
        psi = psi + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); pi = pi + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        t += dt
        if t > avg_from:
            psix = dx(psi); psit = pi - U * psix
            for j, (p, kk) in enumerate(zip(pieces, kern)):
                Q = np.real(p['q'] * np.exp(-1j * w * t))
                F[j] += Q * np.sum(kk * psix) * h; P[j] += Q * np.sum(kk * psit) * h
            cnt += 1
    return F / cnt, P / cnt


def simulate_absorbed(pieces, kappa, c=C, w=W, L=(-30.0, 30.0), h=0.01, T=120.0, avg_from=70.0, width=0.03):
    """The two one-way amplitudes with the counter-streaming part absorbed (outflow from x = 0), upwind differences
    (second order), RK4, open ends. Pieces drive both amplitudes equally and feel psi_x = (u_+ - u_-)/2."""
    x = np.arange(L[0], L[1], h); n = len(x)
    dt = 0.3 * h / c
    kL = np.where(x > 0, kappa, 0.0); kR = np.where(x < 0, kappa, 0.0)
    kern = [np.exp(-0.5 * ((x - p['x']) / width) ** 2) / (np.sqrt(2 * np.pi) * width) for p in pieces]

    def d_left(f):   # derivative for a left-mover (information from the right): upwind uses f[i+1]
        g = np.zeros_like(f); g[:-2] = (-3 * f[:-2] + 4 * f[1:-1] - f[2:]) / (2 * h); return -g

    def d_right(f):  # derivative for a right-mover: upwind uses f[i-1]
        g = np.zeros_like(f); g[2:] = (3 * f[2:] - 4 * f[1:-1] + f[:-2]) / (2 * h); return g

    def rhs(t, up, um):
        ramp = min(1.0, t / 20.0)
        S = sum(ramp * np.real(p['q'] * np.exp(-1j * w * t)) * kk for p, kk in zip(pieces, kern))
        # (d_t - c d_x) u_+ = S - c kL u_+  ->  u_+,t = c d_x u_+ + ...;  (d_t + c d_x) u_- = S - c kR u_-
        return -c * d_left(up) + S - c * kL * up, -c * d_right(um) + S - c * kR * um

    up = np.zeros(n); um = np.zeros(n); t = 0.0
    F = np.zeros(len(pieces)); P = np.zeros(len(pieces)); cnt = 0
    while t < T:
        k1 = rhs(t, up, um); k2 = rhs(t + dt / 2, up + dt / 2 * k1[0], um + dt / 2 * k1[1])
        k3 = rhs(t + dt / 2, up + dt / 2 * k2[0], um + dt / 2 * k2[1]); k4 = rhs(t + dt, up + dt * k3[0], um + dt * k3[1])
        up = up + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); um = um + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        t += dt
        if t > avg_from:
            psix = 0.5 * (up - um) / c; psit = 0.5 * (up + um)
            for j, (p, kk) in enumerate(zip(pieces, kern)):
                Q = np.real(p['q'] * np.exp(-1j * w * t))
                F[j] += Q * np.sum(kk * psix) * h; P[j] += Q * np.sum(kk * psit) * h
            cnt += 1
    return F / cnt, P / cnt


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--no-sim', action='store_true', help='skip the time-domain checks')
    args = ap.parse_args(); out = args.output_dir
    if out.exists():
        raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    res = dict(experiment='round 18: the companion as a flowing medium, one dimension, exact', units='c = w = 1')

    # A. waves carried by the stream: speeds, self terms, the receiver's pull and its energy bill
    rows = []
    for M in (0.0, 0.5, 0.9, 1.05, 1.1, 1.2, 1.3, 2 ** 0.5, 1.5, 2.0, 3.0):
        U = M * C
        Ps, Fs, Pm = self_terms_carried(U)
        rec = receiver_carried(U)
        speeds = dict(fast_crest=U + C, slow_crest=U - C, fast_energy=U + C, slow_energy=U - C)
        rows.append(dict(M=M, one_way=bool(U > C), speeds=speeds, self_pseudo_energy=Ps, self_force_downstream=Fs,
                         stream_frame_power=Pm, drag_over_stream_power_per_U=(Fs * U / Pm) if Pm else None,
                         receiver=rec, v_eff_over_U=(rec['v_eff'] / U) if U else None,
                         fluxes=mode_fluxes_carried(U)))
        print(f"A  M = {M:5.3f} ({'one-way' if U > C else 'two-way'}): lone piece feeds pseudo-energy {Ps:.3f}, feels drag {Fs:+.3f} "
              f"(stream-frame power {Pm:.3f}); a locked receiver: force {rec['F_per_amp']:+.3f} per |q||phi|, "
              f"P/|F| = {rec['v_eff']:.3f}" + (f" = {rec['v_eff'] / U:.3f} U" if U else ''), flush=True)
    res['carried'] = rows
    # the window: one-way (U > c) and the galaxies' energy bill P/|F| <= U/2  ->  (U^2 - c^2)/U <= U/2  ->  c >= U/sqrt 2
    res['carried_window'] = dict(one_way='c < U', energy_bill='P/|F| = (U^2 - c^2)/U <= U/2  <=>  c >= U/sqrt(2)',
                                 c_over_U=[1 / np.sqrt(2), 1.0])

    # B. waves absorbed when moving against the stream: coupling in and out, no drag
    rowsB = []
    for kappa in (0.0, 0.5, 1.0, 2.3, 5.0):
        inward = abs(g_absorbed(1.0, 2.0, kappa)) / abs(g_absorbed(2.0, 1.0, kappa))
        across = abs(g_absorbed(-1.5, 1.5, kappa)) / abs(g_absorbed(1.5, -1.5, 0.0))
        rowsB.append(dict(kappa=kappa, leak_inward_over_1=float(inward), across_centre_from_1p5=float(across),
                          self_force=0.0, v_eff=C))
        print(f"B  kappa = {kappa:3.1f}: a wave from x = 2 reaches x = 1 at {inward:.3f} of the outward coupling; "
              f"across the centre from 1.5 to -1.5 at {across:.3f}; lone piece drag 0; P/|F| = c", flush=True)
    res['absorbed'] = rowsB

    if not args.no_sim:
        checks = []
        for U in (0.5, 1.3):
            pieces = [dict(x=0.0, q=1.0 + 0j), dict(x=7.3, q=0.0 + 0j)]
            # the receiver leads the source's wave at its position by a quarter cycle, |q_R| = 1
            phiR = complex(g_carried(7.3, U)); pieces[1]['q'] = np.exp(-1j * np.pi / 2) * phiR / abs(phiR)
            F, P = simulate_carried(U, pieces)
            # exact: source's self terms, the receiver's self terms plus its coherent force from the source, and the
            # source's force from the receiver's wave (zero if the flow is supersonic: the receiver is downstream)
            Ps, Fs, _ = self_terms_carried(U)
            qS, qR = pieces[0]['q'], pieces[1]['q']
            FR = Fs * abs(qR) ** 2 + 0.5 * np.real(np.conj(qR) * qS * dg_carried(7.3, U))
            PR = Ps * abs(qR) ** 2 + 0.5 * W * np.imag(np.conj(qR) * qS * g_carried(7.3, U))
            FS = Fs * abs(qS) ** 2 + 0.5 * np.real(np.conj(qS) * qR * dg_carried(-7.3, U))
            PS = Ps * abs(qS) ** 2 + 0.5 * W * np.imag(np.conj(qS) * qR * g_carried(-7.3, U))
            checks.append(dict(case=f'carried, U = {U}', sim=dict(F=F.tolist(), P=P.tolist()),
                               exact=dict(F=[float(FS), float(FR)], P=[float(PS), float(PR)])))
            print(f"check A, U = {U}: force source/receiver sim {F[0]:+.4f}/{F[1]:+.4f} exact {FS:+.4f}/{FR:+.4f}; "
                  f"power sim {P[0]:+.4f}/{P[1]:+.4f} exact {PS:+.4f}/{PR:+.4f}", flush=True)
        for kappa in (0.0, 2.0):
            pieces = [dict(x=1.0, q=1.0 + 0j), dict(x=6.1, q=0.0 + 0j)]
            phiR = complex(g_absorbed(6.1, 1.0, kappa)); pieces[1]['q'] = np.exp(-1j * np.pi / 2) * phiR / abs(phiR)
            F, P = simulate_absorbed(pieces, kappa)
            qS, qR = pieces[0]['q'], pieces[1]['q']
            k = W / C
            def dgA(xr, xs):
                return (1j * k * np.sign(xr - xs)) * g_absorbed(xr, xs, kappa)
            FR = 0.5 * np.real(np.conj(qR) * qS * dgA(6.1, 1.0))
            PR = abs(qR) ** 2 / (4 * C) + 0.5 * W * np.imag(np.conj(qR) * qS * g_absorbed(6.1, 1.0, kappa))
            FS = 0.5 * np.real(np.conj(qS) * qR * dgA(1.0, 6.1))
            PS = abs(qS) ** 2 / (4 * C) + 0.5 * W * np.imag(np.conj(qS) * qR * g_absorbed(1.0, 6.1, kappa))
            checks.append(dict(case=f'absorbed, kappa = {kappa}', sim=dict(F=F.tolist(), P=P.tolist()),
                               exact=dict(F=[float(FS), float(FR)], P=[float(PS), float(PR)])))
            print(f"check B, kappa = {kappa}: force source/receiver sim {F[0]:+.4f}/{F[1]:+.4f} exact {FS:+.4f}/{FR:+.4f}; "
                  f"power sim {P[0]:+.4f}/{P[1]:+.4f} exact {PS:+.4f}/{PR:+.4f}", flush=True)
        res['time_domain_checks'] = checks
    res['seconds'] = time.monotonic() - t0
    (out / 'flowing_medium_v18.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
