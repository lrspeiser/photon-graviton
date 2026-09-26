"""Round 25: the proposal's transport picture, checked. Three questions:

1. Flicker. Where two coherent companion streams meet head-on, their net flow cancels but their energy does not: it sits
   in a standing pattern. With f = |J| / (u U) the flow efficiency, the pattern's depth m obeys m^2 = 1 - f^2 exactly
   for two opposed streams; for any set of coherent plane waves 2 <dI^2> / I^2 >= 1 - f^2 (equality only for two
   opposed waves); incoherent (hot) glow shows no pattern at all once averaged over its random phases.
2. Universal turning. If the companion's direction were randomised everywhere at the rate 1/tau, tau = u/a
   (the turning length L = u^2/a = 14.8 kpc), the steady moment equations dU/dt + div J = q, u^2 grad(U)/3 = -J/tau
   give U = 3 q / (4 pi u L r) beyond L, times erfc(r / sqrt(4 D t)) for a source of age t (D = u L / 3): the energy
   density rises by 3 r / L over free streaming, and the field stops within sqrt(4 D t) of its source. A random walk
   with the same numbers (Monte Carlo) checks both.
3. The same with tau = v/a (v a local speed), as the proposal also suggests.

    python code/frustration_transport_v25.py --output run-frustration-v25/transport_v25.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.special import erfc

KPC, MYR = 3.0856775814913673e19, 3.15576e13


def flicker(rng):
    """(a) two opposed plane waves: m^2 + f^2 = 1; (b) random coherent wave sets: 2 <dI^2>/I^2 >= 1 - f^2;
    (c) incoherent: the phase-averaged intensity is flat."""
    out = dict(two_opposed=[], many=[])
    x = np.linspace(0, 2 * np.pi, 20001)
    for r in (1.0, 0.5, 0.2, 0.05):                          # intensity ratio I2 / I1
        A1, A2 = 1.0, np.sqrt(r)
        I = np.abs(A1 * np.exp(1j * x) + A2 * np.exp(-1j * x)) ** 2
        m = (I.max() - I.min()) / (I.max() + I.min()); f = (A1 ** 2 - A2 ** 2) / (A1 ** 2 + A2 ** 2)
        out['two_opposed'].append(dict(I2_over_I1=r, f=float(f), m=float(m), m2_plus_f2=float(m ** 2 + f ** 2)))
    worst = np.inf
    for trial in range(400):
        n = rng.integers(2, 9)
        k = rng.normal(size=(n, 3)); k /= np.linalg.norm(k, axis=1)[:, None]
        Iw = rng.exponential(size=n)
        f = np.linalg.norm((Iw[:, None] * k).sum(0)) / Iw.sum()
        pts = rng.uniform(-50, 50, size=(20000, 3))
        psi = (np.sqrt(Iw)[None, :] * np.exp(1j * (pts @ k.T + rng.uniform(0, 2 * np.pi, n)[None, :]))).sum(1)
        II = np.abs(psi) ** 2
        lhs = 2 * II.var() / II.mean() ** 2
        exact = 2 * (1 - np.sum(Iw ** 2) / Iw.sum() ** 2)
        worst = min(worst, exact - (1 - f ** 2))
        if trial < 6:
            out['many'].append(dict(n_waves=int(n), f=float(f), two_var_over_mean2_sampled=float(lhs), exact=float(exact),
                                    one_minus_f2=float(1 - f ** 2)))
    out['many_min_margin'] = float(worst)                   # >= 0: the bound holds in every trial
    # incoherent: average the pattern over random relative phases
    ph = rng.uniform(0, 2 * np.pi, 4000)
    Iavg = np.mean(np.abs(np.exp(1j * x[None, :]) + np.sqrt(0.5) * np.exp(-1j * x[None, :] + 1j * ph[:, None])) ** 2, 0)
    out['incoherent_depth'] = float((Iavg.max() - Iavg.min()) / (Iavg.max() + Iavg.min()))
    return out


def turning(u_kms, a_SI, age_gyr=13.0, v_list=(9.1, 230.0, 1000.0)):
    """The moment equations with universal turning: energy excess over free streaming and the reach in 13 Gyr."""
    u = u_kms * 1e3
    rows = []
    for label, tau in [('u/a', u / a_SI)] + [(f'{v:g} km/s / a', v * 1e3 / a_SI) for v in v_list]:
        L = u * tau / KPC                                     # kpc: the turning length
        D = u * (L * KPC) / 3 / KPC ** 2 * MYR                # kpc^2 / Myr
        reach = np.sqrt(4 * D * age_gyr * 1e3)                # kpc: where erfc falls off
        r = np.array([10.0, 30.0, 100.0, 300.0, 1000.0, 2000.0])
        excess = 1 + 3 * r / L                                # U / U_free (P1, steady)
        aged = excess * erfc(r / reach)                       # with the source's finite age
        rows.append(dict(tau=label, tau_Myr=float(tau / MYR), turning_length_kpc=float(L), reach_kpc=float(reach),
                         r_kpc=r.tolist(), U_over_free_streaming_steady=excess.tolist(), U_over_free_streaming_13Gyr=aged.tolist(),
                         free_streaming_reach_kpc=float(u * age_gyr * 1e3 * MYR / KPC)))
    return rows


def random_walk(u_kms, a_SI, rng, n=4000, age_gyr=13.0):
    """Monte Carlo check of the turning picture: particles leave a steady source at u, turn isotropically after
    exponential free paths of mean L = u^2/a, for the source's age; the time spent in radial shells gives U(r)."""
    u = u_kms * 1e3; L = u * u / a_SI / KPC
    T = u * age_gyr * 1e3 * MYR / KPC / L                     # the age in free paths (u t / L)
    edges = np.geomspace(3.0, 3000.0, 37) / L
    acc = np.zeros(edges.size - 1)
    dt = 0.05                                                 # sampling interval (units L / u)

    def advance(pos, d, left, span):
        """Move every particle for the times `span`, turning at the end of each free path."""
        rem = np.array(span, float)
        while (rem > 1e-12).any():
            step = np.minimum(left, rem)
            pos += d * step[:, None]; left -= step; rem -= step
            turn = left <= 1e-12
            if turn.any():
                nd = rng.normal(size=(turn.sum(), 3)); d[turn] = nd / np.linalg.norm(nd, axis=1)[:, None]
                left[turn] = rng.exponential(size=turn.sum())

    for _ in range(n // 200):
        m = 200
        pos = np.zeros((m, 3)); d = rng.normal(size=(m, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
        left = rng.exponential(size=m)
        advance(pos, d, left, rng.uniform(0, dt, m))          # a random sampling phase for each particle
        samples = [np.linalg.norm(pos, axis=1)]
        for _ in range(int(T / dt) - 1):                      # sample times spread uniformly over the source's age
            advance(pos, d, left, np.full(m, dt))
            samples.append(np.linalg.norm(pos, axis=1))
        acc += np.histogram(np.concatenate(samples), bins=edges)[0]
    # a steady source of age T: U(r) dV = q x (time particles emitted over the last T spend in the shell)
    rc = np.sqrt(edges[1:] * edges[:-1]); vol = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    U = acc * dt / vol / n                                    # per unit emission rate (units: 1/(u L^2))
    free = 1 / (4 * np.pi * rc ** 2)                          # free streaming, per unit rate, in the same units
    return dict(turning_length_kpc=float(L), age_in_free_paths=float(T), r_kpc=(rc * L).tolist(),
                U_over_free_streaming=(U / free).tolist(),
                P1_prediction=((1 + 3 * rc) * erfc(rc * L / np.sqrt(4 * (u * L * KPC / 3 / KPC ** 2 * MYR) * age_gyr * 1e3))).tolist())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    rng = np.random.default_rng(25)
    law = dict(a_SI=6.297889390049439e-11, u_kms=169.4433686899687)     # the adopted law (round 12)
    out = dict(source='code/frustration_transport_v25.py', law=law, flicker=flicker(rng), turning=turning(law['u_kms'], law['a_SI']),
               random_walk=random_walk(law['u_kms'], law['a_SI'], rng))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=1) + '\n')
    fl = out['flicker']
    print('two opposed streams: m^2 + f^2 =', [round(r['m2_plus_f2'], 6) for r in fl['two_opposed']])
    print('many coherent waves: smallest margin of 2<dI^2>/I^2 - (1 - f^2) =', round(fl['many_min_margin'], 4))
    print('incoherent pattern depth:', round(fl['incoherent_depth'], 4))
    for t in out['turning']:
        print(f"tau = {t['tau']:18s} {t['tau_Myr']:8.1f} Myr  L = {t['turning_length_kpc']:7.1f} kpc  reach (13 Gyr) = {t['reach_kpc']:7.0f} kpc "
              f"(free streaming {t['free_streaming_reach_kpc']:.0f})  U/U_free at 10/100/1000 kpc: "
              + ' / '.join(f'{x:.3g}' for x in np.array(t['U_over_free_streaming_13Gyr'])[[0, 2, 4]]))
    rw = out['random_walk']
    for r, x, p in list(zip(rw['r_kpc'], rw['U_over_free_streaming'], rw['P1_prediction']))[::5]:
        print(f'  random walk r = {r:7.1f} kpc: U/U_free = {x:8.3g}   (P1 with finite age: {p:8.3g})')
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
