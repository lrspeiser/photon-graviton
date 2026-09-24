"""Round 15: does random motion open radiation from dark states when nothing but the wave couples the emitters?

The independent calculation (independent-r15/) postulates the coupling: relative motion detunes a quiet internal
oscillation from radiating ones, delta = chi n.(v1 - v2), the rate of change of the distance between particles.
Here nothing is postulated. N ordinary oscillators (energy |a_j|^2, radiation rate 1) sit in a ball of radius Rb and
talk only through the scalar wave, G(R) = e^{ikR}/(4 pi R) (wavelength 1, k = 2 pi, envelopes quasi-static):
    da_j/dt = -a_j/2 + (i/2) sum_l C_jl a_l,   C_jl = cos(kR)/(k sqrt(R^2 + eps^2)) + i sin(kR)/(kR)
(eps = 0.02: a soft core in the near-field, reactive part only, so the energy law below stays exact). The radiated power
is P = a^H S a with S_jj = 1, S_jl = sin(kR_jl)/(kR_jl), and d(a^H a)/dt = -P exactly. S depends only on the distances
between emitters, so a rigid rotation or a uniform boost cannot change anything; free and colliding motion can.
'dissipative' drops the reactive part (C = i S): the cleanest dark states. Motion: free (Gaussian speeds, rms sigma per
component, bouncing off the ball's wall), colliding (speeds redrawn at rate nu), rigid rotation about z (angular speed
sqrt(5) sigma / Rb), uniform boost (all at sqrt(3) sigma). q = k sigma (Doppler rate in units of the radiation rate).

Parts:
  spectrum  the decay rates of the collective modes at rest (the eigenvalues of -(1/2)(1 - iC)): how dark, how gapped.
  window    settle at rest for 400 (the bright content radiates away), then run 200 more at rest or moving, from the
            same state: extra energy radiated / energy radiated at rest.
  steady    start from a random state and move for 3000: the leak rate P/E averaged over the second half, against q
            (does it grow as q^2, as the heat weight needs?), and with collisions.

    python code/wave_dark_v15.py --output run-reservoir-force-v15/wave_dark_v15.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from numba import njit

K = 2 * np.pi
EPS = 0.02


@njit(cache=True)
def coupling(x, dissip_only, eps):
    n = x.shape[0]
    C = np.zeros((n, n), dtype=np.complex128)
    S = np.eye(n).astype(np.complex128)
    for j in range(n):
        for l in range(j + 1, n):
            d0 = x[j, 0] - x[l, 0]; d1 = x[j, 1] - x[l, 1]; d2 = x[j, 2] - x[l, 2]
            R = np.sqrt(d0 * d0 + d1 * d1 + d2 * d2); kr = 2 * np.pi * R
            s = np.sin(kr) / kr if kr > 1e-9 else 1.0
            c = 0.0 if dissip_only else np.cos(kr) / (2 * np.pi * np.sqrt(R * R + eps * eps))
            C[j, l] = c + 1j * s; C[l, j] = C[j, l]; S[j, l] = s; S[l, j] = s
    return C, S


@njit(cache=True)
def rhs(a, C):
    return -0.5 * a + 0.5j * (C @ a)


@njit(cache=True)
def run(x0, a0, kind, sigma, nu, Rb, T, dt, seed, dissip_only, Om, eps, nrec):
    """kind: 0 rest, 1 free, 2 colliding, 3 rigid rotation, 4 uniform boost. The state is renormalised every step;
    returns (final state, positions, ln(E_end/E_start), leak rate P/E sampled nrec times, energy radiated / E_start)."""
    np.random.seed(seed)
    n = x0.shape[0]
    x = x0.copy(); a = a0.copy()
    a = a / np.sqrt(np.real(np.conj(a) @ a))
    v = np.zeros((n, 3))
    if kind == 1 or kind == 2:
        for j in range(n):
            for c in range(3): v[j, c] = sigma * np.random.normal()
    if kind == 4:
        for j in range(n): v[j, 0] = sigma * np.sqrt(3.0)
    steps = int(T / dt); every = max(steps // nrec, 1)
    rec = np.zeros(steps // every + 1); lnE = 0.0; ir = 0; Erad = 0.0
    for s in range(steps):
        C, S = coupling(x, dissip_only, eps)
        P = np.real(np.conj(a) @ (S @ a))
        Erad += dt * P * np.exp(lnE)
        if s % every == 0 and ir < rec.shape[0]:
            rec[ir] = P; ir += 1
        k1 = rhs(a, C); k2 = rhs(a + 0.5 * dt * k1, C); k3 = rhs(a + 0.5 * dt * k2, C); k4 = rhs(a + dt * k3, C)
        a = a + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        nn = np.real(np.conj(a) @ a); lnE += np.log(nn); a = a / np.sqrt(nn)
        if kind == 1 or kind == 2:
            x += v * dt
            for j in range(n):
                r = np.sqrt(x[j, 0] ** 2 + x[j, 1] ** 2 + x[j, 2] ** 2)
                if r > Rb:
                    nr0 = x[j, 0] / r; nr1 = x[j, 1] / r; nr2 = x[j, 2] / r
                    vn = v[j, 0] * nr0 + v[j, 1] * nr1 + v[j, 2] * nr2
                    v[j, 0] -= 2 * vn * nr0; v[j, 1] -= 2 * vn * nr1; v[j, 2] -= 2 * vn * nr2
                    f = (2 * Rb - r) / r
                    x[j, 0] *= f; x[j, 1] *= f; x[j, 2] *= f
            if kind == 2:
                for j in range(n):
                    if np.random.random() < nu * dt:
                        for c in range(3): v[j, c] = sigma * np.random.normal()
        elif kind == 3:
            c_ = np.cos(Om * dt); s_ = np.sin(Om * dt)
            for j in range(n):
                x0_ = x[j, 0]; x1_ = x[j, 1]
                x[j, 0] = c_ * x0_ - s_ * x1_; x[j, 1] = s_ * x0_ + c_ * x1_
        elif kind == 4:
            x += v * dt
    return a, x, lnE, rec[:ir], Erad


def ball(n, R, rng, dmin=0.04):
    pts = []
    while len(pts) < n:
        p = rng.uniform(-1, 1, 3) * R
        if np.linalg.norm(p) < R and all(np.linalg.norm(p - q) > dmin for q in pts):
            pts.append(p)
    return np.array(pts)


KINDS = dict(rest=0, free=1, collisional=2, rotation=3, boost=4)


def one(job):
    part, cp, seed, kind, q, nu, N, Rb = job
    dis = cp == 'dissipative'
    rng = np.random.default_rng(seed)
    x0 = ball(N, Rb, rng); a0 = (rng.normal(size=N) + 1j * rng.normal(size=N)).astype(np.complex128)
    sig = q / K; Om = np.sqrt(5) * sig / Rb
    if part == 'window':
        a, x, lnE, _, _ = run(x0, a0, 0, 0.0, 0.0, Rb, 400.0, 0.01, seed, dis, 0.0, EPS, 10)
        _, _, lnE1, rec, Erad = run(x, a, KINDS[kind], sig, nu, Rb, 200.0, 0.01, 7, dis, Om, EPS, 200)
        return dict(part=part, coupling=cp, seed=seed, kind=kind, q=q, nu=nu, radiated=float(Erad), settle_lnE=float(lnE))
    _, _, lnE, rec, _ = run(x0, a0, KINDS[kind], sig, nu, Rb, 3000.0, 0.01, 7, dis, Om, EPS, 400)
    m = len(rec)
    return dict(part=part, coupling=cp, seed=seed, kind=kind, q=q, nu=nu, leak_rate=float(np.mean(rec[m // 2:])),
                leak_rate_last_quarter=float(np.mean(rec[3 * m // 4:])), lnE=float(lnE))


def spectrum(N, Rb, seeds):
    out = []
    for cp in ('full', 'dissipative'):
        for seed in seeds:
            rng = np.random.default_rng(seed); x0 = ball(N, Rb, rng)
            C, S = coupling(x0, cp == 'dissipative', EPS)
            rates = np.sort(-2 * np.linalg.eigvals(0.5j * C - 0.5 * np.eye(N)).real)
            out.append(dict(coupling=cp, seed=seed, rates=rates.tolist(),
                            decades=[int(np.sum((rates >= 10.0 ** e) & (rates < 10.0 ** (e + 1)))) for e in range(-11, 2)]))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    N, Rb = 40, 0.25
    jobs = []
    for cp in ('full', 'dissipative'):
        for seed in (1, 2, 3, 4):
            jobs.append(('window', cp, seed, 'rest', 0.0, 0.0, N, Rb))
            for q in (0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128):
                jobs.append(('window', cp, seed, 'free', q, 0.0, N, Rb))
            for nu in (1.0, 10.0):
                jobs.append(('window', cp, seed, 'collisional', 0.016, nu, N, Rb))
            jobs.append(('window', cp, seed, 'rotation', 0.016, 0.0, N, Rb))
            jobs.append(('window', cp, seed, 'boost', 0.016, 0.0, N, Rb))
        for seed in (1, 2, 3):
            for q in (0.0, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128):
                jobs.append(('steady', cp, seed, 'free' if q else 'rest', q, 0.0, N, Rb))
            for nu in (1.0, 10.0, 100.0):
                jobs.append(('steady', cp, seed, 'collisional', 0.032, nu, N, Rb))
            jobs.append(('steady', cp, seed, 'rotation', 0.032, 0.0, N, Rb))
    run(np.zeros((2, 3)) + np.array([[0, 0, 0], [0.1, 0, 0]]), np.ones(2, np.complex128), 0, 0.0, 0.0, 1.0, 0.02, 0.01, 1, False, 0.0, EPS, 1)
    with Pool(args.processes) as pool:
        res = []
        for r in pool.imap_unordered(one, jobs):
            res.append(r)
            print(f"[{time.monotonic() - t0:5.0f} s] {r['part']:6s} {r['coupling']:11s} seed {r['seed']} {r['kind']:11s} q {r['q']:.3f} nu {r['nu']:5.1f}: "
                  + (f"radiated {r['radiated']:.4e}" if r['part'] == 'window' else f"leak {r['leak_rate']:.4e}"), flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 15: motion-opened radiation from dark states of wave-coupled emitters',
                                           N=N, Rb=Rb, eps=EPS, spectrum=spectrum(N, Rb, (1, 2, 3, 4)), runs=res,
                                           seconds=time.monotonic() - t0), indent=1) + '\n')
    print(f'wrote {args.output} ({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
