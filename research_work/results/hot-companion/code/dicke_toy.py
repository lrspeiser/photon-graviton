#!/usr/bin/env python3
"""Why colliding gas does not feed the heat term: collisional (Dicke) narrowing.

    python dicke_toy.py --output-dir ../run-dicke

Round 2 derived the heat term from Doppler scrambling: random motion shifts each emitter's
companion phase by k_c * x(t) (companion wavenumber times the emitter's displacement along
the line of sight), so hot emitters add as a scalar sum instead of cancelling like vectors.

That argument has a known loophole in laboratory physics. If an emitter changes direction
many times before it has travelled one companion wavelength, its displacement grows only
diffusively and its Doppler phase barely moves. Collisions switch the Doppler effect off.
This is Dicke narrowing (R. H. Dicke, Phys. Rev. 89, 472, 1953); the confined-emitter version
is the recoil-free (Mossbauer) line. It separates matter into two kinds with the same
temperature:

  * collisionless (stars, galaxies): free paths of kiloparsecs  -> Doppler-scrambled -> hot
  * collisional (gas, plasma):        free paths far below that -> phase stays put   -> cold

The script checks three things numerically, all in units k_c = 1, sigma = 1:

  1. Spectral narrowing: the companion phase of a colliding emitter decorrelates on the time
     1/(D k_c^2), D = sigma^2 tau_c, instead of the Doppler time 1/(k_c sigma).
  2. Steady state with phase locking at rate Gamma (emitters re-lock to the field they
     share): the phase variance is  k_c^2 sigma^2 tau_c / (Gamma (1 + Gamma tau_c)).
     Free streaming (tau_c -> infinity) gives k_c^2 sigma^2 / Gamma^2, proportional to
     sigma^2: the form of the heat weight k = 3 sigma^2/u^2, with u = sqrt(3) Gamma/k_c.
     Collisions multiply it by Gamma tau_c / (1 + Gamma tau_c).
  3. The summed field (the round-2 step-2 set-up): at one and the same sigma, free-streaming
     emitters give the scalar sum; colliding emitters give the squared vector sum, like cold
     matter.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov


def exact_step(tau, gamma, kc, sigma, dt):
    """Exact one-step map for the linear SDE  dv = -v/tau dt + sigma sqrt(2/tau) dW,
    dphi = kc v dt - gamma phi dt.

    gamma > 0: both modes decay, so the step noise is Q = P - F P F^T with P the stationary
    covariance.  gamma = 0 (no locking): closed forms for the integrated Ornstein-Uhlenbeck
    process (Gillespie 1996, Phys. Rev. E 54, 2084), with phi = kc x."""
    if gamma > 0:
        A = np.array([[-1.0 / tau, 0.0], [kc, -gamma]])
        BB = np.array([[2.0 * sigma ** 2 / tau, 0.0], [0.0, 0.0]])
        P = solve_continuous_lyapunov(A, -BB)
        F = expm(A * dt)
        Q = P - F @ P @ F.T
    else:
        e1, e2 = np.exp(-dt / tau), np.exp(-2 * dt / tau)
        F = np.array([[e1, 0.0], [kc * tau * (1 - e1), 1.0]])
        vv = sigma ** 2 * (1 - e2)
        xx = sigma ** 2 * tau ** 2 * (2 * dt / tau - 3 + 4 * e1 - e2)
        xv = sigma ** 2 * tau * (1 - e1) ** 2
        Q = np.array([[vv, kc * xv], [kc * xv, kc ** 2 * xx]])
    Q = 0.5 * (Q + Q.T)
    w, V = np.linalg.eigh(Q)
    return F, V * np.sqrt(np.clip(w, 0, None))


def simulate(tau, gamma, n, steps, dt, rng, kc=1.0, sigma=1.0, burn=None):
    F, Lq = exact_step(tau, gamma, kc, sigma, dt)
    v = rng.normal(0, sigma, n); phi = np.zeros(n)
    X = np.vstack([v, phi])
    out = np.empty((steps, n))
    burn = burn if burn is not None else 0
    for s in range(burn + steps):
        X = F @ X + Lq @ rng.normal(size=(2, n))
        if s >= burn:
            out[s - burn] = X[1]
    return out


def coherence_time(tau, rng, n=4000, dt=0.02, steps=3000):
    """1/e decay time of |<exp(i(phi(t) - phi(0)))>| with no locking (gamma = 0)."""
    ph = simulate(tau, 0.0, n, steps, dt, rng)
    c = np.abs(np.mean(np.exp(1j * (ph - 0.0)), axis=1))
    i = np.flatnonzero(c < np.exp(-1))
    return float((i[0] + 1) * dt) if len(i) else float('inf')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    out = ap.parse_args().output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True)
    rng = np.random.default_rng(20260923)

    # 1. spectral narrowing: coherence time of one emitter's companion phase
    narrowing = []
    for tau in (1e6, 1.0, 0.1, 0.01):
        D = tau  # sigma^2 tau with sigma = 1
        steps = 3000 if tau >= 0.1 else 6000
        dt = 0.02 if tau >= 0.1 else 0.05
        t_num = coherence_time(tau, rng, dt=dt, steps=steps)
        # theory: ballistic Gaussian decay exp(-t^2/2) -> 1/e at sqrt(2);
        # diffusive exp(-D t) for t >> tau (exact OU: exp(-D (t - tau (1 - e^{-t/tau}))))
        if tau > 100:
            t_th = np.sqrt(2.0)
        else:
            from scipy.optimize import brentq
            t_th = brentq(lambda t: D * (t - tau * (1 - np.exp(-t / tau))) - 1.0, 1e-9, 1e6)
        narrowing.append(dict(tau_c=tau, kc_times_free_path=tau, coherence_time_numeric=t_num,
                              coherence_time_theory=float(t_th)))

    # 2. steady-state phase variance with locking
    gamma = 0.05
    locking = []
    for tau in (1e6, 200.0, 1.0, 0.1, 0.01, 0.001):
        ph = simulate(tau, gamma, n=2000, steps=4000, dt=2.0, rng=rng, burn=400)
        var_num = float(np.var(ph))
        var_th = tau / (gamma * (1 + gamma * tau))
        locking.append(dict(tau_c=tau, var_numeric=var_num, var_theory=float(var_th),
                            suppression_vs_free_streaming=float(gamma * tau / (1 + gamma * tau)),
                            incoherent_fraction_theory=float(1 - np.exp(-var_th))))

    # 3. the summed field, as in round-2 step 2, at a single sigma
    n = 200
    amps = rng.uniform(0.5, 1.5, n)
    dirs = rng.normal(size=(n, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    vec2 = float(np.sum((amps[:, None] * dirs).sum(0) ** 2)); scal = float(np.sum(amps ** 2))
    field = []
    for label, tau, sig in (('cold (at rest)', 1e6, 0.0), ('hot, free-streaming (orbits)', 200.0, 1.0),
                            ('hot, colliding, k_c l = 0.1', 0.1, 1.0), ('hot, colliding, k_c l = 0.01', 0.01, 1.0),
                            ('hot, colliding, k_c l = 0.001', 0.001, 1.0)):
        vals = []
        for seed in range(6):
            r2 = np.random.default_rng(1000 + seed)
            if sig == 0:
                ph = np.zeros((10, n))
            else:
                ph = simulate(tau, gamma, n, steps=5000, dt=2.0, rng=r2, burn=400, sigma=sig)
            E = (amps[None, :, None] * dirs[None, :, :] * np.exp(1j * ph)[:, :, None]).sum(1)
            vals.append(float(np.mean(np.sum(np.abs(E) ** 2, axis=1))))
        I = float(np.mean(vals))
        var_th = 0.0 if sig == 0 else tau / (gamma * (1 + gamma * tau))
        c2 = np.exp(-var_th)
        field.append(dict(emitters=label, kc_times_free_path=None if sig == 0 else tau,
                          intensity=I, intensity_over_scalar_sum=I / scal,
                          intensity_over_squared_vector_sum=I / vec2,
                          predicted=float(c2 * vec2 + (1 - c2) * scal),
                          seeds_spread=float(np.std(vals) / np.sqrt(len(vals)))))

    payload = dict(experiment='Dicke narrowing of the companion phase (round 3)',
                   units='k_c = 1, sigma = 1 (so the Doppler rate k_c sigma = 1); locking rate Gamma = 0.05',
                   narrowing=narrowing, locking=locking,
                   field=dict(scalar_sum=scal, squared_vector_sum=vec2, runs=field),
                   heat_weight_mapping='free streaming: incoherent fraction ~ k_c^2 sigma^2 / Gamma^2 = 3 sigma^2/u^2 with u = sqrt(3) Gamma / k_c; '
                                       'colliding: times Gamma tau_c / (1 + Gamma tau_c)')
    (out / 'dicke.json').write_text(json.dumps(payload, indent=2) + '\n')

    print('1. coherence time of the companion phase (units 1/(k_c sigma)); free streaming ~ 1.41')
    for r in narrowing:
        print(f"   k_c x free path {r['kc_times_free_path']:>9g}: numeric {r['coherence_time_numeric']:9.2f}  theory {r['coherence_time_theory']:9.2f}")
    print('2. steady-state phase variance with locking Gamma = 0.05')
    for r in locking:
        print(f"   tau_c {r['tau_c']:>9g}: numeric {r['var_numeric']:10.4f}  theory {r['var_theory']:10.4f}  suppression {r['suppression_vs_free_streaming']:.4f}")
    print(f'3. summed field: scalar sum {scal:.2f}, squared vector sum {vec2:.2f}')
    for r in field:
        print(f"   {r['emitters']:32s} I = {r['intensity']:8.2f}  (predicted {r['predicted']:8.2f})  I/scalar {r['intensity_over_scalar_sum']:.3f}  I/vector {r['intensity_over_squared_vector_sum']:.3f}")


if __name__ == '__main__':
    main()
