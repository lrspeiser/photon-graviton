"""Round 13: three checks on the joint experiment (code/coherent_force_v13.py), with the same cloud and test bodies.

1. offsets: the locking offset delta over a full turn, at rest. The claim tested: in a locked steady state every
   emitter, test bodies included, sits at the same phase to its own local wave, so the test bodies feed the wave
   (and are pulled) exactly when the cloud's sources feed one another, which is when the cloud outshines
   independent emitters.
2. doubling: the absorbing offset (delta = +pi/2, a dark cloud at rest), the speed spread doubled step by step,
   free and colliding (the same collision rule as the main run), in the main cloud and in one twice as small
   (8 times denser). How does the power released by random motion grow with sigma?
3. amplitude: emitters whose strength can change as well as their timing (Stuart-Landau oscillators),
       dA_j/dt = (mu + i dw_j - (1 + i c)|A_j|^2) A_j + K e^{i beta} E_j,
   with beta = -pi/2 + arctan(c), so that a lone emitter in a weak wave still locks a quarter cycle ahead (it
   feeds and is pulled), at three coupling strengths K and four shears c. Does a cold cloud of them go dark?

    python code/joint_checks_v13.py --output run-coherent-force-v13/joint_checks_v13.json [--processes 4]
"""
from __future__ import annotations
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_force_v13 import K0, ball, sphere_points, simulate   # noqa: E402

BASE = dict(N=100, Rb=1.0, Gamma=0.1, T=900.0, burn=350.0, rp=6.0, P=32, w=0.01)


def simulate_amplitude(cfg):
    """One cloud at rest of Stuart-Landau emitters (mu = 1), test bodies of the same kind."""
    N, Rb, mu, K, c = cfg['N'], cfg['Rb'], 1.0, cfg['K'], cfg['c']
    beta = -np.pi / 2 + np.arctan(c)
    rng = np.random.default_rng(cfg['seed'])
    x = ball(N, Rb, rng)
    A = np.exp(1j * rng.uniform(0, 2 * np.pi, N)) * (1 + 0.1 * rng.normal(size=N))
    dw = rng.uniform(-cfg['w'], cfg['w'], N)
    xp = sphere_points(cfg['P'], cfg['rp']); rhat = xp / cfg['rp']
    Ap = np.exp(1j * rng.uniform(0, 2 * np.pi, cfg['P']))
    eK = K * np.exp(1j * beta)
    D = x[:, None, :] - x[None, :, :]; R = np.linalg.norm(D, axis=2); np.fill_diagonal(R, 1.0)
    G = np.exp(1j * K0 * R) / (4 * np.pi * R); np.fill_diagonal(G, 0.0)
    Rz = R.copy(); np.fill_diagonal(Rz, 0.0); S = np.sinc(K0 * Rz / np.pi)
    Dp = xp[:, None, :] - x[None, :, :]; Rp = np.linalg.norm(Dp, axis=2)
    Gp = np.exp(1j * K0 * Rp) / (4 * np.pi * Rp)
    dGp = (Gp * (1j * K0 - 1 / Rp))[:, :, None] * (Dp / Rp[:, :, None])
    inc = (np.abs(Gp) ** 2).sum(1)
    dt = cfg['dt']; n = int(cfg['T'] / dt); nb = int(cfg['burn'] / dt)
    F_acc = np.zeros(cfg['P']); rad = amp2 = 0.0; ph = np.zeros(cfg['P'], complex); cnt = 0
    for s in range(n):
        E = G @ A; Ep = Gp @ A
        if s >= nb:
            gE = (dGp * A[None, :, None]).sum(1)
            F_acc += -(0.5 * np.real(np.conj(Ap)[:, None] * gE) * rhat).sum(1)
            rad += np.real(np.conj(A) @ (S @ A)) / N; amp2 += np.mean(np.abs(A) ** 2)
            ph += np.exp(1j * (np.angle(Ap) - np.angle(Ep))); cnt += 1
        A = A + dt * ((mu + 1j * dw - (1 + 1j * c) * np.abs(A) ** 2) * A + eK * E)
        Ap = Ap + dt * ((mu - (1 + 1j * c) * np.abs(Ap) ** 2) * Ap + eK * Ep)
    pull = F_acc / cnt
    indep = np.mean(0.5 * K0 * (np.sqrt(np.pi) / 2) * np.sqrt(inc))
    ph /= cnt
    return dict(cfg=cfg, beta_over_pi=beta / np.pi, radiated_over_independent=float(rad / cnt), mean_A2=float(amp2 / cnt),
                pull_over_independent=float(pull.mean() / indep), test_body_locking=float(np.mean(np.abs(ph))),
                test_body_offset_over_pi=float(np.angle(np.mean(ph)) / np.pi))


def run(cfg):
    return simulate_amplitude(cfg) if cfg['part'] == 'amplitude' else simulate(cfg)


def configs():
    runs = []
    for d in np.arange(-1.0, 1.0, 0.125):
        for sd in (1, 2):
            runs.append(dict(BASE, part='offsets', tag='offsets', seed=sd, kind='rest', sigma=0.0, nu=0.0, delta=float(d * np.pi)))
    for Rb in (1.0, 0.5):
        for sd in (1, 2):
            runs.append(dict(BASE, part='doubling', tag='doubling', Rb=Rb, seed=sd, kind='rest', sigma=0.0, nu=0.0, delta=np.pi / 2, ratio=0.0))
            for q in (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
                sig = q * BASE['Gamma'] / K0
                for kind, nu in (('free', 0.0), ('collisional', 20 * q * BASE['Gamma'])):
                    runs.append(dict(BASE, part='doubling', tag='doubling', Rb=Rb, seed=sd, kind=kind, sigma=sig, nu=nu, delta=np.pi / 2, ratio=q))
    for c in (0.0, -1.0, -3.0, 3.0):
        for K in (0.1, 0.3, 1.0):
            runs.append(dict(N=100, Rb=1.0, rp=6.0, P=32, w=0.01, part='amplitude', seed=1, c=c, K=K, T=300.0, burn=150.0, dt=0.02))
    return runs


def cost(c):
    if c['part'] == 'amplitude': return c['T'] / c['dt']
    return c['T'] / min([0.1, 0.05 / c['Gamma']] + ([0.2 / c['nu']] if c['nu'] > 0 else []) + ([0.1 / (K0 * c['sigma'])] if c['sigma'] > 0 else []))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4)
    args = ap.parse_args(); t0 = time.monotonic()
    runs = sorted(configs(), key=lambda c: -cost(c))
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(run, runs):
            res.append(r); c = r['cfg']
            if c['part'] == 'amplitude':
                print(f"[{time.monotonic() - t0:5.0f} s] amplitude c {c['c']:+.0f} K {c['K']:.1f}: radiated/indep {r['radiated_over_independent']:.2f},"
                      f" pull/indep {r['pull_over_independent']:+.3f}, test-body locking {r['test_body_locking']:.2f}", flush=True)
            else:
                print(f"[{time.monotonic() - t0:5.0f} s] {c['part']} Rb {c['Rb']} {c['kind']:11s} ratio {c.get('ratio', 0):5.3f} delta {c['delta'] / np.pi:+.3f}pi"
                      f" seed {c['seed']}: radiated/indep {r['radiated_over_independent']:.3f}, pull/indep {r['pull'] / r['independent_pull']:+.3f},"
                      f" feeding: sources {r['feeding_sources']:+.3f}, test bodies {r['feeding_test_bodies']:+.3f}", flush=True)
    for r in res:
        if 'delta' in r['cfg']: r['cfg']['delta_over_pi'] = r['cfg']['delta'] / np.pi
    args.output.write_text(json.dumps(dict(experiment='round 13: checks on the joint experiment (offsets, doubling, free amplitudes)',
                                           runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
