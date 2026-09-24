"""Round 13: the local-force model and the coherence model in one experiment.

Round 10 derived the pull from one local coupling: a self-sustained emitter of the companion wave that runs a quarter
cycle ahead of the wave passing it is pulled toward the wave's source, with a force (q0 k / 2) x the local amplitude
(code/first_principles_v10.py). Round 3's Dicke toy showed the other half: random free motion scrambles an emitter's
phase, and frequent collisions keep it in step (code/dicke_toy.py). The heat term of the law says the two combine:
freely moving random matter pulls harder (k = 3 sigma^2 / u^2, so in the heat-dominated regime the extra pull grows as
sigma), and colliding matter does not.

This experiment joins them with one set of rules and changes only how the sources move.
  * Matter: N point emitters, uniform in a ball of radius Rb (in companion wavelengths). Each has a fixed amplitude, its
    own pitch (a small spread w of natural frequencies) and a phase that locks, at rate Gamma, to the phase of the wave
    around it plus a fixed offset delta:  d theta_j/dt = dw_j + Gamma sin(arg E_j + delta - theta_j),
    E_j = sum_{l != j} G(R_jl) e^{i theta_l},  G(R) = e^{i k0 R} / (4 pi R)  (the scalar wave's outgoing field).
    delta = -pi/2 is round 10's "a quarter cycle ahead": the offset that feeds the passing wave and is pulled.
  * Test bodies: the same kind of emitter, with the same rule and rate, held at P points on a sphere of radius rp around
    the cloud. Each feels only the local wave: F = (1/2) Re(e^{-i theta_p} grad E(x_p)) (the time average of q grad psi).
    Nothing about the law's formula enters.
  * Motion (the only thing that changes): at rest; ordered (rigid rotation, line-of-sight spread sigma in the equator);
    random and free (Maxwellian, 1D spread sigma, bouncing off the ball's surface); random with collisions (the same
    Maxwellian, each velocity redrawn at rate nu, so the free path is short compared with a wavelength).
The coupling is quasi-static (retardation of the envelope neglected: Gamma * (rp + Rb) / c < 1 rad); the carrier phase
k0 R, which carries every Doppler effect, is exact. Units: companion wavelength 1, wave speed 1, period 1.

Reported per run: the mean pull on the test bodies (toward the cloud), the cloud's radiated power relative to
independent emitters (sum_jl Re(a_j a_l*) sinc(k0 R_jl) / N), the locking order <cos(theta - target)>, and the mean
intensity at the test bodies relative to independent emitters. feeding_sources and feeding_test_bodies, the mean
sin(arg E - theta) (positive: feeding the local wave), and field_sources and field_test_bodies, the mean |E|, were
added for code/joint_checks_v13.py and code/strength_offset_v13.py; the main run's file predates them (the dynamics
are unchanged).

    python code/coherent_force_v13.py --output run-coherent-force-v13/coherent_force_v13.json [--processes 4] [--quick]
"""
from __future__ import annotations
import argparse, json, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

K0 = 2 * np.pi


def ball(n, R, rng):
    p = rng.uniform(-1, 1, (8 * n, 3)); p = p[np.linalg.norm(p, axis=1) < 1][:n]
    return p * R


def sphere_points(n, r):
    i = np.arange(n) + 0.5; phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return r * np.vstack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]).T


def simulate(cfg):
    """One run. cfg: N, Rb, kind, sigma, nu, w, Gamma, delta, T, burn, rp, P, seed.

    Optional (code/strength_offset_v13.py): E_s. Then the offset depends on the strength of the local wave,
    delta(|E|) = delta + (delta_strong - delta) x^n / (1 + x^n), x = |E| / E_s (delta_strong = +pi/2, n = 4 by
    default), for sources and test bodies alike: delta in weak waves, delta_strong in strong ones.

    Optional (code/singer_absorber_v14.py): c_abs = c0 > 0. Then every emitter also absorbs: besides its own
    locked oscillation a_j = e^{i theta_j} it carries a passive part driven by the wave around it, i c0 E_j (a
    quarter cycle behind: it takes energy from the wave), so its total source is s_j = a_j + i c0 E_j with
    E_j = sum_{l != j} G_jl s_l, solved each step as s = (1 - i c0 G)^{-1} a. Test bodies do the same in the cloud's
    wave. c0 <= 4 pi / k0 = 2 keeps the passive part passive (c0 = 1: the largest absorption).
    Optional (round 14): supply = P_s and loss = c_L. Then the offset is set by each emitter's power balance, as a
    synchronous machine's load angle is: it has a fixed supply P_s and loses c_L |E|^2 to the wave around it, so it
    must feed the wave the difference, |E| sin(arg E - theta) = P_s - c_L |E|^2 (units of a unit-strength emitter):
    delta(|E|) = -arcsin(clip(P_s/|E| - c_L |E|, -1, 1)). A quarter cycle ahead (feeding) in weak waves, behind
    (absorbing) in strong ones, with the switch at |E| = sqrt(P_s / c_L); the same for sources and test bodies.
    With gamma_abs as well, the passive part is resonant, with line width gamma_abs: it follows its drive only as
    fast as db_j/dt = -gamma_abs (b_j - i c0 E_j) lets it (implicit steps), so a wave whose phase at j drifts
    faster than gamma_abs (a Doppler shift) is absorbed less. Test bodies' passive parts do the same.
    """
    N, Rb, kind = cfg['N'], cfg['Rb'], cfg['kind']
    sigma, nu, w, Gamma, delta = cfg['sigma'], cfg['nu'], cfg['w'], cfg['Gamma'], cfg['delta']
    Es = cfg.get('E_s'); d_strong = cfg.get('delta_strong', np.pi / 2); n_s = cfg.get('n_s', 4)
    c0 = cfg.get('c_abs', 0.0); gam_a = cfg.get('gamma_abs')
    t_stop = cfg.get('stop_at')                       # optional: the sources stop moving at this time (round 14 check)

    Ps, cL = cfg.get('supply'), cfg.get('loss', 0.0)

    def offset(A):
        if Ps is not None:
            return -np.arcsin(np.clip(Ps / np.maximum(A, 1e-12) - cL * A, -1.0, 1.0))
        if not Es: return delta
        y = (A / Es) ** n_s
        return delta + (d_strong - delta) * y / (1 + y)
    rng = np.random.default_rng(cfg['seed'])
    x = ball(N, Rb, rng); th = rng.uniform(0, 2 * np.pi, N); dw = rng.uniform(-w, w, N)
    v = rng.normal(0, sigma, (N, 3)) if kind in ('free', 'collisional') else np.zeros((N, 3))
    Om = np.sqrt(5) * sigma / Rb if kind == 'ordered' else 0.0         # equatorial line-of-sight spread = sigma
    xp = sphere_points(cfg['P'], cfg['rp']); rhat = xp / cfg['rp']
    polar = np.abs(rhat[:, 2]) > 0.7
    thp = rng.uniform(0, 2 * np.pi, cfg['P'])
    rates = [0.1, 0.05 / max(Gamma, 1e-9)]
    if nu > 0: rates.append(0.2 / nu)
    if sigma > 0: rates.append(0.1 / (K0 * sigma))
    dt = min(rates)
    n = int(cfg['T'] / dt); nb = int(cfg['burn'] / dt)
    F_acc = np.zeros(cfg['P']); I_acc = np.zeros(cfg['P']); inc_acc = np.zeros(cfg['P'])
    lock = rad = feed = feed_p = fld = fld_p = absb = 0.0; cnt = 0
    for s in range(n):
        D = x[:, None, :] - x[None, :, :]; R = np.linalg.norm(D, axis=2); np.fill_diagonal(R, 1.0)
        G = np.exp(1j * K0 * R) / (4 * np.pi * R); np.fill_diagonal(G, 0.0)
        a = np.exp(1j * th)
        if c0 and gam_a:                                                   # resonant passive parts (implicit step)
            if s == 0:
                b = np.linalg.solve(np.eye(N) - 1j * c0 * G, a) - a
            h = dt * gam_a
            b = np.linalg.solve((1 + h) * np.eye(N) - 1j * h * c0 * G, b + 1j * h * c0 * (G @ a))
            src = a + b
        else:
            src = np.linalg.solve(np.eye(N) - 1j * c0 * G, a) if c0 else a     # the sources, with their passive parts
        E = G @ src
        tgt = np.angle(E) + offset(np.abs(E))
        Dp = xp[:, None, :] - x[None, :, :]; Rp = np.linalg.norm(Dp, axis=2)
        Gp = np.exp(1j * K0 * Rp) / (4 * np.pi * Rp)
        Ep = Gp @ src
        if c0 and gam_a:
            bp = (1j * c0 * Ep) if s == 0 else (bp + dt * gam_a * 1j * c0 * Ep) / (1 + dt * gam_a)
        if s >= nb:
            gE = ((Gp * (1j * K0 - 1 / Rp))[:, :, None] * (Dp / Rp[:, :, None]) * src[None, :, None]).sum(1)
            sp = np.exp(1j * thp) + ((bp if gam_a else 1j * c0 * Ep) if c0 else 0.0)   # a test body's own source
            F = 0.5 * np.real(np.conj(sp)[:, None] * gE)
            F_acc += -(F * rhat).sum(1)                                   # toward the cloud: positive
            I_acc += np.abs(Ep) ** 2; inc_acc += (np.abs(Gp) ** 2).sum(1)
            lock += np.mean(np.cos(th - tgt))
            feed += np.mean(np.sin(np.angle(E) - th)); feed_p += np.mean(np.sin(np.angle(Ep) - thp))   # > 0: feeding the local wave
            fld += np.mean(np.abs(E)); fld_p += np.mean(np.abs(Ep))
            if c0: absb += (-np.mean(np.imag(np.conj(src - a) * E)) if gam_a else c0 * np.mean(np.abs(E) ** 2)) / (K0 / (4 * np.pi))   # taken from the wave by the passive parts
            if s % 10 == 0:
                Rz = R.copy(); np.fill_diagonal(Rz, 0.0)
                rad += 10 * np.real(np.conj(src) @ (np.sinc(K0 * Rz / np.pi) @ src)) / N
            cnt += 1
        th = th + dt * (dw + Gamma * np.sin(tgt - th))
        thp = thp + dt * Gamma * np.sin(np.angle(Ep) + offset(np.abs(Ep)) - thp)
        if t_stop is not None and s * dt >= t_stop:
            continue
        if kind == 'ordered':
            c, s_ = np.cos(Om * dt), np.sin(Om * dt)
            x = np.c_[c * x[:, 0] - s_ * x[:, 1], s_ * x[:, 0] + c * x[:, 1], x[:, 2]]
        elif kind in ('free', 'collisional'):
            x = x + v * dt
            r = np.linalg.norm(x, axis=1); out = r > Rb
            if out.any():
                nrm = x[out] / r[out, None]
                v[out] -= 2 * (v[out] * nrm).sum(1)[:, None] * nrm
                x[out] = nrm * (2 * Rb - r[out])[:, None]
            if kind == 'collisional':
                hit = rng.random(N) < nu * dt
                v[hit] = rng.normal(0, sigma, (int(hit.sum()), 3))
    pull = F_acc / cnt
    inc = inc_acc / cnt
    return dict(cfg=cfg, dt=dt, pull=float(pull.mean()), pull_equator=float(pull[~polar].mean()), pull_polar=float(pull[polar].mean()),
                pull_se=float(pull.std(ddof=1) / np.sqrt(len(pull))),
                intensity_over_independent=float(np.mean(I_acc / cnt / inc)), radiated_over_independent=float(rad / cnt),
                locking=float(lock / cnt), feeding_sources=float(feed / cnt), feeding_test_bodies=float(feed_p / cnt),
                field_sources=float(fld / cnt), field_test_bodies=float(fld_p / cnt), absorbed_over_independent=float(absb / cnt),
                independent_pull=float(np.mean(0.5 * K0 * (np.sqrt(np.pi) / 2) * np.sqrt(inc))))


def validate():
    """One static source and one test body locked a quarter cycle ahead: the pull is (k0/2) |E| (round 10)."""
    out = []
    for rp in (0.3, 1.0, 3.0):
        x = np.zeros((1, 3)); xp = np.array([[rp, 0, 0]])
        Dp = xp[:, None, :] - x[None, :, :]; Rp = np.linalg.norm(Dp, axis=2)
        Gp = np.exp(1j * K0 * Rp) / (4 * np.pi * Rp); Ep = Gp @ np.ones(1)
        thp = np.angle(Ep)[0] - np.pi / 2
        gE = ((Gp * (1j * K0 - 1 / Rp))[:, :, None] * (Dp / Rp[:, :, None])).sum(1)
        F = 0.5 * np.real(np.exp(-1j * thp) * gE[0])
        out.append(dict(r=rp, pull=float(-F[0]), expected=float(0.5 * K0 * abs(Ep[0]))))
    return out


def configs(quick=False):
    base = dict(N=100, Rb=1.0, Gamma=0.1, T=500.0 if quick else 900.0, burn=250.0 if quick else 350.0, rp=6.0, P=32)
    seeds = (1,) if quick else (1, 2, 3)
    ratios = (0.3, 1.0, 3.0, 10.0)
    runs = []
    def add(tag, **kw):
        for sd in seeds:
            runs.append(dict(base, tag=tag, seed=sd, **kw))
    # A: round 10's rule (a quarter cycle ahead), emitters able to lock together (pitch spread 0.1 Gamma)
    for tag, delta, w in (('A', -np.pi / 2, 0.01), ('B', -np.pi / 2, 1.0), ('D', np.pi / 2, 0.01)):
        add(tag, kind='rest', sigma=0.0, nu=0.0, delta=delta, w=w)
        rr = ratios if tag != 'D' else (0.1, 0.3, 1.0, 3.0, 10.0)
        for q in rr:
            sig = q * base['Gamma'] / K0
            for kind, nu in (('ordered', 0.0), ('free', 0.0), ('collisional', 20 * q * base['Gamma'])):
                add(tag, kind=kind, sigma=sig, nu=nu, delta=delta, w=w, ratio=q)
    # C: the locking offset, at rest
    for d in (-0.5, -0.25, 0.0, 0.25, 0.5, 1.0):
        runs.append(dict(base, tag='C', seed=1, kind='rest', sigma=0.0, nu=0.0, delta=d * np.pi, w=0.01))
    # E: a larger, denser cloud with round 10's rule (the size check)
    for kind, nu in (('rest', 0.0), ('free', 0.0), ('collisional', 20 * 1.0 * base['Gamma'])):
        runs.append(dict(base, tag='E', seed=1, N=200, Rb=1.26, kind=kind, sigma=0.0 if kind == 'rest' else base['Gamma'] / K0,
                         nu=nu, delta=-np.pi / 2, w=0.01, ratio=1.0))
    return runs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4); ap.add_argument('--quick', action='store_true')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    val = validate()
    print('validation (one source, one test body a quarter cycle ahead):',
          ', '.join(f"r {v['r']}: pull {v['pull']:.5f} vs (k0/2)|E| {v['expected']:.5f}" for v in val), flush=True)
    runs = configs(args.quick)
    runs.sort(key=lambda c: -(c['nu'] * 5 + c.get('ratio', 0)))            # the slow collisional runs first
    with Pool(args.processes) as pool:
        res = []
        for i, r in enumerate(pool.imap_unordered(simulate, runs)):
            res.append(r); c = r['cfg']
            print(f"[{time.monotonic() - t0:5.0f} s] {c['tag']} {c['kind']:11s} ratio {c.get('ratio', 0):5.2f} seed {c['seed']}: pull {r['pull']:.4e}"
                  f" (eq {r['pull_equator']:.3e}, pol {r['pull_polar']:.3e}); radiated/indep {r['radiated_over_independent']:.3f}; lock {r['locking']:.3f}", flush=True)
    for r in res:
        r['cfg']['delta_over_pi'] = r['cfg']['delta'] / np.pi
    args.output.write_text(json.dumps(dict(experiment='round 13: the local-force and coherence models in one experiment',
                                           validation=val, runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
