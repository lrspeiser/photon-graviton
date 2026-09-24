"""Round 15: motion-opened radiation from a finite internal store, driving test bodies through the local force.

The independent calculation "Motion-opened radiation from ordinary matter's internal oscillations" (reproduced here
bit for bit, independent-r15/) gives each piece of matter a quiet (dark) internal oscillation D, holding its internal
energy, and three radiating (bright) ones B_m. Random relative motion shifts their frequencies apart, delta_m = chi w_m,
which mixes D into B linearly:
    dD/dt   = -gamma_0 D - i sum_m delta_m B_m
    dB_m/dt = -gamma   B_m - i delta_m D
so the extra radiation is 3 chi^2 sigma^2 / (gamma_0 gamma) times the quiet leak: the law's heat weight k = 3 sigma^2/u^2
with u^2 = gamma_0 gamma / chi^2, from linear mixing followed by quadratic energy. Collisions (delta redrawn at rate
nu) leave gamma / (gamma + nu) of it. Their README asks for the next step, done here: "use emitted waves from these
evolving internal modes to drive test emitters, derive forces from a common local interaction, and measure force,
phase, energy and momentum separately. Do not convert sqrt(radiated power) into a force."

Set-up (companion wavelength 1, wave speed 1, k = omega = 2 pi; the internal rates in units of gamma = 1; envelopes
quasi-static, as in code/coherent_force_v13.py):
  * Sources: N pieces of matter at fixed points, uniform in a ball of radius Rb, each with its own D, B (the equations
    above, integrated exactly as written, RK4), a random phase of D and its own detunings: at rest (delta = 0), moving
    freely (delta fixed, Gaussian with rms q gamma per component) or colliding (delta an Ornstein-Uhlenbeck process
    with the same rms, redrawn at rate nu). Rigid rotation gives delta = 0 by the coupling's construction (it depends
    on the rate of change of the distances between particles), so it is the rest case here; code/wave_dark_v15.py
    derives that from the wave itself. Optionally the pieces also move in bulk (bulk=True), freely inside the ball
    with the same q in Doppler units (k sigma = q gamma), which re-times their waves at the test bodies.
  * Each piece radiates D's quiet leak as a monopole, s = c0 sqrt(2 gamma_0) D, and each B_m as a dipole along axis m,
    p_m = c1 sqrt(2 gamma) B_m, with c0^2 = 8 pi/(omega k), c1^2 = 24 pi/(omega k^3), so the power each radiates into
    the field equals the energy its modes lose (checked in validate()).
  * Test bodies: P points on a sphere of radius rp. Each carries its own oscillation a_p and feels only the local wave:
        force F_p = (1/2) Re(conj(a_p) grad E(x_p))            (the time average of q grad psi, as in rounds 10-14)
        power it feeds into the wave  P_p = (omega/2) Im(conj(a_p) E(x_p))
    Three kinds of test body share every wave:
      'pll'       round 10/14's emitter: fixed amplitude, phase locked a quarter cycle ahead of the local wave at rate
                  Gamma: d theta/dt = Gamma sin(arg E - pi/2 - theta).
      'bloch'     a self-sustained, inverted emitter (code/receivers_v15.py): a pumped ensemble whose coherence s is kept
                  up by its own collective emission; the local wave drives it only through its Bloch equations,
                     ds/dt = -g_perp s + (Gc/2) w s - (i/2) b E w,   dw/dt = W(1 - w) - g_par(1 + w) - 2 Gc |s|^2
                                                                              + i b (E conj(s) - conj(E) s),
                  and nothing sets its timing: the quarter-cycle lead, if it appears, comes out of these equations.
      'amplifier' an inverted emitter below its own threshold (linear stimulated emission): a_p = -i alpha E.
    pull_p = -F_p . rhat_p (toward the cloud: positive). Nothing about the law's formula enters.
Reported per run: the mean pull of each kind, the mean sin(arg E - arg a) (+1: a quarter cycle ahead, feeding), the mean
power fed by each kind and pull/(fed/c) (momentum against energy), the pieces' radiated power (from their modes) and
stored energy, the mean |E|^2 at the test bodies, and the flux of the pieces' wave through the test-body sphere
(from the fields) against the power their modes lose.

    python code/reservoir_force_v15.py --output run-reservoir-force-v15/reservoir_force_v15.json [--processes 4] [--quick]
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')                  # one thread per run; runs go in parallel
import argparse, json, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

K = 2 * np.pi
OMEGA = K                                   # wave speed 1
C0 = np.sqrt(8 * np.pi / (OMEGA * K))       # monopole: radiated power = 2 gamma_0 |D|^2
C1 = np.sqrt(24 * np.pi / (OMEGA * K ** 3))  # dipole:   radiated power = 2 gamma |B_m|^2
BLOCH = dict(Gc=1.0, W=0.4, g_par=0.01)      # the inverted receiver's constants (code/receivers_v15.py)


def ball(n, R, rng):
    p = rng.uniform(-1, 1, (8 * n, 3)); p = p[np.linalg.norm(p, axis=1) < 1][:n]
    return p * R


def sphere_points(n, r):
    i = np.arange(n) + 0.5; phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return r * np.vstack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]).T


def propagators(xp, x):
    """Field and gradient maps from the pieces' monopoles s (N) and dipoles p (N, 3) to points xp (P, 3):
    E = Mm @ s + Md @ p.ravel(),  grad E[:, c] = Gm[:, c] @ s + Gd[:, c] @ p.ravel()."""
    Dv = xp[:, None, :] - x[None, :, :]; R = np.linalg.norm(Dv, axis=2); Rh = Dv / R[:, :, None]
    G = np.exp(1j * K * R) / (4 * np.pi * R); G1 = (1j * K - 1 / R) * G; G2 = (1 / R ** 2 + (1j * K - 1 / R) ** 2) * G
    P, N = R.shape
    Mm = G
    Md = (-Rh * G1[:, :, None]).reshape(P, 3 * N)                       # a dipole p at x: field p . grad_x G = -(p.Rh) G'
    Gm = np.transpose(G1[:, :, None] * Rh, (0, 2, 1))                   # (P, 3, N)
    eye = np.eye(3)
    RR = Rh[:, :, :, None] * Rh[:, :, None, :]                          # (P, N, 3, 3)
    T = -(G1 / R)[:, :, None, None] * (eye[None, None] - RR) - G2[:, :, None, None] * RR   # d/dxp_c of -(p.Rh) G'
    Gd = np.transpose(T, (0, 3, 1, 2)).reshape(P, 3, 3 * N)             # [p, c, (j, m)]
    return Mm, Md, Gm, Gd


def modes_rhs(D, B, dl, g0, g):
    return -g0 * D - 1j * np.einsum('jm,jm->j', dl, B), -g * B - 1j * dl * D[:, None]


def simulate(cfg):
    """One run. cfg: N, Rb, rp, P, gamma0, kind ('rest' | 'free' | 'collisional'), q, nu, Gamma (the pll's locking rate),
    T, burn, dt, seed; optional bulk (bool), alpha (amplifier gain), b (bloch coupling), flux_points, stop_at (the motion
    stops at this time: detunings and bulk speeds set to zero)."""
    N, Rb, P = cfg['N'], cfg['Rb'], cfg['P']
    g0, g = cfg['gamma0'], 1.0
    kind, q, nu = cfg['kind'], cfg['q'], cfg.get('nu', 0.0)
    Gam, dt = cfg['Gamma'], cfg['dt']
    alpha, bcoup = cfg.get('alpha', 1.0), cfg.get('b', 1.0)
    rng = np.random.default_rng(cfg['seed'])
    x = ball(N, Rb, rng)
    D = np.exp(1j * rng.uniform(0, 2 * np.pi, N)); B = np.zeros((N, 3), complex)
    dl = np.zeros((N, 3))
    if kind in ('free', 'collisional'):
        dl = q * g * rng.normal(size=(N, 3))
    bulk = cfg.get('bulk', False)
    vb = rng.normal(0, q * g / K, (N, 3)) if bulk else None
    xp = sphere_points(P, cfg['rp']); rhat = xp / cfg['rp']
    Mm, Md, Gm, Gd = propagators(xp, x)
    nf = cfg.get('flux_points', 400); xf = sphere_points(nf, cfg['rp'])
    Fm, Fd, FGm, FGd = propagators(xf, x)
    # the test bodies
    th = rng.uniform(0, 2 * np.pi, P)                                  # pll phases
    w0 = 2 * (BLOCH['W'] + BLOCH['g_par']) / 2 / BLOCH['Gc']
    s_b = np.sqrt((BLOCH['W'] * (1 - w0) - BLOCH['g_par'] * (1 + w0)) / (2 * BLOCH['Gc'])) * np.exp(1j * rng.uniform(0, 2 * np.pi, P))
    w_b = np.full(P, w0)
    gperp = (BLOCH['W'] + BLOCH['g_par']) / 2
    rho = np.exp(-nu * dt); kick = q * g * np.sqrt(max(0.0, 1 - rho ** 2))
    n = int(round(cfg['T'] / dt)); nb = int(round(cfg['burn'] / dt))
    acc = {k: 0.0 for k in ('pull_pll', 'pull_bloch', 'pull_amp', 'fed_pll', 'fed_bloch', 'fed_amp', 'lead_pll', 'lead_bloch',
                            'lead_amp', 'I', 'radiated', 'stored', 'bright_frac', 'absE')}
    pull_p = np.zeros(P); flux = []; cnt = 0
    E0 = float(np.sum(np.abs(D) ** 2))

    def field(D, B, M1, M2):
        return M1 @ (C0 * np.sqrt(2 * g0) * D) + M2 @ (C1 * np.sqrt(2 * g) * B).ravel()

    def grad(D, B):
        sm = C0 * np.sqrt(2 * g0) * D; pd = (C1 * np.sqrt(2 * g) * B).ravel()
        return np.einsum('pcj,j->pc', Gm, sm) + np.einsum('pcj,j->pc', Gd, pd)

    def bloch_rhs(s, w, E):
        Om = bcoup * E
        ds = (-gperp + 0.5 * BLOCH['Gc'] * w) * s - 0.5j * Om * w
        dw = BLOCH['W'] * (1 - w) - BLOCH['g_par'] * (1 + w) - 2 * BLOCH['Gc'] * np.abs(s) ** 2 + np.real(1j * (Om * np.conj(s) - np.conj(Om) * s))
        return ds, dw

    t_stop = cfg.get('stop_at')
    for st in range(n):
        if t_stop is not None and st * dt >= t_stop and kind != 'rest':
            dl = np.zeros((N, 3)); kind = 'stopped'
            if bulk: vb = np.zeros((N, 3))
        # the pieces' internal modes (RK4, detunings frozen over the step)
        k1 = modes_rhs(D, B, dl, g0, g); k2 = modes_rhs(D + 0.5 * dt * k1[0], B + 0.5 * dt * k1[1], dl, g0, g)
        k3 = modes_rhs(D + 0.5 * dt * k2[0], B + 0.5 * dt * k2[1], dl, g0, g); k4 = modes_rhs(D + dt * k3[0], B + dt * k3[1], dl, g0, g)
        D = D + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6; B = B + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
        if kind == 'collisional':
            dl = rho * dl + kick * rng.normal(size=(N, 3))
        if bulk:
            x = x + vb * dt
            r = np.linalg.norm(x, axis=1); out = r > Rb
            if out.any():
                nrm = x[out] / r[out, None]
                vb[out] -= 2 * (vb[out] * nrm).sum(1)[:, None] * nrm
                x[out] = nrm * (2 * Rb - r[out])[:, None]
            if kind == 'collisional':
                hit = rng.random(N) < nu * dt
                vb[hit] = rng.normal(0, q * g / K, (int(hit.sum()), 3))
            Mm, Md, Gm, Gd = propagators(xp, x)
        E = field(D, B, Mm, Md)
        # the test bodies respond to the local wave only
        th = th + dt * Gam * np.sin(np.angle(E) - np.pi / 2 - th)
        k1 = bloch_rhs(s_b, w_b, E); k2 = bloch_rhs(s_b + 0.5 * dt * k1[0], w_b + 0.5 * dt * k1[1], E)
        k3 = bloch_rhs(s_b + 0.5 * dt * k2[0], w_b + 0.5 * dt * k2[1], E); k4 = bloch_rhs(s_b + dt * k3[0], w_b + dt * k3[1], E)
        s_b = s_b + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6; w_b = w_b + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
        if st >= nb:
            gE = grad(D, B)
            for name, a in (('pll', np.exp(1j * th)), ('bloch', s_b), ('amp', -1j * alpha * E)):
                F = 0.5 * np.real(np.conj(a)[:, None] * gE)
                pl = -(F * rhat).sum(1)
                acc['pull_' + name] += pl.mean()
                acc['fed_' + name] += (0.5 * OMEGA * np.imag(np.conj(a) * E)).mean()
                acc['lead_' + name] += np.mean(np.sin(np.angle(E) - np.angle(a)))
                if name == 'pll': pull_p += pl
            acc['I'] += np.mean(np.abs(E) ** 2); acc['absE'] += np.mean(np.abs(E))
            rad = 2 * g0 * np.abs(D) ** 2 + 2 * g * (np.abs(B) ** 2).sum(1)
            acc['radiated'] += rad.sum(); acc['stored'] += (np.abs(D) ** 2 + (np.abs(B) ** 2).sum(1)).sum()
            acc['bright_frac'] += (2 * g * (np.abs(B) ** 2).sum()) / rad.sum()
            if (st - nb) % max(1, (n - nb) // 8) == 0 and not bulk:
                Ef = field(D, B, Fm, Fd)
                sm = C0 * np.sqrt(2 * g0) * D; pd = (C1 * np.sqrt(2 * g) * B).ravel()
                gEf = np.einsum('pcj,j->pc', FGm, sm) + np.einsum('pcj,j->pc', FGd, pd)
                dEr = (gEf * (xf / cfg['rp'])).sum(1)
                phi = 0.5 * OMEGA * np.mean(np.imag(np.conj(Ef) * dEr)) * 4 * np.pi * cfg['rp'] ** 2
                flux.append((float(phi), float(rad.sum())))
            cnt += 1
    res = {k: float(v / cnt) for k, v in acc.items()}
    res.update(cfg=cfg, pull_pll_se=float((pull_p / cnt).std(ddof=1) / np.sqrt(P)),
               momentum_over_energy_pll=res['pull_pll'] / res['fed_pll'] if res['fed_pll'] else None,
               momentum_over_energy_bloch=res['pull_bloch'] / res['fed_bloch'] if res['fed_bloch'] else None,
               momentum_over_energy_amp=res['pull_amp'] / res['fed_amp'] if res['fed_amp'] else None,
               flux_over_mode_loss=[f / m for f, m in flux], store_start=E0, store_end=float(np.sum(np.abs(D) ** 2 + (np.abs(B) ** 2).sum(1))),
               detuning_rms=float(np.sqrt(np.mean(dl ** 2))))
    return res


def validate():
    """(1) One piece: the flux of its monopole and of each dipole through a sphere equals the power its modes lose.
    (2) A pll test body a quarter cycle ahead in one piece's wave: pull = fed / c; (3) the bloch test body's timing."""
    out = {}
    x = np.zeros((1, 3))
    xf = sphere_points(4000, 7.0)
    Fm, Fd, FGm, FGd = propagators(xf, x)
    for lab, D, B, g0 in (('monopole', np.array([1.0 + 0j]), np.zeros((1, 3), complex), 0.3),
                          ('dipole_x', np.array([0j]), np.array([[1.0, 0, 0]], complex), 0.3),
                          ('dipoles_xyz', np.array([0j]), np.array([[0.3, 1j, -0.5]], complex), 0.3)):
        s = C0 * np.sqrt(2 * g0) * D; p = (C1 * np.sqrt(2) * B).ravel()
        E = Fm @ s + Fd @ p
        gE = np.einsum('pcj,j->pc', FGm, s) + np.einsum('pcj,j->pc', FGd, p)
        flux = 0.5 * OMEGA * np.mean(np.imag(np.conj(E) * (gE * xf / 7.0).sum(1))) * 4 * np.pi * 49
        loss = float(2 * g0 * abs(D[0]) ** 2 + 2 * (np.abs(B) ** 2).sum())
        out[lab] = dict(flux=float(flux), mode_loss=loss)
    xp = np.array([[6.0, 0, 0]]); Mm, Md, Gm, Gd = propagators(xp, x)
    s = C0 * np.sqrt(2 * 0.3) * np.array([1.0 + 0j]); p = (C1 * np.sqrt(2) * np.array([[0.2, 0.4j, 0]])).ravel()
    E = Mm @ s + Md @ p; gE = np.einsum('pcj,j->pc', Gm, s) + np.einsum('pcj,j->pc', Gd, p)
    a = np.exp(1j * (np.angle(E) - np.pi / 2))
    pull = float(-0.5 * np.real(np.conj(a) * gE[:, 0])[0]); fed = float((0.5 * OMEGA * np.imag(np.conj(a) * E))[0])
    out['pll_quarter_ahead'] = dict(pull=pull, fed_over_c=fed, ratio=pull / fed)
    return out


def configs(quick=False):
    base = dict(N=100, Rb=1.0, rp=6.0, P=64, gamma0=1e-6, Gamma=0.1, b=50.0, T=600.0 if quick else 1200.0,
                burn=300.0 if quick else 600.0, dt=0.05, nu=0.0)
    seeds = (1,) if quick else (1, 2, 3, 4, 5)
    runs = []
    for sd in seeds:
        runs.append(dict(base, tag='rest', kind='rest', q=0.0, seed=sd))
        for q in (0.00025, 0.0005, 0.001, 0.002, 0.004, 0.008, 0.016):
            runs.append(dict(base, tag='free', kind='free', q=q, seed=sd))
        for nu in (0.3, 1.0, 3.0, 10.0, 30.0):
            runs.append(dict(base, tag='collisional', kind='collisional', q=0.002, nu=nu, seed=sd))
        for nu in (1.0, 10.0):
            runs.append(dict(base, tag='collisional_fastlock', kind='collisional', q=0.002, nu=nu, seed=sd, Gamma=10.0, dt=0.01))
        runs.append(dict(base, tag='free_fastlock', kind='free', q=0.002, seed=sd, Gamma=10.0, dt=0.01))
        runs.append(dict(base, tag='free_bulk', kind='free', q=0.002, bulk=True, seed=sd))
        runs.append(dict(base, tag='stop', kind='free', q=0.004, stop_at=900.0, burn=950.0, seed=sd))
    return runs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4); ap.add_argument('--quick', action='store_true')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    val = validate()
    print('validation:', json.dumps(val), flush=True)
    runs = configs(args.quick)
    runs.sort(key=lambda c: -(c['dt'] < 0.05) - c.get('bulk', False))
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(simulate, runs):
            res.append(r); c = r['cfg']
            print(f"[{time.monotonic() - t0:5.0f} s] {c['tag']:21s} q {c['q']:.5f} nu {c['nu']:5.1f} seed {c['seed']}: "
                  f"pull pll {r['pull_pll']:.4e} (lead {r['lead_pll']:.3f}) bloch {r['pull_bloch']:.4e} (lead {r['lead_bloch']:.3f}) "
                  f"amp {r['pull_amp']:.4e}; radiated {r['radiated']:.3e}, bright {r['bright_frac']:.3f}", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 15: motion-opened radiation from a finite internal store, through the local force',
                                           validation=val, runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
