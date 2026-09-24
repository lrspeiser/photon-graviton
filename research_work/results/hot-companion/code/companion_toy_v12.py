"""Round 12: a dynamical toy of the companion. Which local rule gives one stream, no whirlpools and no lost energy?

Round 11 (README §21.1) derived the ordered companion's rule from three requirements on its steady energy flow J,
fed by matter at a rate q per unit area:
  1. no energy lost:   div J = q (power out of any region = power fed in);
  2. one stream:       n u = |J|, f = |J| / (n u) = 1 (all the energy at a point moves one way, at the full speed u);
  3. no whirlpools:    curl J = 0 (graded as 'whirl': the share of J's spatial variation that is rotation,
                       |curl J| / (|dJx/dx| + |dJx/dy| + |dJy/dx| + |dJy/dy|), against the grid's floor for J_N itself).
1 and 3 fix J uniquely: the Newtonian pattern J_N = (1/2 pi) int q(x') (x - x') / |x - x'|^2 d^2x' in two dimensions
(Gauss). 3 was assumed. This toy builds the companion from local rules for its quanta, on a 2D grid, and grades
each rule on the three requirements:
  free       quanta fly straight from where they are made, in every direction (waves passing through each other);
  scatter    ... and scatter isotropically off each other (sigma = 10 per unit length: diffusion);
  annihilate ... and head-on pairs annihilate (counter-streams are removed);
  align      ... and each quantum turns towards the local mean flow (a Vicsek/BGK rule that keeps energy);
  guided     quanta move along the local Newtonian field line, away from the matter, the way Alfven waves are guided
             along magnetic field lines in a plasma.
The first four are solved as a kinetic equation on 64 directions (explicit upwind, run to a steady state); the guided
rule is one stream advected along the field direction e = J_N / |J_N|. Sources: a uniform disk, two equal blobs, an
unequal pair (4:1) and three blobs, each normalized to unit total power; the box is [-1, 1]^2 with nothing coming in.

    python code/companion_toy_v12.py --output run-companion-toy-v12/companion_toy_v12.json
    python code/companion_toy_v12.py --convergence --output run-companion-toy-v12/convergence_v12.json
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np

U = 1.0


def grid(N):
    h = 2.0 / N
    x = -1 + h * (np.arange(N) + 0.5)
    X, Y = np.meshgrid(x, x, indexing='ij')
    return X, Y, h


def gauss(X, Y, x0, y0, s):
    return np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2 * s ** 2))


def sources(kind, X, Y, h):
    r = np.hypot(X, Y)
    if kind == 'uniform disk':
        q = 0.5 * (1 - np.tanh((r - 0.3) / 0.01))
    elif kind == 'two equal blobs':
        q = gauss(X, Y, -0.3, 0, 0.06) + gauss(X, Y, 0.3, 0, 0.06)
    elif kind == 'unequal pair 4:1':
        q = 4 * gauss(X, Y, -0.2, 0, 0.06) + gauss(X, Y, 0.4, 0, 0.06)
    elif kind == 'three blobs':
        q = sum(gauss(X, Y, 0.3 * np.cos(a), 0.3 * np.sin(a), 0.06) for a in (np.pi / 2, np.pi / 2 + 2 * np.pi / 3, np.pi / 2 + 4 * np.pi / 3))
    else:
        raise ValueError(kind)
    return q / (q.sum() * h * h)


def newtonian_flux(q, h):
    """J_N = (1/2 pi) int q(x') (x - x')/|x - x'|^2 d^2x', by FFT on a zero-padded grid (isolated sources)."""
    N = q.shape[0]
    k = h * (np.arange(2 * N) - N)
    KX, KY = np.meshgrid(k, k, indexing='ij')
    r2 = KX ** 2 + KY ** 2; r2[N, N] = np.inf
    pad = np.zeros((2 * N, 2 * N)); pad[:N, :N] = q
    out = []
    for K in (KX / (2 * np.pi * r2), KY / (2 * np.pi * r2)):
        c = np.real(np.fft.ifft2(np.fft.fft2(pad) * np.fft.fft2(np.fft.ifftshift(K)))) * h * h
        out.append(c[:N, :N])
    return out[0], out[1]


def advect_x(I, c, h):
    """Upwind d(c I)/dx for a constant speed c per direction (shape (M,1,1)); nothing enters at the edges."""
    pos = c > 0
    d = np.zeros_like(I)
    Ip = I[pos[:, 0, 0]]; d[pos[:, 0, 0]] = (Ip - np.concatenate([np.zeros_like(Ip[:, :1]), Ip[:, :-1]], 1)) / h * c[pos[:, 0, 0]]
    In = I[~pos[:, 0, 0]]; d[~pos[:, 0, 0]] = (np.concatenate([In[:, 1:], np.zeros_like(In[:, :1])], 1) - In) / h * c[~pos[:, 0, 0]]
    return d


def kinetic(q, h, rule, M=64, sigma=10.0, kappa=200.0, nu=200.0, kappa_vm=40.0, t_end=None, tol=2e-5):
    N = q.shape[0]
    th = 2 * np.pi * (np.arange(M) + 0.5) / M
    cx, cy = (U * np.cos(th))[:, None, None], (U * np.sin(th))[:, None, None]
    I = np.zeros((M, N, N))
    dt = 0.4 * h / U
    t_end = t_end or (40.0 if rule == 'scatter' else 6.0)
    t = 0.0; last = None; loss_rate = 0.0
    while t < t_end:
        dI = -advect_x(I, cx, h) - np.swapaxes(advect_x(np.swapaxes(I, 1, 2), cy, h), 1, 2) + q[None] / M
        I = I + dt * dI
        if rule == 'scatter':
            n = I.sum(0)
            I = I + (1 - np.exp(-sigma * U * dt)) * (n[None] / M - I)
        elif rule == 'annihilate':
            a, b = I[:M // 2], I[M // 2:]
            d = np.minimum(a, b) * (1 - np.exp(-kappa * dt * 0.5 * (a + b)))
            I = np.concatenate([a - d, b - d]); loss_rate = 2 * d.sum() * h * h / dt
        elif rule == 'align':
            n = I.sum(0); jx = (I * cx).sum(0); jy = (I * cy).sum(0)
            phi = np.arctan2(jy, jx)
            K = np.exp(kappa_vm * np.cos(th[:, None, None] - phi[None]))
            K /= K.sum(0, keepdims=True)
            ok = (np.hypot(jx, jy) > 1e-9 * (n + 1e-30))[None]
            I = np.where(ok, I + (1 - np.exp(-nu * dt)) * (n[None] * K - I), I)
        t += dt
        if int(round(t / dt)) % 200 == 0:
            n = I.sum(0)
            if last is not None and np.abs(n - last).max() < tol * n.max():
                break
            last = n.copy()
    n = I.sum(0); jx = (I * cx).sum(0); jy = (I * cy).sum(0)
    out = h * ((np.clip(cx, 0, None)[:, 0, 0] * I[:, -1, :].sum(1)).sum() + (np.clip(-cx, 0, None)[:, 0, 0] * I[:, 0, :].sum(1)).sum()
               + (np.clip(cy, 0, None)[:, 0, 0] * I[:, :, -1].sum(1)).sum() + (np.clip(-cy, 0, None)[:, 0, 0] * I[:, :, 0].sum(1)).sum())
    return n, jx, jy, dict(power_out=float(out), power_lost=float(loss_rate), t=t)


def guided(q, h, jnx, jny, tol=1e-6):
    """One stream carried along the field direction e = J_N/|J_N| at speed U: dn/dt + div(n U e) = q."""
    N = q.shape[0]
    m = np.hypot(jnx, jny); ex = np.where(m > 0, jnx / np.maximum(m, 1e-300), 0.); ey = np.where(m > 0, jny / np.maximum(m, 1e-300), 0.)
    vx = U * 0.5 * (ex[1:, :] + ex[:-1, :]); vy = U * 0.5 * (ey[:, 1:] + ey[:, :-1])     # face speeds
    n = np.zeros((N, N)); dt = 0.4 * h / U; t = 0.0; last = None
    while t < 20.0:
        Fx = np.zeros((N + 1, N)); Fy = np.zeros((N, N + 1))
        Fx[1:-1] = np.where(vx > 0, vx * n[:-1], vx * n[1:])
        Fx[-1] = np.clip(U * ex[-1], 0, None) * n[-1]; Fx[0] = np.clip(U * ex[0], None, 0) * n[0]
        Fy[:, 1:-1] = np.where(vy > 0, vy * n[:, :-1], vy * n[:, 1:])
        Fy[:, -1] = np.clip(U * ey[:, -1], 0, None) * n[:, -1]; Fy[:, 0] = np.clip(U * ey[:, 0], None, 0) * n[:, 0]
        n = n + dt * (q - (Fx[1:] - Fx[:-1]) / h - (Fy[:, 1:] - Fy[:, :-1]) / h)
        t += dt
        if int(round(t / dt)) % 200 == 0:
            if last is not None and np.abs(n - last).max() < tol * n.max():
                break
            last = n.copy()
    out = h * (Fx[-1].sum() - Fx[0].sum() + Fy[:, -1].sum() - Fy[:, 0].sum())
    return n, n * U * ex, n * U * ey, dict(power_out=float(out), power_lost=0.0, t=t)


def grade(n, jx, jy, jnx, jny, q, h, X, Y, info):
    """The three requirements, inside |x|, |y| < 0.8 (away from the box's edge), weighted by the Newtonian flux."""
    inner = (np.abs(X) < 0.8) & (np.abs(Y) < 0.8)
    mN = np.hypot(jnx, jny); mJ = np.hypot(jx, jy)
    w = np.where(inner, mN, 0.)
    sel = inner & (mN > 0.02 * mN[inner].max())
    ws = np.where(sel, mN, 0.)
    f = np.where(n > 0, mJ / (U * np.maximum(n, 1e-300)), 0.)
    dlog = np.where(sel, np.log10(np.maximum(mJ, 1e-300) / np.maximum(mN, 1e-300)), 0.)
    cosd = np.where(sel, (jx * jnx + jy * jny) / np.maximum(mJ * mN, 1e-300), 1.)
    dxx, dxy = np.gradient(jx, h, axis=0), np.gradient(jx, h, axis=1)
    dyx, dyy = np.gradient(jy, h, axis=0), np.gradient(jy, h, axis=1)
    curl = np.abs(dyx - dxy)[inner].sum() / (np.abs(dxx) + np.abs(dxy) + np.abs(dyx) + np.abs(dyy))[inner].sum()
    return dict(power_out_over_in=info['power_out'], power_lost=info['power_lost'],
                one_stream_f=float((f * w).sum() / w.sum()), f_below_0p9=float(((f < 0.9) & inner).sum() / inner.sum()),
                flux_vs_newton_rms_dex=float(np.sqrt((ws * dlog ** 2).sum() / ws.sum())),
                direction_miss_deg=float(np.degrees(np.arccos(np.clip((ws * cosd).sum() / ws.sum(), -1, 1)))),
                whirl=float(curl), t_steady=info['t'])


def convergence(out):
    """Which departures are the grid's and which are the rule's: the guided and free rules (exact Newtonian flux in
    the continuum) and the alignment rule on the two equal blobs, at three resolutions."""
    rows = []
    for N, M in ((64, 32), (128, 64), (256, 64)):
        X, Y, h = grid(N)
        q = sources('two equal blobs', X, Y, h); jnx, jny = newtonian_flux(q, h)
        for rule in ('guided', 'free', 'align'):
            if rule != 'guided' and N == 256:
                continue
            n, jx, jy, info = guided(q, h, jnx, jny) if rule == 'guided' else kinetic(q, h, rule, M=M)
            g = grade(n, jx, jy, jnx, jny, q, h, X, Y, info)
            rows.append(dict(N=N, M=M, rule=rule, **g))
            print(f"N {N:3d} M {M:2d} {rule:7s} f {g['one_stream_f']:.3f}; flux vs Newton {g['flux_vs_newton_rms_dex']:.3f} dex, "
                  f"{g['direction_miss_deg']:.1f} deg; whirl {g['whirl']:.3f}", flush=True)
    out.write_text(json.dumps(dict(experiment='round 12: the toy at three resolutions (two equal blobs)', rows=rows), indent=1) + '\n')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--N', type=int, default=128); ap.add_argument('--M', type=int, default=64)
    ap.add_argument('--convergence', action='store_true', help='only the resolution study, written to --output')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.convergence:
        convergence(args.output)
        return
    X, Y, h = grid(args.N)
    res = dict(experiment='round 12: the dynamical toy of the companion, local rules graded on the three requirements',
               grid=dict(N=args.N, M_directions=args.M, box=[-1, 1], graded_region='|x|, |y| < 0.8'), cases=[])
    t0 = time.monotonic()
    for kind in ('uniform disk', 'two equal blobs', 'unequal pair 4:1', 'three blobs'):
        q = sources(kind, X, Y, h)
        jnx, jny = newtonian_flux(q, h)
        ref = grade(np.hypot(jnx, jny) / U, jnx, jny, jnx, jny, q, h, X, Y, dict(power_out=None, power_lost=0., t=0.))
        res['cases'].append(dict(sources=kind, rule='(the Newtonian pattern itself, a discretization check)', **ref))
        print(f"{kind:17s} (Newton)   whirl {ref['whirl']:.3f} (the grid's floor)", flush=True)
        for rule in ('free', 'scatter', 'annihilate', 'align', 'guided'):
            if rule == 'guided':
                n, jx, jy, info = guided(q, h, jnx, jny)
            else:
                n, jx, jy, info = kinetic(q, h, rule, M=args.M)
            g = grade(n, jx, jy, jnx, jny, q, h, X, Y, info)
            res['cases'].append(dict(sources=kind, rule=rule, **g))
            print(f"{kind:17s} {rule:10s} power out/in {g['power_out_over_in']:.3f} (lost {g['power_lost']:.3f}); one stream f {g['one_stream_f']:.3f}"
                  f" (f<0.9 on {100 * g['f_below_0p9']:.0f}%); flux vs Newton {g['flux_vs_newton_rms_dex']:.3f} dex, {g['direction_miss_deg']:.1f} deg;"
                  f" whirl {g['whirl']:.3f} [t {g['t_steady']:.1f}; {time.monotonic() - t0:.0f} s]", flush=True)
    res['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(res, indent=1) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
