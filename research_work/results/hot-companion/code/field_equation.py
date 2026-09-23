#!/usr/bin/env python3
"""A field equation for hot-companion gravity, solved in 3D, and three checks.

    python field_equation.py --output-dir ../run-field

Field equation (proposed here; quasi-linear, same structure as Milgrom's QUMOND, which it
contains as its cold limit):
    (1)  laplacian Phi_N = 4 pi G rho                         ordinary (Newtonian) potential
    (2)  S_hot = G * integral k rho / d^2                      incoherent companion intensity
    (3)  h = g_N + exp(-|g_N|/g_d) sqrt(a (|g_N| + S_hot)) g_N/|g_N|
    (4)  laplacian Phi = - div h,   g = - grad Phi              the pull on matter and light
(4) makes the pull conservative: it is the curl-free part of h. In spherical symmetry
g = h exactly, which is the algebraic law tested on data.

Checks:
  A. disk galaxy with a hot bulge: in-plane pull from (4) versus the algebraic law (3)
  B. the same with the heat switched off: the cold limit (QUMOND with our derived nu)
  C. momentum: the net force on an isolated pair of blobs, one hot and one cold
Units: kpc, km/s, Msun.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.special import erfcx

G = 4.30091727003628e-6
K = 1e6 / 3.0856775814913673e19
A = 6.5797e-11 / K                      # (km/s)^2/kpc
GD = 4.281 * A
U = 874.1


def grid(n, box):
    dx = box / n
    x = (np.arange(n) - n / 2 + 0.5) * dx
    return x, dx


class Conv:
    """Isolated-boundary convolution with a radial kernel via zero-padded FFT (float32)."""
    def __init__(self, n, dx, kernel_fn, origin_value):
        m = 2 * n
        k = (np.arange(m) - (np.arange(m) >= n) * m) * dx
        X, Y, Z = np.meshgrid(k, k, k, indexing='ij')
        r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
        with np.errstate(divide='ignore'):
            ker = np.where(r > 0, kernel_fn(np.maximum(r, 1e-30)), origin_value)
        self.n = n
        self.K = np.fft.rfftn(ker.astype(np.float32))
        del X, Y, Z, r, ker

    def __call__(self, f):
        n = self.n; m = 2 * n
        pad = np.zeros((m, m, m), np.float32); pad[:n, :n, :n] = f
        return np.fft.irfftn(np.fft.rfftn(pad) * self.K, s=(m, m, m), axes=(0, 1, 2))[:n, :n, :n]


def cube_average(fn, samples=400000, seed=1):
    rng = np.random.default_rng(seed); p = rng.uniform(-.5, .5, (samples, 3))
    return float(np.mean(fn(np.linalg.norm(p, axis=1))))


def grad(f, dx):
    return np.array(np.gradient(f, dx))


def div(v, dx):
    return sum(np.gradient(v[i], dx, axis=i) for i in range(3))


def solve(rho, krho, dx, inv_r, inv_r2, heat=True):
    dV = dx ** 3
    phiN = -G * inv_r(rho * dV)
    gN = -grad(phiN, dx)
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    S = G * inv_r2(krho * dV) if heat else 0.
    h = gN + (np.exp(-mag / GD) * np.sqrt(A * (mag + S)) / mag) * gN
    phi = inv_r(div(h, dx) * dV) / (4 * np.pi)
    g = -grad(phi, dx)
    return gN, h, g


def dW_dS(q, S):
    """dW/dS = int_0^q exp(-q'/g_d) sqrt(a/(q'+S)) dq', in closed form (stable via erfcx)."""
    x1 = np.sqrt(np.maximum(S, 0) / GD); x2 = np.sqrt((q + np.maximum(S, 0)) / GD)
    return np.sqrt(A * np.pi * GD) * (erfcx(x1) - np.exp(-q / GD) * erfcx(x2))


def reaction_force(gN, S, k, dx, inv_r2):
    """Extra force per unit mass on hot matter from the action: (k/8pi) grad Psi,
    Psi = int dW/dS(x) / |x - x'|^2 d^3x. Restores momentum conservation (Noether)."""
    q = np.sqrt(np.sum(gN ** 2, axis=0))
    Psi = inv_r2(dW_dS(q, S) * dx ** 3)
    return (k / (8 * np.pi)) * grad(Psi, dx)


def exp_disk(X, Y, Z, M, Rd, hz):
    R = np.sqrt(X ** 2 + Y ** 2)
    return M / (4 * np.pi * Rd ** 2 * hz) * np.exp(-R / Rd) / np.cosh(Z / hz) ** 2


def hernquist(X, Y, Z, M, a, c=(0, 0, 0)):
    r = np.sqrt((X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2) + 1e-9
    return M * a / (2 * np.pi * r * (r + a) ** 3)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--n', type=int, default=160); ap.add_argument('--box', type=float, default=80.)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    x, dx = grid(args.n, args.box)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    c1 = cube_average(lambda r: 1 / r) / dx          # cell average of 1/r
    c2 = cube_average(lambda r: 1 / r ** 2) / dx ** 2
    inv_r = Conv(args.n, dx, lambda r: 1 / r, c1)
    inv_r2 = Conv(args.n, dx, lambda r: 1 / r ** 2, c2)
    res = dict(grid=dict(n=args.n, box_kpc=args.box, dx_kpc=dx))

    # ---- A/B: disk + gas disk + hot bulge
    sig_b = 150.; kb = 3 * sig_b ** 2 / U ** 2
    rho_d = exp_disk(X, Y, Z, 4e10, 3.0, 0.4) + exp_disk(X, Y, Z, 1e10, 7.0, 0.3)
    rho_b = hernquist(X, Y, Z, 1.2e10, 0.7)
    rho = rho_d + rho_b; krho = kb * rho_b
    i0 = args.n // 2; j0 = args.n // 2          # plane z ~ 0 (two cells straddle it; average them)
    Rs = np.array([1, 2, 3, 5, 8, 12, 16, 20, 25, 30])
    for tag, heat in (('hot_bulge', True), ('cold_limit', False)):
        gN, h, g = solve(rho, krho, dx, inv_r, inv_r2, heat)
        rows = []
        for R in Rs:
            ix = np.argmin(abs(x - R))
            # radial components along +x at z = +-dx/2, averaged
            pick = lambda v: 0.5 * (v[0][ix, j0, i0 - 1] + v[0][ix, j0, i0])
            gn, hh, gg = -pick(gN), -pick(h), -pick(g)
            rows.append(dict(R_kpc=float(x[ix]), gN=float(gn), algebraic=float(hh), field=float(gg),
                             field_over_algebraic=float(gg / hh)))
        res[tag] = rows
        # conservativeness: the curl of g must vanish (it is a gradient); report curl of h
        curl_h = np.sqrt(sum(c ** 2 for c in (
            np.gradient(h[2], dx, axis=1) - np.gradient(h[1], dx, axis=2),
            np.gradient(h[0], dx, axis=2) - np.gradient(h[2], dx, axis=0),
            np.gradient(h[1], dx, axis=0) - np.gradient(h[0], dx, axis=1))))
        res[tag + '_max_curl_h_over_div_h'] = float(np.max(curl_h[20:-20, 20:-20, 20:-20]) /
                                                    np.max(np.abs(div(h, dx))[20:-20, 20:-20, 20:-20]))

    # ---- C: momentum, one hot and one cold blob 30 kpc apart, in growing boxes.
    # Truncating the extra field at the box edge leaves a residual that must shrink as the
    # box grows; a real violation does not shrink.
    mom = []
    for n_c, box_c in ((160, 160.), (200, 240.), (240, 360.)):
        xc, dxc = grid(n_c, box_c); Xc, Yc, Zc = np.meshgrid(xc, xc, xc, indexing='ij')
        ir = Conv(n_c, dxc, lambda r: 1 / r, cube_average(lambda r: 1 / r) / dxc)
        ir2 = Conv(n_c, dxc, lambda r: 1 / r ** 2, cube_average(lambda r: 1 / r ** 2) / dxc ** 2)
        rho1 = hernquist(Xc, Yc, Zc, 5e10, 2.0, (-15, 0, 0)); rho2 = hernquist(Xc, Yc, Zc, 2e10, 1.5, (15, 0, 0))
        dV = dxc ** 3; row = dict(box_kpc=box_c)
        for tag, sig1 in (('cold', 0.), ('hot800', 800.)):
            k1 = 3 * sig1 ** 2 / U ** 2
            gN, h, g = solve(rho1 + rho2, k1 * rho1, dxc, ir, ir2, heat=sig1 > 0)
            F1 = np.array([np.sum(rho1 * g[i]) * dV for i in range(3)])
            F2 = np.array([np.sum(rho2 * g[i]) * dV for i in range(3)])
            row[tag + '_net_over_pair'] = float(np.linalg.norm(F1 + F2) / np.linalg.norm(F1))
            if sig1 > 0:
                S = G * ir2(k1 * rho1 * dV)
                Fx = reaction_force(gN, S, 1.0, dxc, ir2)
                R1 = np.array([np.sum(k1 * rho1 * Fx[i]) * dV for i in range(3)])
                row['hot800_with_action_reaction_net_over_pair'] = float(np.linalg.norm(F1 + F2 + R1) / np.linalg.norm(F1))
                row['reaction_over_pair'] = float(np.linalg.norm(R1) / np.linalg.norm(F1))
        mom.append(row)
        del Xc, Yc, Zc, ir, ir2
    res['momentum'] = mom
    res['seconds'] = time.monotonic() - t0
    (out / 'field.json').write_text(json.dumps(res, indent=2) + '\n')

    print(f'grid {args.n}^3, box {args.box} kpc, dx {dx:.2f} kpc; bulge sigma {sig_b} km/s (k = {kb:.3f})')
    for tag in ('hot_bulge', 'cold_limit'):
        print(f'\n{tag}: in-plane pull, field equation / algebraic law   (max |curl h|/|div h| = {res[tag + "_max_curl_h_over_div_h"]:.3f})')
        for r in res[tag]:
            print(f"   R={r['R_kpc']:5.1f} kpc  gN {r['gN']:8.1f}  algebraic {r['algebraic']:8.1f}  field {r['field']:8.1f}  ratio {r['field_over_algebraic']:.3f}")
    print('\nmomentum: net force on an isolated pair / force on one body (must shrink with box size)')
    for v in mom:
        print(f"   box {v['box_kpc']:4.0f} kpc: cold {v['cold_net_over_pair']:.4f} | hot, law as written "
              f"{v['hot800_net_over_pair']:.4f} | hot, with the action's reaction force "
              f"{v['hot800_with_action_reaction_net_over_pair']:.4f}")
    print(f"\n{res['seconds']:.0f} s")


if __name__ == '__main__':
    main()
