"""Round 22, step 1: in which frame does the collision's glow travel?

Round 16's crossing heat (code/crossing_heat_v16.py) puts the glow a star emitted a time t ago on a sphere of radius
u t centred on the star now: the glow keeps its star's motion (the companion's memory). Round 21
(code/hot_mode_speed_v21.py) let that sphere grow at a faster v_h instead: it fits the Bullet Cluster, but only if the
collision puts about 3.3 times more power into the glow than the heat rule gives, because a glow spreading faster with
the same power is thinner.

Here the glow keeps only a fraction mu of its star's motion relative to the other system:
  mu = 1   round 16: the sphere is centred on the star;
  mu = 0   the glow stays in the other system's frame, like a boat's wake: the sphere is centred where it was emitted;
  between  the centre falls behind the star by (1 - mu) x (the path the star has travelled relative to the other system
           since the emission).
The glow still travels at u in its own frame, so each time slice carries exactly round 16's energy: nothing is added,
only where it sits changes. Seen from the star, the glow spreads backwards at up to (1 - mu) v + u, which is 600-900
km/s for mu = 0.8-0.85 and the Bullet's 3,000-3,900 km/s.

Kernels. Round 16's S-kernel k(d)/d^2 carries 4 pi k(t) u dt of integrated weight in each time slice [t, t + dt]
(t = d/u). Here each slice deposits the same weight on a sphere of radius u t about the displaced centre c(t) (a
Gaussian shell of width sigma, normalised on the grid), and the matching pull kernel points from the point back to c(t).
The time dependence k(t) = p(t) (v^2 - 2 u v cos theta)/u^2 and the straight-pass history are round 16's. The motion is
taken along the line joining the two systems' galaxies.

    python code/crossing_frame_v22.py --mu 1,0.9,0.85,0.8,0.7,0.5,0 --output run-crossing-frame-v22/crossing_frame_v22.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import bullet_v3 as B                             # noqa: E402
import bullet_v4 as V                             # noqa: E402
import crossing_heat_v16 as CH                    # noqa: E402
import collisions_v4 as C4                        # noqa: E402
import collisions_v10 as C10                      # noqa: E402

G = B.G
MYR = CH.MYR_PER_KPC_PER_KMS                      # 977.8 Myr per (kpc per km/s)
SETTINGS = dict(mu=1.0, sigma_cells=0.5, dt_max=5.0)
_ORIG_MAKE, _ORIG_MAP = CH.make_kernel, CH.kappa_map_cross


class FrameKernel:
    """Stands in for round 16's k(d): keeps the history so that kappa_map_frame can build the displaced shells."""

    def __init__(self, u, I_other, I_own_at_shell, D_now, v_out, v_in, reach, kmax_d, b=0.0):
        self.u, self.I_other, self.I_own, self.D_now = u, I_other, I_own_at_shell, D_now
        self.v_out, self.v_in, self.reach, self.b = v_out, v_in, reach, b
        self._iso = _ORIG_MAKE(u, I_other, I_own_at_shell, D_now, v_out, v_in, reach, kmax_d, b)

    def __call__(self, d):                        # round 16's k(d), used only for its bookkeeping printout
        return self._iso(d)

    def k_of_t(self, t):
        u = self.u; r_o, I_o = self.I_other
        D, v, cos, inside = CH.history(t, self.D_now, self.v_out, self.v_in, self.reach, self.b)
        Io = np.exp(np.interp(np.log(np.maximum(D, r_o[0])), np.log(r_o), np.log(I_o)))
        p = np.where(inside, Io / (Io + self.I_own), 0.0)
        return p * np.maximum(v ** 2 - 2 * u * v * cos, 0.0) / u ** 2

    def path_of_t(self, t):
        """Path length travelled relative to the other system in the last t Myr (kpc)."""
        vo, vi = self.v_out / MYR, self.v_in / MYR
        tp = np.sqrt(max(self.D_now ** 2 - self.b ** 2, 0.0)) / vo
        return np.where(t <= tp, vo * t, vo * tp + vi * (t - tp))

    def t_max(self):
        """Last time (Myr ago) at which the other system's flow reached the star."""
        vo, vi = self.v_out / MYR, self.v_in / MYR
        tp = np.sqrt(max(self.D_now ** 2 - self.b ** 2, 0.0)) / vo
        return tp + np.sqrt(max(self.reach ** 2 - self.b ** 2, 0.0)) / vi


def frame_kernels(n, dx, kernels, direction, mu, sigma_cells=0.5, dt_max=5.0):
    """Real-space S-kernels and pull kernels (float32, on the zero-padded FFT grid) of several shells of one system, the
    glow keeping a fraction mu of the system's motion relative to the other one. `direction`: unit 3-vector of that
    motion (receding now). Slices whose sphere would reach beyond the kernel grid (more than n dx from the star, well
    outside the map) are left out rather than wrapped around."""
    m = 2 * n
    kk = ((np.arange(m) - (np.arange(m) >= n) * m) * dx).astype(np.float32)
    X, Y, Z = kk[:, None, None], kk[None, :, None], kk[None, None, :]
    u = kernels[0].u; ukpc = u / MYR
    sig = np.float32(sigma_cells * dx)
    vmax = max(kernels[0].v_out, kernels[0].v_in) / MYR
    dt = min(dt_max, 0.4 * dx / ((1 - mu) * vmax + ukpc))
    tmax = max(kf.t_max() for kf in kernels)
    edges = np.arange(0.0, tmax + dt, dt)
    tm = 0.5 * (edges[:-1] + edges[1:])
    ks = [kf.k_of_t(tm) for kf in kernels]
    path = kernels[0].path_of_t(tm)
    accS = [np.zeros((m, m, m), np.float32) for _ in kernels]
    accV = [[np.zeros((m, m, m), np.float32) for _ in range(3)] for _ in kernels]
    tmp = np.empty((m, m, m), np.float32)
    dvol = dx ** 3
    used = dropped = 0.0
    for i, t in enumerate(tm):
        w = [4 * np.pi * float(k[i]) * ukpc * dt for k in ks]
        if max(w) <= 0:
            continue
        c = -(1 - mu) * float(path[i]) * np.asarray(direction, float)
        R = ukpc * t
        if np.linalg.norm(c) + R + 3 * sig > n * dx:
            dropped += w[0]
            continue
        used += w[0]
        dxg = X - np.float32(c[0]); dyg = Y - np.float32(c[1]); dzg = Z - np.float32(c[2])
        rho = np.sqrt(dxg * dxg + dyg * dyg + dzg * dzg)
        prof = np.exp(-0.5 * ((rho - np.float32(R)) / sig) ** 2).astype(np.float32)
        norm = float(prof.sum(dtype=np.float64)) * dvol
        np.maximum(rho, np.float32(1e-3), out=rho)
        np.divide(prof, rho, out=rho)                      # rho now holds prof / distance
        pu = [rho * dxg, rho * dyg, rho * dzg]             # prof times the unit vector from the centre
        for s, ws in enumerate(w):
            if ws <= 0:
                continue
            a = np.float32(ws / norm)
            np.multiply(prof, a, out=tmp); accS[s] += tmp
            for j in range(3):
                np.multiply(pu[j], a, out=tmp); accV[s][j] -= tmp
        del dxg, dyg, dzg, rho, prof, pu
    info = dict(dt_myr=float(dt), t_max_myr=float(tmax), slices=int(len(tm)),
                weight_dropped_fraction=float(dropped / max(used + dropped, 1e-30)))
    return accS, accV, info


def _conv_fft(K, f, n):
    m = 2 * n
    pad = np.zeros((m, m, m), np.float32); pad[:n, :n, :n] = f
    return np.fft.irfftn(np.fft.rfftn(pad) * K, s=(m, m, m), axes=(0, 1, 2))[:n, :n, :n].astype(np.float32)


def kappa_map_frame(current, ghost_gas, ghost_stars, pos, consts, crossing=(), n=192, dx=15.0, centre=(360., 50.),
                    fresh_kpc=30.0, heat=True, memory=True, report=None):
    """crossing_heat_v16.kappa_map_cross with the crossing heat on displaced shells (FrameKernel items)."""
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']; gd = lam * a
    x = (np.arange(n) - n / 2 + 0.5) * dx + centre[0]
    y = (np.arange(n) - n / 2 + 0.5) * dx + centre[1]
    z = (np.arange(n) - n / 2 + 0.5) * dx
    dV = dx ** 3
    rho_now = B.build_density(current, pos, x, y, z, 'gas') + B.build_density(current, pos, x, y, z, 'st')
    rho_gas_now = B.build_density(current, pos, x, y, z, 'gas')
    rho_ghost_gas = B.build_density(ghost_gas, pos, x, y, z, 'gas') if memory else rho_gas_now
    rho_ghost_st = np.zeros_like(rho_now); krho = np.zeros_like(rho_now)
    import law as L
    for comp, prof in ghost_stars:
        rho_c = B.build_density({'c': comp}, pos, x, y, z, 'c')
        cx, cy = pos[comp['centre']]
        r3 = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
        rho_ghost_st += rho_c
        if heat:
            krho += L.k_from_sig2(np.interp(r3, prof[0], prof[1]) ** 2, u).astype(np.float32) * rho_c
        del r3, rho_c
    inv_r = B.Conv(n, dx, lambda r: 1.0 / r, B.cube_average(lambda r: 1 / r) / dx)
    gN = -np.array(np.gradient(-G * inv_r(rho_now * dV), dx))
    F = -np.array(np.gradient(-G * inv_r((rho_ghost_gas + rho_ghost_st) * dV), dx))
    g_hot = -np.array(np.gradient(-G * inv_r(krho * dV), dx)) if heat else 0.0
    del inv_r
    if memory and fresh_kpc > 0:
        Kin = V.vec_conv(n, dx, rmax=fresh_kpc)
        F = F + G * dV * (V.apply_vec(Kin, rho_gas_now, n) - V.apply_vec(Kin, rho_ghost_gas, n))
        del Kin
    inv_r2 = B.Conv(n, dx, lambda r: 1.0 / r ** 2, B.cube_average(lambda r: 1 / r ** 2) / dx ** 2)
    S = G * inv_r2(krho * dV) if heat else 0.0
    del inv_r2
    S_cross_total = 0.0
    info = []
    # group the crossing items by system (the centre their stars belong to)
    groups = {}
    for comps_, shell, kf in crossing:
        groups.setdefault(comps_[0]['centre'], []).append((comps_, shell, kf))
    names = list(groups)
    for who in names:
        items = groups[who]
        other = [nm for nm in names if nm != who][0]
        d = np.array([*(np.asarray(pos[who]) - np.asarray(pos[other])), 0.0]); d = d / np.linalg.norm(d)
        accS, accV, kinfo = frame_kernels(n, dx, [kf for _, _, kf in items], d, SETTINGS['mu'], SETTINGS['sigma_cells'], SETTINGS['dt_max'])
        info.append(dict(system=who, **kinfo))
        for s_i, (comps_, (r0, r1), kf) in enumerate(items):
            rho_c = np.zeros_like(rho_now)
            for comp in comps_:
                rc = B.build_density({'c': comp}, pos, x, y, z, 'c')
                cx, cy = pos[comp['centre']]
                r3 = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
                rho_c += (rc * ((r3 >= r0) & (r3 < r1))).astype(np.float32); del r3, rc
            Sx = G * dV * _conv_fft(np.fft.rfftn(accS[s_i]), rho_c, n); accS[s_i] = None
            gx = np.array([G * dV * _conv_fft(np.fft.rfftn(accV[s_i][j]), rho_c, n) for j in range(3)]); accV[s_i] = None
            S = S + Sx; g_hot = g_hot + gx
            S_cross_total = S_cross_total + Sx
            del rho_c
        del accS, accV
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    Fmag = np.sqrt(np.sum(F ** 2, axis=0)) + 1e-30
    hmag = np.sqrt(np.sum(g_hot ** 2, axis=0)) if heat else 0.0
    extra = np.exp(-mag / gd) * np.sqrt(a * (Fmag + S))
    h = gN + extra * (F + g_hot) / (Fmag + hmag + 1e-30)
    divh = sum(np.gradient(h[i], dx, axis=i) for i in range(3))
    Sig_eff = (-divh / (4 * np.pi * G)).sum(axis=2) * dx
    Sig_b = rho_now.sum(axis=2) * dx
    if report is not None and crossing:
        k0 = n // 2
        report['S_cross_over_S_mid_plane'] = float(np.sum(S_cross_total[:, :, k0]) / max(np.sum(S[:, :, k0]), 1e-30))
        report['kernels'] = info
        report['_kappa_inputs'] = (x, y)
    return x, y, Sig_eff, Sig_b


def signed_offsets(x, y, kap, pos, fs):
    """Lensing peaks' offsets from their galaxies along the collision axis, positive toward the system's own gas."""
    pk = C4.refined_peaks(x, y, kap, pos, smooth_kpc=40.0 * fs)
    out = {}
    for w in ('main', 'sub'):
        p = min(pk, key=lambda d: d[f'dist_{w}_bcg'])
        g = np.asarray(pos[f'{w}_plasma']) - np.asarray(pos[f'{w}_bcg'])
        out[w] = float((np.array([p['x'], p['y']]) - np.asarray(pos[f'{w}_bcg'])) @ (g / np.linalg.norm(g)))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--mu', type=str, default='1,0.9,0.85,0.8,0.7,0.5,0')
    ap.add_argument('--grid', choices=['test', 'coarse', 'fine'], default='coarse')
    ap.add_argument('--sigma-cells', type=float, default=0.5)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    n, dx = dict(test=(64, 45.0), coarse=(128, 22.5), fine=(192, 15.0))[args.grid]
    cross = dict(v_out=3900.0, v_in=3000.0, b=150.0)
    SETTINGS['sigma_cells'] = args.sigma_cells
    out = dict(experiment='round 22: the frame the collision glow travels in (memory fraction mu), energy booked', law=law['name'],
               grid=dict(n=n, dx_kpc=dx, sigma_cells=args.sigma_cells), crossing=cross, runs=[])
    # the capture of the lensing map for the signed offsets
    captured = {}

    def capture_map(*a, **k):
        res = kappa_map_frame(*a, **k)
        captured['pos'] = a[3]; captured['x'], captured['y'], captured['Se'] = res[0], res[1], res[2]
        return res
    CH.make_kernel = FrameKernel
    CH.kappa_map_cross = capture_map
    try:
        for mu in [float(s) for s in args.mu.split(',')]:
            SETTINGS['mu'] = mu
            r = CH.run_bullet_cross(law, dict(f), cross=cross, n=n, dx=dx)
            import bullet_static_v11 as BS
            kap = captured['Se'] / BS.static_sigma_crit()
            off = signed_offsets(captured['x'], captured['y'], kap, captured['pos'], f['size'])
            row = dict(mu=mu, peak_toward_gas_main=off['main'], peak_toward_gas_sub=off['sub'], m250_sub=r['m250_sub'], m250_main=r['m250_main'], peak_main=r['peak_main'], peak_sub=r['peak_sub'],
                       gas_main=r['gas_main'], gas_sub=r['gas_sub'], kappa_sub=r['kappa_sub'], kappa_main=r['kappa_main'],
                       S_cross_share=r['report'].get('S_cross_over_S_mid_plane'), kernels=r['report'].get('kernels'),
                       targets=r['targets'], peak_limits=r['peak_limits'])
            out['runs'].append(row)
            print(f"[{time.monotonic() - t0:5.0f} s] mu {mu:4.2f}: sub {r['m250_sub']:.3e} (target {r['targets']['m250_sub'][0]:.2e}-"
                  f"{r['targets']['m250_sub'][1]:.2e}), main {r['m250_main']:.3e} (target {r['targets']['m250_main'][0]:.2e}-"
                  f"{r['targets']['m250_main'][1]:.2e}); peaks {r['peak_main']:.0f}, {r['peak_sub']:.0f} kpc (toward gas {off['main']:+.0f}, {off['sub']:+.0f}); gas {r['gas_main']:.3f}, "
                  f"{r['gas_sub']:.3f}; crossing share {row['S_cross_share']:.2f}", flush=True)
            args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    finally:
        CH.make_kernel, CH.kappa_map_cross = _ORIG_MAKE, _ORIG_MAP
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
