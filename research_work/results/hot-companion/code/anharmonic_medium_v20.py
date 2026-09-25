"""Round 20, step 2: does a stiffening medium, once energized, reach farther or only change its range?

The review of 25 September (energy-shift-v20/README.md, "Next discriminating experiment", item 3) proposes adding one
energy-bounded nonlinearity to the positive-energy (gapped) medium, a local quartic term beta X^4 / 4, and asking
whether a finite excitation density changes only the medium's range (a Yukawa length) or its long-wavelength response
(the force's exponent). Its stop/redirect rule: stop if the nonlinearity only changes a Yukawa length.

The medium, on a cubic lattice (spacing 1, wave speed 1):

    E = sum_n [ (m0^2 / 2) X_n^2 + (beta / 4) X_n^4 ] + (1/2) sum_<nm> (X_n - X_m)^2

Two static sources coupled linearly to the medium (energy -q X at their sites) interact through the medium's static
response chi(r) (in a thermal state, chi = <X_0 X_r> / T): energy -q1 q2 chi(r), force q1 q2 chi'(r). The exchange
kernel between two pieces at a transition frequency below the gap has the same structure with m^2 replaced by
m^2 - omega_a^2 (a Yukawa length that only grows shorter as m grows).

Parts:
  mc        exact classical statistics (Metropolis, checkerboard) at several excitation levels T and stiffness beta:
            the wall-to-wall correlator gives the range; the point-to-point correlator gives the force law; both are
            compared with the self-consistent harmonic estimate m_eff^2 = m0^2 + 3 beta <X^2>.
  hot       a medium energized non-uniformly (a hot region around a source, or heating falling as 1/r), in the local
            self-consistent harmonic approximation: does a probe outside feel a slower-falling force?
  critical  the gap closed (m0 = 0): the static nonlinear response to a point source, -lap phi + beta phi^3 = source,
            solved exactly in radius: how the force on a test body depends on distance and on the source's strength.

    python code/anharmonic_medium_v20.py --part mc --output run-anharmonic-medium-v20/mc.json
    python code/anharmonic_medium_v20.py --part hot --output run-anharmonic-medium-v20/hot.json
    python code/anharmonic_medium_v20.py --part critical --output run-anharmonic-medium-v20/critical.json
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

M0SQ = 0.1                       # bare gap^2 (lattice units): bare range 1/kappa0 = 3.2 spacings


# ----------------------------------------------------------------------------------------------- lattice sums
def g_local(msq, L=64):
    """<X^2> / T for the Gaussian lattice medium of mass^2 msq: (1/N) sum_k 1 / (msq + sum_a 2(1 - cos k_a))."""
    k = 2 * np.pi * np.arange(L) / L
    c = 2 * (1 - np.cos(k))
    s = c[:, None, None] + c[None, :, None] + c[None, None, :]
    return float(np.mean(1.0 / (msq + s)))


def scha_msq(T, beta, L=64, iters=200):
    """Self-consistent harmonic mass: m^2 = m0^2 + 3 beta <X^2>, <X^2> = T g_local(m^2)."""
    m = M0SQ
    for _ in range(iters):
        new = M0SQ + 3 * beta * T * g_local(m, L)
        if abs(new - m) < 1e-10:
            break
        m = 0.5 * m + 0.5 * new
    return m


def kappa_of(msq):
    """Decay rate of the zero-transverse-momentum (wall) correlator on the lattice: 2(cosh k - 1) = m^2."""
    return float(np.arccosh(1 + msq / 2))


# ----------------------------------------------------------------------------------------------- Monte Carlo
@njit(cache=True)
def sweep(X, T, beta, msq, step, rng_u, rng_a):
    L = X.shape[0]
    acc = 0
    idx = 0
    for parity in range(2):
        for i in range(L):
            for j in range(L):
                for k in range((i + j + parity) % 2, L, 2):
                    x = X[i, j, k]
                    h = (X[(i + 1) % L, j, k] + X[(i - 1) % L, j, k] + X[i, (j + 1) % L, k] + X[i, (j - 1) % L, k]
                         + X[i, j, (k + 1) % L] + X[i, j, (k - 1) % L])
                    y = x + step * (2 * rng_u[idx] - 1)
                    dE = ((0.5 * msq + 3.0) * (y * y - x * x) + 0.25 * beta * (y ** 4 - x ** 4) - (y - x) * h)
                    if dE <= 0 or rng_a[idx] < np.exp(-dE / T):
                        X[i, j, k] = y
                        acc += 1
                    idx += 1
    return acc


def run_mc(job):
    T, beta, L, n_therm, n_meas, every, seed = job
    rng = np.random.default_rng(seed)
    msq_guess = scha_msq(T, beta, L=L)
    X = rng.normal(0, np.sqrt(T * g_local(msq_guess, L)), (L, L, L))
    step = 2.0 * np.sqrt(T / (msq_guess + 6 + 3 * beta * T * g_local(msq_guess, L)))
    N = L ** 3
    t0 = time.monotonic()
    for s in range(n_therm):
        a = sweep(X, T, beta, M0SQ, step, rng.random(N), rng.random(N)) / N
        if s % 20 == 19:                              # tune the step toward ~50% acceptance during thermalisation
            step *= 1.1 if a > 0.55 else (0.9 if a < 0.45 else 1.0)
    rmax = L // 2
    wall = np.zeros(rmax + 1); point = np.zeros(rmax + 1); x2 = 0.0; x4 = 0.0; nm = 0; accs = 0.0
    walls_series = []
    for s in range(n_meas):
        accs += sweep(X, T, beta, M0SQ, step, rng.random(N), rng.random(N)) / N
        if s % every:
            continue
        nm += 1
        x2 += float(np.mean(X ** 2)); x4 += float(np.mean(X ** 4))
        w_this = np.zeros(rmax + 1)
        for ax in range(3):
            phi = X.mean(axis=tuple(a for a in range(3) if a != ax))       # plane averages along axis ax
            for r in range(rmax + 1):
                w_this[r] += np.mean(phi * np.roll(phi, r))
                point[r] += np.mean(X * np.roll(X, r, axis=ax))
        wall += w_this / 3
        walls_series.append((w_this / 3).tolist())
    wall /= nm; point /= (3 * nm); x2 /= nm; x4 /= nm
    # errors on the wall correlator from blocks of 10 measurements
    ws = np.array(walls_series); nb = max(len(ws) // 10, 2)
    blocks = np.array([b.mean(0) for b in np.array_split(ws, nb)])
    wall_err = blocks.std(0, ddof=1) / np.sqrt(nb)
    return dict(T=T, beta=beta, L=L, sweeps=n_therm + n_meas, measurements=nm, acceptance=accs / n_meas, step=step,
                x2=x2, x4=x4, wall=wall.tolist(), wall_err=wall_err.tolist(), point=point.tolist(),
                seconds=time.monotonic() - t0)


def analyse_mc(row):
    L = row['L']; T = row['T']; beta = row['beta']
    w = np.array(row['wall']); we = np.array(row['wall_err'])
    # effective mass from the cosh form of the wall correlator on the periodic lattice: C(t) ~ cosh(k (t - L/2))
    kap = []
    for t in range(1, L // 2 - 1):
        r = (w[t - 1] + w[t + 1]) / (2 * w[t])
        kap.append(float(np.arccosh(r)) if r > 1 and w[t + 1] > 3 * we[t + 1] else np.nan)
    kap = np.array(kap)
    good = np.isfinite(kap[2:12])
    kappa_mc = float(np.nanmedian(kap[2:12][good])) if good.any() else float('nan')
    msq_mc = 2 * (np.cosh(kappa_mc) - 1)
    msq_scha = scha_msq(T, beta, L=L)
    x2_scha = T * g_local(msq_scha, L)
    # the point-to-point force law: chi(r) = <X_0 X_r> / T, force ~ chi'(r); local exponent of |chi'|
    chi = np.array(row['point']) / T
    rr = np.arange(len(chi))
    F = -np.gradient(chi, rr)
    p = []
    for r in range(2, 11):
        if F[r - 1] > 0 and F[r + 1] > 0:
            p.append(float(-(np.log(F[r + 1]) - np.log(F[r - 1])) / (np.log(r + 1) - np.log(r - 1))))
        else:
            p.append(float('nan'))
    # the same exponent for a continuum Yukawa force of the measured range: 2 + (k r)^2 / (1 + k r)
    p_yuk = [2 + (kappa_mc * r) ** 2 / (1 + kappa_mc * r) for r in range(2, 11)]
    return dict(T=T, beta=beta, acceptance=row['acceptance'], x2_mc=row['x2'], x2_scha=x2_scha,
                kappa_mc=kappa_mc, kappa_scha=kappa_of(msq_scha), kappa_bare=kappa_of(M0SQ),
                range_mc=1 / kappa_mc, range_bare=1 / kappa_of(M0SQ),
                msq_mc=float(msq_mc), msq_scha=msq_scha, force_exponent_r2_to_10=p, yukawa_exponent_same_range=p_yuk,
                chi_0=float(chi[0]))


# ----------------------------------------------------------------------------------------------- a hot region
def static_green(msq_field, src):
    """Solve (-lap + m^2(x)) G = delta_src on the periodic lattice by conjugate gradients."""
    from scipy.sparse.linalg import cg, LinearOperator
    L = msq_field.shape[0]; N = L ** 3
    def mv(v):
        X = v.reshape(L, L, L)
        out = (6 + msq_field) * X
        for ax in range(3):
            out -= np.roll(X, 1, ax) + np.roll(X, -1, ax)
        return out.ravel()
    A = LinearOperator((N, N), matvec=mv, dtype=float)
    b = np.zeros(N); b[np.ravel_multi_index(src, (L, L, L))] = 1.0
    x, info = cg(A, b, rtol=1e-11, maxiter=20000)
    return x.reshape(L, L, L), info


def run_hot(L=96, beta=1.0):
    """Probe outside a heated region: local self-consistent harmonic masses, then the static Green's function."""
    c = L // 2
    ii, jj, kk = np.indices((L, L, L))
    r = np.sqrt((ii - c) ** 2 + (jj - c) ** 2 + (kk - c) ** 2)
    table_T = np.geomspace(1e-4, 30, 120)
    table_m = np.array([scha_msq(t, beta, L=32) for t in table_T])      # m^2 as a function of the local T (LDA)
    msq_of_T = lambda T: np.interp(np.log(np.maximum(T, 1e-4)), np.log(table_T), table_m)
    cases = {
        'cold medium (bare range)': np.zeros_like(r),
        'uniformly energized, T = 1': np.ones_like(r),
        'hot region of radius 6 (T = 3 inside, 0 outside)': np.where(r <= 6, 3.0, 0.0),
        'heating falling as 1/r (T = 3 at r = 1)': 3.0 / np.maximum(r, 1.0),
        'heating falling as 1/r^2 (T = 3 at r = 1)': 3.0 / np.maximum(r, 1.0) ** 2,
    }
    out = {}
    radii = np.arange(2, 31)
    for tag, Tf in cases.items():
        msq = msq_of_T(Tf) if Tf.any() else np.full_like(r, M0SQ)
        G, info = static_green(msq, (c, c, c))
        prof = np.array([np.mean([G[c + d, c, c], G[c - d, c, c], G[c, c + d, c], G[c, c - d, c], G[c, c, c + d], G[c, c, c - d]])
                         for d in range(0, 34)])
        F = -np.gradient(prof)
        p = [float(-(np.log(F[d + 1]) - np.log(F[d - 1])) / (np.log(d + 1) - np.log(d - 1))) for d in radii]
        out[tag] = dict(cg_info=int(info), msq_at=dict(r1=float(msq[c + 1, c, c]), r6=float(msq[c + 6, c, c]), r20=float(msq[c + 20, c, c])),
                        green=prof.tolist(), force_exponent=dict(zip([int(x) for x in radii], p)),
                        min_force_exponent=float(np.nanmin(p)))
        print(f"{tag:52s} m^2 at r=1/6/20: {msq[c+1,c,c]:.3f}/{msq[c+6,c,c]:.3f}/{msq[c+20,c,c]:.3f}; force exponent at r = 4, 10, 20: "
              f"{p[2]:.2f}, {p[8]:.2f}, {p[18]:.2f}; smallest {np.nanmin(p):.2f}", flush=True)
    return dict(L=L, beta=beta, m0sq=M0SQ, cases=out)


# ----------------------------------------------------------------------------------------------- the gap closed
def run_critical(beta=1.0):
    """-phi'' - (2/r) phi' + m^2 phi + beta phi^3 = 0 for r > r0, with -4 pi r0^2 phi'(r0) = s and phi -> 0 far away.
    Solved as a boundary-value problem in u = ln r. The force on a weak test body (energy -q phi) is q |phi'(r)|."""
    from scipy.integrate import solve_bvp
    r0, r1 = 1.0, 1.0e5
    out = {}
    for msq in (0.0, 1e-4, 1e-2):
        rows = {}
        for s in np.geomspace(1e-3, 1e5, 33):
            # y0 = phi, y1 = r dphi/dr, as functions of u = ln r:  y0' = y1,  y1' = y1 + r^2 (m^2 phi + beta phi^3) - y1 ... see below
            def f(u, y):
                rr = np.exp(u)
                # d/du (r phi') = r (phi' + r phi'') = r^2 phi'' + r phi'; with phi'' = -(2/r) phi' + m^2 phi + beta phi^3:
                # d/du (r phi') = -r phi' + r^2 (m^2 phi + beta phi^3)
                return np.vstack([y[1], -y[1] + rr ** 2 * (msq * y[0] + beta * y[0] ** 3)])
            def bc(ya, yb):
                # at r0: r phi' = -s / (4 pi r0);  far away: the linear decaying solution, r phi' = -(1 + m r) phi
                return np.array([ya[1] + s / (4 * np.pi * r0), yb[1] + (1 + np.sqrt(msq) * r1) * yb[0]])
            u = np.linspace(np.log(r0), np.log(r1), 4000)
            phi_lin = s / (4 * np.pi * np.exp(u)) * np.exp(-np.sqrt(msq) * (np.exp(u) - r0))
            y_init = np.vstack([np.minimum(phi_lin, 1.0 / (np.sqrt(beta) * np.exp(u)) + 1e-12), -np.minimum(phi_lin, 1.0 / (np.sqrt(beta) * np.exp(u)))])
            sol = solve_bvp(f, bc, u, y_init, tol=1e-7, max_nodes=200000)
            rr = np.exp(sol.x); phi = sol.y[0]; dphi = sol.y[1] / rr
            probe = np.geomspace(3, 3e4, 9)
            Fp = np.interp(np.log(probe), sol.x, np.abs(dphi))
            rows[float(s)] = dict(success=bool(sol.success), force=Fp.tolist(), phi=np.interp(np.log(probe), sol.x, phi).tolist())
        ss = np.array(sorted(rows)); Fm = np.array([rows[s]['force'] for s in ss])      # [n_s, n_probe]
        probe = np.geomspace(3, 3e4, 9)
        p_r = -np.diff(np.log(Fm), axis=1) / np.diff(np.log(probe))[None, :]              # distance exponent
        p_s = np.diff(np.log(Fm), axis=0) / np.diff(np.log(ss))[:, None]                  # source-strength exponent
        out[f'm2={msq:g}'] = dict(sources=ss.tolist(), probe_r=probe.tolist(), force=Fm.tolist(),
                                  distance_exponent=p_r.tolist(), source_exponent=p_s.tolist(),
                                  all_solved=bool(all(rows[s]['success'] for s in ss)))
        print(f"m^2 = {msq:g}: all solved {out[f'm2={msq:g}']['all_solved']}; source exponent at r = 30 for s = 1e-3, 1, 1e3, 1e5: "
              + ', '.join(f'{p_s[i, 2]:.2f}' for i in (0, 12, 24, 31)) + '; distance exponent (r 30-100) for the same: '
              + ', '.join(f'{p_r[i, 2]:.2f}' for i in (0, 12, 24, 32)), flush=True)
    return dict(beta=beta, r0=r0, cases=out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', choices=['mc', 'hot', 'critical'], required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--L', type=int, default=32)
    ap.add_argument('--sweeps', type=int, default=24000)
    ap.add_argument('--procs', type=int, default=3)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    if args.part == 'mc':
        jobs = []
        for beta in (0.0, 0.3, 1.0):
            for T in ((1.0,) if beta == 0 else (0.03, 0.1, 0.3, 1.0, 3.0)):
                jobs.append((T, beta, args.L, 2000, args.sweeps, 10, 1000 + len(jobs)))
        with Pool(args.procs) as pool:
            raw = pool.map(run_mc, jobs)
        rows = [analyse_mc(r) for r in raw]
        for a in rows:
            print(f"beta {a['beta']:.1f}  T {a['T']:5.2f}  <X^2> MC {a['x2_mc']:.4f} (SCHA {a['x2_scha']:.4f})  range MC {a['range_mc']:.2f} "
                  f"(SCHA {1 / a['kappa_scha']:.2f}; bare {a['range_bare']:.2f})  force exponent r=3,6,9: "
                  + ', '.join(f"{a['force_exponent_r2_to_10'][i]:.2f}" for i in (1, 4, 7)), flush=True)
        out = dict(part='mc', m0sq=M0SQ, L=args.L, rows=rows, raw=raw)
    elif args.part == 'hot':
        out = dict(part='hot', **run_hot())
    else:
        out = dict(part='critical', **run_critical())
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
