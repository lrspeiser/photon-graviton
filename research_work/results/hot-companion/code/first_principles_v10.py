"""Round 10: where "the companion pulls with its amplitude" can come from (proposal 1, first step).

The model: a scalar wave field psi (the companion) obeying
    (1/u^2) psi_tt - lap psi = sum_i q_i(t) delta(x - x_i),
with point emitters whose coupling q_i(t) psi(x_i) is the only interaction. A body at x feels
F = q(t) grad psi_other(x, t). Each body is a self-sustained emitter, q(t) = q0 sin(w t + phase),
whose phase locks to the wave passing it, with an offset Delta from the local phase of psi.

A. Exact time-averaged force on a locked emitter in the wave of one point source Q:
       <F> = -(q0 A_psi / 2) [ k sin(Delta) + cos(Delta) / R ]  r_hat,      A_psi = Q / (4 pi R)
   Delta = +pi/2 (the emitter a quarter cycle ahead of the passing wave): a pull toward the source
   equal to (q0 k / 2) x the local amplitude, at every distance; the emitter then feeds the
   passing wave coherently at power P = F v_phase (stimulated emission; its recoil is the pull).
   Delta = 0: a 1/R^2 pull, Newton-like.  Delta = -pi/2: a push (a passive, absorbing response).
B. The same from a 3D finite-difference simulation of the wave equation alone (the force law is
   not programmed in; the emitter only reads the local phase of the field it sits in).
C. Many sources with independent phases (one frequency band): the ensemble-averaged pull on a
   locked emitter, against two rules:
       round 3 ("cold" rule):  (q0 k / 2) sqrt(| sum_j A_j^2 r_j |)           along the net flux
       derived (Gaussian):     (q0 k / 2) (sqrt(pi)/2) sum_j A_j^2 r_j / sqrt(sum_j A_j^2)
   They agree for one source and for any set of sources seen from one side; inside extended
   systems the derived pull is smaller by the anisotropy factor sqrt(|g_N| / S_N), where
   S_N = G int rho / d^2 is the scalar (unsigned) sum of the Newtonian pulls.
D. Bookkeeping: the pull is paid for by coherent emission, P = m g_extra v_phase per body. With
   v_phase = u this would amplify a galaxy's outgoing companion by a factor comparable to one,
   against the tight v^4 = G M a relation, unless v_phase << u.
E. The anisotropy factor for a uniform sphere, the Milky Way model at several radii (with its
   effect on the circular speed, as a local estimate), and a cluster-like distribution.

    python code/first_principles_v10.py --output-dir run-first-principles-v10
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))


# --------------------------------------------------------------------------- A. analytic force
def force_one_source(R, k, Delta, q0=1.0, Q=1.0):
    """Exact <F_r> on an emitter locked at offset Delta to the local phase of psi = Q sin(w(t-R/u))/(4 pi R)."""
    A = Q / (4 * np.pi * R)
    return -(q0 * A / 2) * (k * np.sin(Delta) + np.cos(Delta) / R), A


def power_one_source(R, k, Delta, q0=1.0, Q=1.0, w=1.0):
    """<q psi_t,ext>: the power the emitter feeds into the passing wave (positive = feeds it)."""
    A = Q / (4 * np.pi * R)
    return (q0 * A * w / 2) * np.sin(Delta)


# --------------------------------------------------------------------------- B. 3D simulation
def fdtd(N=192, lam=16.0, dt=0.45, sponge=24, periods_avg=6, ramp_periods=4, sigma_src=1.0, verbose=True):
    """Leapfrog FDTD of (1/u^2) psi_tt - lap psi = s, u = h = 1, with a smoothed oscillating source at the
    centre and damping layers at the edges. Returns probe time series of psi and its radial gradient."""
    u = 1.0; w = 2 * np.pi * u / lam; T = 2 * np.pi / w
    c = N // 2
    x = np.arange(N) - c
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r2 = (X ** 2 + Y ** 2 + Z ** 2).astype(float)
    src = np.exp(-r2 / (2 * sigma_src ** 2)); src /= src.sum()                  # sum * h^3 = 1
    d = np.minimum(np.arange(N), N - 1 - np.arange(N)).astype(float)
    prof = np.clip((sponge - d) / sponge, 0, None) ** 2 * 0.25                    # damping rate
    gam = prof[:, None, None] + prof[None, :, None] + prof[None, None, :]
    del X, Y, Z, r2
    a1 = (1 - gam * dt / 2) / (1 + gam * dt / 2); a2 = 1 / (1 + gam * dt / 2)
    Rx = [4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 56, 64]
    Rd = [6, 12, 24, 40]                                                          # diagonal (1,1,0) offsets
    T_ramp = ramp_periods * T
    t_end = T_ramp + (max(Rx) + 8) / u + periods_avg * T + 2 * T
    nsteps = int(np.ceil(t_end / dt))
    n_avg = int(round(periods_avg * T / dt))
    rec = dict(t=[], psi_x={R: [] for R in Rx}, gr_x={R: [] for R in Rx}, psi_d={R: [] for R in Rd}, gr_d={R: [] for R in Rd})
    prev = np.zeros((N, N, N)); cur = np.zeros((N, N, N))
    t0 = time.monotonic()
    for n in range(nsteps):
        t = n * dt
        ramp = 0.5 * (1 - np.cos(np.pi * min(t, T_ramp) / T_ramp))
        lap = -6 * cur
        lap[1:] += cur[:-1]; lap[:-1] += cur[1:]
        lap[:, 1:] += cur[:, :-1]; lap[:, :-1] += cur[:, 1:]
        lap[:, :, 1:] += cur[:, :, :-1]; lap[:, :, :-1] += cur[:, :, 1:]
        nxt = a2 * (2 * cur + dt * dt * u * u * (lap + ramp * np.sin(w * t) * src)) - a1 * prev
        prev, cur = cur, nxt
        if n >= nsteps - n_avg:
            tt = (n + 1) * dt
            rec['t'].append(tt)
            for R in Rx:
                rec['psi_x'][R].append(cur[c + R, c, c])
                rec['gr_x'][R].append(0.5 * (cur[c + R + 1, c, c] - cur[c + R - 1, c, c]))
            for R in Rd:
                i = c + R
                gx = 0.5 * (cur[i + 1, i, c] - cur[i - 1, i, c]); gy = 0.5 * (cur[i, i + 1, c] - cur[i, i - 1, c])
                rec['psi_d'][R].append(cur[i, i, c]); rec['gr_d'][R].append((gx + gy) / np.sqrt(2))
        if verbose and n % 200 == 0:
            print(f'    step {n}/{nsteps}  ({time.monotonic() - t0:.0f} s)', flush=True)
    # discrete wavenumber along x for this w and dt (leapfrog + 7-point Laplacian)
    kx = 2 * np.arcsin(np.sin(w * dt / 2) / dt)
    return dict(w=w, kx=kx, lam=lam, dt=dt, N=N, Rx=Rx, Rd=Rd, rec=rec, seconds=time.monotonic() - t0)


def locked_force(t, psi, gr, w, Delta, q0=1.0):
    """Fit psi = A sin(w t + phi); the locked emitter q = q0 sin(w t + phi + Delta); return <q gr>, A."""
    t = np.asarray(t); psi = np.asarray(psi); gr = np.asarray(gr)
    M = np.vstack([np.sin(w * t), np.cos(w * t), np.ones_like(t)]).T
    cs, cc, _ = np.linalg.lstsq(M, psi, rcond=None)[0]
    A = np.hypot(cs, cc); phi = np.arctan2(cc, cs)
    q = q0 * np.sin(w * t + phi + Delta)
    return float(np.mean(q * gr)), float(A)


# --------------------------------------------------------------------------- C. many sources
def locked_force_many(src_pos, Qs, probe, k, Delta, alphas, chunk=200):
    """Exact <F> on an emitter at `probe` locked to the total field of point sources with phases alphas
    (shape [nreal, nsrc]); returns [nreal, 3]."""
    d = probe[None, :] - src_pos                        # from source to probe
    r = np.linalg.norm(d, axis=1); rhat = d / r[:, None]
    out = []
    for i in range(0, alphas.shape[0], chunk):
        beta = alphas[i:i + chunk] - k * r[None, :]     # phase of each wave at the probe
        Zc = (Qs / (4 * np.pi * r))[None, :] * np.exp(1j * beta)
        Phi = np.angle(Zc.sum(1))[:, None]
        # <F> = -(q0/(8 pi)) sum_j Q_j [ (k/r_j) sin(Phi + Delta - beta_j) + (1/r_j^2) cos(Phi + Delta - beta_j) ] rhat_j
        x = Phi + Delta - beta
        coef = -(1 / (8 * np.pi)) * Qs[None, :] * ((k / r)[None, :] * np.sin(x) + (1 / r ** 2)[None, :] * np.cos(x))
        out.append(coef @ rhat)
    return np.vstack(out)


def rules(src_pos, Qs, probe, k):
    d = probe[None, :] - src_pos; r = np.linalg.norm(d, axis=1); rhat = d / r[:, None]
    A2 = (Qs / (4 * np.pi * r)) ** 2
    vec = (A2[:, None] * rhat).sum(0)                    # net flux vector (direction: away from the sources)
    scal = A2.sum()
    nv = np.linalg.norm(vec)
    round3 = -(k / 2) * np.sqrt(nv) * vec / nv if nv > 1e-12 * scal else np.zeros(3)
    derived = -(k / 2) * (np.sqrt(np.pi) / 2) * vec / np.sqrt(scal)
    return round3, derived, float(nv / scal)


def many_source_tests(rng, k=1.0, nreal=2000):
    out = []
    def run(name, pos, Qs, probe):
        al = rng.uniform(0, 2 * np.pi, (nreal, len(Qs)))
        F = locked_force_many(pos, Qs, probe, k, np.pi / 2, al)
        Fm = F.mean(0); Fse = F.std(0) / np.sqrt(nreal)
        r3, dv, aniso = rules(pos, Qs, probe, k)
        A1 = np.abs(Qs / (4 * np.pi * np.linalg.norm(probe[None, :] - pos, axis=1))).max()
        if np.linalg.norm(r3) == 0:                      # no net flow: report the mean pull against the largest single amplitude
            out.append(dict(config=name, n_sources=len(Qs), anisotropy_net_over_scalar=aniso,
                            mean_force_magnitude_over_one_source_pull=float(np.linalg.norm(Fm) / (k * A1 / 2)),
                            se_over_one_source_pull=float(np.linalg.norm(Fse) / (k * A1 / 2)),
                            round3_rule=0.0, derived_rule=0.0, ratio_to_round3=None, ratio_to_derived=None))
            return
        e = r3 / np.linalg.norm(r3)                      # compare along the net-flux direction
        out.append(dict(config=name, n_sources=len(Qs), anisotropy_net_over_scalar=aniso,
                        mean_force_along_net=float(Fm @ e), se=float(np.linalg.norm(Fse)),
                        round3_rule=float(r3 @ e), derived_rule=float(dv @ e),
                        ratio_to_round3=float((Fm @ e) / (r3 @ e)), ratio_to_derived=float((Fm @ e) / (dv @ e))))
    # one source (exact: no fluctuation)
    run('one source', np.array([[0.0, 0.0, 0.0]]), np.array([1.0]), np.array([60.0, 0, 0]))
    # two equal sources on opposite sides
    run('two equal, opposite sides', np.array([[-60.0, 0, 0], [60.0, 0, 0]]), np.array([1.0, 1.0]), np.zeros(3))
    # two sources on the same side
    run('two, same side', np.array([[0.0, 0, 0], [10.0, 0, 0]]), np.array([1.0, 0.7]), np.array([80.0, 0, 0]))
    # a clump of 3000 sources seen from outside (a "point mass" made of many parts)
    p = rng.normal(0, 5.0, (3000, 3)); run('3000-part clump seen from outside', p, np.ones(3000) / np.sqrt(3000), np.array([200.0, 0, 0]))
    # probe inside a uniform sphere at half the radius
    n = 6000; R = 400.0
    p = rng.uniform(-1, 1, (4 * n, 3)); p = p[np.linalg.norm(p, axis=1) < 1][:n] * R
    probe = np.array([R / 2, 0, 0]); p = p[np.linalg.norm(p - probe, axis=1) > 3.0]
    run('inside a uniform sphere (r = R/2)', p, np.ones(len(p)) / np.sqrt(len(p)), probe)
    # probe in the plane of an exponential disk (scale 1000, thickness 60), at 1, 2 and 4 scale lengths
    n = 8000; Rd = 1000.0
    Rr = rng.gamma(2.0, Rd, n); ph = rng.uniform(0, 2 * np.pi, n); zz = rng.laplace(0, 60.0, n)
    p = np.vstack([Rr * np.cos(ph), Rr * np.sin(ph), zz]).T
    for m in (1.0, 2.0, 4.0):
        probe = np.array([m * Rd, 0, 0]); keep = np.linalg.norm(p - probe, axis=1) > 3.0
        run(f'in an exponential disk, R = {m:g} scale lengths', p[keep], np.ones(keep.sum()) / np.sqrt(n), probe)
    return out


# --------------------------------------------------------------------------- E. anisotropy factor
def uniform_sphere_factor(x):
    """|g_N|/S_N at radius x R inside a uniform sphere: S_N = G rho int s(Omega) dOmega."""
    mu = np.linspace(-1, 1, 20001)
    s = -x * mu + np.sqrt(1 - x ** 2 * (1 - mu ** 2))
    I = 2 * np.pi * np.trapezoid(s, mu)                 # in units of R
    S = 3 / (4 * np.pi) * I                             # G M / R^2 units
    g = x
    return g / S


def mw_scalar_sum(rho_fn, R0, z0=0.0, s_in=1.0):
    """S_N / G = int rho / d^2 dV for an axisymmetric density rho_fn(R, z) at (R0, z0) (kpc, Msun)."""
    def ring(Rp, zp):
        a = R0 ** 2 + Rp ** 2 + (zp - z0) ** 2; b = 2 * R0 * Rp
        return 2 * np.pi * Rp / np.sqrt(np.clip(a * a - b * b, 1e-300, None))
    s_in = min(s_in, 0.5 * R0)
    # inner: polar coordinates about (R0, z0); dR dz = s ds dchi
    s = np.geomspace(1e-6, s_in, 700); chi = np.linspace(0, 2 * np.pi, 721)[:-1] + np.pi / 720
    ss, cc = np.meshgrid(s, chi, indexing='ij')
    Rp = R0 + ss * np.cos(cc); zp = z0 + ss * np.sin(cc)
    f = rho_fn(Rp, zp) * ring(Rp, zp) * ss
    inner = np.trapezoid(f.sum(1) * (2 * np.pi / chi.size), s)
    # outer: sinh-spaced grid in R and z, excluding the inner disc
    Re = 0.02 * np.sinh(np.linspace(0, np.arcsinh(300 / 0.02), 1601)); zs = 0.005 * np.sinh(np.linspace(0, np.arcsinh(150 / 0.005), 801))
    ze = np.concatenate([-zs[::-1], zs[1:]])
    Rc = 0.5 * (Re[1:] + Re[:-1]); zc = 0.5 * (ze[1:] + ze[:-1]); dR = np.diff(Re); dz = np.diff(ze)
    RR, ZZ = np.meshgrid(Rc, zc, indexing='ij')
    mask = np.hypot(RR - R0, ZZ - z0) > s_in
    outer = np.sum(np.where(mask, rho_fn(RR, ZZ) * ring(RR, ZZ), 0.0) * dR[:, None] * dz[None, :])
    return float(inner + outer)


def milky_way_factor(Rs=(4.0, 8.2, 12.0, 20.0, 30.0)):
    """The anisotropy factor in the Milky Way model the suite uses (McMillan 2017 visible components),
    and a local estimate of what it does to the circular speed with the round-3 constants."""
    import common as C
    import t_milky_way as TMW
    import milky_way_v7 as MW
    import mw_model as M
    from law_config import load_law
    ctx = C.Context(tier='quick', verbose=False)
    comps, grid, F = TMW.setup(ctx)
    fac = MW.models(comps)['M17']
    rho = lambda R, z: sum(f * comps[n].rho(R, z) for n, f in fac.items())
    law = load_law('round9')
    a, gd = law['a_code'], law['lam'] * law['a_code']      # (km/s)^2/kpc
    G = M.G
    # Newtonian pull in the plane from the cached component fields
    gR = sum(f * F[n][0] for n, f in fac.items())
    iz = np.argmin(np.abs(grid.z))
    rows = []
    for R0 in Rs:
        gN = float(abs(np.interp(R0, grid.R, gR[:, iz])))
        S = G * mw_scalar_sum(rho, R0)
        f2 = gN / S
        rel = np.exp(-gN / gd)
        ex3 = rel * np.sqrt(a * gN)                        # round 3, cold, local algebraic estimate
        exd = rel * np.sqrt(a) * gN / np.sqrt(S)           # derived rule
        v3 = np.sqrt(R0 * (gN + ex3)); vd = np.sqrt(R0 * (gN + exd))
        rows.append(dict(R_kpc=R0, g_N=gN, S_N=S, ratio_gN_over_SN=f2, anisotropy_factor=float(np.sqrt(f2)),
                         v_newton=float(np.sqrt(R0 * gN)), v_round3_local=float(v3), v_derived_local=float(vd)))
    return rows


def cluster_factor():
    """A cluster-like mixture: gas in a beta-model (r_c = 200 kpc) and galaxies in an NFW-like profile
    (r_s = 300 kpc) with 5% of the gas mass but heat weight k = 3 sigma^2/u^2 for sigma = 1000 km/s.
    Returns the derived/round-3 ratio of the extra pull at several radii (all radial, so the net
    vectors add), using spherical shell weights for the scalar sums."""
    u = 197.41; k = 3 * 1000.0 ** 2 / u ** 2
    r = np.geomspace(1.0, 5000.0, 4000)
    rho_g = (1 + (r / 200.0) ** 2) ** -1.5 * (r < 3000)
    rho_s = 1 / ((r / 300.0) * (1 + r / 300.0) ** 2) * (r < 3000)
    dV = 4 * np.pi * r ** 2 * np.gradient(r)
    Mg = rho_g * dV; Ms = rho_s * dV
    Ms *= 0.05 * Mg.sum() / Ms.sum()
    rows = []
    for x in (100.0, 250.0, 500.0, 1000.0, 1500.0):
        inside = r < x
        gN = (Mg[inside].sum() + Ms[inside].sum()) / x ** 2      # G = 1
        gh = k * Ms[inside].sum() / x ** 2
        # shell weight: int over a shell of radius r' of 1/d^2 dOmega / (4 pi) = (1/(2 r r')) ln |(r + r')/(r - r')|
        w = (1 / (2 * x * r)) * np.log(np.abs((x + r) / np.clip(np.abs(x - r), 1e-9, None)))
        SN = np.sum((Mg + Ms) * w); SH = k * np.sum(Ms * w)
        round3 = np.sqrt(gN + SH) * (gN + gh) / (gN + gh)          # direction factor = 1 (aligned)
        derived = (gN + gh) / np.sqrt(SN + SH)
        rows.append(dict(r_kpc=x, derived_over_round3=float(derived / round3),
                         vector_over_scalar=float((gN + gh) / (SN + SH))))
    return rows


def xcop_test():
    """The derived combination rule on the 12 X-COP clusters (hydrostatic masses at six radii), with a and
    g_d held at their SPARC values and u refitted, against the round-3 rule refitted the same way."""
    from scipy.optimize import minimize_scalar
    import common as C
    import run_v3 as R3
    import law as L
    from law_config import load_law
    ctx = C.Context(tier='quick', verbose=False)
    cls = ctx.xcop()
    law = load_law('round9'); a, lam = law['a_code'], law['lam']
    G = L.G

    def pred_derived(c, u):
        k = L.heat_weight(np.sqrt(c['sig2_star_hse']), u)
        gN = G * c['Mb'] / c['Rk'] ** 2
        S = G * (c['W'] @ (k * c['dms'])) / c['Rk'] ** 2
        SN = G * (c['W'] @ (c['dms'] + c['dmg'])) / c['Rk'] ** 2
        ghot = G * np.array([np.sum((k * c['dms'])[c['s'] < R]) for R in c['Rk']]) / c['Rk'] ** 2
        g = gN + L.released(gN, a, lam) * np.sqrt(a) * (gN + ghot) / np.sqrt(SN + S)
        return g * c['Rk'] ** 2 / G

    def fit(pred):
        f = lambda lu: np.mean(R3.resid(cls, lambda c: pred(c, 10 ** lu)) ** 2)
        r = minimize_scalar(f, bounds=(0.5, 3.5), method='bounded')
        u = 10 ** r.x
        res = R3.resid(cls, lambda c: pred(c, u))                  # [cluster, radius]
        return u, res

    out = {}
    for name, pred in (('round3', lambda c, u: R3.cluster_M3(c, a, u, lam)), ('derived', pred_derived)):
        u, res = fit(pred)
        out[name] = dict(u_best=float(u), rms=float(np.sqrt(np.mean(res ** 2))),
                         mean_residual_by_radius=[float(x) for x in res.mean(0)],
                         inner_minus_outer=float(res[:, :2].mean() - res[:, -2:].mean()))
    # the derived rule at the round-3 companion speed, for scale
    res = R3.resid(cls, lambda c: pred_derived(c, law['u_kms']))
    out['derived_at_u_197'] = dict(rms=float(np.sqrt(np.mean(res ** 2))), mean_residual=float(res.mean()))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--grid', type=int, default=192)
    args = ap.parse_args(); out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic(); rng = np.random.default_rng(20260924)
    res = dict(experiment='round 10: the pull of a locked emitter, from one local coupling')

    # A
    k = 1.0
    A_rows = []
    for R in (0.5, 1, 2, 5, 10, 30, 100):
        row = dict(kR=R)
        for name, D in (('quarter_ahead', np.pi / 2), ('in_step', 0.0), ('quarter_behind', -np.pi / 2)):
            F, A = force_one_source(R, k, D)
            row[name + '_force_over_amplitude'] = float(-F / A)          # positive = pull toward the source
        row['quarter_ahead_power_over_force_times_u'] = float(power_one_source(R, k, np.pi / 2) / (-force_one_source(R, k, np.pi / 2)[0] * 1.0))
        A_rows.append(row)
    res['A_analytic'] = dict(formula='<F_r> = -(q0 A/2) [k sin Delta + cos Delta / R], A = Q/(4 pi R)', rows=A_rows)
    print('A. analytic: pull / local amplitude (q0 = 1, k = 1): quarter ahead = k/2 = 0.5 at every R')
    for r in A_rows:
        print(f"   kR {r['kR']:6.1f}: quarter ahead {r['quarter_ahead_force_over_amplitude']:.4f}   in step {r['in_step_force_over_amplitude']:.4f}   quarter behind {r['quarter_behind_force_over_amplitude']:.4f}")

    # B
    print('B. 3D wave simulation ...', flush=True)
    sim = fdtd(N=args.grid)
    w, kx = sim['w'], sim['kx']
    B_rows = []
    for tag, Rs, P, Gr in (('x axis', sim['Rx'], sim['rec']['psi_x'], sim['rec']['gr_x']), ('diagonal', sim['Rd'], sim['rec']['psi_d'], sim['rec']['gr_d'])):
        for R in Rs:
            Rr = R * (np.sqrt(2) if tag == 'diagonal' else 1.0)
            Fq, A = locked_force(sim['rec']['t'], P[R], Gr[R], w, np.pi / 2)
            F0, _ = locked_force(sim['rec']['t'], P[R], Gr[R], w, 0.0)
            Fb, _ = locked_force(sim['rec']['t'], P[R], Gr[R], w, -np.pi / 2)
            th_q = -(A / 2) * kx; th_0 = -(A / 2) / Rr
            B_rows.append(dict(direction=tag, R_cells=float(Rr), kR=float(kx * Rr), amplitude=A, amplitude_times_R=A * Rr,
                               pull_quarter_ahead=-Fq, theory_quarter_ahead=-th_q, ratio_quarter_ahead=Fq / th_q,
                               pull_in_step=-F0, theory_in_step=-th_0, ratio_in_step=F0 / th_0, push_quarter_behind=Fb))
    res['B_simulation'] = dict(grid=sim['N'], wavelength_cells=sim['lam'], dt=sim['dt'], k_discrete_x=kx, seconds=sim['seconds'], rows=B_rows)
    print('   R (cells)  kR    A*R      pull(ahead)/theory   pull(in step)/theory   push(behind)')
    for r in B_rows:
        print(f"   {r['direction'][:4]} {r['R_cells']:6.1f} {r['kR']:5.1f} {r['amplitude_times_R']:.5f}   {r['ratio_quarter_ahead']:.4f}              {r['ratio_in_step']:.4f}              {r['push_quarter_behind']:+.2e}")

    # C
    print('C. many sources with independent phases', flush=True)
    C_rows = many_source_tests(rng)
    res['C_many_sources'] = C_rows
    for r in C_rows:
        if r['ratio_to_round3'] is None:
            print(f"   {r['config']:44s} net flow zero; mean pull / one source's pull {r['mean_force_magnitude_over_one_source_pull']:.4f} (+- {r['se_over_one_source_pull']:.4f})")
        else:
            print(f"   {r['config']:44s} net/scalar {r['anisotropy_net_over_scalar']:.3f}  pull/round3 {r['ratio_to_round3']:.3f}  pull/derived {r['ratio_to_derived']:.3f}  (+- {r['se'] / abs(r['derived_rule']):.3f})")

    # D. bookkeeping, with the round-9 constants
    a_SI, u_SI = 6.561190024871433e-11, 197.41e3
    D_rows = []
    for gN in (1e-12, 1e-11, 1e-10, 3e-10):
        gx = np.exp(-gN / 2.262e-10) * np.sqrt(a_SI * gN)
        D_rows.append(dict(g_N=gN, g_extra=float(gx), stimulated_over_spontaneous_if_vphase_eq_u=float(2 * gx / a_SI)))
    res['D_bookkeeping'] = dict(note='power per kg fed coherently = g_extra * v_phase; spontaneous = l = a u / 2; ratio = 2 g_extra v_phase / (a u)',
                                rows=D_rows)
    print('D. coherent power needed for the pull, over the spontaneous feed l (if v_phase = u):',
          ', '.join(f"g_N {r['g_N']:.0e}: {r['stimulated_over_spontaneous_if_vphase_eq_u']:.2f}" for r in D_rows))

    # E
    E_sphere = [dict(x=x, ratio=float(uniform_sphere_factor(x)), anisotropy_factor=float(np.sqrt(uniform_sphere_factor(x)))) for x in (0.25, 0.5, 0.75, 0.99)]
    print('E. anisotropy factor sqrt(|g_N|/S_N): uniform sphere', ', '.join(f"r/R {r['x']}: {r['anisotropy_factor']:.3f}" for r in E_sphere))
    E_mw = milky_way_factor()
    for r in E_mw:
        print(f"   Milky Way R {r['R_kpc']:5.1f} kpc: |g_N|/S_N {r['ratio_gN_over_SN']:.3f}  factor {r['anisotropy_factor']:.3f}  v (local estimate) round 3 {r['v_round3_local']:.1f}, derived {r['v_derived_local']:.1f} km/s (Newton {r['v_newton']:.1f})")
    E_cl = cluster_factor()
    print('   cluster-like mix: derived/round-3 extra pull', ', '.join(f"{r['r_kpc']:.0f} kpc: {r['derived_over_round3']:.2f}" for r in E_cl))
    res['E_anisotropy'] = dict(uniform_sphere=E_sphere, milky_way=E_mw, cluster_like=E_cl)

    # F
    F = xcop_test()
    res['F_xcop'] = F
    for n in ('round3', 'derived'):
        r = F[n]
        print(f"F. X-COP, {n:8s}: best u {r['u_best']:6.1f} km/s  rms {r['rms']:.3f}  mean residual by radius " + ' '.join(f'{x:+.2f}' for x in r['mean_residual_by_radius']))
    print(f"   derived rule at u = 197: rms {F['derived_at_u_197']['rms']:.3f}, mean residual {F['derived_at_u_197']['mean_residual']:+.2f} (ln M_obs/M_pred)")
    res['seconds'] = time.monotonic() - t0
    (out / 'first_principles_v10.json').write_text(json.dumps(res, indent=1))
    print(f"wrote {out / 'first_principles_v10.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
