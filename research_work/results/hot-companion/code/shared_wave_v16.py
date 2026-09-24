"""Round 16, step 3a: the pieces' shared wave with its energy budget closed exactly.

An independent audit of round 15 (README §25.4) found that code/reservoir_force_v15.py lets each piece of matter lose
energy at its own rate while the pieces' waves are added coherently. Interference changes the power the waves carry
(0.60-1.40 times the modes' loss, depending on the arrangement), and that change was not fed back into the pieces.
The audit's repair: derive the damping from the same coupling that radiates,
    dz/dt = -i H z - (1/2) W^dagger W z,        so  d|z|^2/dt = -P_out  exactly,
with W the map from the internal amplitudes to the outgoing waves. In its reduced test (20 pieces, no test bodies,
the motion-to-detuning rule still assumed) the extra radiation per doubling of the detuning fell from 4.13, 4.05,
4.00, 3.91 to 3.80, 3.28, 2.02, 1.64. This script builds that repair with the full wave coupling and asks when the
sigma^2 survives.

The coupling (companion wavelength 1, wave speed 1, k = omega = 2 pi; G(R) = e^{ikR}/(4 pi R), the exact Green's
function, near and far field): piece j at x_j has four internal amplitudes z_j = (D_j, B_jx, B_jy, B_jz). D radiates
as a monopole q = C0 D, the B's as a dipole p = C1 B, with C0^2 = 2 gamma_0 8 pi/(omega k), C1^2 = 2 gamma 24 pi/
(omega k^3), so an isolated piece radiates 2 gamma_0 |D|^2 + 2 gamma |B|^2 (round 15's normalization). Eliminating the
wave gives one complex symmetric matrix M = (omega/2) C Gfun C:
    M[jD, lD]   = (omega/2) C0^2  G(R)            R = x_j - x_l
    M[jD, lBm]  = (omega/2) C0 C1 (-d_m G)(R)
    M[jBa, lD]  = (omega/2) C1 C0 (d_a G)(R)
    M[jBa, lBb] = (omega/2) C1^2  (-d_a d_b G)(R)
and on the diagonal only the finite radiative part, i 2 gamma_0 and i 2 gamma (the divergent reactive self-part only
renames the frequency). Gamma = Im M is exactly W^dagger W: z^dagger Gamma z is the power through a distant sphere
(checked in validate()). J = Re M is the reactive near field, which moves energy between pieces but radiates none.
A passive piece obeys
    dz/dt = -i H z + (i/2) M z = -i (H - J/2) z - (Gamma/2) z,
H holding the motion's mixing of D into the B's (the upload's postulate, delta = chi w, kept as an input) and any
spread of the pieces' own frequencies. The independent-damping version of round 15 replaces M by its diagonal.

Stage A (this script): the source alone. Pieces start with unit D at random phases and B = 0; for free motion the
evolution is the exact matrix exponential, for collisions (delta an Ornstein-Uhlenbeck process of rate nu) a split
step: the collective part exactly, the mixing as a unitary rotation. Reported at fixed times: the power radiated per
unit of energy still stored, P/E, and its ratio to the same arrangement at rest; the extra (ratio - 1) per doubling
of the detuning. Clouds: dense (radius 0.5 and 1 wavelength) to dilute (3 wavelengths), and the pieces' own
frequencies either identical or spread (rms Delta0).

    python code/shared_wave_v16.py --output run-shared-wave-v16/shared_wave_v16.json [--processes 4] [--quick]
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.linalg import expm

K = 2 * np.pi
OMEGA = K
IM_D2G0 = K ** 3 / (12 * np.pi)          # -d_a d_b Im G at R = 0 is IM_D2G0 delta_ab


def cmono(g0):
    return np.sqrt(2 * g0 * 8 * np.pi / (OMEGA * K))


def cdip(g):
    return np.sqrt(2 * g * 24 * np.pi / (OMEGA * K ** 3))


def ball_min_sep(n, R, sep, rng, tries=200000):
    pts = []
    for _ in range(tries):
        p = rng.uniform(-R, R, 3)
        if p @ p > R * R:
            continue
        if all(np.sum((p - q) ** 2) >= sep * sep for q in pts):
            pts.append(p)
            if len(pts) == n:
                return np.array(pts)
    raise RuntimeError('could not place the pieces')


def green_derivs(R):
    """G, grad G, Hessian G, third derivatives of G = e^{ikr}/(4 pi r) at separation vectors R (..., 3), r > 0."""
    r = np.linalg.norm(R, axis=-1); n = R / r[..., None]
    g = np.exp(1j * K * r) / (4 * np.pi * r)
    a = 1j * K - 1 / r
    g1 = a * g
    g2 = (a * a + 1 / r ** 2) * g
    g3 = (2 * a / r ** 2 - 2 / r ** 3) * g + (a * a + 1 / r ** 2) * a * g
    I3 = np.eye(3)
    nn = n[..., :, None] * n[..., None, :]
    G1 = g1[..., None] * n
    G2 = g2[..., None, None] * nn + (g1 / r)[..., None, None] * (I3 - nn)
    nnn = nn[..., :, :, None] * n[..., None, None, :]
    sym = (I3[:, :, None] * n[..., None, None, :] + I3[:, None, :] * n[..., None, :, None]
           + I3[None, :, :] * n[..., :, None, None])
    G3 = (g3 - 3 * g2 / r + 3 * g1 / r ** 2)[..., None, None, None] * nnn + (g2 / r - g1 / r ** 2)[..., None, None, None] * sym
    return g, G1, G2, G3


def coupling(x, g0, g):
    """The complex symmetric coupling M (4N x 4N) and its derivative dM (3, 4N, 4N) with respect to the first piece's
    position (the field point), including the finite radiative self-parts. Block order per piece: D, Bx, By, Bz."""
    N = len(x); c0, c1 = cmono(g0), cdip(g)
    R = x[:, None, :] - x[None, :, :]
    off = ~np.eye(N, dtype=bool)
    Rs = np.where(off[..., None], R, 1.0)                          # placeholder on the diagonal
    G0, G1, G2, G3 = green_derivs(Rs)
    pre = OMEGA / 2
    M = np.zeros((N, 4, N, 4), complex); dM = np.zeros((3, N, 4, N, 4), complex)
    M[:, 0, :, 0] = pre * c0 * c0 * G0
    M[:, 0, :, 1:] = pre * c0 * c1 * (-G1)
    M[:, 1:, :, 0] = np.transpose(pre * c1 * c0 * G1, (0, 2, 1))
    M[:, 1:, :, 1:] = np.transpose(pre * c1 * c1 * (-G2), (0, 2, 1, 3))
    for c in range(3):
        dM[c, :, 0, :, 0] = pre * c0 * c0 * G1[..., c]
        dM[c, :, 0, :, 1:] = pre * c0 * c1 * (-G2[..., c, :])
        dM[c, :, 1:, :, 0] = np.transpose(pre * c1 * c0 * G2[..., c, :], (0, 2, 1))
        dM[c, :, 1:, :, 1:] = np.transpose(pre * c1 * c1 * (-G3[..., c, :, :]), (0, 2, 1, 3))
    idx = np.arange(N)
    M[idx, :, idx, :] = 0.0; dM[:, idx, :, idx, :] = 0.0
    M[idx, 0, idx, 0] = 2j * g0
    for m in range(1, 4):
        M[idx, m, idx, m] = 2j * g
        # the radiative self-force of a piece whose monopole and dipole are coherent: d_c of the D-B blocks at R = 0
        dM[m - 1, idx, 0, idx, m] = pre * c0 * c1 * 1j * IM_D2G0
        dM[m - 1, idx, m, idx, 0] = -pre * c1 * c0 * 1j * IM_D2G0
    return M.reshape(4 * N, 4 * N), dM.reshape(3, 4 * N, 4 * N)


def far_field(x, z, g0, g, ndir=4000):
    """Power and momentum radiated to infinity per unit time, from the far-field pattern on a Fibonacci sphere."""
    i = np.arange(ndir) + 0.5; ph = np.arccos(1 - 2 * i / ndir); th = np.pi * (1 + 5 ** 0.5) * i
    nh = np.vstack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)]).T
    zz = z.reshape(-1, 4); q = cmono(g0) * zz[:, 0]; p = cdip(g) * zz[:, 1:]
    A = (np.exp(-1j * K * nh @ x.T) * (q[None, :] - 1j * K * (nh @ p.T))).sum(1)
    dP = OMEGA * K / (32 * np.pi ** 2) * np.abs(A) ** 2
    return float(dP.mean() * 4 * np.pi), (dP[:, None] * nh).mean(0) * 4 * np.pi


def forces(dM, z):
    """Time-averaged force on each piece: F_j = (1/omega) Re sum_l z_j^dagger dM_jl z_l."""
    N = len(z) // 4
    v = np.einsum('cab,b->ca', dM, z)
    return (np.real(np.conj(z)[None, :] * v).reshape(3, N, 4).sum(2) / OMEGA).T


def validate():
    rng = np.random.default_rng(5)
    out = {}
    g0, g = 0.3, 1.0
    for lab, n, R in (('one piece', 1, 0.0), ('two pieces a tenth of a wavelength apart', 2, None), ('12 pieces, radius 0.6', 12, 0.6)):
        if n == 1:
            x = np.zeros((1, 3))
        elif R is None:
            x = np.array([[0, 0, 0], [0.1, 0, 0]], float)
        else:
            x = ball_min_sep(n, R, 0.15, rng)
        M, dM = coupling(x, g0, g)
        rows = []
        for _ in range(3):
            z = rng.normal(size=4 * n) + 1j * rng.normal(size=4 * n)
            P = float(np.real(np.conj(z) @ (M.imag @ z)))
            Pf, mom = far_field(x, z, g0, g)
            F = forces(dM, z).sum(0)
            rows.append(dict(power_from_M=P, power_far_field=Pf, total_force=F.tolist(), momentum_radiated=mom.tolist(),
                             momentum_balance=float(np.linalg.norm(F + mom) / max(np.linalg.norm(mom), 1e-30))))
        out[lab] = rows
    # the audit's two-emitter example: in phase and in opposite phase, a tenth of a wavelength apart, monopoles only
    x = np.array([[0, 0, 0], [0.1, 0, 0]], float); M, _ = coupling(x, 1.0, 1.0)
    for lab, s in (('same phase', 1.0), ('opposite phase', -1.0)):
        z = np.zeros(8, complex); z[0] = 1; z[4] = s
        out['two monopoles, ' + lab] = float(np.real(np.conj(z) @ (M.imag @ z)) / (2 * 1.0 * 2))
    return out


def variant(M, lab):
    """'shared': the full coupling; 'dissipative': Gamma only, no reactive near field (the audit's reduced repair);
    'independent': round 15's own-rate damping."""
    if lab == 'independent':
        return np.diag(np.diag(M))
    if lab == 'dissipative':
        return 1j * M.imag
    return M


def evolve_free(M, H, z0, times, lab='shared'):
    """Exact solution for constant H: z(t) = exp(A t) z0, A = -i H + (i/2) M."""
    A = -1j * H + 0.5j * variant(M, lab)
    w, V = np.linalg.eig(A)
    c = np.linalg.solve(V, z0)
    return [V @ (np.exp(w * t) * c) for t in times]


def mixing_H(dl, Delta):
    N = len(dl); H = np.zeros((N, 4, N, 4))
    idx = np.arange(N)
    H[idx, 0, idx, 1:] = dl; H[idx, 1:, idx, 0] = dl
    for m in range(4):
        H[idx, m, idx, m] = Delta
    return H.reshape(4 * N, 4 * N)


def run(cfg):
    rng = np.random.default_rng(cfg['seed'])
    N, Rb, g0, g = cfg['N'], cfg['Rb'], cfg['gamma0'], 1.0
    x = ball_min_sep(N, Rb, cfg.get('sep', 0.15), rng)
    M, _ = coupling(x, g0, g)
    Gam = M.imag
    z0 = np.zeros(4 * N, complex); z0[0::4] = np.exp(1j * rng.uniform(0, 2 * np.pi, N))
    e = rng.normal(size=(N, 3))                                    # the same directions for every q: paired runs
    Delta = cfg.get('Delta0', 0.0) * rng.normal(size=N)
    times = cfg['times']
    out = dict(cfg=cfg, rows=[])
    for q in cfg['qs']:
        dl = q * g * e
        H = mixing_H(dl, Delta)
        row = dict(q=q, k=float(3 * q * q * g / g0))
        for lab in ('shared', 'dissipative', 'independent'):
            zs = evolve_free(M, H, z0, times, lab)
            Gx = variant(M, lab).imag
            Pt = [float(np.real(np.conj(z) @ (Gx @ z))) for z in zs]
            Et = [float(np.real(np.conj(z) @ z)) for z in zs]
            row[lab] = dict(P=Pt, E=Et, P_over_E=[p / e_ for p, e_ in zip(Pt, Et)])
        out['rows'].append(row)
    if cfg.get('collisions'):
        out['collisions'] = [collide(M, x, z0, e, Delta, cfg, q, nu, rng_seed=cfg['seed'] * 100 + i)
                             for i, (q, nu) in enumerate(cfg['collisions'])]
    if cfg.get('check_energy'):
        # the budget: stored energy lost = integral of z^dagger Gamma z, by quadrature on a fine grid
        H = mixing_H(cfg['qs'][-1] * e, Delta)
        tt = np.linspace(0, times[-1], 4001)
        zs = evolve_free(M, H, z0, tt)
        P = np.array([np.real(np.conj(z) @ (Gam @ z)) for z in zs]); E = np.array([np.real(np.conj(z) @ z) for z in zs])
        rad = np.trapezoid(P, tt) if hasattr(np, 'trapezoid') else np.trapz(P, tt)
        Pf, _ = far_field(x, zs[len(zs) // 2], g0, g)
        out['energy_check'] = dict(stored_lost=float(E[0] - E[-1]), radiated=float(rad),
                                   relative_residual=float(abs(E[0] - E[-1] - rad) / (E[0] - E[-1])),
                                   flux_far_field_over_M=float(Pf / P[len(P) // 2]))
    return out


def collide(M, x, z0, e, Delta, cfg, q, nu, rng_seed):
    """Collisions: the detunings an Ornstein-Uhlenbeck process (rms q, rate nu). Split step: exp(A0 dt/2) U_mix exp(A0
    dt/2), A0 the collective part with the pieces' own frequencies, U_mix the unitary mixing rotation of each piece."""
    rng = np.random.default_rng(rng_seed)
    N = len(x); g = 1.0
    dt = min(0.05, 0.1 / nu) if nu > 0 else 0.05
    A0 = 0.5j * M - 1j * np.diag(np.repeat(Delta, 4))
    half = expm(A0 * dt / 2)
    Gam = M.imag
    dl = q * g * rng.normal(size=(N, 3))
    rho = np.exp(-nu * dt); kick = q * g * np.sqrt(1 - rho * rho)
    z = z0.copy(); t = 0.0; res = []; T = cfg['times'][-1]
    marks = list(cfg['times'])
    while marks:
        if t >= marks[0] - 1e-9:
            P = float(np.real(np.conj(z) @ (Gam @ z))); E = float(np.real(np.conj(z) @ z))
            res.append(dict(t=marks[0], P=P, E=E, P_over_E=P / E)); marks.pop(0); continue
        z = half @ z
        zz = z.reshape(N, 4)
        a = np.linalg.norm(dl, axis=1); u = dl / np.where(a > 0, a, 1)[:, None]
        Bpar = (zz[:, 1:] * u).sum(1)
        c, s = np.cos(a * dt), np.sin(a * dt)
        D_new = c * zz[:, 0] - 1j * s * Bpar
        Bpar_new = c * Bpar - 1j * s * zz[:, 0]
        zz[:, 1:] += (Bpar_new - Bpar)[:, None] * u; zz[:, 0] = D_new
        z = half @ zz.ravel()
        dl = rho * dl + kick * rng.normal(size=(N, 3)); t += dt
    return dict(q=q, nu=nu, k=float(3 * q * q / cfg['gamma0']), marks=res)


def configs(quick=False):
    qs = [0.0, 0.00025, 0.0005, 0.001, 0.002, 0.004, 0.008, 0.016]
    times = [50.0, 200.0, 1000.0, 3000.0]
    seeds = (1,) if quick else (1, 2, 3, 4, 5)
    runs = []
    for sd in seeds:
        for N, Rb in ((20, 0.5), (20, 1.0), (20, 3.0), (100, 1.0)):
            for D0 in ((0.0, 0.3) if (N, Rb) != (100, 1.0) else (0.0,)):
                cfg = dict(N=N, Rb=Rb, gamma0=1e-6, qs=qs, times=times, seed=sd, Delta0=D0,
                           tag=f'N{N}_R{Rb}_D{D0}')
                if N == 20 and Rb == 1.0:
                    cfg['check_energy'] = True
                if N == 20 and not quick:
                    cfg['collisions'] = [(0.002, 0.0), (0.002, 1.0), (0.002, 10.0)]
                runs.append(cfg)
    return runs


def summarize(res):
    """Mean over seeds, per cloud, variant and time: the extra radiation per unit stored energy, (P/E)_q/(P/E)_0 - 1,
    its ratio per doubling of q (mean of the per-seed ratios, and its standard error), and the cold P/E in units of the
    isolated piece's 2 gamma_0."""
    out = {}
    tags = sorted({r['cfg']['tag'] for r in res})
    for tag in tags:
        rs = [r for r in res if r['cfg']['tag'] == tag]
        c = rs[0]['cfg']; qs = c['qs'][1:]
        out[tag] = {}
        for lab in ('shared', 'dissipative', 'independent'):
            for ti, t in enumerate(c['times']):
                ex = np.array([[row[lab]['P_over_E'][ti] / r['rows'][0][lab]['P_over_E'][ti] - 1 for row in r['rows'][1:]] for r in rs])
                rat = ex[:, 1:] / ex[:, :-1]
                cold = np.array([r['rows'][0][lab]['P_over_E'][ti] / (2 * c['gamma0']) for r in rs])
                out[tag][f'{lab} t={t:g}'] = dict(q=qs, extra=ex.mean(0).tolist(), per_doubling=rat.mean(0).tolist(),
                                                   per_doubling_se=(rat.std(0, ddof=1) / np.sqrt(len(rs))).tolist() if len(rs) > 1 else None,
                                                   cold_over_isolated=float(cold.mean()), seeds=len(rs))
        if 'collisions' in rs[0]:
            rows = []
            for i, (q, nu) in enumerate(c['collisions']):
                base = np.array([[m['P_over_E'] for m in r['collisions'][i]['marks']] for r in rs])
                rows.append(dict(q=q, nu=nu, P_over_E=base.mean(0).tolist()))
            out[tag]['collisions'] = rows
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4); ap.add_argument('--quick', action='store_true')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    val = validate()
    print('validation:', json.dumps(val)[:1500], flush=True)
    runs = configs(args.quick)
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(run, runs):
            res.append(r); c = r['cfg']
            base = r['rows'][0]
            line = []
            for row in r['rows'][1:]:
                sh = row['shared']['P_over_E'][1] / base['shared']['P_over_E'][1] - 1
                ind = row['independent']['P_over_E'][1] / base['independent']['P_over_E'][1] - 1
                line.append(f"{row['q']:.5f}: {sh:.3g}/{ind:.3g}")
            print(f"[{time.monotonic() - t0:5.0f} s] {c['tag']} seed {c['seed']}: extra at t=200 shared/independent  " + '  '.join(line), flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 16 step 3a: the shared wave with its energy budget closed',
                                           validation=val, summary=summarize(res), runs=res, seconds=time.monotonic() - t0),
                                      indent=1, default=float) + '\n')
    print('wrote', args.output, f'({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
