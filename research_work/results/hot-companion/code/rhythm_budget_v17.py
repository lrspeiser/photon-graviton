"""Round 17, step A: what puts a warm source out of tune? The rhythm budget of code/one_matter_v16.py's matter.

Round 16 (README §26.4) found that in one kind of matter, every piece both sender and receiver with every wave fed
back, a warming source glows as the heat rule says and its wave reaches distant matter as the square root of the glow,
but the distant pieces fall out of step, because warming spreads the source's own rhythms (4.7e-4 at rest, 3.9e-3 at
k = 16). A cold source put out of tune by hand loses them the same way. An independent review asked for the next step
to identify exactly which term of the energy-balanced, all-mutual equations shifts the warm pieces' rhythms.

The reduction. Each piece's radiators B relax at gamma + gi = 5, a thousand times faster than the rhythms change, so
they can be eliminated exactly in the linear response (the adiabatic limit of the same equations):
    0 = (i/2) [M_Bs s + M_BB B - 2 D s] - gi B       (D: the 3N x N block of the pieces' mixings delta_j)
    B = T s,   T = -(i/2) Y^-1 X,   Y = (i/2) M_BB - gi I,   X = M_Bs - 2 D
and the quiet oscillations obey, with every piece at its steady inversion w0 and amplitude a,
    ds_j/dt = -(i/2) w0 sum_l C_jl s_l + (real gain and loss),     C = M_ss - (i/2) X^T Y^-1 X   (complex symmetric)
so the rhythm of piece j is
    dphi_j/dt = omega_j - (w0/2) sum_{l != j} |C_jl| cos(phi_l - phi_j + arg C_jl),     omega_j = -(w0/2) Re C_jj.
omega_j is each piece's own rhythm offset: for an isolated piece C_jj = 2 i gamma_0 + 2 i |delta|^2/(gamma + gi) is
purely imaginary (no offset); in a cloud the neighbours' radiators make its own radiators' response complex, and a
warm piece's offset grows as |delta_j|^2. The off-diagonal C_jl are how pieces tug each other: directly (M_ss) and
through each other's radiators.

This script (1) computes omega_j and C for the main run's arrangements, cold and warm; (2) runs the full model again
recording each piece's rhythm, to test the budget against it; (3) runs the reduced rhythm model with parts switched off,
to find which part keeps a warm source from settling and its receivers from keeping step.

(4) splits C exactly by order in the mixing, C = C0 + C1 + C2 (Y does not depend on delta), for the mixing round 16
used (parity 'even') and for the one a velocity must give ('odd', see code/rhythm_protect_v17.py).

    python code/rhythm_budget_v17.py --output run-rhythm-budget-v17/rhythm_budget_v17.json [--processes 4]
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import one_matter_v16 as om                                        # noqa: E402
from shared_wave_v16 import ball_min_sep, coupling               # noqa: E402

MAIN = dict(Ns=48, Nr=8, Rb=3.0, rp=6.0, gamma0=0.01, gi=4.0, piece=dict(G=4.0, W0=1.6, g_par=0.01))


def arrangement(cfg, seed, q):
    """The main run's arrangement for this seed: positions, the pieces' mixings, the coupling (same draws as
    one_matter_v16.simulate)."""
    rng = np.random.default_rng(seed)
    Ns, Nr = cfg['Ns'], cfg['Nr']
    xs = ball_min_sep(Ns, cfg['Rb'], 0.15, rng); xr = om.sphere_points(Nr, cfg['rp']); x = np.vstack([xs, xr])
    N = Ns + Nr
    rng.uniform(0, 2 * np.pi, N)                                    # the initial phases (drawn in simulate)
    e_dir = rng.normal(size=(Ns, 3))
    dl = np.zeros((N, 3)); dl[:Ns] = q * e_dir
    M, dM = coupling(x, cfg['gamma0'], 1.0)
    return x, dl, M, dM


def effective(M, dl, gi, parts='all'):
    """C = M_ss - (i/2) X^T Y^-1 X and the radiators' response T (B = T s). parts: 'all'; 'no_heat' (delta = 0 in the
    coupling, the cold structure); 'heat_only_diag' (heat kept only in each piece's own offset)."""
    N = len(dl)
    idx_s = np.arange(N) * 4; idx_B = (np.arange(N)[:, None] * 4 + np.arange(1, 4)[None, :]).ravel()
    Mss = M[np.ix_(idx_s, idx_s)]; MsB = M[np.ix_(idx_s, idx_B)]; MBs = M[np.ix_(idx_B, idx_s)]; MBB = M[np.ix_(idx_B, idx_B)]
    D = np.zeros((3 * N, N))
    for j in range(N):
        D[3 * j:3 * j + 3, j] = dl[j]
    Y = 0.5j * MBB - gi * np.eye(3 * N)
    X = MBs - 2 * D
    Yi = np.linalg.inv(Y)
    C = Mss - 0.5j * X.T @ Yi @ X
    T = -0.5j * Yi @ X
    return C, T, dict(Mss=Mss, MsB=MsB, MBs=MBs, MBB=MBB, D=D, Yi=Yi)


def orders(cfg, seed, q, parity='even'):
    """C split exactly by order in the mixing (Y does not depend on delta): C = C0 + C1 + C2. C0 is the cold coupling
    (the pieces' passive radiators included), C1 the part linear in delta (a piece's quiet oscillation reaching another
    piece's quiet oscillation through one set of stirred radiators), C2 the part quadratic in delta (through two).
    parity 'odd': the mixing a velocity must give (imaginary and antisymmetric; code/rhythm_protect_v17.py), which
    makes C1 antisymmetric, so that its diagonal (a piece's own offset linear in delta) vanishes identically."""
    x, dl, M, dM = arrangement(cfg, seed, q)
    Ns = cfg['Ns']; N = len(dl)
    p = dict(om.PIECE, **cfg['piece']); w0, a = steady_amplitude(p)
    _, _, parts = effective(M, dl, cfg['gi'])
    Mss, MsB, MBs, D, Yi = parts['Mss'], parts['MsB'], parts['MBs'], parts['D'], parts['Yi']
    KsB, KBs = (-2 * D.T, -2 * D) if parity == 'even' else (-2j * D.T, 2j * D)
    C0 = Mss - 0.5j * MsB @ Yi @ MBs
    C1 = -0.5j * (MsB @ Yi @ KBs + KsB @ Yi @ MBs)
    C2 = -0.5j * KsB @ Yi @ KBs
    off = ~np.eye(N, dtype=bool)
    om_ = lambda C: -0.5 * w0 * np.real(np.diag(C))[:Ns]
    tug = lambda C: float(np.sqrt(np.mean(np.sum(np.abs((C * off)[:Ns, :Ns]) ** 2, 1))) * w0 / 2)
    grad = lambda C: float(np.sum(((C * off)[:Ns, :Ns]).imag ** 2) / max(np.sum(np.abs((C * off)[:Ns, :Ns]) ** 2), 1e-300))
    return dict(seed=seed, q=q, parity=parity, k=float(3 * q * q / (cfg['gamma0'] * (1 + cfg['gi']) ** 2)),
                offset_spread_cold=float(om_(C0).std()), offset_spread_linear=float(om_(C1).std()),
                offset_spread_quadratic=float(om_(C2).std()), offset_spread_all=float(om_(C0 + C1 + C2).std()),
                offset_mean_all=float(om_(C0 + C1 + C2).mean()),
                tug_cold=tug(C0), tug_linear=tug(C1), tug_quadratic=tug(C2), tug_all=tug(C0 + C1 + C2),
                gradient_share_linear=grad(C1), gradient_share_quadratic=grad(C2))


def offsets(C, w0):
    return -0.5 * w0 * np.real(np.diag(C))


def steady_amplitude(piece):
    lam = (piece['W0'] + piece['g_par']) / 2; w0 = 2 * lam / piece['G']
    a = np.sqrt((piece['W0'] * (1 - w0) - piece['g_par'] * (1 + w0)) / (2 * piece['G']))
    return w0, a


def budget(cfg, seed, q):
    """Part (1): the predicted own-rhythm offsets omega_j and the tugging couplings, split by origin."""
    x, dl, M, dM = arrangement(cfg, seed, q)
    Ns = cfg['Ns']
    p = dict(om.PIECE, **cfg['piece']); w0, a = steady_amplitude(p)
    C, T, parts = effective(M, dl, cfg['gi'])
    C0, _, _ = effective(M, np.zeros_like(dl), cfg['gi'])
    om_all = offsets(C, w0)[:Ns]
    off = ~np.eye(len(dl), dtype=bool)
    src = np.zeros(len(dl), bool); src[:Ns] = True
    pair = off & src[:, None] & src[None, :]
    heat = C - C0
    return dict(seed=seed, q=q, k=float(3 * q * q / (cfg['gamma0'] * (1 + cfg['gi']) ** 2)), w0=float(w0), a=float(a),
                offset_mean=float(om_all.mean()), offset_spread=float(om_all.std()),
                tug_cold_rms=float(np.sqrt(np.mean(np.abs(C0[pair]) ** 2)) * w0 / 2),
                tug_heat_rms=float(np.sqrt(np.mean(np.abs(heat[pair]) ** 2)) * w0 / 2),
                tug_cold_sum=float(np.sqrt(np.mean(np.sum(np.abs(C0[:Ns][:, :Ns] * off[:Ns][:, :Ns]) ** 2, 1))) * w0 / 2),
                tug_heat_sum=float(np.sqrt(np.mean(np.sum(np.abs(heat[:Ns][:, :Ns] * off[:Ns][:, :Ns]) ** 2, 1))) * w0 / 2),
                heat_radiative_share=float(np.sum(np.abs(heat[pair].imag) ** 2) / np.sum(np.abs(heat[pair]) ** 2)) if q > 0 else None)


def full_run(cfg, seed, q, T_run=8000.0, dt=0.02, burn=4000.0, protect=None):
    """Part (2): the full model again, recording each source piece's own rhythm (phase advance over [burn, T]) and each
    budget term's time average: the tug of the other pieces' quiet oscillations, of their radiators, and of its own
    radiators through its mixing (the last contains the own-rhythm offset)."""
    x, dl, M, dM = arrangement(cfg, seed, q)
    Ns, Nr = cfg['Ns'], cfg['Nr']; N = Ns + Nr
    p = dict(om.PIECE, **cfg['piece'])
    mat = om.Matter(x, cfg['gamma0'], 1.0, p, np.zeros(N), gi=cfg['gi'])
    rng = np.random.default_rng(seed)
    ball_min_sep(Ns, cfg['Rb'], 0.15, rng)                           # replay the draws: positions, then phases
    ph0 = rng.uniform(0, 2 * np.pi, N)
    w0, a = steady_amplitude(p)
    z = np.zeros((N, 4), complex); z[:, 0] = a * np.exp(1j * ph0); z = z.ravel()
    w = np.full(N, w0); n0 = 1e5; n = np.full(N, n0)
    Kb = om.mixing_block(dl)
    idx_s = np.arange(N) * 4
    offd = ~np.eye(N, dtype=bool)
    Mss_off = M[np.ix_(idx_s, idx_s)] * offd
    idx_B = (np.arange(N)[:, None] * 4 + np.arange(1, 4)[None, :]).ravel()
    MsB = M[np.ix_(idx_s, idx_B)]
    own = np.zeros((N, 3 * N));
    for j in range(N):
        own[j, 3 * j:3 * j + 3] = 1
    MsB_other = MsB * (1 - own)
    nstep = int(round(T_run / dt)); nb = int(round(burn / dt))
    ph = np.zeros(N); s_prev = None
    acc = dict(qq=np.zeros(N), rb=np.zeros(N), own=np.zeros(N)); cnt = 0
    for st in range(nstep):
        f = lambda z_, w_, n_: mat.rhs(z_, w_, n_, Kb, n0)
        k1 = f(z, w, n); k2 = f(z + .5 * dt * k1[0], w + .5 * dt * k1[1], n + .5 * dt * k1[2])
        k3 = f(z + .5 * dt * k2[0], w + .5 * dt * k2[1], n + .5 * dt * k2[2]); k4 = f(z + dt * k3[0], w + dt * k3[1], n + dt * k3[2])
        z = z + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6; w = w + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
        n = n + dt * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
        if st >= nb and st % 10 == 0:
            zz = z.reshape(N, 4); s = zz[:, 0]; B = zz[:, 1:].ravel()
            if s_prev is not None:
                ph += np.angle(s * np.conj(s_prev))
            s_prev = s.copy()
            acc['qq'] += -0.5 * w * np.real((Mss_off @ s) / s)
            acc['rb'] += -0.5 * w * np.real((MsB_other @ B) / s)
            acc['own'] += w * np.real(np.einsum('jm,jm->j', dl, zz[:, 1:]) / s)
            cnt += 1
    span = (cnt - 1) * 10 * dt
    rhythm = ph / span
    terms = {k: v / cnt for k, v in acc.items()}
    return dict(rhythm=rhythm[:Ns].tolist(), qq=terms['qq'][:Ns].tolist(), rb=terms['rb'][:Ns].tolist(), own=terms['own'][:Ns].tolist(),
                rhythm_receivers=rhythm[Ns:].tolist())


def reduced_run(cfg, seed, q, variant, T_run=16000.0, dt=0.5, burn=8000.0):
    """Part (3): the reduced rhythm model. variant: 'all'; 'no_offsets' (omega_j removed, tugs kept); 'offsets_only_heat'
    (the cold tugs with the warm offsets); 'no_heat_tugs' (warm offsets, cold tugs). Returns the sources' rhythm spread
    and the receivers' lead on the sources' wave (their field through the same coupling, radiators included)."""
    x, dl, M, dM = arrangement(cfg, seed, q)
    Ns, Nr = cfg['Ns'], cfg['Nr']; N = Ns + Nr
    p = dict(om.PIECE, **cfg['piece']); w0, a = steady_amplitude(p)
    C, T, parts = effective(M, dl, cfg['gi'])
    C0, T0, _ = effective(M, np.zeros_like(dl), cfg['gi'])
    om_w = offsets(C, w0); om_c = offsets(C0, w0)
    if variant == 'all':
        Cu, omv = C, om_w
    elif variant == 'no_offsets':
        Cu, omv = C, om_c
    elif variant == 'cold_tugs_warm_offsets':
        Cu, omv = C0, om_w
    elif variant == 'cold':
        Cu, omv = C0, om_c
    else:
        raise ValueError(variant)
    Coff = Cu * (~np.eye(N, dtype=bool))
    rng = np.random.default_rng(seed + 1000)
    phi = rng.uniform(0, 2 * np.pi, N)
    idx_s = np.arange(N) * 4
    # the sources' wave at each receiver's quiet channel, radiators included: E = (M_ss + M_sB T) s over sources
    Mfull_sources = (parts['Mss'] + parts['MsB'] @ T)[Ns:, :Ns]
    nstep = int(T_run / dt); nb = int(burn / dt)
    ph_acc = np.zeros(N); lead = 0.0; cnt = 0; prev = phi.copy()

    def rhs(ph_):
        e = np.exp(1j * ph_)
        return omv - 0.5 * w0 * np.real((Coff @ e) * np.conj(e))

    for st in range(nstep):
        k1 = rhs(phi); k2 = rhs(phi + .5 * dt * k1); k3 = rhs(phi + .5 * dt * k2); k4 = rhs(phi + dt * k3)
        phi = phi + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if st >= nb:
            E = Mfull_sources @ np.exp(1j * phi[:Ns])
            lead += np.mean(np.sin(np.angle(E) - phi[Ns:])); cnt += 1
            if st == nb:
                ph_start = phi.copy()
    rhythm = (phi - ph_start) / ((nstep - 1 - nb) * dt)
    return dict(variant=variant, spread=float(rhythm[:Ns].std()), lead=float(lead / cnt))


def job(args):
    kind, cfg, seed, q, extra = args
    t0 = time.time()
    if kind == 'budget':
        r = budget(cfg, seed, q)
    elif kind == 'orders':
        r = orders(cfg, seed, q, extra)
    elif kind == 'full':
        r = full_run(cfg, seed, q)
    else:
        r = reduced_run(cfg, seed, q, extra)
    r.update(kind=kind, seed=seed, q=q, extra=extra, seconds=time.time() - t0)
    return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4); ap.add_argument('--quick', action='store_true')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    qs = (0.0, 0.2, 0.41, 0.82, 1.15)
    seeds = (1,) if args.quick else (1, 2, 3)
    jobs = [('budget', MAIN, sd, q, None) for sd in seeds for q in qs]
    jobs += [('orders', MAIN, sd, q, par) for sd in seeds for q in qs[1:] for par in ('even', 'odd')]
    jobs += [('full', MAIN, 1, q, None) for q in (0.0, 0.41, 0.82)]
    jobs += [('reduced', MAIN, sd, q, v) for sd in seeds for q in (0.0, 0.41, 0.82) for v in
             (('all', 'no_offsets', 'cold_tugs_warm_offsets') if q > 0 else ('all',))]
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(job, jobs):
            res.append(r)
            if r['kind'] == 'budget':
                print(f"[{time.monotonic() - t0:5.0f} s] budget seed {r['seed']} q {r['q']:.2f} (k {r['k']:.2f}): own offsets mean "
                      f"{r['offset_mean']:+.2e} spread {r['offset_spread']:.2e}; tugs cold {r['tug_cold_rms']:.2e}, heat {r['tug_heat_rms']:.2e}", flush=True)
            elif r['kind'] == 'orders':
                print(f"[{time.monotonic() - t0:5.0f} s] orders seed {r['seed']} q {r['q']:.2f} {r['parity']}: own offsets cold {r['offset_spread_cold']:.2e} "
                      f"linear {r['offset_spread_linear']:.2e} quadratic {r['offset_spread_quadratic']:.2e}; tugs cold {r['tug_cold']:.2e} "
                      f"linear {r['tug_linear']:.2e} quadratic {r['tug_quadratic']:.2e}", flush=True)
            elif r['kind'] == 'reduced':
                print(f"[{time.monotonic() - t0:5.0f} s] reduced seed {r['seed']} q {r['q']:.2f} {r['variant']:>24}: rhythm spread {r['spread']:.2e}, lead {r['lead']:+.2f}", flush=True)
            else:
                rh = np.array(r['rhythm']); tot = np.array(r['qq']) + np.array(r['rb']) + np.array(r['own'])
                print(f"[{time.monotonic() - t0:5.0f} s] full q {r['q']:.2f}: rhythm spread {rh.std():.2e}; budget sum vs rhythm corr "
                      f"{np.corrcoef(rh, tot)[0, 1]:.3f}; term spreads qq {np.std(r['qq']):.2e} rb {np.std(r['rb']):.2e} own {np.std(r['own']):.2e}", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 17 step A: the rhythm budget of warm and cold sources', runs=res,
                                           seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
