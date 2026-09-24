"""Round 17, step B (reduced model): which internal structure keeps a warm source in tune?

code/rhythm_budget_v17.py showed that the reduced rhythm model (each piece's fast internal modes eliminated exactly in
linear response) reproduces the full model: a cold source settles and its receivers keep step, a warm one does not. Here
the same reduction is applied to a family of internal structures for the "loud" part of a piece, all energy-consistent
(Hermitian couplings, positive losses, every loss an explicit internal channel), all the same for every piece, and
compared at equal released glow.

A structure is a list of internal families (each a triplet, one mode per direction) plus an optional chamber:
  bright family   radiates as a dipole into the shared wave (rate gamma = 1 from M) and loses `loss` inside; its own
                  frequency sits at Delta; the motion mixes the quiet oscillation s into it with weight c.
  dark family     the same, but it does not couple to the wave at all (its losses are all inside).
  chamber (valve) a strongly damped internal triplet R (rate gR) that the motion opens from s (weight c) and that passes
                  on to one bright family by a fixed internal coupling gB: a reservoir-mediated, dissipative route.
The motion's mixing has a time-reversal parity:
  even  round 16's (and the independent calculation's) choice: a real symmetric coupling, K = -2 c delta both ways;
  odd   what a velocity must give: a velocity changes sign when time runs backwards, and a coupling linear in it
        between two modes that do not must be imaginary and antisymmetric (a Coriolis- or Fizeau-like drag):
        K_Fs = +2i c delta, K_sF = -2i c delta^T. It is Hermitian, so it exchanges energy without creating any.
For one isolated piece the two parities give the same glow and no rhythm shift; they differ in how a piece's radiators
and its neighbours' waves feed back on its rhythm (the part of the coupling linear in delta becomes antisymmetric, and a
piece's own rhythm offset from it vanishes identically).

Round 16 and the independent review's structures are special cases:
  single     one bright family, c = 1, Delta = 0, loss gi = 4.
  paired     two bright families at +Delta and -Delta with c = sqrt(1/2) each (the review's proposal).
  valve      the chamber (gR = 20, gB = 5, c = 1) feeding one bright family that the motion does not mix directly.

For each structure: the fast modes F obey 0 = P F + Q s, the quiet oscillations are driven by M_ss s + H F, so
C = M_ss - H P^-1 Q, each piece's own rhythm offset is -(w0/2) Re C_jj and its tug on the others is C_jl. The heat is set
by the isolated piece's released glow, k = extra/cold, matched between structures. The receivers keep step on the
sources' part of their drive, C[receivers, sources] e^{i phi_sources}.

    python code/rhythm_protect_v17.py --output run-rhythm-budget-v17/rhythm_protect_v17.json [--processes 4]
    python code/rhythm_protect_v17.py --set compact --output run-rhythm-budget-v17/rhythm_protect_compact_v17.json
    python code/rhythm_protect_v17.py --set trace --output run-rhythm-budget-v17/rhythm_trace_v17.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rhythm_budget_v17 import MAIN, arrangement, steady_amplitude          # noqa: E402
import one_matter_v16 as om                                              # noqa: E402


def blocks(M, N):
    idx_s = np.arange(N) * 4; idx_B = (np.arange(N)[:, None] * 4 + np.arange(1, 4)[None, :]).ravel()
    return (M[np.ix_(idx_s, idx_s)], M[np.ix_(idx_s, idx_B)], M[np.ix_(idx_B, idx_s)], M[np.ix_(idx_B, idx_B)])


def mixing_matrix(dl):
    N = len(dl); D = np.zeros((3 * N, N))
    for j in range(N):
        D[3 * j:3 * j + 3, j] = dl[j]
    return D


def motion_blocks(D, c, parity):
    """The motion's mixing between s and one fast triplet, weight c: (K_Fs, K_sF) in the convention dF/dt = (i/2) K_Fs s,
    ds/dt = -(i/2) w K_sF F (round 16's)."""
    if parity == 'even':
        return -2 * c * D, -2 * c * D.T
    if parity == 'odd':
        return 2j * c * D, -2j * c * D.T
    raise ValueError(parity)


def reduce(M, dl, st):
    """C (N x N), the fast modes' response Fs = -P^-1 Q, and the radiators' total dipole per unit s, R_B Fs."""
    N = len(dl)
    Mss, MsB, MBs, MBB = blocks(M, N)
    D = mixing_matrix(dl)
    I3 = np.eye(3 * N)
    fams = st['families']; val = st.get('valve'); par = st.get('parity', 'even')
    nf = len(fams) + (1 if val else 0)
    P = np.zeros((nf * 3 * N, nf * 3 * N), complex); Q = np.zeros((nf * 3 * N, N), complex)
    H = np.zeros((N, nf * 3 * N), complex); RB = np.zeros((3 * N, nf * 3 * N))
    sl = lambda a: slice(a * 3 * N, (a + 1) * 3 * N)
    for a, fa in enumerate(fams):
        bright_a = fa['kind'] == 'bright'; sa = fa.get('sign', 1.0)
        Kfs, Ksf = motion_blocks(D, fa['c'], par)
        P[sl(a), sl(a)] += -(fa['loss'] + 1j * fa['Delta']) * I3
        if bright_a:
            for b, fb in enumerate(fams):
                if fb['kind'] == 'bright':
                    P[sl(a), sl(b)] += 0.5j * sa * fb.get('sign', 1.0) * MBB
            Q[sl(a)] += 0.5j * sa * MBs; H[:, sl(a)] += sa * MsB; RB[:, sl(a)] = sa * I3
        Q[sl(a)] += 0.5j * Kfs; H[:, sl(a)] += Ksf
    if val:
        r = len(fams); t = val['into']
        Krs, Ksr = motion_blocks(D, val['c'], par)
        P[sl(r), sl(r)] = -val['gR'] * I3
        P[sl(r), sl(t)] += -1j * val['gB'] * I3; P[sl(t), sl(r)] += -1j * val['gB'] * I3
        Q[sl(r)] += 0.5j * Krs; H[:, sl(r)] += Ksr
    Fs = -np.linalg.solve(P, Q)
    C = Mss + H @ Fs
    return C, Fs, RB @ Fs, dict(Mss=Mss, MsB=MsB)


def k_isolated(st, q, g0):
    """The released glow of one isolated piece moving with rms delta q per component: extra radiated power over the cold
    2 gamma_0 |s|^2, from the same reduction."""
    from shared_wave_v16 import coupling
    M, _ = coupling(np.zeros((1, 3)), g0, 1.0)
    dl = np.array([[q, q, q]])                                         # |delta|^2 = 3 q^2, the rms of the ensemble
    C, Fs, Bsum, _ = reduce(M, dl, st)
    extra = 2 * 1.0 * np.sum(np.abs(Bsum[:, 0]) ** 2)
    return extra / (2 * g0)


def q_for_k(st, k, g0):
    return brentq(lambda q: k_isolated(st, q, g0) - k, 1e-6, 200.0)


def phase_run(C, Ns, w0, seed, T_run=16000.0, dt=0.5, burn=8000.0):
    """The reduced rhythm model: dphi_j/dt = -(w0/2) Re[(C e^{i phi})_j e^{-i phi_j}], own offsets included on the
    diagonal. Returns the sources' rhythm spread and the receivers' lead on the sources' part of their drive."""
    N = len(C)
    omv = -0.5 * w0 * np.real(np.diag(C))
    Coff = C * (~np.eye(N, dtype=bool))
    Esrc = C[Ns:, :Ns]
    rng = np.random.default_rng(seed + 1000)
    phi = rng.uniform(0, 2 * np.pi, N)
    nstep = int(T_run / dt); nb = int(burn / dt)
    lead = 0.0; cnt = 0

    def rhs(ph_):
        e = np.exp(1j * ph_)
        return omv - 0.5 * w0 * np.real((Coff @ e) * np.conj(e))

    for s_ in range(nstep):
        k1 = rhs(phi); k2 = rhs(phi + .5 * dt * k1); k3 = rhs(phi + .5 * dt * k2); k4 = rhs(phi + dt * k3)
        phi = phi + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if s_ == nb:
            ph_start = phi.copy()
        if s_ >= nb:
            E = Esrc @ np.exp(1j * phi[:Ns])
            lead += np.mean(np.sin(np.angle(E) - phi[Ns:])); cnt += 1
    rhythm = (phi - ph_start) / ((nstep - 1 - nb) * dt)
    return float(rhythm[:Ns].std()), float(lead / cnt)


def reduced_run(cfg, seed, st, k, diagnostic=None):
    q = q_for_k(st, k, cfg['gamma0']) if k > 0 else 0.0
    x, dl, M, dM = arrangement(cfg, seed, q)
    Ns = cfg['Ns']; N = len(dl)
    p = dict(om.PIECE, **cfg['piece']); w0, a = steady_amplitude(p)
    C, Fs, Bsum, parts = reduce(M, dl, st)
    if diagnostic:                                                   # parts of the heat switched off by hand
        C0, _, _, _ = reduce(M, np.zeros_like(dl), st)
        heat = C - C0; dg = np.diag(np.diag(heat)); off = heat - dg
        if diagnostic == 'heat_tugs_gradient_part_only':
            C = C0 + dg + 1j * off.imag
        elif diagnostic == 'heat_tugs_nongradient_part_only':
            C = C0 + dg + off.real
        elif diagnostic == 'no_heat_tugs_among_sources':
            mask = np.zeros((N, N), bool); mask[:Ns, :Ns] = True
            C = C0 + dg + off * ~mask
        elif diagnostic == 'no_heat_offsets':
            C = C0 + 1j * dg.imag + off
        else:
            raise ValueError(diagnostic)
    sp, ld = phase_run(C, Ns, w0, seed)
    off_src = -0.5 * w0 * np.real(np.diag(C))[:Ns]
    Coff = C * (~np.eye(N, dtype=bool))
    return dict(q=float(q), k=float(k), spread=sp, lead=ld, offset_spread=float(off_src.std()),
                tug_sum=float(np.sqrt(np.mean(np.sum(np.abs(Coff[:Ns, :Ns]) ** 2, 1))) * w0 / 2),
                tug_gradient_share=float(np.sum(Coff[:Ns, :Ns].imag ** 2) / np.sum(np.abs(Coff[:Ns, :Ns]) ** 2)),
                wave_rms=float(np.sqrt(np.mean(np.sum(np.abs(C[Ns:, :Ns]) ** 2, 1)))))


def fam(c=1.0, Delta=0.0, loss=4.0, kind='bright', sign=1.0):
    return dict(kind=kind, c=c, Delta=Delta, loss=loss, sign=sign)


h = float(np.sqrt(0.5))
STRUCTURES = {
    'single': dict(parity='even', families=[fam()]),
    'single, odd': dict(parity='odd', families=[fam()]),
    'paired, Delta = gamma + gi': dict(parity='even', families=[fam(h, 5.0), fam(h, -5.0)]),
    'paired, Delta = gamma + gi, odd': dict(parity='odd', families=[fam(h, 5.0), fam(h, -5.0)]),
    'paired, Delta = (gamma + gi)/2': dict(parity='even', families=[fam(h, 2.5), fam(h, -2.5)]),
    'paired, Delta = gamma + gi, 1% asymmetry': dict(parity='even', families=[fam(h * 1.01, 5.0), fam(h * 0.99, -5.0)]),
    'paired, Delta = 4 (gamma + gi)': dict(parity='even', families=[fam(h, 20.0), fam(h, -20.0)]),
    'paired, Delta = 4 (gamma + gi), odd': dict(parity='odd', families=[fam(h, 20.0), fam(h, -20.0)]),
    'bright + dark, Delta = gamma + gi': dict(parity='even', families=[fam(h, 5.0), fam(h, -5.0, loss=5.0, kind='dark')]),
    'bright + dark, Delta = gamma + gi, odd': dict(parity='odd', families=[fam(h, 5.0), fam(h, -5.0, loss=5.0, kind='dark')]),
    'valve, gR = 20, gB = 5': dict(parity='even', families=[fam(0.0)], valve=dict(gR=20.0, gB=5.0, c=1.0, into=0)),
    'valve, gR = 20, gB = 5, odd': dict(parity='odd', families=[fam(0.0)], valve=dict(gR=20.0, gB=5.0, c=1.0, into=0)),
}
# round 16's low-f matter (radiators sending a twentieth of their energy into the companion), with the velocity's parity;
# used by code/one_matter_v17.py, not in this script's comparison
EXTRA = {'single, odd, f = 0.05': dict(parity='odd', families=[fam(loss=19.0)])}
DIAGNOSTICS = ('heat_tugs_gradient_part_only', 'heat_tugs_nongradient_part_only', 'no_heat_tugs_among_sources', 'no_heat_offsets')


COMPACT = dict(MAIN, Rb=0.6)          # the same 48 pieces in a ball smaller than a wavelength


def trace(args):
    """A run of the reduced model 24,000 time units long, reported every 2,000: the sources' mean rhythm and spread,
    the receivers' mean rhythm, their keeping step, and the sources' synchrony."""
    name, seed, k = args
    st = STRUCTURES[name]
    q = q_for_k(st, k, MAIN['gamma0']) if k > 0 else 0.0
    x, dl, M, dM = arrangement(MAIN, seed, q)
    Ns = MAIN['Ns']; N = len(dl)
    p = dict(om.PIECE, **MAIN['piece']); w0, a = steady_amplitude(p)
    C, Fs, Bsum, parts = reduce(M, dl, st)
    omv = -0.5 * w0 * np.real(np.diag(C)); Coff = C * (~np.eye(N, dtype=bool)); Esrc = C[Ns:, :Ns]
    rng = np.random.default_rng(seed + 1000); phi = rng.uniform(0, 2 * np.pi, N)
    dt = 0.5; win = 2000.0; out = []; mark = phi.copy(); leads = []

    def rhs(ph_):
        e = np.exp(1j * ph_)
        return omv - 0.5 * w0 * np.real((Coff @ e) * np.conj(e))

    for s_ in range(1, int(24000.0 / dt) + 1):
        k1 = rhs(phi); k2 = rhs(phi + .5 * dt * k1); k3 = rhs(phi + .5 * dt * k2); k4 = rhs(phi + dt * k3)
        phi = phi + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        E = Esrc @ np.exp(1j * phi[:Ns]); leads.append(np.mean(np.sin(np.angle(E) - phi[Ns:])))
        if s_ % int(win / dt) == 0:
            r = (phi - mark) / win
            out.append(dict(t=s_ * dt, sources_rhythm=float(r[:Ns].mean()), sources_spread=float(r[:Ns].std()),
                            receivers_rhythm=float(r[Ns:].mean()), receivers_spread=float(r[Ns:].std()),
                            keeping_step=float(np.mean(leads)), sync=float(abs(np.mean(np.exp(1j * phi[:Ns]))))))
            mark = phi.copy(); leads = []
    return dict(structure=name, seed=seed, k=k, q=float(q), windows=out)


def job(args):
    name, seed, k, diag, cfg = (tuple(args) + (MAIN,))[:5]
    t0 = time.time()
    r = reduced_run(cfg, seed, STRUCTURES[name], k, diag)
    r.update(structure=name, seed=seed, diagnostic=diag, seconds=time.time() - t0)
    return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4)
    ap.add_argument('--set', default='main', help="'main' (radius 3 wavelengths/2 pi), 'compact' (0.6) or 'trace' (long "
                    "time traces: the valve with the odd parity at k = 8, and a cold source)")
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    if args.set == 'trace':
        jobs = [('valve, gR = 20, gB = 5, odd', 1, 8.0), ('valve, gR = 20, gB = 5, odd', 2, 8.0), ('single', 1, 0.0)]
        with Pool(args.processes) as pool:
            res = list(pool.imap(trace, jobs))
        for r in res:
            print(f"{r['structure']} seed {r['seed']} k {r['k']}:")
            for w in r['windows']:
                print(f"   t {w['t']:6.0f}: sources {w['sources_rhythm']:+.2e} (spread {w['sources_spread']:.1e}), receivers "
                      f"{w['receivers_rhythm']:+.2e} (spread {w['receivers_spread']:.1e}), keeping step {w['keeping_step']:+.2f}")
        args.output.write_text(json.dumps(dict(experiment='round 17 step B: time traces (reduced model)', runs=res), indent=1) + '\n')
        print('wrote', args.output)
        return
    if args.set == 'compact':
        names = ['single', 'single, odd', 'bright + dark, Delta = gamma + gi', 'bright + dark, Delta = gamma + gi, odd',
                 'paired, Delta = 4 (gamma + gi), odd', 'valve, gR = 20, gB = 5', 'valve, gR = 20, gB = 5, odd']
        jobs = [(n, sd, k, None, COMPACT) for n in names for sd in (1, 2, 3) for k in (0.0, 2.0, 8.0)]
        jobs += [('single', sd, k, dg, COMPACT) for sd in (1, 2, 3) for k in (2.0, 8.0) for dg in DIAGNOSTICS[:2]]
    else:
        jobs = [(name, sd, k, None) for name in STRUCTURES for sd in (1, 2, 3) for k in (0.0, 2.0, 8.0, 16.0)]
        jobs += [(name, sd, k, dg) for name in ('single', 'single, odd') for dg in DIAGNOSTICS for sd in (1, 2, 3) for k in (2.0, 8.0)]
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(job, jobs):
            res.append(r)
            print(f"[{time.monotonic() - t0:5.0f} s] {r['structure']:>40} {r['diagnostic'] or '':>32} seed {r['seed']} k {r['k']:5.1f} (q {r['q']:.3f}): "
                  f"spread {r['spread']:.2e} (own offsets {r['offset_spread']:.2e}, tugs {r['tug_sum']:.2e}, gradient share {r['tug_gradient_share']:.2f}), "
                  f"lead {r['lead']:+.2f}", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 17 step B: internal structures that keep a warm source in tune (reduced model)',
                                           set=args.set, source_radius=(COMPACT if args.set == 'compact' else MAIN)['Rb'],
                                           structures=STRUCTURES, diagnostics=DIAGNOSTICS, runs=res,
                                           seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
