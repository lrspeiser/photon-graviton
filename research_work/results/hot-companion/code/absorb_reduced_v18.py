"""Round 18: which absorbing stream keeps warm matter in tune? A scan in round 17's reduced rhythm model.

The full model (code/absorbing_stream_v18.py) with a stream that absorbs only the inward-moving part of a wave removes
the push but leaves a warm source's pieces out of tune: neighbours at the same radius still hear each other both ways,
because the wave between them travels sideways. A stream whose absorption grows with the wave's speed relative to the
stream absorbs sideways waves too. Two media, in the eikonal form of their local transport equations (stream along
e = x/|x|, a wave travelling along n):
    'inward':  rate kappa max(0, -n . e)      ->  T_jl = exp(-kappa (r_l - smallest radius reached on the path))
    'kinetic': rate kappa (1 - n . e)/2       ->  T_jl = exp(-kappa (d_jl - (r_j - r_l))/2)
(the second integrates exactly, since the integral of n . e along a straight path is the change in radius: outward
radial paths are untouched, sideways ones lose kappa d/2, inward ones kappa d). The reduced model (fast internal modes
eliminated exactly, phases evolved; code/one_way_v17.py) runs with the coupling multiplied pair by pair by T, for
round 16's matter at rest and at k = 8, three arrangements, and reports the sources' rhythm spread and how well the
receivers on the four shells keep step. Round 17's imposed radius mask (eps = 0) and the two-way wave are the
references.
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
import rhythm_protect_v17 as rp                                      # noqa: E402
import one_matter_v16 as om                                          # noqa: E402
from rhythm_budget_v17 import MAIN, steady_amplitude                 # noqa: E402
from one_way_v17 import geometry, directional, RADII                 # noqa: E402
from absorbing_stream_v18 import absorption_factors                  # noqa: E402


def kinetic_factors(x, kappa):
    r = np.linalg.norm(x, axis=1)
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    T = np.exp(-kappa * np.clip(d - (r[:, None] - r[None, :]), 0, None) / 2)
    np.fill_diagonal(T, 1.0)
    return T


def run(args):
    medium, kappa, seed, k = args
    t0 = time.time()
    st = rp.STRUCTURES['single']
    q = rp.q_for_k(st, k, MAIN['gamma0']) if k > 0 else 0.0
    Ns = 48
    x, dl, M = geometry(seed, q)
    if medium == 'two-way':
        pass
    elif medium == 'radius mask':
        M = directional(M, x, 0.0)
    elif medium == 'inward':
        T, _ = absorption_factors(x, kappa); M = M * np.kron(T, np.ones((4, 4)))
    elif medium == 'kinetic':
        T = kinetic_factors(x, kappa); M = M * np.kron(T, np.ones((4, 4)))
    N = len(x)
    p = dict(om.PIECE, **MAIN['piece']); w0, a = steady_amplitude(p)
    C, Fs, Bsum, parts = rp.reduce(M, dl, st)
    omv = -0.5 * w0 * np.real(np.diag(C)); Coff = C * (~np.eye(N, dtype=bool)); Esrc = C[Ns:, :Ns]
    rng = np.random.default_rng(seed + 1000); phi = rng.uniform(0, 2 * np.pi, N)
    dt = 0.5; T_run = 16000.0; nb = int(8000.0 / dt); nstep = int(T_run / dt)

    def rhs(ph_):
        e = np.exp(1j * ph_)
        return omv - 0.5 * w0 * np.real((Coff @ e) * np.conj(e))

    lead = np.zeros(N - Ns); amp = np.zeros(N - Ns); cnt = 0
    for s_ in range(nstep):
        k1 = rhs(phi); k2 = rhs(phi + .5 * dt * k1); k3 = rhs(phi + .5 * dt * k2); k4 = rhs(phi + dt * k3)
        phi = phi + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if s_ == nb:
            ph0 = phi.copy()
        if s_ >= nb:
            E = Esrc @ np.exp(1j * phi[:Ns]); lead += np.sin(np.angle(E) - phi[Ns:]); amp += np.abs(E); cnt += 1
    rh = (phi - ph0) / ((nstep - 1 - nb) * dt)
    lead /= cnt; amp /= cnt
    shell = lambda v: [float(v[4 * i:4 * i + 4].mean()) for i in range(len(RADII))]
    # how two-way the sources' coupling still is: the largest |M_lj|/|M_jl| asymmetry summarized by the mean over
    # neighbouring pairs of min(|M_jl|, |M_lj|)/max(...) for the monopole blocks
    Mm = np.abs(M[0::4, 0::4][:Ns, :Ns]); iu = np.triu_indices(Ns, 1)
    a1, a2 = Mm[iu], Mm.T[iu]
    two_way = float(np.sum(np.minimum(a1, a2)) / np.sum(np.maximum(a1, a2)))
    return dict(medium=medium, kappa=kappa, seed=seed, k=k, spread=float(rh[:Ns].std()), common=float(rh[:Ns].mean()),
                lead_by_radius=shell(lead), wave_by_radius=shell(amp), two_way_share=two_way, seconds=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    media = [('two-way', 0.0), ('radius mask', 0.0)] + [('inward', kp) for kp in (2.3, 5.0, 10.0, 20.0)] + \
            [('kinetic', kp) for kp in (1.0, 2.3, 5.0, 10.0)]
    jobs = [(m, kp, sd, k) for m, kp in media for sd in (1, 2, 3) for k in (0.0, 8.0)]
    t0 = time.monotonic(); res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(run, jobs):
            res.append(r)
            print(f"[{time.monotonic() - t0:5.0f} s] {r['medium']:>11} kappa {r['kappa']:4.1f} seed {r['seed']} k {r['k']:3.0f}: spread {r['spread']:.1e}, "
                  f"keeping step {', '.join('%.2f' % v for v in r['lead_by_radius'])}, two-way share {r['two_way_share']:.2f}", flush=True)
    table = []
    for m, kp in media:
        row = dict(medium=m, kappa=kp)
        for k in (0.0, 8.0):
            rs = [r for r in res if r['medium'] == m and r['kappa'] == kp and r['k'] == k]
            row[f'k{k:g}'] = dict(spread=float(np.mean([r['spread'] for r in rs])),
                                  lead_by_radius=np.mean([r['lead_by_radius'] for r in rs], 0).tolist(),
                                  wave_by_radius=np.mean([r['wave_by_radius'] for r in rs], 0).tolist(),
                                  two_way_share=float(np.mean([r['two_way_share'] for r in rs])))
        c, w = row['k0'], row['k8']
        row['R_by_radius'] = [(lw / lc) if lc > 0 else None for lw, lc in zip(w['lead_by_radius'], c['lead_by_radius'])]
        table.append(row)
        print(f"{m:>11} kappa {kp:4.1f}: spread cold {c['spread']:.1e} warm {w['spread']:.1e}; keeping step warm "
              f"{', '.join('%.2f' % v for v in w['lead_by_radius'])} (cold {', '.join('%.2f' % v for v in c['lead_by_radius'])}); "
              f"two-way share {c['two_way_share']:.2f}", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 18: absorbing streams in the reduced rhythm model', runs=res,
                                           table=table, seconds=time.monotonic() - t0), indent=1) + '\n')


if __name__ == '__main__':
    main()
