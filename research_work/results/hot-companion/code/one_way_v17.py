"""Round 17, a hypothesis test (reduced model): does a wave that travels only outward keep warm matter in tune?

code/rhythm_protect_v17.py found that what puts a warm source out of tune is the two-way exchange between its pieces:
each piece's wave comes back to it off its neighbours (its own rhythm offset) and neighbours tug each other back and
forth (the non-gradient part of the tugs). No internal structure of a piece removes that: at a given glow, reciprocity
fixes its size. What would remove it is a wave that cannot come back: if the companion's wave is carried outward by the
companion's own outflow faster than it can travel against it (the review's "one medium for crest and transport
speeds"), a piece hears only matter nearer the source's centre, never its own echo.

This script repeats the reduced rhythm model (code/rhythm_protect_v17.py: fast internal modes eliminated exactly,
phases evolved) with the coupling M made one-way before the reduction: the part carrying a wave from piece l to piece
j is kept when l is nearer the centre than j, and multiplied by a leak factor eps otherwise (eps = 1 is the ordinary
two-way wave, eps = 0 strictly outward). A 'core' variant keeps the wave two-way inside a central core (where an
outflow starts from rest) and one-way outside. Every piece's own self-coupling is kept. Receivers of the same matter
sit on four shells (radius 6, 9, 13.5 and 20, four each), all outside the source.

It is a test of the hypothesis only: the one-way coupling here is made by hand, not derived from a model of the
companion's flow, and the reduced model tracks rhythms, not energy or forces. In the reduced model the pull on a locked
receiver is its lead times the wave's amplitude, so R(sigma, r) = lead(sigma, r)/lead(0, r).

    python code/one_way_v17.py --output run-rhythm-budget-v17/one_way_v17.json [--processes 3]
    python code/one_way_v17.py --set mass --output run-rhythm-budget-v17/one_way_mass_v17.json
    python code/one_way_v17.py --set order --output run-rhythm-budget-v17/one_way_order_v17.json
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
from shared_wave_v16 import ball_min_sep, coupling                   # noqa: E402

RADII = (6.0, 9.0, 13.5, 20.0)


def shells(radii, per=4):
    from one_matter_v16 import sphere_points
    out = []
    for i, r in enumerate(radii):
        a = 0.7 * i; ca, sa = np.cos(a), np.sin(a)
        rot = np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]]) @ np.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]])
        out.append(r * sphere_points(per, 1.0) @ rot.T)
    return np.vstack(out)


def geometry(seed, q, Ns=48, Rb=3.0):
    """The main run's arrangement (same draws as code/one_matter_v17.py), receivers on the four shells."""
    rng = np.random.default_rng(seed)
    xs = ball_min_sep(Ns, Rb, 0.15, rng); x = np.vstack([xs, shells(RADII)]); N = len(x)
    rng.uniform(0, 2 * np.pi, N)                                     # the phases drawn in the full model (kept in step)
    e_dir = rng.normal(size=(Ns, 3)); dl = np.zeros((N, 3)); dl[:Ns] = q * e_dir
    M, _ = coupling(x, MAIN['gamma0'], 1.0)
    return x, dl, M


def directional(M, x, eps=0.0, core=0.0, order=None):
    """M with the inward part multiplied by eps: j receives from l fully when l is nearer the centre (or both lie inside
    the core), and eps times otherwise. order: a rank for each piece to use in place of its distance from the centre
    (e.g. a random one-way order among the sources, the receivers after them)."""
    r = np.linalg.norm(x, axis=1) if order is None else np.asarray(order, float); N = len(x)
    w = np.where(r[:, None] > r[None, :], 1.0, eps)
    if core > 0:
        inside = r < core
        w[np.ix_(inside, inside)] = 1.0
    np.fill_diagonal(w, 1.0)
    return M * np.kron(w, np.ones((4, 4)))


def run(args):
    name, seed, k, eps, core = args[:5]
    Ns = args[5] if len(args) > 5 else 48
    random_order = len(args) > 6 and args[6] == 'random order'
    t0 = time.time()
    st = rp.STRUCTURES[name]
    q = rp.q_for_k(st, k, MAIN['gamma0']) if k > 0 else 0.0
    x, dl, M = geometry(seed, q, Ns=Ns, Rb=3.0 * (Ns / 48) ** (1 / 3))
    order = None
    if random_order:                                                # a random one-way order among the sources
        rr = np.linalg.norm(x, axis=1)
        order = np.concatenate([np.random.default_rng(seed + 77).permutation(Ns).astype(float), 1e3 + rr[Ns:]])
    M = directional(M, x, eps, core, order)
    N = len(x)
    p = dict(om.PIECE, **MAIN['piece']); w0, a = steady_amplitude(p)
    C, Fs, Bsum, parts = rp.reduce(M, dl, st)
    omv = -0.5 * w0 * np.real(np.diag(C)); Coff = C * (~np.eye(N, dtype=bool)); Esrc = C[Ns:, :Ns]
    rng = np.random.default_rng(seed + 1000); phi = rng.uniform(0, 2 * np.pi, N)
    dt = 0.5; T = 16000.0; nb = int(8000.0 / dt); nstep = int(T / dt)

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
    return dict(structure=name, seed=seed, k=k, q=float(q), eps=eps, core=core, Ns=Ns, order='random' if random_order else 'radius',
                spread=float(rh[:Ns].std()),
                common_rhythm=float(rh[:Ns].mean()), receivers_rhythm=shell(rh[Ns:]), own_offset_spread=float(omv[:Ns].std()),
                lead_by_radius=shell(lead), wave_by_radius=shell(amp), seconds=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=3)
    ap.add_argument('--set', default='leak', help="'leak': how much inward wave the protection tolerates; 'mass': "
                    "sources of 24, 48 and 96 pieces at the same density, does the wave at the receivers grow as the "
                    "square root of the number of pieces (power in proportion to mass) or faster")
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    if args.set == 'order':
        jobs = [('single', sd, k, 0.0, 0.0, 48, 'random order') for sd in (1, 2, 3) for k in (0.0, 8.0, 16.0)]
    elif args.set == 'mass':
        jobs = [('single', sd, k, eps, 0.0, Ns) for eps in (1.0, 0.0) for Ns in (24, 48, 96) for sd in (1, 2, 3) for k in (0.0, 8.0)]
    else:
        variants = [(1.0, 0.0), (0.5, 0.0), (0.3, 0.0), (0.1, 0.0), (0.0, 0.0), (0.0, 1.5)]
        jobs = [(n, sd, k, eps, core) for (eps, core) in variants for n in ('single', 'single, odd') for sd in (1, 2, 3)
                for k in (0.0, 2.0, 8.0, 16.0)]
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(run, jobs):
            res.append(r)
            print(f"[{time.monotonic() - t0:5.0f} s] eps {r['eps']:.1f} core {r['core']:.1f} {r['structure']:>12} Ns {r['Ns']} seed {r['seed']} k {r['k']:4.0f}: "
                  f"sources' spread {r['spread']:.1e}, common rhythm {r['common_rhythm']:+.1e}; keeping step at r = "
                  f"{', '.join('%.1f' % v for v in RADII)}: {' '.join('%+.2f' % v for v in r['lead_by_radius'])}; wave "
                  f"{' '.join('%.2e' % v for v in r['wave_by_radius'])}", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 17: a wave that travels only outward (reduced model, hypothesis test)',
                                           set=args.set, radii=RADII, runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
