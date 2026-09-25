"""Round 18, step 1 in three dimensions: a one-way wave from a local property of the stream, with no push.

code/flowing_medium_v18.py found that waves carried by the companion's own stream are one-way but push every emitter
downstream (wave drag), while waves that travel through matter's frame at their own speed, with the stream absorbing
the part that moves against it, are one-way with no push. This builds the second medium into round 16's full model
of one kind of matter and repeats round 17's distance test with it, in place of the inward couplings deleted by hand.

The medium, in the eikonal (ray) form of its local transport equation: a wave travelling in direction n through a
point x where the stream moves along e(x) = x/|x| (the source's outflow) loses amplitude at the rate
kappa max(0, -n . e) per unit length. Along the straight path from piece l to piece j this integrates exactly to
    T_jl = exp(-kappa L_jl),   L_jl = r_l - (the smallest radius reached before the path turns outward or ends),
so kappa is the attenuation per unit of distance travelled inward. Outward paths are untouched; a wave from the far
side of a source loses exp(-kappa (r_l - impact parameter)). Every block of the coupling between the two pieces
(near and far field, monopole and dipole) is multiplied by T_jl, the force carries the derivative of T_jl with respect
to the receiving piece's position, and each piece's own terms are unchanged: emission is symmetric, so there is no
drag. The energy the pieces give the wave, Im(z^dagger M z), is split into what reaches infinity (the far-field
pattern of every piece, each direction attenuated by its own inward travel) and what the stream absorbs.
"""
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def path_inward(xl, xj):
    """L (distance travelled inward along the straight path from xl to xj) and dL/dxj."""
    d_vec = xj - xl; d = np.linalg.norm(d_vec)
    if d == 0:
        return 0.0, np.zeros(3)
    n = d_vec / d; b = float(n @ xl); rl = float(np.linalg.norm(xl))
    if b >= 0:
        return 0.0, np.zeros(3)
    if -b >= d:                                                       # the whole path moves inward: ends at xj
        rj = float(np.linalg.norm(xj))
        return rl - rj, -(xj / max(rj, 1e-12))
    foot = xl - b * n; p = float(np.linalg.norm(foot))                # closest approach to the centre
    dp = -b * foot / (d * max(p, 1e-12))
    return rl - p, -dp


def absorption_factors(x, kappa):
    """T (N, N): T[j, l] multiplies the wave from l at j; dT (3, N, N): its derivative with respect to x_j."""
    N = len(x); T = np.ones((N, N)); dT = np.zeros((3, N, N))
    if kappa == 0:
        return T, dT
    for j in range(N):
        for l in range(N):
            if j == l:
                continue
            L, dL = path_inward(x[l], x[j])
            if L > 0:
                T[j, l] = np.exp(-kappa * L)
                dT[:, j, l] = -kappa * T[j, l] * dL
    return T, dT


def far_field_absorbed(x, z, g0, g, kappa, ndir=2000):
    """Power reaching infinity when every piece's emission in direction n is attenuated by its inward travel along the
    ray, and the unattenuated (static) value for comparison."""
    from shared_wave_v16 import cmono, cdip, K, OMEGA
    i = np.arange(ndir) + 0.5; ph = np.arccos(1 - 2 * i / ndir); th = np.pi * (1 + 5 ** 0.5) * i
    nh = np.vstack([np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)]).T
    zz = z.reshape(-1, 4); q = cmono(g0) * zz[:, 0]; p = cdip(g) * zz[:, 1:]
    b = nh @ x.T                                                      # (ndir, N)
    r = np.linalg.norm(x, axis=1)[None, :]
    imp = np.sqrt(np.clip(r ** 2 - b ** 2, 0, None))
    att = np.where(b < 0, np.exp(-kappa * (r - imp)), 1.0)
    phase = np.exp(-1j * K * b)
    src = q[None, :] - 1j * K * (nh @ p.T)
    A = (phase * att * src).sum(1); A0 = (phase * src).sum(1)
    dP = OMEGA * K / (32 * np.pi ** 2) * np.abs(A) ** 2; dP0 = OMEGA * K / (32 * np.pi ** 2) * np.abs(A0) ** 2
    return float(dP.mean() * 4 * np.pi), float(dP0.mean() * 4 * np.pi), (dP[:, None] * nh).mean(0) * 4 * np.pi


def passivity(M):
    """The smallest and largest eigenvalue of the coupling's dissipative part (M - M^dagger)/2i: a passive medium,
    which can only take energy from the pieces, has none negative."""
    H = (M - M.conj().T) / 2j
    ev = np.linalg.eigvalsh(0.5 * (H + H.conj().T))
    return float(ev.min()), float(ev.max())


def configs(which):
    from one_matter_v17 import case, shells
    runs = []
    if which == 'distance':
        rec = shells((6.0, 9.0, 13.5, 20.0))
        for sd in (1, 2):
            for kappa in (2.3, 5.0):
                kw = dict(receivers=rec, T=16000.0, burn=8000.0, momentum=False, seed=sd, absorb=kappa)
                runs.append(case('single', 'cold', **kw))
                runs.append(case('single', 'free', k=8.0, **kw))
                runs.append(case('single', 'free', k=2.0, **kw))
                runs.append(case('single', 'collisional', k=8.0, nu=50.0, **kw))
        return runs
    if which == 'strong':
        # the strongest inward absorption of the reduced scan (code/absorb_reduced_v18.py), in the full model
        rec = shells((6.0, 9.0, 13.5, 20.0))
        for sd in (1, 2):
            kw = dict(receivers=rec, T=16000.0, burn=8000.0, momentum=False, seed=sd, absorb=20.0)
            runs.append(case('single', 'cold', **kw))
            runs.append(case('single', 'free', k=8.0, **kw))
            runs.append(case('single', 'free', k=2.0, **kw))
            runs.append(case('single', 'collisional', k=8.0, nu=50.0, **kw))
        return runs
    if which == 'quick':
        rec = shells((6.0, 13.5))
        kw = dict(receivers=rec, T=2000.0, burn=1000.0, momentum=False, seed=1, absorb=2.3)
        return [case('single', 'cold', **kw), case('single', 'free', k=8.0, **kw)]
    raise ValueError(which)


def main():
    from one_matter_v17 import simulate
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=4); ap.add_argument('--set', default='distance')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic(); runs = configs(args.set); res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(simulate, runs):
            res.append(r); c = r['cfg']; ab = r.get('absorb_booking', {})
            print(f"[{time.monotonic() - t0:5.0f} s] kappa {c['absorb']:3.1f} {c['structure']:>8} {c['tag']:12s} k {c['k']:4.1f} nu {c.get('nu', 0):4.1f} "
                  f"seed {c['seed']}: output {r['src_out']:.3e}, pull {r['pull']:+.3e}, keeping step {r['lead_src']:+.2f}, rhythm spread "
                  f"{r['rhythm_spread_sources']:.2e}, absorbed {ab.get('absorbed_fraction', float('nan')):.2f} of the power, "
                  f"E res {r['energy']['residual']:.1e}", flush=True)
            args.output.write_text(json.dumps(dict(experiment=f'round 18: one kind of matter, a stream that absorbs counter-moving waves, set {args.set}',
                                                   runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
