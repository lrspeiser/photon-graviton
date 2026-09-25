"""Round 20, step 3: in the energy-shift model, how does the pull on a separate probe depend on the source's mass,
the distance and the number of excitations?

The supplied check (energy-shift-v20/checks.py) diagonalised eight identical two-level pieces sharing one medium, with
the medium's exchange kernel J(r) = -C exp(-kappa r)/r, and found attraction through pair correlations at every
partial filling. It measured the internal forces of one fixed cloud, not how an outside body is pulled by a source of
varying mass, which the review lists as the decisive next measurement (its protocol, item 6).

Here the same Hamiltonian and constants,
    H = omega_a n + sum_{i<j} J_ij (s_i^+ s_j^- + s_j^+ s_i^-),   J = -C exp(-kappa r)/r,   C = 0.04/(4 pi),
    kappa = 1.3537825795626388,
describe a source cluster of Ns pieces at unit mean density (Ns = 3 ... 13; its mass) and one probe piece at distance
D from the cluster's centre. Each total excitation number n is diagonalised exactly; the lowest state of each sector is
used, as in the supplied check (how it would form is not shown). The pull on the probe is -<dH/dR_probe> (exact
Hellmann-Feynman, checked by a finite difference of the energy).

Measured: the probe's pull against D (the distance exponent) and against Ns at fixed D (the mass exponent), for a
fixed number of excitations and for fixed fractions of the pieces excited, over several random clusters.

    python code/finite_population_probe_v20.py --output run-finite-population-v20/finite_population_probe_v20.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, time
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

C = 0.04 / (4 * np.pi)
KAPPA = 1.3537825795626388
OMEGA = 1.2


def cluster(Ns, seed, dmin=0.55):
    """Ns positions at unit mean density inside a ball (radius (3 Ns / 4 pi)^(1/3)), no two closer than dmin."""
    rng = np.random.default_rng(seed)
    R = (3 * Ns / (4 * np.pi)) ** (1 / 3)
    pts = []
    while len(pts) < Ns:
        p = rng.uniform(-R, R, 3)
        if np.linalg.norm(p) <= R and all(np.linalg.norm(p - q) >= dmin for q in pts):
            pts.append(p)
    pts = np.array(pts)
    return pts - pts.mean(0)


def sector_ground(pos, n):
    """Lowest state of the n-excitation sector: energy, state, basis and the pair list."""
    N = len(pos)
    states = [sum(1 << i for i in c) for c in combinations(range(N), n)]
    idx = {s: a for a, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    pairs = []
    for i in range(N):
        for j in range(i + 1, N):
            d = pos[i] - pos[j]; r = float(np.linalg.norm(d))
            J = -C * np.exp(-KAPPA * r) / r
            rows, cols = [], []
            for a, s in enumerate(states):
                if ((s >> i) & 1) != ((s >> j) & 1):
                    b = idx[s ^ (1 << i) ^ (1 << j)]
                    H[a, b] += J
                    rows.append(a); cols.append(b)
            pairs.append((i, j, np.array(rows, int), np.array(cols, int)))
    vals, vecs = eigh(H)
    return vals, vecs[:, 0], pairs


def probe_pull(pos, n, p):
    """Pull on piece p toward the cluster centre (positive = attraction), its pair correlations, and the energy gap."""
    vals, gs, pairs = sector_ground(pos, n)
    F = np.zeros(3)
    corr_p = []
    for i, j, rows, cols in pairs:
        if p not in (i, j) or not len(rows):
            continue
        q = j if i == p else i
        c = float(np.dot(gs[rows], gs[cols]))            # <s_i^+ s_j^- + s_j^+ s_i^->
        d = pos[p] - pos[q]; r = float(np.linalg.norm(d))
        dJdr = C * np.exp(-KAPPA * r) * (1 / r ** 2 + KAPPA / r)
        F += -dJdr * c * d / r
        corr_p.append(c)
    centre = np.delete(pos, p, axis=0).mean(0)
    u = (centre - pos[p]) / np.linalg.norm(centre - pos[p])
    gap = float(vals[1] - vals[0]) if len(vals) > 1 else float('nan')
    return float(F @ u), float(np.sum(corr_p)), gap, float(vals[0])


def job(args):
    Ns, seed, D, n = args
    src = cluster(Ns, seed)
    direction = np.array([0.0, 0.0, 1.0])
    pos = np.vstack([src, D * direction])
    Fr, csum, gap, E = probe_pull(pos, n, Ns)
    return dict(Ns=Ns, seed=seed, D=D, n=n, pull=Fr, probe_corr_sum=csum, gap=gap, energy=E)


def fd_check(Ns=6, seed=1, D=2.5, n=3, h=1e-6):
    src = cluster(Ns, seed)
    pos = np.vstack([src, [0, 0, D]])
    Fr, _, _, _ = probe_pull(pos, n, Ns)
    pp = pos.copy(); pm = pos.copy(); pp[Ns, 2] += h; pm[Ns, 2] -= h
    Ep = sector_ground(pp, n)[0][0]; Em = sector_ground(pm, n)[0][0]
    Fz = -(Ep - Em) / (2 * h)
    centre = src.mean(0); u = (centre - pos[Ns]) / np.linalg.norm(centre - pos[Ns])
    return dict(hellmann_feynman=Fr, finite_difference_along_z=float(Fz * u[2]), difference=float(abs(Fr - Fz * u[2])))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--procs', type=int, default=3); ap.add_argument('--seeds', type=int, default=6)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    sizes = list(range(3, 14))
    Ds = [2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
    jobs = []
    for Ns in sizes:
        N = Ns + 1
        fillings = sorted({1, 2, max(1, round(N / 4)), max(1, round(N / 2))})
        for seed in range(1, args.seeds + 1):
            for D in Ds:
                for n in fillings:
                    jobs.append((Ns, seed, D, n))
    with Pool(args.procs) as pool:
        rows = pool.map(job, jobs, chunksize=4)
    check = fd_check()
    out = dict(experiment='round 20: pull on a probe in the energy-shift (exchange) model, against source size, distance and filling',
               constants=dict(C=C, kappa=KAPPA, omega=OMEGA), fd_check=check, rows=rows)

    # summaries: median over clusters of the pull; exponents by finite differences of the medians
    def med(Ns, D, n):
        v = [r['pull'] for r in rows if r['Ns'] == Ns and r['D'] == D and r['n'] == n]
        return float(np.median(v)) if v else float('nan')

    summ = {}
    for label, nfun in (('one excitation', lambda N: 1), ('two excitations', lambda N: 2),
                        ('a quarter excited', lambda N: max(1, round(N / 4))), ('half excited', lambda N: max(1, round(N / 2)))):
        tab = {D: [med(Ns, D, nfun(Ns + 1)) for Ns in sizes] for D in Ds}
        mass_exp = {}
        for D in Ds:
            y = np.array(tab[D]); x = np.array(sizes, float)
            ok = y > 0
            mass_exp[D] = float(np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)[0]) if ok.sum() >= 3 else float('nan')
        dist_exp = {}
        for Ns in (5, 9, 13):
            y = np.array([med(Ns, D, nfun(Ns + 1)) for D in Ds]); x = np.array(Ds)
            ok = y > 0
            dist_exp[Ns] = [float(-(np.log(y[i + 1]) - np.log(y[i])) / (np.log(x[i + 1]) - np.log(x[i]))) if ok[i] and ok[i + 1] else float('nan')
                            for i in range(len(Ds) - 1)]
        attractive = float(np.mean([r['pull'] > 0 for r in rows if r['n'] == nfun(r['Ns'] + 1)]))
        summ[label] = dict(median_pull={str(D): v for D, v in tab.items()}, mass_exponent_by_D={str(D): v for D, v in mass_exp.items()},
                           distance_exponent_by_Ns={str(k): v for k, v in dist_exp.items()}, fraction_attractive=attractive)
        print(f"{label:18s} attractive in {100 * attractive:.0f}% of cases; mass exponent (log-fit over Ns = 3-13) at D = 3, 4, 6: "
              + ', '.join(f'{mass_exp[D]:.2f}' for D in (3.0, 4.0, 6.0)) + '; distance exponent (Ns = 9) for D 2.5 -> 6: '
              + ', '.join(f'{v:.1f}' for v in dist_exp[9]), flush=True)
    out['summary'] = summ
    out['yukawa_first_order_exponent'] = {str(D): 2 + (KAPPA * D) ** 2 / (1 + KAPPA * D) for D in Ds}
    out['seconds'] = time.monotonic() - t0
    print('Hellmann-Feynman vs finite difference:', check)
    print('first-order Yukawa exponent 2 + (kD)^2/(1+kD) at D = 2.5 ... 6:', ', '.join(f'{v:.1f}' for v in out['yukawa_first_order_exponent'].values()))
    args.output.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {args.output} ({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()
