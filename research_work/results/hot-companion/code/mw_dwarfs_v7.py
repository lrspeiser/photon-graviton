#!/usr/bin/env python3
"""Round 7: the Milky Way's dwarf spheroidal galaxies under our law, with the Galaxy's pull.

    python mw_dwarfs_v7.py --output-dir ../run-mw-dwarfs-v7

Dwarf spheroidals are the faintest galaxies we know: a few hundred thousand to twenty million
suns spread over a few hundred parsecs, orbiting the Milky Way 76-254 kpc out. Their own pull is
10^-12 to 10^-14 m/s^2, deep in the regime where galaxies stop following Newton, and the
Galaxy's pull on them (a few x 10^-13 m/s^2) is often larger than their own. That makes them a
classic test: a law that depends on the TOTAL pull (as ours and MOND do) predicts that a dwarf
bathed in the Galaxy's field is less boosted than an isolated one (the 'external field effect').

Our law in a dwarf at distance D from the Galaxy (all pulls far below g_d, so the release factor
is 1 and only the square root and the direction matter):
    g_N   = g_N,MW(D) e  +  g_i(r) (-n)           e: toward the Galaxy;  n: outward in the dwarf
    g_hot = S_MW(D) e    +  k_d g_i(r) (-n)        the Galaxy's hot bulge and stellar halo seen from D
    S     = S_MW(D)      +  k_d g_i(r)
    h     = g_N + sqrt(a (|g_N| + S)) (g_N + g_hot) / (|g_N| + |g_hot|)
The pull that holds the dwarf's stars, averaged over directions, follows exactly from Gauss's law
(the flux of the true pull through a sphere equals the flux of h): gbar(r) = <-n . h>. The
stars' spread of speeds then follows from the isotropic Jeans equation for a Plummer sphere
(scale = the measured half-light radius), weighted by light over the whole dwarf.

The same code gives QUMOND with the 'simple' function (a0 = 1.2 x 10^-10 m/s^2) and Newton.
Inputs: ../data/mw_dwarfs.json (distances, luminosities, half-light radii, measured speed spreads,
with sources). Stellar mass-to-light ratio: V-band 2 (old, metal-poor stars), varied 1-3.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import mw_model as MW
import law as L                                    # the heat weight (round 14: its exponent)

G, K_SI = MW.G, MW.K_SI
A0_SI = 1.2e-10


def plummer(M, b):
    rho = lambda r: 3 * M / (4 * np.pi * b ** 3) * (1 + (r / b) ** 2) ** -2.5
    gi = lambda r: G * M * r / (r ** 2 + b ** 2) ** 1.5
    return rho, gi


def mean_pull(r, gi, env, law, consts, k_d=0.0, n_mu=256):
    """Direction-averaged inward pull <-n . h> at radius r (Gauss's law makes this exact)."""
    mu, w = np.polynomial.legendre.leggauss(n_mu)          # mu = n . e
    ge, Se = env['gN'], env['S']
    g_i = gi(r)[:, None]
    # components along n (radial) and perpendicular, with e = mu n + sqrt(1-mu^2) t
    gN_par = ge * mu[None, :] - g_i                        # along n
    gN_perp = ge * np.sqrt(1 - mu ** 2)[None, :] + 0 * g_i
    g = np.sqrt(gN_par ** 2 + gN_perp ** 2) + 1e-300
    if law == 'newton':
        extra_par = 0 * g
    elif law == 'mond':
        a0 = A0_SI / K_SI
        nu = 0.5 + np.sqrt(0.25 + a0 / g)
        extra_par = (nu - 1) * gN_par
    else:
        a, gd = consts['a_code'], consts['lam'] * consts['a_code']
        S = Se + k_d * g_i
        F = np.exp(-g / gd) * np.sqrt(a * (g + S))
        gh_par = Se * mu[None, :] - k_d * g_i
        gh_perp = Se * np.sqrt(1 - mu ** 2)[None, :] + 0 * g_i
        gh = np.sqrt(gh_par ** 2 + gh_perp ** 2)
        extra_par = F * (gN_par + gh_par) / (g + gh)
    # inward pull = -(component along n)
    return gi(r) + 0.5 * ((-extra_par) * w[None, :]).sum(1)


def sigma_los(M, b, env, law, consts, k_d=0.0):
    rho, gi = plummer(M, b)
    r = np.geomspace(1e-3 * b, 300 * b, 3000)
    gbar = mean_pull(r, gi, env, law, consts, k_d)
    integrand = rho(r) * gbar
    # rho sigma_r^2 (r) = int_r^inf rho gbar dr'
    cum = np.concatenate([np.cumsum((0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r))[::-1])[::-1], [0.0]])
    sr2 = cum / rho(r)
    dV = 4 * np.pi * r ** 2 * np.gradient(r)
    glob = np.sqrt(np.sum(rho(r) * sr2 * dV) / np.sum(rho(r) * dV))
    # inside the half-light radius (projected): Sigma sigma_los^2 = 2 int_R^inf rho sr2 r dr / sqrt(r^2 - R^2)
    Rp = np.linspace(0.01 * b, b, 60)
    num = []; den = []
    for R in Rp:
        m = r > R * 1.000001
        x = r[m]; f = rho(x) * sr2[m] * x / np.sqrt(x ** 2 - R ** 2); fr = rho(x) * x / np.sqrt(x ** 2 - R ** 2)
        num.append(2 * np.trapezoid(f, x)); den.append(2 * np.trapezoid(fr, x))
    num, den = np.array(num), np.array(den)
    inner = np.sqrt(np.trapezoid(num * Rp, Rp) / np.trapezoid(den * Rp, Rp))
    return float(glob), float(inner), float(np.interp(b, r, gbar) * K_SI), float(gi(np.array([b]))[0] * K_SI)


def galaxy_env(D, mw, consts):
    """The Galaxy's Newtonian pull, heat and heat-weighted pull at Galactocentric distance D (kpc)."""
    s, dm, dmk = mw['s'], mw['dm'], mw['dmk']
    gN = G * np.sum(dm[s < D]) / D ** 2 + 0.0
    x = s / D
    w = np.where(x < 1, np.arctanh(np.clip(x, 1e-12, 1 - 1e-12)) / np.clip(x, 1e-12, None),
                 0.5 * np.log((x + 1) / np.maximum(x - 1, 1e-12)) / x)
    S = G * np.sum(w * dmk) / D ** 2
    return dict(gN=gN, S=S)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((HERE.parent / 'run-v3/results.json').read_text())['constants']
    data = json.loads((HERE.parent / 'data/mw_dwarfs.json').read_text())
    mwj = json.loads((HERE.parent / 'run-milky-way-v7/milky_way_v7.json').read_text())['galaxy_mass_profile']
    mw = dict(s=np.array(mwj['s_kpc']), dm=np.array(mwj['dm_Msun']), dmk=np.array(mwj['dmk_Msun']))
    u = consts['u_kms']
    rows = []
    print(f"{'dwarf':12s} {'D':>5s} {'L_V':>9s} {'r_h':>6s} {'obs':>11s} | {'ours':>17s} | {'MOND':>17s} | {'Newton':>6s}   (M/L_V 1 / 2 / 3; global light-weighted)", flush=True)
    for d in data['dwarfs']:
        env = galaxy_env(d['D_gc_kpc'], mw, consts)
        b = d['r_h_pc'] / 1000.0
        res = {}
        for law in ('ours', 'mond', 'newton'):
            res[law] = {}
            for ml in (1.0, 2.0, 3.0):
                M = ml * d['L_V']
                k_d = L.k_from_sig2(d['sigma_obs'] ** 2, u)
                glob, inner, gbar_b, gi_b = sigma_los(M, b, env, law, consts, k_d)
                res[law][f'ML{ml:g}'] = dict(sigma_global=glob, sigma_inside_rh=inner, pull_at_rh_SI=gbar_b, newton_pull_at_rh_SI=gi_b)
            # isolated (no Galaxy) for reference
            glob_iso, _, _, _ = sigma_los(2.0 * d['L_V'], b, dict(gN=0.0, S=0.0), law, consts, L.k_from_sig2(d['sigma_obs'] ** 2, u))
            res[law]['isolated_ML2'] = glob_iso
        rows.append(dict(name=d['name'], env_SI=dict(gN=env['gN'] * K_SI, S=env['S'] * K_SI), internal_gN_at_rh_SI=res['newton']['ML2']['newton_pull_at_rh_SI'],
                         observed=dict(sigma=d['sigma_obs'], err=d['sigma_err']), predictions=res, source=d.get('source', '')))
        f = lambda law: '/'.join(f"{res[law][f'ML{m:g}']['sigma_global']:.1f}" for m in (1.0, 2.0, 3.0))
        print(f"{d['name']:12s} {d['D_gc_kpc']:5.0f} {d['L_V']:9.2e} {d['r_h_pc']:6.0f} {d['sigma_obs']:5.1f}+-{d['sigma_err']:3.1f} | {f('ours'):>17s} | {f('mond'):>17s} | {res['newton']['ML2']['sigma_global']:6.1f}"
              f"   isolated M/L 2: ours {res['ours']['isolated_ML2']:.1f}, MOND {res['mond']['isolated_ML2']:.1f};  Galaxy's pull {env['gN'] * K_SI:.2e}, own pull {res['newton']['ML2']['newton_pull_at_rh_SI']:.2e}", flush=True)
    # summary: mean log ratio predicted/observed at M/L 2, and chi2
    summ = {}
    for law in ('ours', 'mond', 'newton'):
        lr = [np.log10(r['predictions'][law]['ML2']['sigma_global'] / r['observed']['sigma']) for r in rows]
        chi = [((r['predictions'][law]['ML2']['sigma_global'] - r['observed']['sigma']) / r['observed']['err']) ** 2 for r in rows]
        summ[law] = dict(mean_log_ratio=float(np.mean(lr)), rms_log_ratio=float(np.sqrt(np.mean(np.square(lr)))), chi2=float(np.sum(chi)), n=len(rows))
        print(f"{law:7s}: mean log10(predicted/observed) {np.mean(lr):+.3f}, rms {np.sqrt(np.mean(np.square(lr))):.3f}, chi2 {np.sum(chi):.1f} for {len(rows)} dwarfs (M/L_V 2)", flush=True)
    # what would have to change: weaken the Galaxy's hold (its pull and heat entering each dwarf's law) by a factor c
    efe = {}
    for c in (1.0, 0.5, 0.3, 0.1, 0.03, 0.0):
        sig = []
        for d in data['dwarfs']:
            env = galaxy_env(d['D_gc_kpc'], mw, consts)
            env = dict(gN=c * env['gN'], S=c * env['S'])
            sig.append(sigma_los(2.0 * d['L_V'], d['r_h_pc'] / 1000.0, env, 'ours', consts, L.k_from_sig2(d['sigma_obs'] ** 2, u))[0])
        chi = float(sum(((s_ - d['sigma_obs']) / d['sigma_err']) ** 2 for s_, d in zip(sig, data['dwarfs'])))
        efe[str(c)] = dict(sigma=dict(zip([d['name'] for d in data['dwarfs']], sig)), chi2=chi,
                           mean_log_ratio=float(np.mean([np.log10(s_ / d['sigma_obs']) for s_, d in zip(sig, data['dwarfs'])])))
        print(f"Galaxy's hold x{c}: " + ', '.join(f"{d['name']} {s_:.1f}" for s_, d in zip(sig, data['dwarfs'])) + f"; chi2 {chi:.0f}", flush=True)
    (out / 'mw_dwarfs_v7.json').write_text(json.dumps(dict(hold_scan=efe,
        experiment='Milky Way dwarf spheroidals: stars speed spread under our law with the Galaxy pull, vs MOND and Newton (round 7)',
        constants=consts, a0_mond_SI=A0_SI, dwarfs=rows, summary=summ, seconds=time.monotonic() - t0), indent=2) + '\n')


if __name__ == '__main__':
    main()
