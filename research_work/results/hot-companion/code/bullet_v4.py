#!/usr/bin/env python3
"""Round 4: the Bullet Cluster's subcluster, with the companion's memory.

    python bullet_v4.py --output-dir ../run-bullet-v4

Round 3 left the subcluster's lensing 2.5 sigma low (0.07 against 0.20 +- 0.05). Round 4
adds no new rule. It follows two consequences of the rules already in place:

1. **The companion is slow.** It streams at u = 197 km/s (fitted on clusters), so the field
   now arriving at 100 kpc from a galaxy left it about 500 million years ago. The collision
   happened about 150 million years ago, so most of the field around each cluster was
   emitted BEFORE the collision.
2. **The companion keeps its emitter's velocity.** This is required, not assumed: galaxies
   move at hundreds of km/s, faster than u, yet all obey the same rotation law. A field that
   moved relative to a fixed medium would make a galaxy's gravity depend on its speed.

So the companion now surrounding each cluster is that of the cluster as it was before the
collision, still travelling with the pre-collision motion. The collisionless galaxies kept
that motion, so the old field rides with them. Only inside a fresh sphere of radius
u x (time since the gas was stopped), about 30 kpc, has the companion been rebuilt around
the matter where it is now. The ordinary pull and the release factor follow the matter where
it is now.

Before the collision each cluster was a settled cluster:
* **Main:** its present gas and stars, re-centred on its galaxies.
* **Subcluster:** its cool core (the bullet), its gas atmosphere, since stripped, and its
  galaxies, many since stripped, all centred on its galaxies.
  * Baryons: 1/8 of the main cluster's (merger reconstructions give 1:6 to 1:10).
  * Stars: a typical cluster fraction, 7% of the gas (X-COP clusters: 4-7% at R500, more in
    the core).
  * Star speeds: computed by our law from that cluster's own visible matter.

Nothing about the lensing is fitted. The constants are round 3's.
"""
from __future__ import annotations
import argparse, copy, json, time
from pathlib import Path
import numpy as np
import bullet_v3 as B
import law as L

G = B.G


def shells_on(r, dr, c):
    prof = B.rho_beta if c['kind'] == 'beta' else B.rho_nfw
    return c['M'] / B.mass3d(prof, (c['scale'],), c['rt']) * prof(r, c['scale'], c['rt']) * 4 * np.pi * r ** 2 * dr


def own_sigma_multi(gas_list, star_list, consts, n=500, iters=400, rmax=3000.0):
    """Stars' dispersion (all stellar components together) from the isotropic Jeans equation in
    the gravity our law makes from this cluster's own gas and stars, iterated to a fixed point."""
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']
    r = np.geomspace(1.0, rmax, n); dr = np.gradient(r)
    dmg = sum(shells_on(r, dr, c) for c in gas_list); dms = sum(shells_on(r, dr, c) for c in star_list)
    gN = G * np.cumsum(dmg + dms) / r ** 2
    W = L.shell_weights(r, r); rho_s = dms / (4 * np.pi * r ** 2 * dr)
    def jeans(g):
        tail = np.cumsum((rho_s * g * dr)[::-1])[::-1]
        return np.where(rho_s > 1e-30, tail / np.maximum(rho_s, 1e-300), 0.0)
    sig2 = jeans(gN + np.sqrt(a * gN))
    for _ in range(iters):
        S = G * (W @ (L.k_from_sig2(sig2, u) * dms)) / r ** 2
        g = gN + np.exp(-gN / (lam * a)) * np.sqrt(a * (gN + S))
        new = jeans(g)
        if np.max(np.abs(new - sig2) * dms) < 1e-7 * np.max(new * dms): sig2 = new; break
        sig2 = 0.5 * sig2 + 0.5 * new
    return r, np.sqrt(np.maximum(sig2, 0.0)), g


def vec_conv(n, dx, rmax=None, rmin=None):
    """Three convolutions giving G * sum m (x' - x)/|x' - x|^3, optionally restricted in distance."""
    m = 2 * n
    kk = ((np.arange(m) - (np.arange(m) >= n) * m) * dx).astype(np.float32)
    X = kk[:, None, None]; Y = kk[None, :, None]; Z = kk[None, None, :]
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2, dtype=np.float32)
    mask = np.ones_like(r, dtype=bool)
    if rmax is not None: mask &= r < rmax
    if rmin is not None: mask &= r >= rmin
    mask &= r > 0
    out = []
    for comp in (X, Y, Z):
        ker = np.where(mask, -comp / np.maximum(r, 1e-6) ** 3, 0.0).astype(np.float32)   # pull toward the source
        out.append(np.fft.rfftn(ker)); del ker
    del r, mask
    return out


def apply_vec(Ks, f, n):
    m = 2 * n
    pad = np.zeros((m, m, m), np.float32); pad[:n, :n, :n] = f
    P = np.fft.rfftn(pad); del pad
    return np.array([np.fft.irfftn(P * K, s=(m, m, m), axes=(0, 1, 2))[:n, :n, :n].astype(np.float32) for K in Ks])


def kappa_map_v4(current, ghost_gas, ghost_stars, pos, consts, n=192, dx=15.0, centre=(360., 50.),
                 fresh_kpc=30.0, heat=True, memory=True):
    """Lensing map. `current`: today's gas and stars (ordinary pull, release factor, baryons).
    `ghost_gas`: pre-collision gas components, where the old companion is centred.
    `ghost_stars`: list of (component, (r, sigma)) pre-collision stars with their speeds."""
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']; gd = lam * a
    x = (np.arange(n) - n / 2 + 0.5) * dx + centre[0]
    y = (np.arange(n) - n / 2 + 0.5) * dx + centre[1]
    z = (np.arange(n) - n / 2 + 0.5) * dx
    dV = dx ** 3
    rho_now = B.build_density(current, pos, x, y, z, 'gas') + B.build_density(current, pos, x, y, z, 'st')
    rho_gas_now = B.build_density(current, pos, x, y, z, 'gas')
    rho_ghost_gas = B.build_density(ghost_gas, pos, x, y, z, 'gas') if memory else rho_gas_now
    rho_ghost_st = np.zeros_like(rho_now); krho = np.zeros_like(rho_now)
    for comp, prof in ghost_stars:
        cname = {'c': comp}
        rho_c = B.build_density(cname, pos, x, y, z, 'c')
        cx, cy = pos[comp['centre']]
        r3 = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
        rho_ghost_st += rho_c
        if heat:
            krho += L.k_from_sig2(np.interp(r3, prof[0], prof[1]) ** 2, u).astype(np.float32) * rho_c
        del r3, rho_c
    inv_r = B.Conv(n, dx, lambda r: 1.0 / r, B.cube_average(lambda r: 1 / r) / dx)
    gN = -np.array(np.gradient(-G * inv_r(rho_now * dV), dx))
    # coherent companion flow: stars (unchanged motion) + gas, old part from the ghost
    F = -np.array(np.gradient(-G * inv_r((rho_ghost_gas + rho_ghost_st) * dV), dx))
    g_hot = -np.array(np.gradient(-G * inv_r(krho * dV), dx)) if heat else 0.0
    del inv_r
    if memory and fresh_kpc > 0:
        # inside the fresh sphere around the matter where it is now, the companion was emitted
        # after the collision: swap the ghost gas's near-field for today's gas's near-field
        Kin = vec_conv(n, dx, rmax=fresh_kpc)
        F = F + G * dV * (apply_vec(Kin, rho_gas_now, n) - apply_vec(Kin, rho_ghost_gas, n))
        del Kin
    inv_r2 = B.Conv(n, dx, lambda r: 1.0 / r ** 2, B.cube_average(lambda r: 1 / r ** 2) / dx ** 2)
    S = G * inv_r2(krho * dV) if heat else 0.0
    del inv_r2
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    Fmag = np.sqrt(np.sum(F ** 2, axis=0)) + 1e-30
    hmag = np.sqrt(np.sum(g_hot ** 2, axis=0)) if heat else 0.0
    extra = np.exp(-mag / gd) * np.sqrt(a * (Fmag + S))
    h = gN + extra * (F + g_hot) / (Fmag + hmag + 1e-30)
    divh = sum(np.gradient(h[i], dx, axis=i) for i in range(3))
    Sig_eff = (-divh / (4 * np.pi * G)).sum(axis=2) * dx
    Sig_b = rho_now.sum(axis=2) * dx
    return x, y, Sig_eff, Sig_b


def pre_collision_models(comps, ratio=8.0, star_frac=0.07, atm_core=150.0):
    """Settled clusters before the collision, from standard cluster ingredients."""
    Mmain_b = comps['gas_main']['M'] + comps['st_main']['M']
    Msub_b = Mmain_b / ratio
    gas_pre = Msub_b / (1 + star_frac); stars_pre = Msub_b - gas_pre
    atm = gas_pre - comps['gas_sub']['M']                       # the atmosphere the bullet lost
    sats = max(stars_pre - comps['st_sub']['M'], 0.0)           # galaxies beyond today's core
    main_gas = dict(comps['gas_main']); main_gas['M'] = comps['gas_main']['M'] - atm; main_gas['centre'] = 'main_bcg'
    ghost_gas = dict(gas_main_ghost=main_gas,
                     gas_bullet_ghost=dict(comps['gas_sub'], centre='sub_bcg'),
                     gas_atm_ghost=dict(kind='beta', M=atm, scale=atm_core, rt=1000.0, centre='sub_bcg'))
    sub_sats = dict(kind='nfw', M=sats, scale=atm_core, rt=1000.0, centre='sub_bcg')
    return ghost_gas, sub_sats, dict(sub_baryons=Msub_b, sub_gas=gas_pre, sub_atmosphere=atm, sub_stars=stars_pre, sub_satellites=sats)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--n', type=int, default=192); ap.add_argument('--dx', type=float, default=15.0)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((Path(__file__).resolve().parent.parent / 'run-v3/results.json').read_text())['constants']
    pos = B.positions(); comps, _ = B.fit_components(pos); scrit = B.sigma_crit()
    obs = B.OBS['kappa']
    runs = []

    def run(label, ratio=8.0, star_frac=0.07, atm_core=150.0, fresh=30.0, memory=True, ml=2.0):
        # ml: stellar mass-to-light ratio for the published star masses (Clowe et al.: 2, range 0.5-3,
        # upper limits because foreground galaxies were not removed)
        cm = copy.deepcopy(comps)
        for k in ('st_main', 'st_sub'): cm[k]['M'] = comps[k]['M'] * ml / 2.0
        ghost_gas, sub_sats, pre = pre_collision_models(cm, ratio, star_frac, atm_core)
        rm, sm, _ = own_sigma_multi([cm['gas_main']], [cm['st_main']], consts)
        gas_sub_pre = [dict(cm['gas_sub']), ghost_gas['gas_atm_ghost']]
        rs, ss, gs = own_sigma_multi(gas_sub_pre, [cm['st_sub'], sub_sats], consts)
        ghost_stars = [(dict(cm['st_main']), (rm, sm)), (dict(cm['st_sub']), (rs, ss))]
        if sub_sats['M'] > 0: ghost_stars.append((sub_sats, (rs, ss)))
        if not memory:
            ghost_stars = [(dict(cm['st_main']), (rm, sm)), (dict(cm['st_sub']), (rs, ss))]
        x, y, Se, Sb = kappa_map_v4(cm, ghost_gas, ghost_stars, pos, consts, n=args.n, dx=args.dx, fresh_kpc=fresh, memory=memory)
        kap = Se / scrit
        dec = B.clowe_decomposition(x, y, kap, pos)
        pk = B.peaks(x, y, kap, pos)
        near = lambda key: min(pk, key=lambda d: d[f'dist_{key}'])
        mpk, spk = near('main_bcg'), near('sub_bcg')
        row = dict(case=label, ratio=ratio, star_frac=star_frac, atm_core=atm_core, fresh_kpc=fresh, memory=memory, mass_to_light=ml,
                   pre_collision=pre,
                   star_speeds_kms=dict(main={f'{r0}kpc': float(np.interp(r0, rm, sm)) for r0 in (20, 50, 100, 200, 400)},
                                        sub_before={f'{r0}kpc': float(np.interp(r0, rs, ss)) for r0 in (20, 50, 100, 200, 400)}),
                   clowe_style=dec, kappa_min=float(kap.min()),
                   main_peak=dict(dist_bcg=mpk['dist_main_bcg'], dist_gas=mpk['dist_main_plasma'],
                                  fraction_stars_to_gas=B.along_axis_fraction(mpk, pos['main_bcg'], pos['main_plasma'])),
                   sub_peak=dict(dist_bcg=spk['dist_sub_bcg'], dist_gas=spk['dist_sub_plasma'],
                                 fraction_stars_to_gas=B.along_axis_fraction(spk, pos['sub_bcg'], pos['sub_plasma'])),
                   peaks=pk[:6])
        runs.append(row)
        np.savez_compressed(out / f'kappa_{len(runs)}.npz', x=x, y=y, kappa=kap.astype(np.float32))
        print(f"{label}\n   main {dec['main_bcg']:.3f} ({obs['main_bcg'][0]} +- {obs['main_bcg'][1]})  sub {dec['sub_bcg']:.3f} ({obs['sub_bcg'][0]} +- {obs['sub_bcg'][1]})"
              f"  gas {dec['main_plasma']:.3f} / {dec['sub_plasma']:.3f} (0.05 / 0.02)   min kappa {kap.min():.3f}\n"
              f"   main peak {mpk['dist_main_bcg']:.0f} kpc from its galaxies; sub peak {spk['dist_sub_bcg']:.0f} kpc from its galaxies, {spk['dist_sub_plasma']:.0f} from its gas\n"
              f"   subcluster before: baryons {pre['sub_baryons']:.2e}, stars {pre['sub_stars']:.2e}; its star speeds at 50/100/200 kpc: "
              + ', '.join(f"{np.interp(r0, rs, ss):.0f}" for r0 in (50, 100, 200)) + ' km/s', flush=True)

    run('A. no memory (round 3 as it stands)', memory=False)
    run('B. memory: pre-collision clusters, sub baryons 1/8 of main, stars 7% of gas')
    run('C. memory, sub baryons 1/10 of main', ratio=10.0)
    run('D. memory, sub baryons 1/6 of main', ratio=6.0)
    run('E. memory, stars 5% of gas', star_frac=0.05)
    run('F. memory, stars 10% of gas', star_frac=0.10)
    run('G. memory, compact atmosphere (100 kpc core)', atm_core=100.0)
    run('H. memory, gas stopped 300 Myr ago (fresh sphere 60 kpc)', fresh=60.0)
    run('I. memory, star mass-to-light 1.5 for both clusters', ml=1.5)
    run('J. memory, star mass-to-light 1.0 for both clusters', ml=1.0)
    payload = dict(experiment='Bullet Cluster with the companion memory (round 4)', data=B.OBS, constants=consts,
                   sigma_crit=scrit, runs=runs, grid=dict(n=args.n, dx_kpc=args.dx), seconds=time.monotonic() - t0)
    (out / 'bullet_v4.json').write_text(json.dumps(payload, indent=2, default=float) + '\n')
    print(f'{time.monotonic() - t0:.0f} s')


if __name__ == '__main__':
    main()
