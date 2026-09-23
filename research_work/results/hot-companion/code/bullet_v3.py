#!/usr/bin/env python3
"""The Bullet Cluster (1E0657-558) under round-3 hot-companion gravity, on published data.

    python bullet_v3.py --output-dir ../run-bullet-v3

Data (Clowe et al. 2006, ApJ 648, L109, Table 2; 100 kpc apertures, 4.413 kpc/arcsec):
                         RA (J2000)   Dec (J2000)    gas (1e12)   stars (1e12)   kappa
    main cluster BCG     06:58:35.3   -55:56:56.3    5.5 +- 0.6   0.54 +- 0.08   0.36 +- 0.06
    main cluster plasma  06:58:30.2   -55:56:35.9    6.6 +- 0.7   0.23 +- 0.02   0.05 +- 0.06 (excess)
    subcluster BCG       06:58:16.0   -55:56:35.1    2.7 +- 0.3   0.58 +- 0.09   0.20 +- 0.05
    subcluster plasma    06:58:21.2   -55:56:30.0    5.8 +- 0.6   0.12 +- 0.01   0.02 +- 0.06 (excess)
The BCG kappas are after subtracting the other peak's circularly symmetric profile; the
plasma kappas are what is left after subtracting both.  Stellar masses assume M/L_I = 2
(range 0.5-3) and are upper limits (no colour selection).
Galaxy velocity dispersions (Barrena et al. 2002, A&A 386, 816): main cluster 1249 +109/-100
km/s (71 galaxies); subcluster 212 +67/-52 km/s (7 galaxies; described there as the stripped
remnant core of a more massive cluster).

Model: two gas clouds (beta = 2/3) and two stellar concentrations (NFW-shaped), each fitted to
reproduce the four gas and four stellar aperture masses above.  Nothing about the lensing is
fitted: the constants a, lambda, u come from round 3 (galaxies and X-COP clusters).  Gas
feeds only the ordinary (coherent) part; stars feed the heat term with k = 3 sigma^2/u^2.
The lensing map is the projected effective density  -div h / (4 pi G)  of the field equation
(round 2), converted to kappa with sources at z = 1 in the same cosmology used for the
published masses (H0 = 70, Omega_m = 0.3).
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.integrate import quad
from scipy.ndimage import gaussian_filter, maximum_filter

G = 4.30091727003628e-6                   # kpc (km/s)^2 / Msun
K = 1e6 / 3.0856775814913673e19           # (km/s)^2/kpc -> m/s^2
SCALE = 4.413                             # kpc per arcsec (Clowe et al. 2006)
DEC0 = -55.94


def radec_to_kpc(ra, dec):
    h, m, s = [float(x) for x in ra.split(':')]; d, dm, ds = [float(x) for x in dec.split(':')]
    ra_s = (h * 3600 + m * 60 + s) * 15 * np.cos(np.radians(DEC0))          # arcsec
    dec_s = -(abs(d) * 3600 + dm * 60 + ds)
    return ra_s, dec_s


POS_RAW = {'main_bcg': ('06:58:35.3', '-55:56:56.3'), 'main_plasma': ('06:58:30.2', '-55:56:35.9'),
           'sub_bcg': ('06:58:16.0', '-55:56:35.1'), 'sub_plasma': ('06:58:21.2', '-55:56:30.0')}
OBS = dict(gas=dict(main_bcg=5.5, main_plasma=6.6, sub_bcg=2.7, sub_plasma=5.8),
           gas_err=dict(main_bcg=0.6, main_plasma=0.7, sub_bcg=0.3, sub_plasma=0.6),
           stars=dict(main_bcg=0.54, main_plasma=0.23, sub_bcg=0.58, sub_plasma=0.12),
           stars_err=dict(main_bcg=0.08, main_plasma=0.02, sub_bcg=0.09, sub_plasma=0.01),
           kappa=dict(main_bcg=(0.36, 0.06), main_plasma=(0.05, 0.06), sub_bcg=(0.20, 0.05), sub_plasma=(0.02, 0.06)))


def positions():
    ref = radec_to_kpc(*POS_RAW['main_bcg'])
    out = {}
    for k, (ra, dec) in POS_RAW.items():
        x, y = radec_to_kpc(ra, dec)
        out[k] = np.array([-(x - ref[0]) * SCALE, (y - ref[1]) * SCALE])   # x = west, y = north
    return out


def sigma_crit(zl=0.296, zs=1.0, h=0.7, om=0.3):
    """Critical surface density in Msun/kpc^2, flat LCDM."""
    c = 299792.458; H0 = 100 * h
    chi = lambda z: quad(lambda x: 1 / np.sqrt(om * (1 + x) ** 3 + 1 - om), 0, z)[0] * c / H0   # Mpc
    Dl = chi(zl) / (1 + zl); Ds = chi(zs) / (1 + zs); Dls = (chi(zs) - chi(zl)) / (1 + zs)
    c2_4piG = (c ** 2) / (4 * np.pi * G)      # Msun/kpc
    return c2_4piG * Ds / (Dl * Dls) / 1e3     # Dl*Dls/Ds in Mpc -> kpc


# ---------------------------------------------------------------- 3D profiles and projections
def rho_beta(r, rc, rt):
    return np.where(r < rt, 1.0 / (1 + (r / rc) ** 2), 0.0)


def rho_nfw(r, rs, rt):
    x = np.maximum(r, 1e-3) / rs
    return np.where(r < rt, 1.0 / (x * (1 + x) ** 2), 0.0)


def mass3d(rho, p, rt):
    r = np.geomspace(1e-2, rt, 4000)
    return np.trapezoid(4 * np.pi * r ** 2 * rho(r, *p, rt), r)


def surface(rho, p, rt, R):
    """Projected surface density of a unit-normalised profile at radii R."""
    z = np.concatenate([[0], np.geomspace(1e-2, rt, 1500)])
    out = np.empty_like(R, dtype=float)
    for i, RR in enumerate(R):
        out[i] = 2 * np.trapezoid(rho(np.sqrt(RR ** 2 + z ** 2), *p, rt), z)
    return out


def aperture(rho, p, rt, offset, rap=100.0):
    """Projected mass (unit normalisation) inside a circle of radius rap at distance offset."""
    R = np.linspace(0, offset + rap + 1, 700) if offset > 0 else np.linspace(0, rap, 400)
    Sig = surface(rho, p, rt, np.maximum(R, 0.5))
    if offset == 0:
        return np.trapezoid(2 * np.pi * R * Sig, R)
    # arc length of the circle |x - offset| < rap at radius R
    with np.errstate(invalid='ignore'):
        cosang = (R ** 2 + offset ** 2 - rap ** 2) / (2 * R * offset)
    ang = np.where(cosang <= -1, np.pi, np.where(cosang >= 1, 0.0, np.arccos(np.clip(cosang, -1, 1))))
    return np.trapezoid(2 * ang * R * Sig, R)


def fit_components(pos):
    """Fit masses and scale radii of the two gas clouds and two stellar concentrations to the
    eight published aperture masses."""
    keys = ['main_bcg', 'main_plasma', 'sub_bcg', 'sub_plasma']
    dist = lambda a, b: float(np.linalg.norm(pos[a] - pos[b]))
    RT = dict(gas_main=2000., gas_sub=400., st_main=1500., st_sub=600.)

    def model_gas(q):
        Mm, rcm, Ms, rcs = np.exp(q)
        um = Mm / mass3d(rho_beta, (rcm,), RT['gas_main']); us = Ms / mass3d(rho_beta, (rcs,), RT['gas_sub'])
        return {k: um * aperture(rho_beta, (rcm,), RT['gas_main'], dist(k, 'main_plasma')) +
                   us * aperture(rho_beta, (rcs,), RT['gas_sub'], dist(k, 'sub_plasma')) for k in keys}

    def model_st(q):
        Mm, rsm, Ms, rss = np.exp(q)
        um = Mm / mass3d(rho_nfw, (rsm,), RT['st_main']); us = Ms / mass3d(rho_nfw, (rss,), RT['st_sub'])
        return {k: um * aperture(rho_nfw, (rsm,), RT['st_main'], dist(k, 'main_bcg')) +
                   us * aperture(rho_nfw, (rss,), RT['st_sub'], dist(k, 'sub_bcg')) for k in keys}

    rg = least_squares(lambda q: [np.log(model_gas(q)[k] / (OBS['gas'][k] * 1e12)) for k in keys],
                       np.log([1.5e14, 250., 2e13, 60.]), bounds=(np.log([1e12, 5, 1e11, 30]), np.log([1e16, 2000, 1e15, 400])))
    # the bullet's cool core is tens of kpc across in the Chandra image (Markevitch et al. 2002);
    # left free, the four apertures push its core to the 5 kpc floor, which plants an artificial
    # lensing bump on the bullet. A 30 kpc floor still fits all four gas apertures within errors.
    rs = least_squares(lambda q: [np.log(model_st(q)[k] / (OBS['stars'][k] * 1e12)) for k in keys],
                       np.log([3e12, 300., 1e12, 60.]), bounds=(np.log([1e10, 5, 1e10, 5]), np.log([1e15, 3000, 1e15, 3000])))
    mg, ms = model_gas(rg.x), model_st(rs.x)
    comps = dict(gas_main=dict(kind='beta', M=float(np.exp(rg.x[0])), scale=float(np.exp(rg.x[1])), rt=RT['gas_main'], centre='main_plasma'),
                 gas_sub=dict(kind='beta', M=float(np.exp(rg.x[2])), scale=float(np.exp(rg.x[3])), rt=RT['gas_sub'], centre='sub_plasma'),
                 st_main=dict(kind='nfw', M=float(np.exp(rs.x[0])), scale=float(np.exp(rs.x[1])), rt=RT['st_main'], centre='main_bcg'),
                 st_sub=dict(kind='nfw', M=float(np.exp(rs.x[2])), scale=float(np.exp(rs.x[3])), rt=RT['st_sub'], centre='sub_bcg'))
    check = {k: dict(gas_model=mg[k] / 1e12, gas_obs=OBS['gas'][k], stars_model=ms[k] / 1e12, stars_obs=OBS['stars'][k]) for k in keys}
    return comps, check


# ---------------------------------------------------------------- 3D field solution
class Conv:
    """Isolated convolution with a radial kernel (zero-padded FFT, float32)."""
    def __init__(self, n, dx, fn, origin):
        m = 2 * n
        k = ((np.arange(m) - (np.arange(m) >= n) * m) * dx).astype(np.float32)
        r2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2
        r = np.sqrt(r2, dtype=np.float32); del r2
        with np.errstate(divide='ignore'):
            ker = np.where(r > 0, fn(np.maximum(r, 1e-6)), origin).astype(np.float32)
        del r
        self.n = n; self.K = np.fft.rfftn(ker); del ker

    def __call__(self, f):
        n = self.n; m = 2 * n
        pad = np.zeros((m, m, m), np.float32); pad[:n, :n, :n] = f
        return np.fft.irfftn(np.fft.rfftn(pad) * self.K, s=(m, m, m), axes=(0, 1, 2))[:n, :n, :n].astype(np.float32)


def build_density(comps, pos, x, y, z, which):
    rho = np.zeros((len(x), len(y), len(z)), np.float32)
    for name, c in comps.items():
        if which not in name: continue
        cx, cy = pos[c['centre']]
        r = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
        prof = rho_beta if c['kind'] == 'beta' else rho_nfw
        norm = c['M'] / mass3d(prof, (c['scale'],), c['rt'])       # analytic normalisation
        rho += (norm * prof(r, c['scale'], c['rt'])).astype(np.float32)
    return rho


def cube_average(fn, samples=2000000, seed=1):
    """Average of a radial kernel over a unit cube centred on the origin."""
    p = np.random.default_rng(seed).uniform(-.5, .5, (samples, 3))
    return float(np.mean(fn(np.linalg.norm(p, axis=1))))


def own_sigma_profile(gas, stars, consts, extra_gas=None, n=500, iters=300):
    """Stars' velocity dispersion from the isotropic Jeans equation in the gravity our law
    makes from the cluster's own gas and stars (spherical, iterated to a fixed point).
    Nothing observed about the lensing or the galaxy speeds goes in."""
    import law as Lw
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']
    r = np.geomspace(1.0, 3000.0, n); dr = np.gradient(r)
    def shells(c):
        prof = rho_beta if c['kind'] == 'beta' else rho_nfw
        norm = c['M'] / mass3d(prof, (c['scale'],), c['rt'])
        return norm * prof(r, c['scale'], c['rt']) * 4 * np.pi * r ** 2 * dr
    dmg = shells(gas) + (shells(extra_gas) if extra_gas else 0.0)
    dms = shells(stars)
    Mb = np.cumsum(dmg + dms); gN = G * Mb / r ** 2
    W = Lw.shell_weights(r, r)
    rho_s = dms / (4 * np.pi * r ** 2 * dr)
    def jeans(g):
        tail = np.cumsum((rho_s * g * dr)[::-1])[::-1]
        return np.where(rho_s > 1e-30, tail / np.maximum(rho_s, 1e-300), 0.0)
    sig2 = jeans(gN + np.sqrt(a * gN))
    for _ in range(iters):
        k = 3 * sig2 / u ** 2
        S = G * (W @ (k * dms)) / r ** 2
        g = gN + np.exp(-gN / (lam * a)) * np.sqrt(a * (gN + S))
        new = jeans(g)
        if np.max(np.abs(new - sig2) * dms) < 1e-7 * np.max(new * dms): sig2 = new; break
        sig2 = 0.5 * sig2 + 0.5 * new
    return r, np.sqrt(np.maximum(sig2, 0.0))


def kappa_map(comps, pos, consts, sig_main, sig_sub, n=192, dx=15.0, centre=(360., 50.), heat=True, direction='flow', want_force=False):
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']
    gd = lam * a
    x = (np.arange(n) - n / 2 + 0.5) * dx + centre[0]
    y = (np.arange(n) - n / 2 + 0.5) * dx + centre[1]
    z = (np.arange(n) - n / 2 + 0.5) * dx
    rho_g = build_density(comps, pos, x, y, z, 'gas')
    rho_sm = build_density({k: v for k, v in comps.items() if k == 'st_main'}, pos, x, y, z, 'st')
    rho_ss = build_density({k: v for k, v in comps.items() if k == 'st_sub'}, pos, x, y, z, 'st')
    rho_b = rho_g + rho_sm + rho_ss
    def kfield(sig, cen):
        if np.isscalar(sig):
            return 3 * sig ** 2 / u ** 2
        rr, ss = sig                                     # a radial dispersion profile about cen
        r3 = np.sqrt((x[:, None, None] - cen[0]) ** 2 + (y[None, :, None] - cen[1]) ** 2 + z[None, None, :] ** 2)
        return (3 * np.interp(r3, rr, ss) ** 2 / u ** 2).astype(np.float32)
    k_main = kfield(sig_main, pos[comps['st_main']['centre']])
    k_sub = kfield(sig_sub, pos[comps['st_sub']['centre']])
    krho = (k_main * rho_sm + k_sub * rho_ss) if heat else np.zeros_like(rho_b)
    dV = dx ** 3
    # cell-averaged 1/r and 1/r^2 at the origin cell (cube of side dx)
    inv_r = Conv(n, dx, lambda r: 1.0 / r, cube_average(lambda r: 1 / r) / dx)
    phiN = -G * inv_r(rho_b * dV)
    gN = -np.array(np.gradient(phiN, dx)); del phiN
    # heat-weighted flow of the free-streaming matter: the incoherent energy fluxes still add
    # as vectors, so the net companion flow is along g_N + g_hot
    g_hot = -np.array(np.gradient(-G * inv_r(krho * dV), dx)) if (heat and direction in ('flow', 'mix')) else 0.0
    inv_r2 = Conv(n, dx, lambda r: 1.0 / r ** 2, cube_average(lambda r: 1 / r ** 2) / dx ** 2)
    S = G * inv_r2(krho * dV) if heat else 0.0
    del inv_r, inv_r2
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    flow = gN + g_hot
    extra = np.exp(-mag / gd) * np.sqrt(a * (mag + S))
    if direction == 'mix':
        # pull along the net flow, weighted by how much of the flow survives: the
        # directions of the ordinary and heat-fed flows mix in proportion to their sizes
        hmag = np.sqrt(np.sum(g_hot ** 2, axis=0)) if heat else 0.0
        h = gN + extra * flow / (mag + hmag + 1e-30)
    else:
        fmag = np.sqrt(np.sum(flow ** 2, axis=0)) + 1e-30
        h = gN + extra * flow / fmag
    divh = sum(np.gradient(h[i], dx, axis=i) for i in range(3))
    rho_eff = -divh / (4 * np.pi * G)
    Sig_eff = rho_eff.sum(axis=2) * dx
    Sig_b = rho_b.sum(axis=2) * dx; Sig_g = rho_g.sum(axis=2) * dx
    if want_force:
        # conservative pull g = -grad Phi with  lap Phi = -div h ;  Phi = (1/4pi) int div h / |x - x'|
        inv_r = Conv(n, dx, lambda r: 1.0 / r, cube_average(lambda r: 1 / r) / dx)
        phi = inv_r(divh * dV) / (4 * np.pi); del inv_r
        g = -np.array(np.gradient(phi, dx))
        inner = (np.abs(x)[:, None, None] < 1e9)                     # whole box
        Fnet = np.array([np.sum(rho_b * g[i]) * dV for i in range(3)])
        Fabs = float(np.sum(rho_b * np.sqrt(np.sum(g ** 2, axis=0))) * dV)
        Mb = float(rho_b.sum() * dV)
        budget = float(0.5 * a * (rho_b.sum() + krho.sum()) * dV)      # companion momentum flux P/u = (a/2) sum (1+k) m
        force = dict(net_force_over_mass=(Fnet / Mb).tolist(), mean_pull=Fabs / Mb, momentum_budget_over_mass=budget / Mb,
                     net_over_mean_pull=float(np.linalg.norm(Fnet) / Fabs), net_over_budget=float(np.linalg.norm(Fnet) / budget))
        return x, y, Sig_eff, Sig_b, Sig_g, force
    return x, y, Sig_eff, Sig_b, Sig_g


def aperture_mean(x, y, Sig, p, rap=100.0):
    X, Y = np.meshgrid(x, y, indexing='ij')
    m = (X - p[0]) ** 2 + (Y - p[1]) ** 2 <= rap ** 2
    return float(Sig[m].mean())


def clowe_decomposition(x, y, kap, pos, rmax=1200.0):
    """Fit two circularly symmetric (pseudo-isothermal) profiles centred on the BCGs, as Clowe
    et al. did, then report each BCG's own kappa and the residual at the plasma peaks."""
    X, Y = np.meshgrid(x, y, indexing='ij')
    mid = 0.5 * (pos['main_bcg'] + pos['sub_bcg'])
    sel = (X - mid[0]) ** 2 + (Y - mid[1]) ** 2 <= rmax ** 2
    Rm = np.hypot(X - pos['main_bcg'][0], Y - pos['main_bcg'][1])
    Rs = np.hypot(X - pos['sub_bcg'][0], Y - pos['sub_bcg'][1])
    prof = lambda R, k0, rc: k0 / np.sqrt(1 + (R / rc) ** 2)
    def res(q):
        k0m, rcm, k0s, rcs, c0 = q[0], np.exp(q[1]), q[2], np.exp(q[3]), q[4]
        return (prof(Rm, k0m, rcm) + prof(Rs, k0s, rcs) + c0 - kap)[sel]
    best = least_squares(res, [1.0, np.log(100.), 0.5, np.log(100.), 0.0])
    k0m, rcm, k0s, rcs, c0 = best.x[0], np.exp(best.x[1]), best.x[2], np.exp(best.x[3]), best.x[4]
    main_only = kap - prof(Rs, k0s, rcs) - c0
    sub_only = kap - prof(Rm, k0m, rcm) - c0
    resid = kap - prof(Rm, k0m, rcm) - prof(Rs, k0s, rcs) - c0
    return dict(main_bcg=aperture_mean(x, y, main_only, pos['main_bcg']),
                sub_bcg=aperture_mean(x, y, sub_only, pos['sub_bcg']),
                main_plasma=aperture_mean(x, y, resid, pos['main_plasma']),
                sub_plasma=aperture_mean(x, y, resid, pos['sub_plasma']),
                fit=dict(k0_main=float(k0m), rc_main=float(rcm), k0_sub=float(k0s), rc_sub=float(rcs), sheet=float(c0)))


def peaks(x, y, kap, pos, smooth_kpc=40.0):
    dx = x[1] - x[0]
    ks = gaussian_filter(kap, smooth_kpc / dx)
    mx = (ks == maximum_filter(ks, size=7)) & (ks > 0.25 * ks.max())
    out = []
    for i, j in zip(*np.nonzero(mx)):
        p = np.array([x[i], y[j]])
        out.append(dict(x=float(p[0]), y=float(p[1]), kappa=float(ks[i, j]),
                        **{f'dist_{k}': float(np.linalg.norm(p - v)) for k, v in pos.items()}))
    return sorted(out, key=lambda d: -d['kappa'])


def along_axis_fraction(peak, stars, gas):
    """Where the peak sits between the stars (0) and the gas (1), projected on that axis."""
    d = gas - stars
    return float(np.dot(np.array([peak['x'], peak['y']]) - stars, d) / np.dot(d, d))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--n', type=int, default=192); ap.add_argument('--dx', type=float, default=15.0)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((Path(__file__).resolve().parent.parent / 'run-v3/results.json').read_text())['constants']
    pos = positions(); scrit = sigma_crit()
    comps, check = fit_components(pos)
    print('positions (kpc, x = west, y = north):', {k: np.round(v, 1).tolist() for k, v in pos.items()})
    print(f'Sigma_crit (z_s = 1) = {scrit:.3e} Msun/kpc^2')
    for k, v in comps.items(): print(f"  {k:9s} M = {v['M']:.3e}  scale = {v['scale']:.1f} kpc")
    for k, v in check.items(): print(f"  aperture {k:12s} gas {v['gas_model']:.2f} (obs {v['gas_obs']})  stars {v['stars_model']:.2f} (obs {v['stars_obs']})")

    runs = []
    prof_main = own_sigma_profile(comps['gas_main'], comps['st_main'], consts)
    prof_sub = own_sigma_profile(comps['gas_sub'], comps['st_sub'], consts)
    # before the collision the subcluster held the gas that ram pressure has since stripped;
    # merger reconstructions put it at 1:6 to 1:10 of the main cluster (we restore 1:8)
    Mmain_b = comps['gas_main']['M'] + comps['st_main']['M']
    lost = dict(kind='beta', M=Mmain_b / 8 - comps['gas_sub']['M'] - comps['st_sub']['M'], scale=150.0, rt=1000.0)
    prof_sub_pre = own_sigma_profile(comps['gas_sub'], comps['st_sub'], consts, extra_gas=lost)
    speeds = {}
    for nm, (rr, ss) in (('main', prof_main), ('sub_today', prof_sub), ('sub_pre_collision', prof_sub_pre)):
        speeds[nm] = {f'{int(r0)}kpc': float(np.interp(r0, rr, ss)) for r0 in (20, 50, 100, 200, 400)}
        print(f"  own-gravity star speeds, {nm}: " + ', '.join(f"{k} {v:.0f}" for k, v in speeds[nm].items()) + ' km/s')
    cases = [('A. prediction: star speeds from our law and the visible matter (subcluster as it was before the collision)', prof_main, prof_sub_pre, 'mix'),
             ('B. same, subcluster as it is today', prof_main, prof_sub, 'mix'),
             ('C. measured galaxy speeds, constant (main 1249, sub 212)', 1249., 212., 'mix'),
             ('D. speeds implied by the lensing mass (circular; reference only: main 660, sub 500)', 660., 500., 'mix'),
             ('E. no heat at all (the cold, MOND-like limit of our law)', 1249., 212., 'mix'),
             ('F. rounds 1-2 direction rule (pull along ordinary gravity), speeds as A', prof_main, prof_sub_pre, 'newton')]
    for label, sm, ss, direction in cases:
        heat = 'no heat' not in label
        res = kappa_map(comps, pos, consts, sm, ss, n=args.n, dx=args.dx, heat=heat, direction=direction, want_force=label.startswith('A.'))
        x, y, Sig_eff, Sig_b, Sig_g = res[:5]
        force = res[5] if len(res) > 5 else None
        kap = Sig_eff / scrit
        total = {k: aperture_mean(x, y, kap, p) for k, p in pos.items()}
        baryon = {k: aperture_mean(x, y, Sig_b / scrit, p) for k, p in pos.items()}
        grid_masses = {k: aperture_mean(x, y, Sig_b, p) * np.pi * 100.0 ** 2 / 1e12 for k, p in pos.items()}
        gasfrac = {k: aperture_mean(x, y, Sig_g / scrit, p) / total[k] for k, p in pos.items()}
        dec = clowe_decomposition(x, y, kap, pos)
        kap_min = float(kap.min())
        pk = peaks(x, y, kap, pos)
        main_pk = min(pk, key=lambda d: d['dist_main_bcg']) if pk else None
        sub_pk = min(pk, key=lambda d: d['dist_sub_bcg']) if pk else None
        row = dict(case=label, sigma_main=sm if np.isscalar(sm) else {'r_kpc': [20, 50, 100, 200, 400], 'sigma': [float(np.interp(r0, *sm)) for r0 in (20, 50, 100, 200, 400)]},
                   sigma_sub=ss if np.isscalar(ss) else {'r_kpc': [20, 50, 100, 200, 400], 'sigma': [float(np.interp(r0, *ss)) for r0 in (20, 50, 100, 200, 400)]},
                   heat=heat, direction=direction,
                   k_main=(3 * sm ** 2 / consts['u_kms'] ** 2) if np.isscalar(sm) else 'radial profile',
                   k_sub=(3 * ss ** 2 / consts['u_kms'] ** 2) if np.isscalar(ss) else 'radial profile',
                   kappa_total_apertures=total, kappa_baryons_apertures=baryon, gas_share_of_kappa=gasfrac,
                   grid_baryon_aperture_masses_1e12=grid_masses,
                   clowe_style=dec, peaks=pk[:6], kappa_min=kap_min, momentum=force,
                   main_peak_fraction_stars_to_gas=along_axis_fraction(main_pk, pos['main_bcg'], pos['main_plasma']) if main_pk else None,
                   sub_peak_fraction_stars_to_gas=along_axis_fraction(sub_pk, pos['sub_bcg'], pos['sub_plasma']) if sub_pk else None)
        runs.append(row)
        print(f"\n{label}:  k_main {row['k_main'] if np.isscalar(row['k_main']) else 'profile'}  k_sub {row['k_sub'] if np.isscalar(row['k_sub']) else 'profile'}")
        print('   kappa (Clowe-style)  main BCG %.3f (obs 0.36+-0.06)  sub BCG %.3f (obs 0.20+-0.05)  main plasma %.3f (obs 0.05+-0.06)  sub plasma %.3f (obs 0.02+-0.06)'
              % (dec['main_bcg'], dec['sub_bcg'], dec['main_plasma'], dec['sub_plasma']))
        print('   total kappa in apertures', {k: round(v, 3) for k, v in total.items()}, ' baryons only', {k: round(v, 3) for k, v in baryon.items()}, ' min kappa %.3f' % kap_min)
        for p in pk[:4]:
            print(f"   peak kappa {p['kappa']:.3f} at ({p['x']:.0f},{p['y']:.0f}): to main BCG {p['dist_main_bcg']:.0f}, main plasma {p['dist_main_plasma']:.0f}, sub BCG {p['dist_sub_bcg']:.0f}, sub plasma {p['dist_sub_plasma']:.0f} kpc")
        print(f"   peak position between stars (0) and gas (1): main {row['main_peak_fraction_stars_to_gas']}, sub {row['sub_peak_fraction_stars_to_gas']}")
        if force: print('   momentum:', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in force.items()})
        np.savez_compressed(out / f"kappa_{len(runs)}.npz", x=x, y=y, kappa=kap.astype(np.float32), kappa_baryons=(Sig_b / scrit).astype(np.float32))

    payload = dict(experiment='Bullet Cluster under round-3 hot-companion gravity',
                   data=dict(apertures=OBS, source='Clowe et al. 2006 Table 2; Barrena et al. 2002 dispersions'),
                   sigma_crit_Msun_per_kpc2=scrit, positions_kpc={k: v.tolist() for k, v in pos.items()},
                   components=comps, aperture_check=check, constants=consts, own_gravity_star_speeds_kms=speeds, runs=runs,
                   grid=dict(n=args.n, dx_kpc=args.dx), seconds=time.monotonic() - t0)
    (out / 'bullet.json').write_text(json.dumps(payload, indent=2, default=float) + '\n')
    print(f'{time.monotonic() - t0:.0f} s')


if __name__ == '__main__':
    main()
