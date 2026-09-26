"""Round 25: the trapped companion -- consequence A of the self-energy repair (results README 36.8) -- as a switch that
reaches every test of the regression suite. Registered in round25-trapped-suite.md before it was run.

With the repair, a hot-glow packet of free energy eps0 has the energy n eps0 (aligned case), where

    n = 1 - 2 (sqrt(Q + T) - sqrt(T)) / sqrt(a),        n < 0  exactly where  Q > a/4 + sqrt(a T),

Q = |g_N| and T the hot brightness where the packet is. In a static field glow born where n < 0 can never reach a place
where n > 0, so it is heard only where n < 0; a place where n < 0 hears all of it, as in the law.

MODE (set by install(), from the law's 'companion_trapping'):
    None            the law (every hook removed; the suite reproduces its baseline to the last digit)
    'own_heat'      a cluster galaxy's own glow, inside it, carries the heat of its internal motions (SDSS, 10^11.05)
    'cluster_heat'  it carries the cluster's heat (the law's formula taken literally inside the galaxy)

Where each test meets the rule (the registration's rules 2-5). Each hook site reads MODE from this module when it runs:
    SPARC (run.galaxy_g -> galaxy_g), SLACS (lenses_t35.analyse -> slacs_S), the Milky Way (milky_way_v7.Evaluator.run ->
        mw_trap, with mw_model.heat_fields' shell_factor): resolved systems; n from the law's own fields at each place;
        emitters where n < 0 heard only where n < 0
    KiDS and Mistele (kids_heat_v12.kids / mistele -> lens_factor): each lens a Hernquist or exponential profile at its
        SDSS size; its heat times the share of its stars outside its zone
    X-COP (run_v3.cluster_M3 -> cluster_M3), the pre-collision clusters (bullet_v4.own_sigma_multi and
        bullet_main_v5.own_sigma -> jeans_S) and the collision maps (bullet_v4.kappa_map_v4 -> map_trap): each emitter is a
        typical cluster galaxy; its escaping share depends on the cluster's brightness there, solved together with it
    the dwarfs (regression/t_dwarfs.py -> mw_shell_escape): the Galaxy's glow from its shells' escaping shares
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
K_SI = 1e6 / 3.0856775814913673e19
MODES = (None, 'own_heat', 'cluster_heat')
MODE = None
CLUSTER_GALAXY_LOGM = 11.05          # the SDSS bin 10^11.0-11.1: M* of a typical (stellar-mass-weighted) cluster galaxy
LENS_STAR_SHARE = 0.85               # the lens models' stars / baryons (lensing_census_v7.pull_profile's fstar)
STATS = {}                           # what the suite's last calls found (escaping shares, zones); written by run_suite.py
_CACHE = {}


# ------------------------------------------------------------------------------------------------ the rule
def n_factor(Q, T, a):
    """A packet's energy over its free energy (aligned case)."""
    Q, T = np.asarray(Q, float), np.maximum(np.asarray(T, float), 0.0)
    return 1 - 2 * (np.sqrt(Q + T) - np.sqrt(T)) / np.sqrt(a)


def held(Q, T, a):
    """True where n < 0: Q > a/4 + sqrt(a T)."""
    return np.asarray(Q) > a / 4 + np.sqrt(a * np.maximum(np.asarray(T, float), 0.0))


# ------------------------------------------------------------------------------------------------ unresolved galaxies
def sizes():
    if 'sizes' not in _CACHE:
        _CACHE['sizes'] = json.loads((RESULTS / 'data/lens_sizes_sdss_v25.json').read_text())
    return _CACHE['sizes']


def lens_size(sample, logM):
    """('hernquist', circularised half-light radius) for early types, ('disk', scale length) for late types, in kpc
    (the surveys' units), interpolated in log M* between the SDSS bins' centres."""
    rows = sizes()['samples'][sample]
    early = sample in ('red', 'bulge')
    key = 'median_Re_deV_circ_kpc' if early else 'median_Rd_exp_kpc'
    c = np.array([r['logM_lo'] + 0.05 for r in rows]); v = np.array([r[key] for r in rows])
    return ('hernquist' if early else 'disk'), float(np.interp(logM, c, v))


def heat_table():
    if 'heat' not in _CACHE:
        _CACHE['heat'] = json.loads((RESULTS / 'data/lens_heat_sdss_v12.json').read_text())
    return _CACHE['heat']


def k_own_cluster_galaxy(u):
    """The heat of a typical cluster galaxy's internal motions (SDSS red galaxies at 10^11.05; round 12's measure)."""
    import lens_heat_sdss_v12 as LH
    return float(LH.k_eff(heat_table(), 'red', u, logM=CLUSTER_GALAXY_LOGM))


class Galaxy:
    """A galaxy the models do not resolve: stellar mass M (Msun), a Hernquist profile with half-light radius R or an
    exponential disk with scale length R (kpc), its other baryons following its stars (baryons = M x baryon_ratio).
    Q(r) = G M_b(<r) / r^2; its own glow per unit heat, shape(r) = G sum dm <1/d^2> over its shells (the law's two-way
    weights), so T(r) = T_outside + k shape(r)."""
    def __init__(self, M, kind, R, baryon_ratio=1.0, n=600):
        import law as L
        r = R * np.geomspace(1e-3, 3e2, n)
        m = (r / (r + R / 1.815)) ** 2 if kind == 'hernquist' else 1 - (1 + r / R) * np.exp(-r / R)
        self.dm = np.diff(np.r_[0.0, m])                     # shell i lies between r[i-1] and r[i]
        rm = np.sqrt(r * np.r_[r[0] / (r[1] / r[0]), r[:-1]])  # its middle, where it is placed for the brightness
        self.r, self.M, self.kind, self.R = r, M, kind, R
        self.Q = L.G * M * baryon_ratio * m / r ** 2
        self.shape = L.G * M * (L.shell_weight_fn(rm[None, :] / r[:, None]) @ self.dm) / r ** 2

    def trapped_share(self, T_out, k, a, chunk=2000):
        """The share of the stars where n < 0, for outside brightness T_out and own heat k (arrays of one length, or
        scalars)."""
        T_out = np.atleast_1d(np.asarray(T_out, float)); k = np.broadcast_to(np.asarray(k, float), T_out.shape)
        out = np.empty(T_out.size)
        for i in range(0, T_out.size, chunk):
            T = T_out[i:i + chunk, None] + k[i:i + chunk, None] * self.shape[None, :]
            out[i:i + chunk] = (held(self.Q[None, :], T, a) * self.dm[None, :]).sum(1)
        return out


class ClusterTable:
    """Escaping share of the typical cluster galaxy against the cluster's brightness T (code units) and the galaxy's own
    heat k, tabulated in log T and log k and read bilinearly."""
    LT = np.linspace(-3.0, 10.0, 261)
    LK = np.linspace(-3.0, 4.0, 141)

    def __init__(self, a):
        kind, R = lens_size('red', CLUSTER_GALAXY_LOGM)
        self.gal = Galaxy(10 ** CLUSTER_GALAXY_LOGM, kind, R)
        TT, KK = np.meshgrid(10 ** self.LT, 10 ** self.LK, indexing='ij')
        self.esc = 1 - self.gal.trapped_share(TT.ravel(), KK.ravel(), a).reshape(TT.shape)

    def __call__(self, T, k):
        lt = np.clip(np.log10(np.maximum(np.asarray(T, float), 1e-30)), self.LT[0], self.LT[-1])
        if np.ndim(k) == 0:
            lk = float(np.clip(np.log10(max(float(k), 1e-30)), self.LK[0], self.LK[-1]))
            j = min(int((lk - self.LK[0]) / (self.LK[1] - self.LK[0])), self.LK.size - 2)
            w = (lk - self.LK[j]) / (self.LK[1] - self.LK[0])
            col = (1 - w) * self.esc[:, j] + w * self.esc[:, j + 1]
            return np.interp(lt, self.LT, col)
        lk = np.clip(np.log10(np.maximum(np.asarray(k, float), 1e-30)), self.LK[0], self.LK[-1])
        dt, dk = self.LT[1] - self.LT[0], self.LK[1] - self.LK[0]
        i = np.minimum(((lt - self.LT[0]) / dt).astype(int), self.LT.size - 2)
        j = np.minimum(((lk - self.LK[0]) / dk).astype(int), self.LK.size - 2)
        wt = (lt - self.LT[i]) / dt; wk = (lk - self.LK[j]) / dk
        e = self.esc
        return ((1 - wt) * ((1 - wk) * e[i, j] + wk * e[i, j + 1]) + wt * ((1 - wk) * e[i + 1, j] + wk * e[i + 1, j + 1]))


def cluster_table(a):
    key = ('table', round(float(a), 9))
    if key not in _CACHE:
        _CACHE[key] = ClusterTable(a)
    return _CACHE[key]


def _own_heat(u, k_stars):
    return k_own_cluster_galaxy(u) if MODE == 'own_heat' else k_stars


def solve_escape(bright, kd, k_stars, a, u, esc0=None, tol=1e-6, iters=400):
    """The cluster galaxies' escaping share at each emitter, esc = table(T, k_own) with T = bright(kd x esc), iterated
    from the law's brightness (esc = 1) to the first self-consistent state."""
    tab = cluster_table(a); kown = _own_heat(u, k_stars)
    esc = np.ones_like(kd) if esc0 is None else esc0.copy()
    for it in range(iters):
        new = tab(bright(kd * esc), kown)
        d = float(np.max(np.abs(new - esc))) if new.size else 0.0
        esc = new
        if d < tol:
            break
    return esc, it + 1


# ------------------------------------------------------------------------------------------------ hooks: SPARC
def galaxy_g(g, a, u, lam=np.inf):
    """run.galaxy_g with the bulge's glow held where n < 0 (the law's fields at each shell and each radius)."""
    import law as L
    gN = g['gN']
    if u is None or g['sigb'] == 0:
        return L.total(gN, 0., a, lam)
    k = L.heat_weight(g['sigb'], u)
    S = L.scalar_sum(g['r'], g['sfine'], g['dmb'], k)
    s, dmb = g['sfine'], g['dmb']
    sm = np.r_[s[0] / 2, np.sqrt(s[1:] * s[:-1])]                 # shells between the nodes (never on a node)
    hold = held(np.interp(s, g['r'], np.abs(gN)), L.scalar_sum(s, sm, dmb, k), a)
    share = float(np.sum(dmb[hold]) / max(np.sum(dmb), 1e-300))
    if hold.any():
        zone = held(np.abs(gN), S, a)
        S_out = L.scalar_sum(g['r'], s, np.where(hold, 0.0, dmb), k)
        band = float(np.mean(~zone & held(np.abs(gN), S_out, a)))
        S = np.where(zone, S, S_out)
        STATS.setdefault('sparc', {})[g['name']] = dict(bulge_share_held=share, radii_in_zone=float(np.mean(zone)), band=band)
    return L.total(gN, S, a, lam)


# ------------------------------------------------------------------------------------------------ hooks: X-COP
def _xcop_grid(c):
    import law as L
    key = ('xcop', c['name'], L.HOT_GEOMETRY, c['s'].size, float(c['s'][0]), float(c['s'][-1]))
    if key not in _CACHE:
        sT = np.geomspace(c['s'][0], c['s'][-1], 160)
        _CACHE[key] = (sT, L.shell_weights(sT, c['s']))
    return _CACHE[key]


def cluster_M3(c, a, u, lam, stars_weight=1.0):
    """run_v3.cluster_M3 (the graded path: stars' heat from hydrostatic equilibrium, gas cold) with each shell's stars
    in typical cluster galaxies that keep the glow born in their zones."""
    import law as L
    G = L.G
    gN = G * c['Mb'] / c['Rk'] ** 2
    k = stars_weight * L.heat_weight(np.sqrt(c['sig2_star_hse']), u)
    sT, WT = _xcop_grid(c)
    ls, lsT = np.log(c['s']), np.log(sT)
    bright = lambda kde: np.exp(np.interp(ls, lsT, np.log(np.maximum(G * (WT @ kde) / sT ** 2, 1e-30))))
    esc, it = solve_escape(bright, k * c['dms'], k, a, u)
    STATS.setdefault('xcop', {})[c['name']] = dict(u=float(u), escaping_share=float(np.sum(k * c['dms'] * esc) / np.sum(k * c['dms'])),
                                                   iterations=it)
    S = G * (c['W'] @ (k * c['dms'] * esc)) / c['Rk'] ** 2
    return L.total(gN, S, a, lam) * c['Rk'] ** 2 / G


# ------------------------------------------------------------------------------------------------ hooks: KiDS, Mistele
def lens_factor(sample, logM, k, consts, f):
    """The share of an isolated lens's glow that leaves it: its stars outside its zone (own glow with heat k)."""
    a = consts['a_SI'] / K_SI
    kind, R = lens_size(sample, logM)
    fs = 1.0 if f is None else f['size']; fm = 1.0 if f is None else f['stars']
    key = ('lens', kind, round(R * fs, 9), round(10 ** logM * fm, 3))
    if key not in _CACHE:
        _CACHE[key] = Galaxy(LENS_STAR_SHARE * 10 ** logM * fm, kind, R * fs, 1 / LENS_STAR_SHARE)
    gal = _CACHE[key]
    share = float(gal.trapped_share(0.0, k, a)[0])
    zone = gal.r[held(gal.Q, k * gal.shape, a)]
    STATS.setdefault('lenses', {})[f'{sample} {logM:.2f}'] = dict(k=float(k), escaping_share=1 - share, size_kpc=R * fs, profile=kind,
                                                                  zone_outer_kpc=float(zone.max()) if zone.size else 0.0)
    return 1 - share


# ------------------------------------------------------------------------------------------------ hooks: SLACS
def slacs_S(r, gN, S, kdm, a):
    """lenses_t35.analyse: the lens's shells where n < 0 heard only where n < 0."""
    import law as L
    zone = held(np.abs(gN), S, a)
    if not zone.any():
        return S
    S_out = L.G * (L.shell_weights(r, r) @ np.where(zone, 0.0, kdm)) / r ** 2
    STATS['slacs_last'] = dict(zone_outer_kpc=float(r[zone].max()), held_share=float(np.sum(kdm[zone]) / np.sum(kdm)))
    return np.where(zone, S, S_out)


# ------------------------------------------------------------------------------------------------ hooks: the Milky Way
def _shell_escape(hot, grid, Q, T, a, s, nmu=64):
    """For each hot component, the share of each spherical shell's mass (its real density on the shell) where n >= 0."""
    from scipy.interpolate import RegularGridInterpolator as RGI
    fQ = RGI((grid.R, grid.z), Q, bounds_error=False, fill_value=None)
    fT = RGI((grid.R, grid.z), T, bounds_error=False, fill_value=None)
    mu, w = np.polynomial.legendre.leggauss(nmu)
    R = np.clip(s[:, None] * np.sqrt(1 - mu ** 2)[None, :], grid.R[0], grid.R[-1])
    z = np.clip(s[:, None] * mu[None, :], grid.z[0], grid.z[-1])
    pts = np.c_[R.ravel(), z.ravel()]
    hold = held(np.maximum(fQ(pts), 0.0), np.maximum(fT(pts), 0.0), a).reshape(R.shape)
    out = {}
    for h in hot:
        rho = h.c.rho(R, z) * w[None, :]
        tot = rho.sum(1)
        out[h.c.name] = np.where(tot > 0, 1 - (rho * hold).sum(1) / np.where(tot > 0, tot, 1.0), 1.0)
    return out


def mw_trap(hot, grid, consts, gR, gz, S, hR, hz, key=None):
    """milky_way_v7.Evaluator.run: the Galaxy's hot glow held where n < 0 (the law's fields on the (R, z) grid)."""
    import mw_model as M
    a = consts['a_code']
    Q = np.hypot(gR, gz)
    zone = held(Q, S, a)
    if not zone.any():
        return S, hR, hz
    s = np.geomspace(1e-3, 3000.0, 1200)                           # mw_model.heat_fields' shells
    esc = _shell_escape(hot, grid, Q, S, a, s)
    Se, hRe, hze = M.heat_fields(hot, grid, consts['u_kms'], shell_factor=lambda c, s_: np.interp(np.log(s_), np.log(s), esc[c.c.name]))
    jz = grid.nzh
    inplane = grid.R[zone[:, jz]]
    band = ~zone & held(Q, Se, a)
    info = dict(zone_outer_R_in_plane_kpc=float(inplane.max()) if inplane.size else 0.0,
                band_outer_R_in_plane_kpc=float(grid.R[band[:, jz]].max()) if band[:, jz].any() else 0.0,
                held_share={h.c.name: float(1 - np.sum(np.diff(np.r_[0.0, h.m_profile(s)]) * esc[h.c.name]) / h.m_profile(s)[-1]) for h in hot},
                S_at_sun_ratio=float(np.interp(8.2, grid.R, np.where(zone, S, Se)[:, jz] / np.maximum(S[:, jz], 1e-30))))
    STATS.setdefault('milky_way', {})[str(key)] = info
    _CACHE[('mw_escape', str(key))] = (s, esc)
    _CACHE[('mw_law', str(key))] = (consts['a_code'], consts['lam'], consts['u_kms'])
    return np.where(zone, S, Se), np.where(zone, hR, hRe), np.where(zone, hz, hze)


def mw_shell_escape(law, ctx, name, s):
    """The dwarfs' Galaxy (McMillan 2017, 'M17', as the Milky Way test solves it for this law): the escaping share of
    the component called `name` at shell radii s."""
    import milky_way_v7 as MW
    import t_milky_way as TMW
    comps, grid, F = TMW.setup(ctx)
    parts = MW.models(comps)['M17']
    ck = ('mw_escape', str(tuple(sorted(parts.items()))))
    if ck not in _CACHE or _CACHE.get(('mw_law', ck[1])) != (law['a_code'], law['lam'], law['u_kms']):
        MW.Evaluator(comps, grid, F, law).run(parts, 'ours', reach=law['reach_kpc'])
    s0, esc = _CACHE[ck]
    return np.interp(np.log(s), np.log(s0), esc[name])


# ------------------------------------------------------------------------------------------------ hooks: collisions
def jeans_S(r, W, kk, dms, consts):
    """bullet_v4.own_sigma_multi, bullet_main_v5.own_sigma: the pre-collision cluster's brightness from the glow that
    escaped its galaxies (receivers and emitters on the same shells, the law's self-weights)."""
    import law as L
    G = L.G
    key = ('jeans', id(W), W.shape)
    esc, it = solve_escape(lambda kde: G * (W @ kde) / r ** 2, kk * dms, kk, consts['a_code'], consts['u_kms'], esc0=_CACHE.get(key))
    _CACHE[key] = esc
    STATS['jeans_last'] = dict(escaping_share=float(np.sum(kk * dms * esc) / max(np.sum(kk * dms), 1e-300)), iterations=it)
    return G * (W @ (kk * dms * esc)) / r ** 2


def map_trap(krho, rho_st, n, dx, dV, consts):
    """bullet_v4.kappa_map_v4: the stars' hot emission times their galaxies' escaping share, solved together with the
    brightness on the map (every cell of the map lies outside the galaxies' zones)."""
    import bullet_v3 as B
    import law as L
    a, u = consts['a_code'], consts['u_kms']
    tab = cluster_table(a)
    conv = B.Conv(n, dx, lambda r: 1.0 / r ** 2, B.cube_average(lambda r: 1 / r ** 2) / dx ** 2)
    m = krho > 0
    kst = (krho[m] / np.maximum(rho_st[m], 1e-30)).astype(np.float64)
    kown = k_own_cluster_galaxy(u) if MODE == 'own_heat' else kst
    esc = np.ones(krho.shape, np.float32)
    for it in range(200):
        T = L.G * conv(krho * esc * dV)
        new = tab(T[m], kown)
        d = float(np.max(np.abs(new - esc[m]))) if m.any() else 0.0
        esc[m] = new
        if d < 1e-4:
            break
    share = float(np.sum(krho * esc) / max(float(np.sum(krho)), 1e-300))
    STATS.setdefault('maps', []).append(dict(n=int(n), dx=float(dx), escaping_share=share, iterations=it + 1))
    return (krho * esc).astype(np.float32)


# ------------------------------------------------------------------------------------------------ switch
def install(mode=None):
    """Switch the trapped companion on ('own_heat', 'cluster_heat') or off (None / 'none': the law, exactly). Each hook
    site (run.galaxy_g, run_v3.cluster_M3, kids_heat_v12.kids / mistele, lenses_t35.analyse, milky_way_v7.Evaluator.run,
    bullet_v4.own_sigma_multi / kappa_map_v4, bullet_main_v5.own_sigma) reads MODE from this module when it runs, so
    nothing is imported here (the suite's import order, which decides which 'run' module some codes see, is kept)."""
    global MODE
    mode = None if mode in (None, '', 'none') else mode
    if mode not in MODES:
        raise ValueError(f'companion_trapping: {mode}')
    MODE = mode
    STATS.clear()
    for k in [k for k in _CACHE if isinstance(k, tuple) and k[0] in ('jeans', 'mw_escape', 'mw_law')]:
        del _CACHE[k]
