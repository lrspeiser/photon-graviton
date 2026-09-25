"""Hot-companion gravity: the law and its ingredients.

Mechanism (proposed here; originality unverified -- see ../README.md for prior art):

 1. All ordinary matter feeds a companion field.
 2. Ordered (cold) matter feeds it COHERENTLY. Contributions from different places keep a
    fixed relative phase and add as vectors, so they cancel exactly the way Newton's pull
    does. The coherent intensity at a point equals |g_N| there (Gauss's law: established).
 3. Randomly moving (hot) matter feeds it INCOHERENTLY, in proportion to its random kinetic
    energy measured against the companion's own speed u: weight k = 3 sigma^2 / u^2.
    Incoherent contributions add as intensities -- no cancellation:
        S_hot(x) = G * integral k(x') rho(x') / |x - x'|^2 dV'
 4. The companion pulls with its AMPLITUDE, the square root of its total intensity times a
    constant a, along the net (coherent) pull.
 5. Strong fields hold the companion ATTACHED (inactive); it is released where the
    ordinary pull is weak, like a thermally activated escape. The released fraction is
        f = exp( -|g_N| / g_d ),     g_d = lambda * a
    This switch-off is required by the Solar System, where any surviving companion pull
    would be measured, and it also improves the galaxy fits.
 6. Total pull on matter and light (no slip assumed):
        g = g_N + f * sqrt( a * ( |g_N| + S_hot ) )

Three constants: a and lambda (fitted once, on SPARC rotation curves) and u (fitted once,
on X-COP cluster masses). Everything else is predicted.

Units: kpc, km/s, Msun; accelerations in (km/s)^2/kpc unless marked _SI.
"""
from __future__ import annotations
import numpy as np

G = 4.30091727003628e-6                     # kpc (km/s)^2 / Msun
KMS2_PER_KPC = 1e6 / 3.0856775814913673e19  # (km/s)^2/kpc -> m/s^2
C_KMS = 299792.458


def w_inside(x):
    """Shell-averaged 1/d^2 for a shell of radius s = x R inside the field point (x<1),
    in units of 1/R^2: artanh(x)/x. Established geometry."""
    x = np.clip(np.asarray(x, float), 1e-12, 1 - 1e-12)
    return np.arctanh(x) / x


def w_outside(x):
    """Same for a shell outside (x>1): arcoth(x)/x -> 1/x^2 far away."""
    x = np.maximum(np.asarray(x, float), 1 + 1e-12)
    return 0.5 * np.log((x + 1) / (x - 1)) / x


HOT_GEOMETRY = 'two_way'   # round 19: how a receiver hears the hot matter's extra glow. 'two_way' is the law: every
                           # shell counts, inside and outside (w_inside, w_outside). 'one_way': a medium whose wave only
                           # travels outward, so shells outside the receiver are not heard (the hot-shell benchmark);
                           # inner shells keep w_inside. 'one_way_vector': the net outward flux of the inner shells, which
                           # counts each inner shell as if at the centre (Gauss: weight 1). 'stream': a stream that absorbs
                           # waves travelling inward at STREAM_KAPPA per kpc (straight rays; README 29.3): as kappa -> infinity
                           # it hears exactly half of every inner shell and nothing outside. Set by the regression suite from
                           # the law's 'hot_geometry' and 'stream_kappa_per_Mpc' (regression/common.py).
STREAM_KAPPA = 0.0         # per kpc, for HOT_GEOMETRY == 'stream'
_STREAM_CACHE = {}


def shell_weight_fn(x):
    """The shell weight at x = s/R (shell radius over field radius) for the ratio-only HOT_GEOMETRY rules."""
    x = np.asarray(x, float)
    if HOT_GEOMETRY == 'two_way':
        return np.where(x < 1, w_inside(x), w_outside(x))
    if HOT_GEOMETRY == 'one_way':
        return np.where(x < 1, w_inside(x), 0.0)
    if HOT_GEOMETRY == 'one_way_vector':
        return np.where(x < 1, 1.0, 0.0)
    raise ValueError(f'{HOT_GEOMETRY} depends on the radii themselves, not only their ratio')


def stream_weights(R, s, kappa, nq=600, chunk=2_000_000, absorbed=False):
    """Shell weights for a stream that absorbs inward-travelling waves at kappa per unit length (straight rays, the
    eikonal form of round 18): the scalar weight R^2 <T/d^2> (the hot glow's intensity) and the net-flux weight
    R^2 <T (R - s mu)/d^3> (the cold glow's vector sum), averaged over each shell. T = exp(-kappa L), L the distance the
    path from the shell point to the receiver travels inward (code/absorbing_stream_v18.path_inward). The shell average
    is taken over d, the distance to the receiver ((1/2) dmu = d dd/(2 R s)), on a logarithmic grid from |R - s| to
    R + s, which resolves the shells next to the receiver: with kappa = 0 the weights are the law's artanh(x)/x and
    arcoth(x)/x and Gauss's 1 and 0, exactly (code/hot_shell_v19.py checks both limits)."""
    R = np.atleast_1d(np.asarray(R, float)); s = np.atleast_1d(np.asarray(s, float))
    ws = np.empty((R.size, s.size)); wv = np.empty((R.size, s.size))
    t = np.linspace(0.0, 1.0, nq)[None, None, :]
    w = np.full(nq, 1.0 / (nq - 1)); w[[0, -1]] *= 0.5
    step = max(1, chunk // max(1, s.size * nq))
    for i in range(0, R.size, step):
        Rr = R[i:i + step, None, None]; ss = s[None, :, None]
        lo = np.maximum(np.abs(Rr - ss), 1e-9 * Rr); hi = Rr + ss
        d = lo * (hi / lo) ** t; du = np.log(hi / lo)[..., 0]
        mu = np.clip((Rr ** 2 + ss ** 2 - d ** 2) / (2 * Rr * ss), -1, 1)
        xs, zs = ss * np.sqrt(1 - mu ** 2), ss * mu                    # the shell point; the receiver at (0, R)
        dx, dz = -xs, Rr - zs; dd = np.sqrt(dx ** 2 + dz ** 2) + 1e-300
        b0 = (dx * xs + dz * zs) / dd                                  # n . x' (< 0: the path starts inward)
        px, pz = xs - b0 * dx / dd, zs - b0 * dz / dd; p = np.sqrt(px ** 2 + pz ** 2)
        L = np.where(b0 >= 0, 0.0, np.where(-b0 >= dd, ss - Rr, ss - p))
        T = np.exp(-kappa * np.maximum(L, 0.0))
        if absorbed:                                                   # the removed part, 1 - T (for the table)
            T = -np.expm1(-kappa * np.maximum(L, 0.0))
        ws[i:i + step] = (T * w).sum(-1) * du * (Rr / (2 * ss))[..., 0]
        wv[i:i + step] = (T * (Rr ** 2 - ss ** 2 + d ** 2) / (4 * ss * d) * w).sum(-1) * du
    return ws, wv


_DW_TABLE = None


def _absorbed_table():
    """The part of the two-way shell weights the stream removes, dW = W_two-way - W_stream (scalar and net flux), as a
    function of t = ln(s/R) and ln(kappa R): the weights are scale-free, and dW is bounded and continuous across s = R
    (the removed part vanishes where the path to the receiver starts outward), so it interpolates cleanly where the
    weights themselves diverge. Built once from stream_weights (about 10 s) and kept in regression/cache."""
    global _DW_TABLE
    if _DW_TABLE is not None:
        return _DW_TABLE
    from pathlib import Path
    f = Path(__file__).resolve().parent.parent / 'regression' / 'cache' / 'stream_dw_table_v19.npz'
    if f.exists():
        z = np.load(f); _DW_TABLE = (z['t'], z['ll'], z['dws'], z['dwv']); return _DW_TABLE
    v = np.linspace(-np.arcsinh(600.0), np.arcsinh(600.0), 1201); t = 0.02 * np.sinh(v)   # dense at s = R
    ll = np.linspace(np.log(1e-8), np.log(1e6), 141)
    x = np.exp(t); dws = np.empty((t.size, ll.size)); dwv = np.empty((t.size, ll.size))
    for j, lam in enumerate(np.exp(ll)):
        ws, wv = stream_weights(np.array([1.0]), x, lam, nq=800, absorbed=True)
        dws[:, j], dwv[:, j] = ws[0], wv[0]
    f.parent.mkdir(parents=True, exist_ok=True); np.savez(f, t=t, ll=ll, dws=dws, dwv=dwv)
    _DW_TABLE = (t, ll, dws, dwv); return _DW_TABLE


def _stream_interp(R, s):
    """Scalar and net-flux stream weights for field radii R and shells s at kappa = STREAM_KAPPA, from the table."""
    t, ll, dws, dwv = _absorbed_table()
    R = np.atleast_1d(np.asarray(R, float))[:, None]; s = np.atleast_1d(np.asarray(s, float))[None, :]
    x = s / R
    tq = np.clip(np.log(x), t[0], t[-1]); lq = np.clip(np.log(np.maximum(STREAM_KAPPA * R, 1e-300)) * np.ones_like(x), ll[0], ll[-1])
    i = np.clip(np.searchsorted(t, tq) - 1, 0, t.size - 2); j = np.clip(np.searchsorted(ll, lq) - 1, 0, ll.size - 2)
    a = (tq - t[i]) / (t[i + 1] - t[i]); b = (lq - ll[j]) / (ll[j + 1] - ll[j])
    bil = lambda T: (1 - a) * (1 - b) * T[i, j] + a * (1 - b) * T[i + 1, j] + (1 - a) * b * T[i, j + 1] + a * b * T[i + 1, j + 1]
    small = STREAM_KAPPA * R * np.ones_like(x) < 1e-8                 # below the table: two-way to 1e-8
    Ds, Dv = np.where(small, 0.0, bil(dws)), np.where(small, 0.0, bil(dwv))
    w2 = np.where(x < 1, w_inside(x), w_outside(x))
    return w2 - Ds, np.where(x < 1, 1.0, 0.0) - Dv


def _stream_cached(R, s):
    R = np.atleast_1d(np.asarray(R, float)); s = np.atleast_1d(np.asarray(s, float))
    key = (STREAM_KAPPA, R.tobytes(), s.tobytes())
    if key not in _STREAM_CACHE:
        if len(_STREAM_CACHE) > 4000:
            _STREAM_CACHE.clear()
        _STREAM_CACHE[key] = _stream_interp(R, s)
    return _STREAM_CACHE[key]


def stream_cold_weights(R, s):
    """The stream's net-flux weights (the cold glow's vector sum, heard through the absorbing stream)."""
    return _stream_cached(R, s)[1]


def shell_weights(R, s):
    """Matrix W[i,j]: 1/d^2 averaged over shell j (radius s_j) at radius R_i, times R_i^2 (for the current
    HOT_GEOMETRY; the law's two-way sum by default)."""
    if HOT_GEOMETRY == 'stream':
        return _stream_cached(R, s)[0]
    R = np.atleast_1d(R)[:, None]; s = np.asarray(s)[None, :]
    return shell_weight_fn(s / R)


def scalar_sum(R, s, dm, k):
    """S_hot at radii R for spherical shells of mass dm at radii s with heat weight k."""
    R = np.atleast_1d(R)
    return G * (shell_weights(R, s) @ (np.asarray(k) * np.asarray(dm))) / R ** 2


HEAT_P = 2.0   # round 14: the heat weight's exponent, k = 3 (sigma / u)^p; the law has 2. Set by the regression
               # suite from the law's 'heat_exponent' (regression/common.py) before any test runs.


def heat_weight(sigma_kms, u_kms):
    """k = 3 sigma^2 / u^2: random kinetic energy against the companion speed. With HEAT_P = p != 2,
    k = 3 (sigma / u)^p (round 14)."""
    if HEAT_P == 2.0:
        return 3.0 * np.asarray(sigma_kms) ** 2 / u_kms ** 2
    return 3.0 * (np.abs(np.asarray(sigma_kms)) / u_kms) ** HEAT_P


def k_from_sig2(sig2, u_kms, pref=3.0):
    """The heat weight from a mean-square speed: pref sig2 / u^2 (pref = 3 for a 1D dispersion sig2; 1 for a
    3D mean square; 3 - 2 beta for an anisotropic radial one). With HEAT_P = p != 2, the same effective
    1D dispersion sigma_eff^2 = pref sig2 / 3 enters as k = 3 (sigma_eff / u)^p (round 14). At p = 2 the
    arithmetic is the scripts' own, so results are unchanged to the last digit."""
    if HEAT_P == 2.0:
        return pref * sig2 / u_kms ** 2
    return 3.0 * (np.maximum(pref * np.asarray(sig2) / 3.0, 0.0) / u_kms ** 2) ** (HEAT_P / 2.0)


def released(gN, a, lam):
    """Fraction of the companion released from attachment: exp(-|g_N| / (lambda a))."""
    return np.exp(-np.abs(np.asarray(gN)) / (lam * a))


def total(gN, S_hot, a, lam=np.inf):
    """g = g_N + f sqrt(a (|g_N| + S_hot)); lam=inf switches the release factor off."""
    gN = np.asarray(gN)
    f = 1.0 if not np.isfinite(lam) else released(gN, a, lam)
    return gN + f * np.sqrt(a * (np.abs(gN) + S_hot))


def total_heard(gN, heard, S_hot, a, lam=np.inf):
    """Round 19: g = g_N + f sqrt(a (heard + S_hot)), where the cold glow's vector sum reaches the receiver through the
    medium as `heard` (|g_N| in the law). The release f still follows the local |g_N|."""
    gN = np.asarray(gN)
    f = 1.0 if not np.isfinite(lam) else released(gN, a, lam)
    return gN + f * np.sqrt(a * np.maximum(np.asarray(heard) + S_hot, 0.0))


def point_mass(a, u, lam=np.inf):
    """Isolated point mass with internal dispersion sigma: S_hot = k g_N exactly."""
    def f(r, M, sigma_kms):
        gN = G * M / np.asarray(r) ** 2
        return total(gN, heat_weight(sigma_kms, u) * gN, a, lam)
    return f


# ---- context laws (NOT ours; used only as labelled benchmarks) ----------------------
def mond_simple(gN, a):
    """Milgrom's MOND with the 'simple' interpolation function. Context only."""
    return gN / 2 + np.sqrt(gN ** 2 / 4 + gN * a)


def nfw_velocity2(r, v200, c, h=0.7):
    """NFW halo circular velocity^2 (Navarro, Frenk & White 1996). Dark-matter context."""
    r200 = v200 / (10 * h * 0.1)      # kpc for v200 in km/s with H0 = 100h km/s/Mpc
    x = np.asarray(r) / r200
    m = lambda y: np.log(1 + y) - y / (1 + y)
    return v200 ** 2 * m(c * x) / (x * m(c))
