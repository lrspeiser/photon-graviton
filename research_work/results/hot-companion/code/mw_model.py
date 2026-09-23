"""The Milky Way's ordinary matter and an exact axisymmetric solver for our law (round 7).

Units: kpc, km/s, Msun; accelerations in (km/s)^2/kpc unless marked _SI.

Pieces (all established Newtonian machinery; the law enters only in `solve`):
* Baryon components: exponential or sech^2 disks with optional central holes, a flattened
  (oblate) bulge, spherical parts (stellar halo, hot corona).
* Newtonian pull g_N on a cylindrical (R, z) grid:
    - disks: Hankel-transform solution, exact for rho = Sigma(R) zeta(z) (Binney & Tremaine 2.6);
    - bulge: the homoeoid formula for a stratified oblate spheroid (Binney & Tremaine 2.5);
    - spherical parts: enclosed mass;
    - beyond 60 kpc: an even-l multipole expansion (l <= 10) of disks and bulge.
* The law: h = g_N + exp(-|g_N|/g_d) sqrt(a (|g_N| + S)) (g_N + g_hot)/(|g_N| + |g_hot|),
  with S and g_hot from the hot, free-streaming stars (bulge and stellar halo; disks are cold,
  as in the SPARC calibration). Zero beyond the companion's reach u t.
* The pull on stars is the curl-free part of h (lap Phi = -div h). The 'as if' density
  rho_ph = -div(h - g_N) / (4 pi G) is formed on the grid and its pull is summed over rings
  (complete elliptic integrals) at each target point.
"""
from __future__ import annotations
import numpy as np
from scipy.special import j0, j1, ellipk, ellipe, eval_legendre

G = 4.30091727003628e-6                       # kpc (km/s)^2 / Msun
K_SI = 1e6 / 3.0856775814913673e19            # (km/s)^2/kpc -> m/s^2
PC2, PC3 = 1e6, 1e9                           # Msun/pc^2 -> Msun/kpc^2 ; Msun/pc^3 -> Msun/kpc^3


# ----------------------------------------------------------------------------- components
class Disk:
    """rho = Sigma(R) zeta(z), Sigma = S0 exp(-R/Rd - Rm/R); zeta exponential or sech^2."""
    def __init__(self, name, S0, Rd, zd, vert='exp', Rm=0.0, kind='stars', sigma=0.0):
        self.name, self.S0, self.Rd, self.zd, self.vert, self.Rm, self.kind, self.sigma = name, S0, Rd, zd, vert, Rm, kind, sigma

    def Sigma(self, R):
        R = np.maximum(np.asarray(R, float), 1e-9)
        hole = self.Rm / R if self.Rm > 0 else 0.0
        return self.S0 * np.exp(-R / self.Rd - hole)

    def zeta(self, z):
        z = np.abs(np.asarray(z, float))
        if self.vert == 'exp':
            return np.exp(-z / self.zd) / (2 * self.zd)
        return 1.0 / np.cosh(z / (2 * self.zd)) ** 2 / (4 * self.zd)

    def rho(self, R, z):
        return self.Sigma(R) * self.zeta(z)

    def mass(self):
        R = np.linspace(0, 60 * self.Rd + 10 * self.Rm, 200001)
        return float(np.trapezoid(2 * np.pi * R * self.Sigma(R), R))

    def hankel(self, k):
        """Sigma~(k) = int Sigma(R) J0(kR) R dR."""
        if self.Rm == 0:
            return self.S0 * self.Rd ** 2 / (1 + (k * self.Rd) ** 2) ** 1.5
        R = np.linspace(0, 40 * self.Rd, 16001)[1:]; w = np.gradient(R)
        f = self.Sigma(R) * R * w
        out = np.empty_like(k)
        for i in range(0, k.size, 500):
            out[i:i + 500] = j0(np.outer(k[i:i + 500], R)) @ f
        return out

    def Zfun(self, k, z):
        """Z(k,z) = int zeta(z') exp(-k|z - z'|) dz' and dZ/dz, for z >= 0 (k: column, z: row)."""
        k = k[:, None]; z = np.abs(z)[None, :]
        if self.vert == 'exp':
            b = 1.0 / self.zd; d = k - b
            dd = np.where(np.abs(d) < 1e-12, 1e-12, d)
            ebz, ekz = np.exp(-b * z), np.exp(-k * z)
            small = np.abs(dd * z) < 0.5                      # near k = b: cancellation-free forms
            x = np.where(small, dd * z, 0.0)
            with np.errstate(over='ignore', invalid='ignore'):
                mid = np.where(small, ebz * (-np.expm1(-x)) / dd, (ebz - ekz) / dd)            # (e^-bz - e^-kz)/(k-b)
                midz = np.where(small, ebz * (np.exp(-x) + b * np.expm1(-x) / dd), (k * ekz - b * ebz) / dd)
            Z = 0.5 * b * (ekz / (b + k) + mid + ebz / (b + k))
            Zz = 0.5 * b * (-k * ekz / (b + k) + midz - b * ebz / (b + k))
            return Z, Zz
        # sech^2: numerical in the layer, exact exponential tail above it
        zc = 30 * self.zd
        zp = np.linspace(-zc, zc, 6001); wz = np.gradient(zp) * self.zeta(zp)
        Z = np.empty((k.shape[0], z.shape[1])); Zz = np.empty_like(Z)
        above = z[0] > zc
        C = (np.cosh(np.outer(k[:, 0], zp)) * wz[None, :]).sum(1)[:, None]
        Z[:, above] = np.exp(-k * z[:, above]) * C
        Zz[:, above] = -k * Z[:, above]
        zin = z[0, ~above]
        for j0_, zz in zip(np.where(~above)[0], zin):
            e = np.exp(-np.outer(k[:, 0], np.abs(zz - zp)))
            Z[:, j0_] = e @ wz
            Zz[:, j0_] = -(k[:, 0]) * (e @ (wz * np.sign(zz - zp)))
        return Z, Zz


class Spheroid:
    """Oblate stratified bulge: rho(m), m^2 = R^2 + z^2/q^2 (McMillan 2017 form)."""
    def __init__(self, name, rho0, r0, rcut, alpha, q, sigma=0.0):
        self.name, self.rho0, self.r0, self.rcut, self.alpha, self.q, self.sigma = name, rho0, r0, rcut, alpha, q, sigma
        self.kind = 'stars'

    def rho_m(self, m):
        m = np.asarray(m, float)
        return self.rho0 / (1 + m / self.r0) ** self.alpha * np.exp(-(m / self.rcut) ** 2)

    def rho(self, R, z):
        return self.rho_m(np.sqrt(np.asarray(R) ** 2 + (np.asarray(z) / self.q) ** 2))

    def mass(self):
        m = np.geomspace(1e-5, 20 * self.rcut, 20001)
        return float(np.trapezoid(4 * np.pi * self.q * m ** 2 * self.rho_m(m), m))

    def m_profile(self, s):
        """mass inside the homoeoid m < s (used as its spherical stand-in for S and g_hot)"""
        m = np.geomspace(1e-5, s.max(), 20001)
        cum = np.concatenate([[0], np.cumsum(0.5 * np.diff(m) * (4 * np.pi * self.q * (m ** 2 * self.rho_m(m))[1:] + 4 * np.pi * self.q * (m ** 2 * self.rho_m(m))[:-1]))])
        return np.interp(s, m, cum)

    def field(self, R, z):
        """homoeoid theorem: g_R = -2 pi G q R int rho(m^2(u)) du / ((1+u)^2 sqrt(q^2+u)), etc."""
        q = self.q
        u = np.geomspace(1e-7, 1e9, 1500); w = np.gradient(u)
        R = np.asarray(R, float); z = np.asarray(z, float)
        gR = np.empty(R.shape); gz = np.empty(R.shape)
        Rf, zf = R.ravel(), z.ravel(); gRf, gzf = gR.ravel(), gz.ravel()
        for i in range(0, Rf.size, 2000):
            Ri, zi = Rf[i:i + 2000, None], zf[i:i + 2000, None]
            m = np.sqrt(Ri ** 2 / (1 + u) + zi ** 2 / (q ** 2 + u))
            r = self.rho_m(m) * w
            gRf[i:i + 2000] = -2 * np.pi * G * q * Ri[:, 0] * (r / ((1 + u) ** 2 * np.sqrt(q ** 2 + u))).sum(1)
            gzf[i:i + 2000] = -2 * np.pi * G * q * zi[:, 0] * (r / ((1 + u) * (q ** 2 + u) ** 1.5)).sum(1)
        return gRf.reshape(R.shape), gzf.reshape(R.shape)


class Sphere:
    """Spherical component given by a density function rho(r) (stellar halo, hot corona)."""
    def __init__(self, name, rho_fn, rmax, kind='stars', sigma=0.0):
        self.name, self.rho_fn, self.rmax, self.kind, self.sigma = name, rho_fn, rmax, kind, sigma
        self.s = np.geomspace(1e-3, rmax, 4000)
        dm = 4 * np.pi * self.s ** 2 * rho_fn(self.s) * np.gradient(self.s)
        self.dm = dm; self.M = np.cumsum(dm)

    def rho(self, R, z):
        r = np.sqrt(np.asarray(R) ** 2 + np.asarray(z) ** 2)
        return np.where(r < self.rmax, self.rho_fn(np.maximum(r, 1e-3)), 0.0)

    def mass(self):
        return float(self.M[-1])

    def m_profile(self, s):
        return np.interp(s, self.s, self.M, left=0.0)

    def field(self, R, z):
        r = np.sqrt(R ** 2 + z ** 2) + 1e-12
        g = -G * np.interp(r, self.s, self.M, left=0.0) / r ** 2
        return g * R / r, g * z / r


# ----------------------------------------------------------------------------- grid
class Grid:
    """Cylindrical grid, sinh-spaced in R and z (fine near the axis and the plane), z symmetric."""
    def __init__(self, Rmax=3000.0, zmax=3000.0, nR=360, nz=300, AR=1.0, Bz=0.05):
        tR = np.linspace(0, np.arcsinh(Rmax / AR), nR + 1)
        self.Re = AR * np.sinh(tR)
        tz = np.linspace(0, np.arcsinh(zmax / Bz), nz + 1)
        ze = Bz * np.sinh(tz)
        self.ze = np.concatenate([-ze[::-1], ze[1:]])
        self.R = 0.5 * (self.Re[1:] + self.Re[:-1]); self.z = 0.5 * (self.ze[1:] + self.ze[:-1])
        self.dR = np.diff(self.Re); self.dz = np.diff(self.ze)
        self.RR, self.ZZ = np.meshgrid(self.R, self.z, indexing='ij')
        self.V = (np.pi * (self.Re[1:] ** 2 - self.Re[:-1] ** 2))[:, None] * self.dz[None, :]
        self.nzh = nz                          # index of the first z > 0 cell


def newton_field(comps, grid, **kw):
    """g_N (R and z components) of all components on the grid."""
    parts = component_fields(comps, grid, **kw)
    return sum(p[0] for p in parts), sum(p[1] for p in parts)


def component_fields(comps, grid, r_split=60.0, lmax=10, kmax=40.0, dk=0.005):
    """[(g_R, g_z)] of each component separately (g_N is linear in the masses, so models that
    rescale components reuse these)."""
    out = []
    for c in comps:
        out.append(_one_field(c, grid, r_split, lmax, kmax, dk))
    return out


def _one_field(c, grid, r_split, lmax, kmax, dk):
    comps = [c]
    gR = np.zeros(grid.RR.shape); gz = np.zeros(grid.RR.shape)
    r = np.sqrt(grid.RR ** 2 + grid.ZZ ** 2)
    inner = r <= r_split
    # rectangular block that contains the inner sphere (upper half only; mirror afterwards)
    iR = np.where(grid.R <= r_split * 1.0001)[0]
    iz = np.where((grid.z > 0) & (grid.z <= r_split * 1.0001))[0]
    Rb, zb = grid.R[iR], grid.z[iz]
    k = np.arange(dk, kmax + dk / 2, dk); wk = np.full(k.size, dk); wk[-1] *= 0.5
    for c in comps:
        if isinstance(c, Disk):
            St = c.hankel(k)
            big = np.where(np.abs(St) > 1e-10 * abs(St[0]))[0]
            nk = min(k.size, big.max() + 2)          # smooth (holed) gas disks need far fewer k
            kk, St, ww = k[:nk], St[:nk], wk[:nk].copy()
            Z, Zz = c.Zfun(kk, zb)
            A = j1(np.outer(Rb, kk)) * (kk * St * ww)[None, :]
            B = j0(np.outer(Rb, kk)) * (St * ww)[None, :]
            # k -> 0 end: the integrand of g_R ~ k^2 and of g_z ~ k Zz -> 0, so the missing [0, dk) is O(dk^2)
            blkR = -2 * np.pi * G * (A @ Z)
            blkz = 2 * np.pi * G * (B @ Zz)
            sub_R = np.zeros((grid.R.size, iz.size)); sub_z = np.zeros_like(sub_R)
            sub_R[iR] = blkR; sub_z[iR] = blkz
            full_R = np.zeros(grid.RR.shape); full_z = np.zeros(grid.RR.shape)
            full_R[:, iz] = sub_R; full_z[:, iz] = sub_z
            # mirror to z < 0 (g_R even, g_z odd)
            mir = grid.z.size - 1 - iz
            full_R[:, mir] = sub_R; full_z[:, mir] = -sub_z
            gR += np.where(inner, full_R, 0.0); gz += np.where(inner, full_z, 0.0)
        elif isinstance(c, Spheroid):
            m = inner & (grid.ZZ > 0)
            fR, fz = c.field(grid.RR[m], grid.ZZ[m])
            tmpR = np.zeros(grid.RR.shape); tmpz = np.zeros(grid.RR.shape)
            tmpR[m] = fR; tmpz[m] = fz
            tmpR[:, :grid.nzh] = tmpR[:, grid.nzh:][:, ::-1]; tmpz[:, :grid.nzh] = -tmpz[:, grid.nzh:][:, ::-1]
            gR += np.where(inner, tmpR, 0.0); gz += np.where(inner, tmpz, 0.0)
        else:                                  # spherical parts: exact everywhere
            fR, fz = c.field(grid.RR, grid.ZZ)
            gR += fR; gz += fz
    # outer region: multipoles of the non-spherical parts
    flat = [c for c in comps if not isinstance(c, Sphere)]
    if not flat:
        return gR, gz
    Ql = multipoles(flat, lmax)
    outer = ~inner
    rr = r[outer]; mu = grid.ZZ[outer] / rr; st = grid.RR[outer] / rr
    gr = np.zeros(rr.shape); gth = np.zeros(rr.shape)
    for l, Q in Ql.items():
        P = eval_legendre(l, mu)
        dP = l * (mu * P - eval_legendre(l - 1, mu)) / (mu ** 2 - 1 + 1e-300) if l > 0 else 0.0 * mu
        # Phi = -G Q P_l(mu) / r^(l+1):  g_r = -dPhi/dr ,  g_theta = -(1/r) dPhi/dtheta = -G Q sin(theta) P_l'(mu) / r^(l+2)
        gr += -G * (l + 1) * Q * P / rr ** (l + 2)
        gth += -G * Q * st * dP / rr ** (l + 2)
    gR[outer] += gr * st + gth * mu
    gz[outer] += gr * mu - gth * st
    return gR, gz


def multipoles(comps, lmax):
    """Q_l = int rho r^l P_l(cos theta) dV for even l (equatorially symmetric)."""
    Re = 0.05 * np.sinh(np.linspace(0, np.arcsinh(90 / 0.05), 3001))        # midpoint rule on sinh-spaced cells
    ze = 0.002 * np.sinh(np.linspace(0, np.arcsinh(40 / 0.002), 3001))
    R, z = 0.5 * (Re[1:] + Re[:-1]), 0.5 * (ze[1:] + ze[:-1]); dR, dz = np.diff(Re), np.diff(ze)
    RR, ZZ = np.meshgrid(R, z, indexing='ij')
    rho = sum(c.rho(RR, ZZ) for c in comps)
    dV = 2 * (2 * np.pi * RR * dR[:, None] * dz[None, :])          # both halves
    r = np.sqrt(RR ** 2 + ZZ ** 2); mu = ZZ / r
    return {l: float(np.sum(rho * dV * r ** l * eval_legendre(l, mu))) for l in range(0, lmax + 1, 2)}


# ----------------------------------------------------------------------------- rings
def ring_sum(tR, tz, grid, rho, mask=None, chunk=4000, potential=False):
    """Pull (and potential) at target points (tR, tz) of the density rho on the grid, summed over rings."""
    m = rho * grid.V
    sel = np.abs(m) > 0 if mask is None else (mask & (np.abs(m) > 0))
    a = grid.RR[sel]; zp = grid.ZZ[sel]; mm = m[sel]
    tR = np.atleast_1d(tR).astype(float); tz = np.atleast_1d(tz).astype(float)
    gR = np.zeros(tR.size); gz = np.zeros(tR.size); ph = np.zeros(tR.size)
    for i in range(tR.size):
        R, z = tR[i], tz[i]
        for s in range(0, a.size, 200000):
            aa, zz, mmm = a[s:s + 200000], zp[s:s + 200000], mm[s:s + 200000]
            dz = z - zz
            Q = (R + aa) ** 2 + dz ** 2
            D2 = (aa - R) ** 2 + dz ** 2
            ok = D2 > 1e-12
            k2 = np.where(ok, 4 * aa * R / Q, 0.0)
            Kk = ellipk(k2); Ee = ellipe(k2); sq = np.sqrt(Q)
            gz[i] += np.sum(np.where(ok, -2 * G * mmm / np.pi * dz * Ee / (np.where(ok, D2, 1) * sq), 0.0))
            if R > 0:
                gR[i] += np.sum(np.where(ok, -G * mmm / (np.pi * R * sq) * (Kk - (aa ** 2 - R ** 2 + dz ** 2) / np.where(ok, D2, 1) * Ee), 0.0))
            if potential:
                ph[i] += np.sum(np.where(ok, -2 * G * mmm / np.pi * Kk / sq, 0.0))
    return (gR, gz, ph) if potential else (gR, gz)


# ----------------------------------------------------------------------------- the law
def heat_fields(comps, grid, u):
    """S and g_hot of the hot free-streaming stars (spherical stand-ins; k = 3 sigma^2 / u^2)."""
    r = np.sqrt(grid.RR ** 2 + grid.ZZ ** 2) + 1e-9
    s = np.geomspace(1e-3, 3000.0, 1200)
    dmk = np.zeros(s.size)
    for c in comps:
        if getattr(c, 'sigma', 0.0) > 0 and c.kind == 'stars' and not isinstance(c, Disk):
            M = c.m_profile(s)
            dmk += 3 * c.sigma ** 2 / u ** 2 * np.diff(np.concatenate([[0.0], M]))
    if not dmk.any():
        return np.zeros(r.shape), np.zeros(r.shape), np.zeros(r.shape)
    # S = G sum k dm <1/d^2>_shell  (law.shell_weights logic, vectorised in chunks)
    S = np.zeros(r.size); rf = r.ravel()
    for i in range(0, rf.size, 4000):
        x = s[None, :] / rf[i:i + 4000, None]
        w = np.where(x < 1, np.arctanh(np.clip(x, 1e-12, 1 - 1e-12)) / np.clip(x, 1e-12, None),
                     0.5 * np.log((x + 1) / np.maximum(x - 1, 1e-12)) / x)
        S[i:i + 4000] = G * (w @ dmk) / rf[i:i + 4000] ** 2
    Mk = np.interp(r, s, np.cumsum(dmk))
    gh = -G * Mk / r ** 2
    return S.reshape(r.shape), gh * grid.RR / r, gh * grid.ZZ / r


def law_extra(gR, gz, S, hR, hz, consts, law='ours', a0_SI=1.2e-10, reach=None, grid=None):
    """h - g_N for our law (or QUMOND with the simple function, or Newton)."""
    g = np.sqrt(gR ** 2 + gz ** 2) + 1e-30
    if law == 'newton':
        return np.zeros_like(gR), np.zeros_like(gz)
    if law == 'mond':
        a0 = a0_SI / K_SI
        nu = 0.5 + np.sqrt(0.25 + a0 / g)
        F = (nu - 1) * g
        return F * gR / g, F * gz / g
    a, gd = consts['a_code'], consts['lam'] * consts['a_code']
    F = np.exp(-g / gd) * np.sqrt(a * (g + S))
    fR, fz = gR + hR, gz + hz
    den = g + np.sqrt(hR ** 2 + hz ** 2)
    eR, ez = F * fR / den, F * fz / den
    if reach is not None:
        r = np.sqrt(grid.RR ** 2 + grid.ZZ ** 2)
        cut = r < reach
        eR, ez = eR * cut, ez * cut
    return eR, ez


def divergence(grid, AR, Az):
    """(1/R) d(R A_R)/dR + dA_z/dz on the nonuniform grid."""
    t1 = np.gradient(grid.RR * AR, grid.R, axis=0, edge_order=2) / grid.RR
    t2 = np.gradient(Az, grid.z, axis=1, edge_order=2)
    return t1 + t2


def solve(comps, grid, consts, law='ours', reach=None, gN=None, heat=None):
    """Returns the phantom density and the pieces needed to evaluate the pull anywhere."""
    gR, gz = gN if gN is not None else newton_field(comps, grid)
    S, hR, hz = heat if heat is not None else heat_fields(comps, grid, consts['u_kms'])
    eR, ez = law_extra(gR, gz, S, hR, hz, consts, law=law, reach=reach, grid=grid)
    rho_ph = -divergence(grid, eR, ez) / (4 * np.pi * G)
    return dict(gR=gR, gz=gz, S=S, eR=eR, ez=ez, rho_ph=rho_ph)
