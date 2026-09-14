"""Axisymmetric AQUAL field solver for RPG-1.

Solves div[mu(|grad Phi|/a) grad Phi] = 4 pi G rho with mu(x) = x/(1+x); a = 0 gives Newtonian
gravity. Finite volumes on a spherical-polar grid, logarithmic in r, with latitude cells clustered
at the equator and equatorial symmetry. Gauss's law fixes the flux through the inner and outer
spheres from the enclosed mass (monopole). The nonlinearity is handled by a Picard (Kacanov)
iteration relaxed in log K with weight (1+x)/(2+x), which removes the first-order error of the
spherical iteration at every x. Units: kpc, km/s, Msun.
"""
import numpy as np
import scipy.sparse as sp
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import splu

G = 4.30091727003628e-6                 # kpc (km/s)^2 / Msun
C_KMS = 299792.458
KPC_M = 3.0856775814913673e19
SI_TO_CODE = KPC_M/1e6                  # m s^-2 -> (km/s)^2 / kpc


def mu(x):
    return x/(1 + x)


def nu(y):
    """Inverse of x mu(x): g = nu(g_N/a) g_N."""
    return .5 + np.sqrt(.25 + 1/np.maximum(y, 1e-300))


def field_energy_function(y):
    """F(y) with F'(y) = mu(sqrt y), F(0) = 0. Field energy density: a^2 F(|grad Phi|^2/a^2)/(8 pi G)."""
    s = np.sqrt(y)
    return y - 2*s + 2*np.log1p(s)


def _gl(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return .5*(x + 1), .5*w


class Grid:
    """Cells in r (logarithmic) and polar angle theta in [0, pi/2]; the equator is theta = pi/2."""

    def __init__(self, r_in, r_out, nr, nth, beta=7.):
        self.rf = np.geomspace(r_in, r_out, nr + 1)
        self.rc = np.sqrt(self.rf[:-1]*self.rf[1:])
        lat = (np.pi/2)*np.sinh(beta*np.linspace(0, 1, nth + 1))/np.sinh(beta)
        self.thf = np.pi/2 - lat[::-1]
        self.thc = .5*(self.thf[:-1] + self.thf[1:])
        self.nr, self.nth = nr, nth
        cf = np.cos(self.thf)
        self.omega = 2*np.pi*(cf[:-1] - cf[1:])                      # sums to 2 pi over the hemisphere
        self.volume = ((self.rf[1:]**3 - self.rf[:-1]**3)/3)[:, None]*self.omega


def spherical_cell_masses(grid, enclosed):
    """Hemisphere cell masses of a spherical component with enclosed mass M(<r), and M(<r_in)."""
    M = enclosed(grid.rf)
    return np.diff(M)[:, None]*grid.omega/(4*np.pi), float(M[0])


def exponential_disk_cell_masses(grid, sigma, h, nq_r=4, nq_u=6):
    """rho = Sigma(R) exp(-|z|/h)/(2h). The vertical integral uses u = exp(-z/h), so it stays exact in z
    for cells thicker than h, up to the slow variation of Sigma(sqrt(r^2 - z^2)) across the cell."""
    xr, wr = _gl(nq_r)
    xu, wu = _gl(nq_u)
    rf, cf = grid.rf, np.cos(grid.thf)
    out = np.zeros((grid.nr, grid.nth))
    for q in range(nq_r):
        r = rf[:-1] + np.diff(rf)*xr[q]
        z_top, z_bot = r[:, None]*cf[:-1], r[:, None]*cf[1:]          # cell j spans z_bot..z_top at radius r
        u_top, u_bot = np.exp(-z_top/h), np.exp(-z_bot/h)
        s = 0.
        for p in range(nq_u):
            z = -h*np.log(np.maximum(u_top + (u_bot - u_top)*xu[p], 1e-300))
            s = s + wu[p]*sigma(np.sqrt(np.maximum(r[:, None]**2 - z*z, 0.)))
        out += 2*np.pi*(r*np.diff(rf)*wr[q])[:, None]*.5*(u_bot - u_top)*s
    return out, float(sigma(np.zeros(1))[0]*np.pi*rf[0]**2)


def smooth_cell_masses(grid, rho, nq_r=4, nq_z=8):
    """Hemisphere cell masses of a smooth axisymmetric density rho(R, z) by Gauss-Legendre quadrature."""
    xr, wr = _gl(nq_r)
    xz, wz = _gl(nq_z)
    rf, cf = grid.rf, np.cos(grid.thf)
    out = np.zeros((grid.nr, grid.nth))
    for q in range(nq_r):
        r = rf[:-1] + np.diff(rf)*xr[q]
        z_top, z_bot = r[:, None]*cf[:-1], r[:, None]*cf[1:]
        s = 0.
        for p in range(nq_z):
            z = z_bot + (z_top - z_bot)*xz[p]
            s = s + wz[p]*rho(np.sqrt(np.maximum(r[:, None]**2 - z*z, 0.)), z)
        out += 2*np.pi*(r*np.diff(rf)*wr[q])[:, None]*(z_top - z_bot)*s
    return out, float(rho(np.zeros(1), np.zeros(1))[0]*4*np.pi*rf[0]**3/3)


def thin_disk_cell_masses(grid, enclosed_cylindrical):
    """Razor-thin disk with cylindrical enclosed mass M(<R): half of each annulus in the equatorial row."""
    out = np.zeros((grid.nr, grid.nth))
    M = enclosed_cylindrical(grid.rf)
    out[:, -1] = .5*np.diff(M)
    return out, float(M[0])


class Solver:
    """Finite-volume AQUAL solve for given hemisphere cell masses (plus a mass inside the inner sphere)."""

    def __init__(self, grid, cell_mass, inner_mass=0.):
        g = self.grid = grid
        self.m_in = float(inner_mass)
        self.m_tot = 2*float(cell_mass.sum()) + self.m_in
        self.dr, self.dth = np.diff(g.rc), np.diff(g.thc)
        self.cr = g.rf[1:-1, None]**2*g.omega/self.dr[:, None]                         # interior radial faces
        self.ct = (np.pi*np.sin(g.thf[1:-1])*(g.rf[1:]**2 - g.rf[:-1]**2)[:, None]
                   /(g.rc[:, None]*self.dth))                                           # interior polar faces
        b = 4*np.pi*G*cell_mass
        b[0] += G*self.m_in*g.omega
        b[-1] -= G*self.m_tot*g.omega
        self.b_full = b.ravel().copy()
        idx = np.arange(g.nr*g.nth).reshape(g.nr, g.nth)
        self.pin = int(idx[-1, -1])                    # gauge: Phi = 0 in the outermost equatorial cell
        self.b = self.b_full.copy()
        self.b[self.pin] = 0.
        i1, i2, j1, j2 = idx[:-1].ravel(), idx[1:].ravel(), idx[:, :-1].ravel(), idx[:, 1:].ravel()
        rows = np.concatenate([i1, i2, j1, j2, i1, i2, j1, j2])
        cols = np.concatenate([i2, i1, j2, j1, i1, i2, j1, j2])
        self._keep = rows != self.pin
        self._rows = np.append(rows[self._keep], self.pin)
        self._cols = np.append(cols[self._keep], self.pin)
        self.a = 0.

    def _linear(self, Kr, Kt):
        wr, wt = (Kr*self.cr).ravel(), (Kt*self.ct).ravel()
        vals = np.append(np.concatenate([wr, wr, wt, wt, -wr, -wr, -wt, -wt])[self._keep], 1.)
        n = self.grid.nr*self.grid.nth
        return splu(sp.csc_matrix((vals, (self._rows, self._cols)), shape=(n, n))).solve(self.b)

    def boundary_gradient(self, M, r):
        gN = G*M/r**2
        return gN*nu(gN/self.a) if self.a > 0 else gN

    def gradients(self, phi):
        """dPhi/dr on radial faces, (1/r) dPhi/dtheta on polar faces, both at cell centers, and |grad Phi| on
        interior faces."""
        g = self.grid
        P = phi.reshape(g.nr, g.nth)
        gr = np.empty((g.nr + 1, g.nth))
        gr[1:-1] = np.diff(P, axis=0)/self.dr[:, None]
        gr[0] = self.boundary_gradient(self.m_in, g.rf[0])
        gr[-1] = self.boundary_gradient(self.m_tot, g.rf[-1])
        gt = np.zeros((g.nr, g.nth + 1))
        gt[:, 1:-1] = np.diff(P, axis=1)/(g.rc[:, None]*self.dth)
        grc, gtc = .5*(gr[:-1] + gr[1:]), .5*(gt[:, :-1] + gt[:, 1:])
        mag_r = np.hypot(gr[1:-1], .5*(gtc[:-1] + gtc[1:]))
        mag_t = np.hypot(gt[:, 1:-1], .5*(grc[:, :-1] + grc[:, 1:]))
        return gr, gt, grc, gtc, mag_r, mag_t

    def conductances(self, phi):
        *_, mr, mt = self.gradients(phi)
        if self.a <= 0:
            return np.ones_like(mr), np.ones_like(mt)
        return mu(np.maximum(mr/self.a, 1e-12)), mu(np.maximum(mt/self.a, 1e-12))

    def operator(self, phi, Kr, Kt):
        g = self.grid
        P = phi.reshape(g.nr, g.nth)
        fr, ft = Kr*self.cr*np.diff(P, axis=0), Kt*self.ct*np.diff(P, axis=1)
        out = np.zeros_like(P)
        out[:-1] += fr
        out[1:] -= fr
        out[:, :-1] += ft
        out[:, 1:] -= ft
        return out.ravel()

    def solve(self, a=0., tol=1e-11, maxit=200):
        self.a = float(a)
        Kr, Kt = np.ones((self.grid.nr - 1, self.grid.nth)), np.ones((self.grid.nr, self.grid.nth - 1))
        phi = self._linear(Kr, Kt)
        self.history = []
        if self.a > 0:
            for _ in range(maxit):
                *_, mr, mt = self.gradients(phi)
                xr, xt = np.maximum(mr/self.a, 1e-12), np.maximum(mt/self.a, 1e-12)
                br, bt = (1 + xr)/(2 + xr), (1 + xt)/(2 + xt)
                Kr = np.exp(br*np.log(mu(xr)) + (1 - br)*np.log(Kr))
                Kt = np.exp(bt*np.log(mu(xt)) + (1 - bt)*np.log(Kt))
                new = self._linear(Kr, Kt)
                change = float(np.max(np.abs(new - phi))/np.ptp(new))
                phi = new
                self.history.append(change)
                if change < tol:
                    break
        self.phi = phi
        Kr, Kt = self.conductances(phi)
        res = self.operator(phi, Kr, Kt) - self.b_full
        self.residual = float(np.max(np.abs(res))/np.max(np.abs(self.b_full)))
        self.converged = self.a <= 0 or self.history[-1] < tol
        return phi

    def midplane_speed(self, R):
        """Circular speed from the radial gradient on the equatorial row (latitude below 2e-4 rad)."""
        gr = self.gradients(self.phi)[0]
        v2 = self.grid.rf*gr[:, -1]
        return np.sqrt(np.maximum(np.interp(np.log(R), np.log(self.grid.rf), v2), 0.))

    def gradient_at(self, R, z):
        """(dPhi/dR, dPhi/d|z|) at points; the acceleration is minus these."""
        _, _, grc, gtc, _, _ = self.gradients(self.phi)
        pts = np.column_stack([np.log(np.hypot(R, z)), np.arctan2(R, np.abs(z))])
        axes = (np.log(self.grid.rc), self.grid.thc)
        ar = RegularGridInterpolator(axes, grc, bounds_error=False, fill_value=None)(pts)
        at = RegularGridInterpolator(axes, gtc, bounds_error=False, fill_value=None)(pts)
        t = pts[:, 1]
        return np.sin(t)*ar + np.cos(t)*at, np.cos(t)*ar - np.sin(t)*at

    def energy_profile(self):
        """Cumulative field energy inside each radial face: AQUAL F-energy, or |grad Phi|^2/(8 pi G) if a = 0."""
        _, _, grc, gtc, _, _ = self.gradients(self.phi)
        s2 = grc**2 + gtc**2
        dens = (self.a**2*field_energy_function(s2/self.a**2) if self.a > 0 else s2)/(8*np.pi*G)
        return np.concatenate([[0.], np.cumsum((2*dens*self.grid.volume).sum(1))])
