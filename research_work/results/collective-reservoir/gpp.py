"""Spherical Gross-Pitaevskii-Poisson solver for CR-1.

The wavefunction is normalized to mass density, rho = |psi|^2, and obeys
    i psi_t = -(kappa/2) Lap psi + (Phi/kappa) psi,   Lap Phi_self = 4 pi G rho,   kappa = hbar/m.
Code units hbar = m = G = 1 are kappa = G = 1. The physical runs use kpc, km/s and Msun, with time in
kpc/(km/s). psi = u/r on a uniform grid r_j = j dr (j = 1..N), with u(0) = u(R) = 0.

Strang splitting: each potential half-step is an exact phase rotation, which leaves |u| and hence the
self-gravity unchanged within the step; the kinetic step is exact in the discrete sine basis (DST-I).
Imaginary time relaxes to stationary states. Optional pieces:
- an absorbing layer near R, whose removed mass is counted;
- source terms applied in the potential half-steps:
  'mode' scales the whole wavefunction (stimulated growth into the occupied mode);
  'local' adds density s(r) at the local phase.
"""
import numpy as np
from scipy.fft import dst, idst

FOUR_PI = 4*np.pi


class Grid:
    def __init__(self, R, N):
        self.N, self.R = N, float(R)
        self.dr = R/(N + 1)
        self.r = self.dr*np.arange(1, N + 1)
        self.k2 = (np.pi*np.arange(1, N + 1)/R)**2

    def sine_filter(self, u, f):
        return idst(f*dst(u.real, type=1), type=1) + 1j*idst(f*dst(u.imag, type=1), type=1)


def _cumulative(f, dr):
    """Trapezoid cumulative integral from r = 0 (f(0) = 0) to each grid point."""
    return dr*(np.cumsum(f) - .5*f)


class GPP:
    def __init__(self, grid, kappa=1., G=1., phi_ext=None, absorb_from=None, absorb_rate=0.):
        self.grid, self.kappa, self.G = grid, float(kappa), float(G)
        r = grid.r
        self.phi_ext = np.zeros_like(r) if phi_ext is None else np.asarray(phi_ext, float)
        self.W = np.zeros_like(r)
        if absorb_from is not None:
            x = np.clip((r - absorb_from)/(grid.R - absorb_from), 0, 1)
            self.W = absorb_rate*x*x
        self.absorbed = 0.
        self.injected = 0.

    def mass(self, u):
        return FOUR_PI*self.grid.dr*float(np.sum(np.abs(u)**2))

    def self_potential(self, u):
        r, dr = self.grid.r, self.grid.dr
        q = FOUR_PI*np.abs(u)**2                                   # dM/dr
        M = _cumulative(q, dr)
        outer = _cumulative((q/r)[::-1], dr)[::-1]                 # int_r^R q/r' dr'
        return -self.G*(M/r + outer)

    def potential(self, u):
        return self.self_potential(u) + self.phi_ext

    def energies(self, u):
        dr = self.grid.dr
        q = FOUR_PI*np.abs(u)**2
        lap = self.grid.sine_filter(u, self.grid.k2)
        T = .5*self.kappa**2*FOUR_PI*dr*float(np.real(np.sum(np.conj(u)*lap)))
        W = .5*dr*float(np.sum(q*self.self_potential(u)))
        X = dr*float(np.sum(q*self.phi_ext))
        M = self.mass(u)
        return dict(kinetic=T, self_gravity=W, external=X, total=T + W + X, mass=M, mu=(T + 2*W + X)/M)

    def stationary_residual(self, u):
        """||H u - mu u||/(|mu| ||u||) with H u = (kappa^2/2)(-u_rr) + V u: zero for a stationary state."""
        Hu = .5*self.kappa**2*self.grid.sine_filter(u, self.grid.k2) + self.potential(u)*u
        mu = float(np.real(np.vdot(u, Hu))/np.vdot(u, u).real)
        return float(np.sqrt(np.sum(np.abs(Hu - mu*u)**2)/np.sum(np.abs(u)**2))/abs(mu))

    def _source(self, u, dt, source):
        """('rate', Mdot): fixed mass rate into the occupied mode (validation).
        ('mode', q, K) or ('local', q, K): conversion mass-rate density q scaled by K, collected only inside
        r_99, the radius enclosing 99% of the current condensate mass. 'mode' puts the collected mass into
        the occupied mode (uniform scaling); 'local' adds it where it is converted, at the local phase."""
        if source is None:
            return u
        if source[0] == 'rate':
            add = source[1]*dt
            self.injected += add
            return u*np.sqrt(1 + add/self.mass(u))
        kind, q, K = source
        r, dr = self.grid.r, self.grid.dr
        M = _cumulative(FOUR_PI*np.abs(u)**2, dr)
        inside = r <= np.interp(.99*M[-1], M, r)
        if kind == 'mode':
            add = K*FOUR_PI*dr*float(np.sum((r*r*q)[inside]))*dt
            self.injected += add
            return u*np.sqrt(1 + add/self.mass(u))
        add = np.where(inside, r*r*K*q*dt, 0.)
        self.injected += FOUR_PI*dr*float(np.sum(add))
        a = np.abs(u)
        phase = np.where(a > 0, u/np.where(a > 0, a, 1.), 1.)
        return phase*np.sqrt(a*a + add)

    def collection_rate(self, u, q):
        """Mass rate collected inside r_99 for conversion mass-rate density q."""
        r, dr = self.grid.r, self.grid.dr
        M = _cumulative(FOUR_PI*np.abs(u)**2, dr)
        inside = r <= np.interp(.99*M[-1], M, r)
        return FOUR_PI*dr*float(np.sum((r*r*q)[inside]))

    def _absorb(self, u, dt):
        if not self.W.any():
            return u
        f = np.exp(-self.W*dt)
        self.absorbed += FOUR_PI*self.grid.dr*float(np.sum(np.abs(u)**2*(1 - f*f)))
        return u*f

    def step(self, u, dt, source=None):
        u = self._source(u, .5*dt, source)
        u = self._absorb(u*np.exp(-.5j*dt*self.potential(u)/self.kappa), .5*dt)
        u = self.grid.sine_filter(u, np.exp(-.5j*dt*self.kappa*self.grid.k2))
        u = self._absorb(u*np.exp(-.5j*dt*self.potential(u)/self.kappa), .5*dt)
        return self._source(u, .5*dt, source)

    def relax(self, u0, mass, dts=(.05, .01, .002, .0005, .0002), tol=1e-14, maxit=400000):
        """Imaginary-time ground state at fixed mass; dts are in units of 1/(kappa k_scale^2), set by the caller."""
        u = u0.astype(complex)
        u *= np.sqrt(mass/self.mass(u))
        for dt in dts:
            e_old = np.inf
            for it in range(maxit):
                # The potential is frozen within a step, so the fixed point is the ground state of the symmetric
                # operator exp(-V dt/2) exp(-T dt) exp(-V dt/2) (a second-order bias; re-evaluating V between the
                # half-steps leaves a first-order one). It is also shifted by its minimum: the constant only
                # rescales the norm, which is restored every step, and the factors cannot overflow in deep wells.
                V = self.potential(u)
                half = np.exp(-.5*dt*(V - V.min())/self.kappa)
                u = half*self.grid.sine_filter(half*u, np.exp(-.5*dt*self.kappa*self.grid.k2))
                u *= np.sqrt(mass/self.mass(u))
                if it % 200 == 0:
                    e = self.energies(u)['total']
                    if not np.isfinite(e):
                        raise FloatingPointError('imaginary-time relaxation produced a non-finite state')
                    if abs(e - e_old) < tol*abs(e):
                        break
                    e_old = e
        return u.real.astype(complex)


def half_mass_radius(grid, u):
    M = _cumulative(FOUR_PI*np.abs(u)**2, grid.dr)
    return float(np.interp(.5*M[-1], M, grid.r))
