"""CF-1 transport controls (see protocol.md): fixed spherical Newtonian potentials built from ordinary matter plus
counted seeds, an isotropic bath specified at an exterior boundary R_b, and the entry rate, unbound density and
no-capture limit computed from trajectories. Units: kpc, km/s, Msun, and time in kpc/(km/s).

With Phi(R_b) = 0, a companion arriving with speed u has w^2(r) = u^2 + v_esc^2(r) at radius r, where
v_esc^2(r) = -2 Phi(r) is the escape speed relative to R_b. Mass outside R_b is ignored, as declared.
"""
import math
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.optimize import brentq

G = 4.30091727003628e-6


class Potential:
    """Spherical potential from an enclosed-mass profile M(<r) on an increasing grid, referenced to Phi(R_b) = 0."""

    def __init__(self, r, M, R_b, label=''):
        keep = r < R_b                            # strict: R_b is appended below, and a repeated node breaks gradients
        self.r = np.r_[r[keep], R_b]
        self.M = np.r_[M[keep], np.interp(R_b, r, M)]
        self.R_b, self.label = R_b, label
        g = G*self.M/self.r**2
        # Phi(r) = -int_r^{R_b} g dr'
        inward = cumulative_trapezoid(g[::-1], self.r[::-1], initial=0)[::-1]
        self.phi = -np.abs(inward)
        self.lr = np.log(self.r)

    def v_esc2(self, rr):
        return -2*np.interp(np.log(rr), self.lr, self.phi)

    def accel(self, rr):
        rr = max(rr, self.r[0])                   # a plunging trial step can reach r <= 0; clamp inside the grid
        return G*np.interp(math.log(rr), self.lr, self.M)/rr**2

    def enclosed(self, rr):
        return np.interp(np.log(rr), self.lr, self.M)


def entry_analytic(pot, R, u):
    """sigma_enter = pi R^2 (1 + v_esc^2(R)/u^2), valid when no centrifugal barrier outside R blocks entry."""
    return math.pi*R*R*(1 + pot.v_esc2(R)/u**2)


def entry_barrier(pot, R, u, n=4001):
    """sigma_enter = pi L_max^2/u^2 with L_max = min over [R, R_b] of r w(r), the effective-potential limit."""
    rr = np.geomspace(R, pot.R_b, n)
    return math.pi*float(np.min(rr*rr*(u*u + pot.v_esc2(rr))))/u**2


def pericentre_by_orbit(pot, L, u):
    """Integrate an orbit inward from R_b with angular momentum L and speed u; return its closest approach."""
    def rhs(t, y):
        r, vr = y
        return [vr, L*L/r**3 - pot.accel(r)]

    def turn(t, y):
        return y[1]
    turn.terminal, turn.direction = True, 1
    vr0 = -math.sqrt(max(u*u - (L/pot.R_b)**2, 0.))
    sol = solve_ivp(rhs, (0, 50*pot.R_b/u), [pot.R_b*(1 - 1e-12), vr0], method='DOP853', rtol=1e-11, atol=1e-12,
                    events=turn)
    if sol.status != 1:
        raise FloatingPointError(f'orbit did not reach pericentre (L={L:g}): {sol.message}')
    return float(sol.y_events[0][0][0])


def entry_by_orbits(pot, R, u):
    """V1: the angular momentum whose integrated orbit has pericentre exactly R gives sigma = pi L^2/u^2."""
    hi = pot.R_b*u*(1 - 1e-9)
    if pericentre_by_orbit(pot, hi, u) <= R:      # boundary-limited: every orbit entering R_b reaches R
        return math.pi*hi*hi/u**2
    L = brentq(lambda l: pericentre_by_orbit(pot, l, u) - R, 1e-9*hi, hi, xtol=1e-10*hi, rtol=1e-12)
    return math.pi*L*L/u**2


def density_analytic(pot, r, u):
    return math.sqrt(1 + pot.v_esc2(r)/u**2)


def density_by_orbits(pot, r, u, n=4000):
    """V2: steady density of a single-speed isotropic bath at r from the time each orbit spends there.
    Orbits with L below the barrier limit reach r; each crosses the shell twice at radial speed
    sqrt(w^2 - L^2/r^2). Gauss–Legendre in s = L/L_max with the substitution that removes the end singularity."""
    w2 = u*u + pot.v_esc2(r)
    rr = np.geomspace(r, pot.R_b, 4001)
    Lmax = math.sqrt(float(np.min(rr*rr*(u*u + pot.v_esc2(rr)))))
    x, wts = np.polynomial.legendre.leggauss(n)
    smax = min(1., Lmax/(r*math.sqrt(w2)))        # s = L/(r w) = sin(theta) runs over [0, smax]
    thmax = math.asin(smax)
    th = .5*thmax*(x + 1)
    s = np.sin(th)
    # integral of 2L dL / sqrt(w^2 - L^2/r^2) with L = r w s: = 2 r^2 w int s ds/sqrt(1 - s^2) = 2 r^2 w int s dtheta
    integral = 2*r*r*math.sqrt(w2)*.5*thmax*float(np.sum(wts*s))
    # n(r)/n_inf = (pi/u) * 2 * integral / (4 pi r^2)
    return 2*math.pi/u*integral/(4*math.pi*r*r)


def no_capture_retained(pot, u, n=200, seed=1):
    """V3: integrate n incoming orbits with random L from the isotropic bath; count those still inside R_b after
    ten crossing times. Without an interaction every one must leave."""
    rng = np.random.default_rng(seed)
    # isotropic flux is uniform in L^2; the floor excludes near-radial plunges (a 1e-6 share of the flux)
    L = pot.R_b*u*np.sqrt(1e-6 + (1 - 1e-6)*rng.random(n))
    retained = 0
    for l in L:
        def rhs(t, y):
            return [y[1], l*l/y[0]**3 - pot.accel(y[0])]

        def out(t, y):
            return y[0] - pot.R_b
        out.terminal, out.direction = True, 1
        vr0 = -math.sqrt(max(u*u - (l/pot.R_b)**2, 0.))
        sol = solve_ivp(rhs, (0, 200*pot.R_b/u), [pot.R_b*(1 - 1e-9), vr0], method='DOP853', rtol=1e-10, atol=1e-12,
                        events=out)
        retained += int(sol.status != 1)
    return retained/n
