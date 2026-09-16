"""PM-2A stage B: the four field equations, each with its own mu, spherical inverse, relaxation weight
and field functional (protocol-pm2a.md, amended fa8b407).

Every equation is div[mu(|grad U|/a_star) grad U] = 4 pi G rho for its own mu. Changing mu alone is not
enough: the boundary gradient, the log-conductance relaxation weight and the functional all follow from
the same mu, and treating the three implementations as one equation is what the distinction test exists
to catch.

    equation                 mu(x)                  spherical x(y)        w = 1/(1+dln mu/dln x)   F(u)
    Newtonian                1                      y                     1                        u
    Completion I (total)     4x/[1+sqrt(1+4x)]^2    y + sqrt(y)           s/(1+s), s=sqrt(1+4x)    t^3(t+2/3)
    Completion II (aux psi)  x                      sqrt(y)               1/2                      (2/3)u^(3/2)
    simple-mu (comparison)   x/(1+x)                (y+sqrt(y^2+4y))/2    (1+x)/(2+x)              u-2sqrt(u)+2ln(1+sqrt u)

with x = |grad U|/a_star, y = g_N/a_star, u = x^2 and t = 2x/[1+sqrt(1+4x)], which equals sqrt(y) on the
spherical solution. RPG-1's aqual.py is imported and subclassed, never modified: its own archive stays
reproducible.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'radiation-polarized-gravity'))
import aqual as Q  # noqa: E402

G = Q.G
SERIES_X = .1     # below this the simple-mu functional is summed as a series, not evaluated as written


def _t(x):
    """t = 2x/[1 + sqrt(1+4x)] = sqrt(1+4x) - 1 over 2; equals sqrt(y) on the spherical solution."""
    return 2*np.asarray(x, float)/(1 + np.sqrt(1 + 4*np.asarray(x, float)))


def _simple_F(u):
    """u - 2 sqrt(u) + 2 ln(1 + sqrt u). Written as given, its leading terms cancel to O(x^3): at x = 1e-8
    it returns exactly 0.0 against a true 6.667e-25, so below SERIES_X it is summed as its series
    2 x^3 sum_k (-x)^k/(k+3) instead."""
    x = np.sqrt(np.asarray(u, float))
    written = x*x - 2*x + 2*np.log1p(x)
    xs = np.minimum(x, SERIES_X)                        # clipped so the unused branch cannot overflow
    ser = np.zeros_like(xs)
    for k in range(60, -1, -1):                         # Horner from the tail
        ser = ser*(-xs) + 1./(k + 3)
    return np.where(x < SERIES_X, 2*xs**3*ser, written)


def _completion_I_F(u):
    """t^3 (t + 2/3): free of cancellation, and (2/3)x^3 exactly as x -> 0."""
    t = _t(np.sqrt(np.asarray(u, float)))
    return t**3*(t + 2/3)


class Equation:
    """One field equation's constitutive functions. x_of_y is the spherical inverse, evaluated directly
    rather than as y*nu(y), because nu diverges as y -> 0 while x does not."""

    def __init__(self, name, mu, x_of_y, weight, F, nu=None):
        self.name, self.mu, self.x_of_y, self.weight, self.F, self.nu = name, mu, x_of_y, weight, F, nu

    def __repr__(self):
        return f'Equation({self.name})'


NEWTONIAN = Equation(
    'newtonian', lambda x: np.ones_like(np.asarray(x, float)), lambda y: np.asarray(y, float),
    lambda x: np.ones_like(np.asarray(x, float)), lambda u: np.asarray(u, float),
    nu=lambda y: np.ones_like(np.asarray(y, float)))

COMPLETION_I = Equation(
    'completion_I', lambda x: 4*np.asarray(x, float)/(1 + np.sqrt(1 + 4*np.asarray(x, float)))**2,
    lambda y: np.asarray(y, float) + np.sqrt(np.asarray(y, float)),
    lambda x: np.sqrt(1 + 4*np.asarray(x, float))/(1 + np.sqrt(1 + 4*np.asarray(x, float))),
    _completion_I_F, nu=lambda y: 1 + 1/np.sqrt(np.asarray(y, float)))

COMPLETION_II_AUX = Equation(
    'completion_II_aux', lambda x: np.asarray(x, float), lambda y: np.sqrt(np.asarray(y, float)),
    lambda x: np.full_like(np.asarray(x, float), .5), lambda u: (2/3)*np.asarray(u, float)**1.5,
    nu=lambda y: 1/np.sqrt(np.asarray(y, float)))

SIMPLE = Equation(
    'simple_mu', lambda x: np.asarray(x, float)/(1 + np.asarray(x, float)),
    lambda y: .5*(np.asarray(y, float) + np.sqrt(np.asarray(y, float)**2 + 4*np.asarray(y, float))),
    lambda x: (1 + np.asarray(x, float))/(2 + np.asarray(x, float)), _simple_F,
    nu=lambda y: .5 + np.sqrt(.25 + 1/np.maximum(np.asarray(y, float), 1e-300)))

EQUATIONS = {e.name: e for e in (NEWTONIAN, COMPLETION_I, COMPLETION_II_AUX, SIMPLE)}


class Solver(Q.Solver):
    """RPG-1's finite-volume solver driven by a named equation instead of its hardcoded simple mu."""

    def __init__(self, grid, cell_mass, inner_mass=0., equation=NEWTONIAN):
        super().__init__(grid, cell_mass, inner_mass)
        self.eq = equation

    def boundary_gradient(self, M, r):
        """The selected equation's own spherical inverse, from the enclosed mass through y = GM/(r^2 a)."""
        gN = G*M/r**2
        if self.a <= 0:
            return gN
        return self.a*self.eq.x_of_y(gN/self.a)

    def conductances(self, phi):
        *_, mr, mt = self.gradients(phi)
        if self.a <= 0:
            return np.ones_like(mr), np.ones_like(mt)
        return self.eq.mu(np.maximum(mr/self.a, 1e-12)), self.eq.mu(np.maximum(mt/self.a, 1e-12))

    def solve(self, a=0., tol=1e-11, res_tol=1e-8, maxit=400):
        """Picard iteration relaxed in log K with this equation's own weight. Convergence requires BOTH a
        small iteration change and a small equation residual: aqual.py sets its flag from the change
        alone, which can call a run converged while its residual is not."""
        self.a = float(a)
        Kr = np.ones((self.grid.nr - 1, self.grid.nth))
        Kt = np.ones((self.grid.nr, self.grid.nth - 1))
        phi = self._linear(Kr, Kt)
        self.history = []
        if self.a > 0:
            for _ in range(maxit):
                *_, mr, mt = self.gradients(phi)
                xr, xt = np.maximum(mr/self.a, 1e-12), np.maximum(mt/self.a, 1e-12)
                br, bt = self.eq.weight(xr), self.eq.weight(xt)
                Kr = np.exp(br*np.log(self.eq.mu(xr)) + (1 - br)*np.log(Kr))
                Kt = np.exp(bt*np.log(self.eq.mu(xt)) + (1 - bt)*np.log(Kt))
                new = self._linear(Kr, Kt)
                change = float(np.max(np.abs(new - phi))/np.ptp(new))
                phi = new
                self.history.append(change)
                if change < tol:
                    break
        self.phi = phi
        Kr, Kt = self.conductances(phi)
        res = (self.operator(phi, Kr, Kt) - self.b_full).reshape(self.grid.nr, self.grid.nth)
        scale = np.max(np.abs(self.b_full))
        self.residual = float(np.max(np.abs(res))/scale)
        off = np.ones(res.size, bool)
        off[self.pin] = False                       # the gauge cell's equation was replaced by the pin
        self.residual_off_pin = float(np.max(np.abs(res.ravel()[off]))/scale)
        self.iteration_change = float(self.history[-1]) if self.history else 0.
        self.converged = bool(self.iteration_change < tol and self.residual_off_pin < res_tol)
        return phi

    def radial_gradient_at(self, r):
        """Equatorial dU/dr off the grid, read from the interior FACE gradients through a cubic spline in
        ln r.

        This is a different numerical path from the inherited readouts, and a far more accurate one.
        `Q.gradient_at` interpolates cell-centred values formed by averaging the two bounding faces, and
        `Q.midplane_speed` interpolates face quantities linearly; both cost O(h^2) and sit near 1e-4 at the
        archive's resolution while the face gradients themselves are accurate to 1e-7. Reading the faces
        directly keeps the solution's own accuracy, so stage C's galaxy readout uses this.
        """
        rf = self.grid.rf[1:-1]
        faces = self.gradients(self.phi)[0][1:-1][:, -1]
        return CubicSpline(np.log(rf), faces)(np.log(np.asarray(r, float)))

    def field_functional_profile(self):
        """Cumulative a^2 F(|grad U|^2/a^2)/(8 pi G) inside each radial face; |grad U|^2/(8 pi G) if a = 0."""
        _, _, grc, gtc, _, _ = self.gradients(self.phi)
        s2 = grc**2 + gtc**2
        dens = (self.a**2*self.eq.F(s2/self.a**2) if self.a > 0 else s2)/(8*np.pi*G)
        return np.concatenate([[0.], np.cumsum((2*dens*self.grid.volume).sum(1))])


def spherical_solve(equation, enclosed, a_star, r_in, r_out, nr=400, nth=24):
    """Solve one equation for a spherical source given its cumulative mass M(<r)."""
    grid = Q.Grid(r_in, r_out, nr, nth)
    cell, m_in = Q.spherical_cell_masses(grid, enclosed)
    s = Solver(grid, cell, m_in, equation)
    s.solve(a_star if equation is not NEWTONIAN else 0.)
    return s
