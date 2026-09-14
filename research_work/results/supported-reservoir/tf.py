"""Thomas–Fermi (n = 1 polytrope) condensate in hydrostatic equilibrium with a fixed host mass profile (CR-2, see
protocol.md). Units: kpc, km/s, Msun.

With pressure P = (K/2) rho^2, k^2 = 4 pi G/K (a condensate alone ends at R_TF = pi/k), x = k r, psi = mu - Phi =
psi_c u and M_c = psi_c m/(G k), the equilibrium
    dpsi/dr = -G [M_h frac(r) + M_c(r)]/r^2,   dM_c/dr = 4 pi r^2 psi/K
becomes
    du/dx = -(lam frac(x/k) + m)/x^2,   dm/dx = x^2 u,   u(0) = 1,   lam = G k M_h/psi_c.
For a fixed host shape and K the family has one parameter, lam. The condensate ends where u = 0, its density is
rho = psi_c u/K, and its mass is M_c(r) = M_h m(k r)/lam = psi_c m(k r)/(G k).
"""
import bisect
import math
import numpy as np
from scipy.integrate import solve_ivp

G = 4.30091727003628e-6


def k_of(R_tf):
    return math.pi/R_tf


class Family:
    """Solutions for one host shape, the enclosed fraction frac on the increasing grid r, at one k."""

    def __init__(self, r, frac, k, x0=1e-7):
        self.lr, self.fr = list(np.log(r)), list(frac)
        self.k, self.x0 = k, x0
        self.inner_power = math.log(frac[1]/frac[0])/(self.lr[1] - self.lr[0])

    def host(self, rr):
        lr = math.log(rr)
        if lr <= self.lr[0]:
            return self.fr[0]*math.exp(self.inner_power*(lr - self.lr[0]))
        if lr >= self.lr[-1]:
            return 1.
        i = bisect.bisect_right(self.lr, lr)
        w = (lr - self.lr[i - 1])/(self.lr[i] - self.lr[i - 1])
        return self.fr[i - 1] + w*(self.fr[i] - self.fr[i - 1])

    def solve(self, lam):
        """Integrate in s = ln x to the edge u = 0. Returns x_edge, m_edge and m(x) (vectorised)."""
        k, host = self.k, self.host
        x0 = self.x0/max(1., lam)          # a strong host confines the condensate to x ~ 1/lam; start well inside

        def rhs(s, y):
            x = math.exp(s)
            return [-(lam*host(x/k) + y[1])/x, x*x*x*y[0]]

        def edge(s, y):
            return y[0]
        edge.terminal, edge.direction = True, -1
        # relative control on m, which spans many decades; absolute control on u, which is of order one
        sol = solve_ivp(rhs, (math.log(x0), math.log(1e3)), [1., x0**3/3], method='DOP853', rtol=1e-10,
                        atol=[1e-13, 1e-300], events=edge, dense_output=True)
        if sol.status != 1 or not np.all(np.isfinite(sol.y)):
            raise FloatingPointError(f'Thomas-Fermi profile did not reach its edge (lam={lam:g}): {sol.message}')
        x_e = math.exp(float(sol.t_events[0][0]))
        m_e = float(sol.y_events[0][0][1])
        dense = sol.sol

        def m(x):
            x = np.asarray(x, float)
            out = np.full(x.shape, m_e)
            inner = x < x0
            mid = (~inner) & (x < x_e)
            out[inner] = x[inner]**3/3
            if mid.any():
                out[mid] = dense(np.log(x[mid]))[1]
            return out
        return dict(x_edge=x_e, m_edge=m_e, m=m, lam=lam)
