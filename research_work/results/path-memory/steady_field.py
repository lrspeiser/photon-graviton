"""The steady written field of a spherical or projected source (CL-1, protocol-cl1.md).

The footprint is the three-dimensional isotropic Gaussian K_w(d) = exp(-d^2/2w^2), whose restriction to a plane
is RUT-1's footprint. Writing proportional to mass at rate alpha, retained for tau_keep, has the steady state
C = Lambda (rho * K_w) with Lambda = alpha tau_keep, Phi_mem = -C and a_mem = grad C. Everything here is PER UNIT
LAMBDA; the caller multiplies. ell = G/Lambda is a length and w/ell = Lambda w/G the dimensionless writing
strength.

Spherical source. Averaging the footprint over a shell of radius r' seen from radius r gives

    k_w(r, r') = exp[-(r^2 + r'^2)/2w^2] sinh(u)/u,  u = r r'/w^2,

written as exp[-(r - r')^2/2w^2] (1 - e^{-2u})/(2u) so that nothing overflows. Its r-derivative and the shell
average of the footprint's Laplacian are analytic; the derivative of sinh(u)/u is (u cosh u - sinh u)/u^2, in the
same scaled form, summed as a series where the scaled form cancels. The field is C(r) = Lambda int k_w dM(r'),
the extra inward acceleration g_mem = -dC/dr, the equivalent Newtonian mass M_eff = r^2 g_mem/G and the equivalent
density rho_eff = -lap C/(4 pi G). Because int lap K_w d^3x = 0, the written field of a bounded source has no net
equivalent mass -- a theorem this module lets the driver verify rather than assume.

Projected source. Under the light rule L1 the memory field lenses as rho_eff. Since int d^2C/dz^2 dz = 0, the
projected equivalent density is a two-dimensional Laplacian of F = Sigma_b * K_2D, the projected baryons
convolved with the two-dimensional Gaussian, whose azimuthal average is RUT-1's ring kernel
exp[-(R - R')^2/2w^2] I0e(R R'/w^2):

    M_proj,eff(<R) = -sqrt(2 pi) w R F'(R)/(2G),   Sigma_eff(R) = -(sqrt(2 pi) w/(4 pi G)) [F'' + F'/R].

This route shares no kernel with the spherical one, and the driver gates their agreement.

Units: kpc, km/s, Msun (G as in aqual.py). Lambda in (km/s)^2 per Msun.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import brentq
from scipy.special import gammainc, i0e, i1e

G = 4.30091727003628e-6                 # kpc (km/s)^2 / Msun
C_KMS = 299792.458
KPC_M = 3.0856775814913673e19
MSUN_KG, MP_G = 1.98847e30, 1.67262192e-24
SI_TO_CODE = KPC_M/1e6                  # m s^-2 -> (km/s)^2 / kpc
ARCSEC = 206264.806


# ---------------------------------------------------------------- scaled Bessel-like factors
def _bracket(u):
    """(1 - e^{-2u})/(2u) = e^{-u} sinh(u)/u, the scaled shell factor; 1 at u = 0."""
    u = np.asarray(u, float)
    out = np.empty_like(u)
    s = u < 1e-6
    out[s] = 1 - u[s] + (2/3)*u[s]**2
    out[~s] = -np.expm1(-2*u[~s])/(2*u[~s])
    return out


def _sprime(u):
    """e^{-u} d/du[sinh(u)/u] = (1 + e^{-2u})/(2u) - (1 - e^{-2u})/(2u^2).

    Written as given it cancels catastrophically below u ~ 0.02 (two terms of order 1/u whose difference is
    of order u), so there the series u/3 - u^2/3 + u^3/5 - 4u^4/45 + 2u^5/63 - u^6/105 is used; at the
    switch the two agree to about 1e-12.
    """
    u = np.asarray(u, float)
    out = np.empty_like(u)
    s = u < .02
    us = u[s]
    out[s] = us/3 - us**2/3 + us**3/5 - 4*us**4/45 + 2*us**5/63 - us**6/105
    ub = u[~s]
    out[~s] = (1 + np.exp(-2*ub))/(2*ub) + np.expm1(-2*ub)/(2*ub*ub)
    return out


# ---------------------------------------------------------------- the shell kernel (spherical sources)
def shell_kernel(r, rp, w):
    """(1/4 pi) int exp(-|x - x'|^2/2w^2) dOmega' for |x| = r (rows) and |x'| = rp (columns)."""
    r, rp = np.asarray(r, float)[:, None], np.asarray(rp, float)[None, :]
    return np.exp(-(r - rp)**2/(2*w*w))*_bracket(r*rp/w**2)


def shell_kernel_dr(r, rp, w):
    """d/dr of the shell kernel: E [-(r/w^2) S + (rp/w^2) S'] in scaled form."""
    r, rp = np.asarray(r, float)[:, None], np.asarray(rp, float)[None, :]
    u = r*rp/w**2
    return np.exp(-(r - rp)**2/(2*w*w))*(-(r/w**2)*_bracket(u) + (rp/w**2)*_sprime(u))


def shell_kernel_laplacian(r, rp, w):
    """Shell average of lap K_w = K_w (d^2/w^4 - 3/w^2): E [((r^2 + rp^2)/w^4 - 3/w^2) S - (2 r rp/w^4) S']."""
    r, rp = np.asarray(r, float)[:, None], np.asarray(rp, float)[None, :]
    u = r*rp/w**2
    E = np.exp(-(r - rp)**2/(2*w*w))
    return E*(((r*r + rp*rp)/w**4 - 3/w**2)*_bracket(u) - (2*r*rp/w**4)*_sprime(u))


class SphericalSource:
    """A spherical source given by its cumulative mass M(<r) on a grid; dM on the segments is the measure,
    placed at the geometric midpoints, and the unresolved mass inside r[0] sits at r[0]/2."""

    def __init__(self, r, M):
        self.r, self.M = np.asarray(r, float), np.asarray(M, float)
        self.nodes = np.concatenate([[self.r[0]/2], np.sqrt(self.r[1:]*self.r[:-1])])
        self.dM = np.concatenate([[self.M[0]], np.diff(self.M)])

    def written_potential(self, r, w):
        """(rho * K_w)(r) per unit Lambda, in Msun."""
        return shell_kernel(r, self.nodes, w)@self.dM

    def written_gradient(self, r, w):
        """d/dr (rho * K_w) per unit Lambda; the extra INWARD acceleration is -Lambda times this."""
        return shell_kernel_dr(r, self.nodes, w)@self.dM

    def written_laplacian(self, r, w):
        return shell_kernel_laplacian(r, self.nodes, w)@self.dM

    def g_mem(self, r, w):
        """Extra inward acceleration per unit Lambda."""
        return -self.written_gradient(r, w)

    def equivalent_mass(self, r, w):
        """M_eff(r) = r^2 g_mem/G per unit Lambda."""
        return np.asarray(r, float)**2*self.g_mem(r, w)/G

    def equivalent_density(self, r, w):
        """rho_eff = -lap C/(4 pi G) per unit Lambda."""
        return -self.written_laplacian(r, w)/(4*np.pi*G)

    def mass(self, r):
        return np.interp(r, self.r, self.M, left=0., right=self.M[-1])

    def g_newton(self, r):
        r = np.asarray(r, float)
        return G*self.mass(r)/r**2


# ---------------------------------------------------------------- the ring kernel (projected sources; RUT-1's)
def ring_kernel(R, Rp, w):
    """(1/2 pi) int exp(-|X - X'|^2/2w^2) dphi' in the plane: exp[-(R - Rp)^2/2w^2] I0e(R Rp/w^2)."""
    R, Rp = np.asarray(R, float)[:, None], np.asarray(Rp, float)[None, :]
    return np.exp(-(R - Rp)**2/(2*w*w))*i0e(R*Rp/w**2)


def ring_kernel_dR(R, Rp, w):
    """d/dR, using d/du [e^{-u} I0(u)] = e^{-u} [I1(u) - I0(u)]."""
    R, Rp = np.asarray(R, float)[:, None], np.asarray(Rp, float)[None, :]
    u = R*Rp/w**2
    E = np.exp(-(R - Rp)**2/(2*w*w))
    return E*(-(R - Rp)/w**2*i0e(u) + (Rp/w**2)*(i1e(u) - i0e(u)))


def ring_kernel_d2R(R, Rp, w, drop_cross_term=False):
    """d^2/dR^2 = E'' I0e + 2 E' I0e' + E I0e'', with d/du I1e = I0e - I1e/u - I1e (I1e(u)/u -> 1/2 at 0).
    drop_cross_term=True is the negative control of gate G3(d): it omits 2 E' I0e'."""
    R, Rp = np.asarray(R, float)[:, None], np.asarray(Rp, float)[None, :]
    u = R*Rp/w**2
    E = np.exp(-(R - Rp)**2/(2*w*w))
    E1 = -(R - Rp)/w**2*E
    E2 = ((R - Rp)**2/w**4 - 1/w**2)*E
    I0, I1 = i0e(u), i1e(u)
    dI0 = (Rp/w**2)*(I1 - I0)
    with np.errstate(divide='ignore', invalid='ignore'):
        I1_over_u = np.where(u > 1e-12, I1/np.where(u > 1e-12, u, 1.), .5)
    dI1 = I0 - I1_over_u - I1
    d2I0 = (Rp/w**2)**2*(dI1 - (I1 - I0))
    cross = 0. if drop_cross_term else 2*E1*dI0
    return E2*I0 + cross + E*d2I0


class ProjectedSource:
    """An axisymmetric projected source given by its cumulative projected mass M_p(<R) on a grid."""

    def __init__(self, R, Mp):
        self.R, self.Mp = np.asarray(R, float), np.asarray(Mp, float)
        self.nodes = np.concatenate([[self.R[0]/2], np.sqrt(self.R[1:]*self.R[:-1])])
        self.dM = np.concatenate([[self.Mp[0]], np.diff(self.Mp)])

    def convolved(self, R, w):
        """F(R) = (Sigma_b * K_2D)(R) per unit amplitude, in Msun."""
        return ring_kernel(R, self.nodes, w)@self.dM

    def convolved_dR(self, R, w):
        return ring_kernel_dR(R, self.nodes, w)@self.dM

    def convolved_d2R(self, R, w, drop_cross_term=False):
        return ring_kernel_d2R(R, self.nodes, w, drop_cross_term)@self.dM

    def projected_equivalent_mass(self, R, w, kernel_dR=None):
        """M_proj,eff(<R) per unit Lambda = -sqrt(2 pi) w R F'(R)/(2G). kernel_dR overrides the kernel
        derivative (the shell kernel is gate G5's negative control)."""
        R = np.asarray(R, float)
        Fp = self.convolved_dR(R, w) if kernel_dR is None else kernel_dR(R, self.nodes, w)@self.dM
        return -np.sqrt(2*np.pi)*w*R*Fp/(2*G)

    def equivalent_surface_density(self, R, w):
        """Sigma_eff(R) per unit Lambda = -(sqrt(2 pi) w/(4 pi G)) [F'' + F'/R]."""
        R = np.asarray(R, float)
        return -(np.sqrt(2*np.pi)*w/(4*np.pi*G))*(self.convolved_d2R(R, w) + self.convolved_dR(R, w)/R)

    def projected_mass(self, R):
        return np.interp(R, self.R, self.Mp, left=0., right=self.Mp[-1])


# ---------------------------------------------------------------- deflection, Einstein radius, shear
def deflection_from_g(gfun, b):
    """alpha(b) = (4/c^2) int_{-inf}^{inf} g(r) (b/r) dz, radians, with g(r) the inward acceleration (scalar in,
    scalar out); the three-dimensional route, exactly RPG-1's quadrature."""
    return 4/C_KMS**2*quad(lambda t: gfun(b/np.cos(t))*b/np.cos(t), 0, np.pi/2,
                           limit=400, epsabs=0, epsrel=1e-10)[0]


def deflection_from_projected_mass(Mp, b):
    """alpha(b) = 4 G M_proj(<b)/(c^2 b), radians; the two-dimensional route."""
    return 4*G*Mp/(C_KMS**2*b)


def einstein_radius(alpha_fun, Dl, ratio, bmax=3000.):
    """theta_E in arcsec solving b/Dl = (Dls/Ds) alpha(b); 0 for a subcritical lens, inf if the ring lies
    beyond bmax."""
    f = lambda b: ratio*alpha_fun(b) - b/Dl
    if f(1e-3) <= 0:
        return 0.
    if f(bmax) > 0:
        return np.inf
    return brentq(f, 1e-3, bmax, xtol=1e-12)/Dl*ARCSEC


def delta_sigma(Mp, Sigma, R):
    """Sigma_bar(<R) - Sigma(R) = Sigma_crit gamma_t for a circular lens."""
    return Mp/(np.pi*np.asarray(R, float)**2) - Sigma


# ---------------------------------------------------------------- analytic sources
def plummer_mass(M, a, r):
    r = np.asarray(r, float)
    return M*r**3/(r*r + a*a)**1.5


def plummer_projected_mass(M, a, R):
    R = np.asarray(R, float)
    return M*R*R/(R*R + a*a)


def plummer_density(M, a, r):
    r = np.asarray(r, float)
    return 3*M/(4*np.pi*a**3)*(1 + r*r/(a*a))**-2.5


def beta_model_mass(r, rho0, rc, beta, r_trunc=np.inf):
    """Cumulative mass of rho0 [1 + (r/rc)^2]^(-3 beta/2), truncated at r_trunc, on the grid r."""
    r = np.asarray(r, float)
    rho = np.where(r <= r_trunc, rho0*(1 + (r/rc)**2)**(-1.5*beta), 0.)
    return cumulative_trapezoid(4*np.pi*r*r*rho, r, initial=0.) + 4/3*np.pi*r[0]**3*rho[0], rho


def project_density(r, rho, R, n_z=1500):
    """Sigma(R) = 2 int_0^inf rho(sqrt(R^2 + z^2)) dz for a density given on the grid r (zero beyond it)."""
    z = np.concatenate([[0.], np.geomspace(1e-3*r[0], 1.5*r[-1], n_z)])
    s = np.sqrt(np.asarray(R, float)[:, None]**2 + z[None, :]**2)
    return 2*np.trapezoid(np.interp(s, r, rho, right=0.), z, axis=1)


def sersic_projected_fraction(components, fractions, Dl, R):
    """Projected cumulative light fraction of Sersic components: the regularised lower incomplete gamma
    P(2n, b_n (R/R_e)^{1/n}), weighted by each component's light fraction."""
    R = np.asarray(R, float)
    out = np.zeros_like(R)
    for c, p in zip(components, fractions):
        Re = c['R_arcsec']*Dl/ARCSEC
        out += p*gammainc(2*c['n'], c['bn']*(R/Re)**(1/c['n']))
    return out
