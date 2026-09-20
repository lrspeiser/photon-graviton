"""Shared library for NL-1 (the root spectrum) and NK-1 (shell footprints), on CL-2 stage 2's machinery.
Provenance: Gauss-Legendre quadrature and the ring/shell averages are established; the shell footprint and the
root basis are proposed here (see the protocols); PM-1's local root law is the Bekenstein toy function."""
import numpy as np
from scipy.integrate import quad
import cl2_sources as CS
import cl2_response as CR
import cl2s2_lib as L2
import steady_field as SF

G, C = SF.G, SF.C_KMS
WIDTHS = L2.WIDTHS
SHELLS = [(d0, d0*f) for d0 in (1., 2., 5., 10., 20., 50.) for f in (.1, .25, .5)]      # (d0, w) in kpc, 18 columns
_N128 = np.polynomial.legendre.leggauss(128)
_N256 = np.polynomial.legendre.leggauss(256)


# ------------------------------------------------------------------ NK-1: the shell footprint's ring and shell averages
def shell_ring(R, Rp, d0, w, nodes=_N128):
    """Ring average over phi of K(d) = exp(-(d - d0)^2/2w^2) in the plane and its R-derivative, d^2 = R^2 + Rp^2 - 2 R Rp cos phi.
    Returns (K, dK/dR) with shape (len(R), len(Rp))."""
    x, wt = nodes
    phi = np.pi*(x + 1)/2
    wphi = wt*np.pi/2
    R = np.atleast_1d(np.asarray(R, float))[:, None, None]
    Rp = np.asarray(Rp, float)[None, :, None]
    c = np.cos(phi)[None, None, :]
    d = np.sqrt(np.maximum(R*R + Rp*Rp - 2*R*Rp*c, 1e-300))
    K = np.exp(-(d - d0)**2/(2*w*w))
    dK = -(d - d0)/(w*w)*K*(R - Rp*c)/d
    return (K*wphi).sum(-1)/np.pi, (dK*wphi).sum(-1)/np.pi


def shell_sphere(r, rp, d0, w, nodes=_N128):
    """Shell average over the polar angle of K(d) for a spherical shell at rp, and its r-derivative,
    d^2 = r^2 + rp^2 - 2 r rp mu. Returns (K, dK/dr) with shape (len(r), len(rp))."""
    x, wt = nodes
    r = np.atleast_1d(np.asarray(r, float))[:, None, None]
    rp = np.asarray(rp, float)[None, :, None]
    mu = x[None, None, :]
    d = np.sqrt(np.maximum(r*r + rp*rp - 2*r*rp*mu, 1e-300))
    K = np.exp(-(d - d0)**2/(2*w*w))
    dK = -(d - d0)/(w*w)*K*(r - rp*mu)/d
    return (K*wt).sum(-1)/2, (dK*wt).sum(-1)/2


def shell_galaxy_operator(gal, R, shells=SHELLS, n_grid=3000, r_max=300., nodes=_N128):
    """Inward written force per unit amplitude at the observed radii, one column per shell footprint: the razor-thin
    disks through the ring average, the bulge through the shell average (the same construction as stage 2's Gaussian
    columns)."""
    comp = CS.BAR.sparc_components(gal)
    Rg = np.geomspace(1e-3, r_max, n_grid)
    sig = comp['sigma_star'](Rg) + (comp['sigma_gas'](Rg) if comp['sigma_gas'] is not None else 0.)
    Rm = np.sqrt(Rg[1:]*Rg[:-1])
    mass = np.pi*(Rg[1:]**2 - Rg[:-1]**2)*.5*(sig[1:] + sig[:-1])
    cols = []
    for d0, w in shells:
        _, dK = shell_ring(R, Rm, d0, w, nodes)
        col = -(dK*mass[None, :]).sum(1)
        if comp['m_bulge'] is not None:
            Mb = comp['m_bulge'](Rg)
            dM = np.diff(Mb)
            _, dKs = shell_sphere(R, Rm, d0, w, nodes)
            col = col - (dKs*dM[None, :]).sum(1)
        cols.append(col)
    return np.array(cols).T


def shell_kernel_gates():
    Rt = np.array([.5, 2., 8., 20.])
    Rg = np.geomspace(1e-2, 100., 400)
    K0, dK0 = shell_ring(Rt, Rg, 0., 2.)
    Kg, dKg = SF.ring_kernel(Rt, Rg, 2.), SF.ring_kernel_dR(Rt, Rg, 2.)
    K5, dK5 = shell_ring(Rt, Rg, 5., 2.)
    Ks0, dKs0 = shell_sphere(Rt, Rg, 0., 2.)
    Ksg, dKsg = SF.shell_kernel(Rt, Rg, 2.), SF.shell_kernel_dr(Rt, Rg, 2.)      # the library adds the axes itself
    _, dK256 = shell_ring(Rt, Rg, 5., 1., _N256)
    _, dK128 = shell_ring(Rt, Rg, 5., 1., _N128)
    return dict(K1_ring_kernel=float(np.max(np.abs(K0 - Kg))/np.max(np.abs(Kg))), K1_ring_derivative=float(np.max(np.abs(dK0 - dKg))/np.max(np.abs(dKg))),
                K1_control_d0_5=float(np.max(np.abs(K5 - Kg))/np.max(np.abs(Kg))),
                K3_shell_kernel=float(np.max(np.abs(Ks0 - Ksg))/np.max(np.abs(Ksg))), K3_shell_derivative=float(np.max(np.abs(dKs0 - dKsg))/np.max(np.abs(dKsg))),
                K2_quadrature=float(np.max(np.abs(dK256 - dK128))/np.max(np.abs(dK256))))


# ------------------------------------------------------------------ NL-1: the root basis
def root_basis(cols):
    """sqrt of the nonnegative part of each column; negative entries (outward written force) are zero for the root law."""
    cols = np.asarray(cols, float)
    neg = cols < 0
    return np.sqrt(np.where(neg, 0., cols)), dict(negative_entries=int(neg.sum()), entries=int(cols.size),
                                                   largest_negative_relative=float(np.max(np.abs(cols[neg]))/np.max(np.abs(cols))) if neg.any() else 0.)


def cluster_root_columns(cl, widths, min_width=1.):
    """Pressure columns for the root basis of one cluster: the local member sqrt(g_N) and sqrt(g_j) for widths >= min_width
    (narrower written forces are not resolved on this grid and are excluded from the root family)."""
    r, src = cl['r'], cl['source']
    gN = G*cl['M_b']/r**2
    basis = [np.sqrt(np.maximum(gN, 0))]
    labels = ['local']
    neg = 0
    for w in widths:
        if w < min_width:
            basis.append(np.zeros_like(r))
        else:
            g = src.g_mem(r, w)
            neg += int((g < 0).sum())
            basis.append(np.sqrt(np.maximum(g, 0)))
        labels.append(float(w))
    ops = np.array([cl['pressure_of'](b) for b in basis]).T
    return ops, labels, neg


class LensRoot(L2.LensRows):
    """A lens with the root basis: the local member sqrt(m starforce) and sqrt(m memforce_j) as force rows; the bend of
    each basis force by the three-dimensional route (the projected route applies to fields of the source, not to roots)."""
    def __init__(self, name, scenario=L2.LENS_GEOMETRY, widths=WIDTHS):
        super().__init__(name, scenario, widths)
        L = self.L
        r = L.model.r
        self.basis = np.vstack([np.sqrt(np.maximum(L.starforce, 0))] + [np.sqrt(np.maximum(row, 0)) for row in L.memforce])
        self.basis_bend = np.array([SF.deflection_from_g(lambda x, row=row: np.interp(x, r, row, left=row[0], right=0.), L.bE) for row in self.basis])

    def coefficients(self, beta):
        key = ('root', round(float(beta), 12))
        if key not in self.cache:
            self.L.model.forces = np.vstack([self.L.starforce, self.basis])
            self.cache[key] = self.L.model.coefficients(beta)
        return self.cache[key]

    def rows(self, imf, beta):
        m = self.L.pop[imf]/1e11
        co = self.coefficients(beta)
        A = (np.sqrt(m)*co[1:]).T                 # the extra force scales as sqrt(m): sqrt(m g) = sqrt(m) sqrt(g)
        y = self.L.y**2 - m*co[0]
        from scipy.linalg import solve_triangular
        return dict(Ak=solve_triangular(self.Lc, A, lower=True), yk=solve_triangular(self.Lc, y, lower=True),
                    ae=np.sqrt(m)*self.basis_bend, ye=self.L.need - m*self.L.starbend, se=L2.EINSTEIN_LIMIT*self.L.need, m=m, co=co, star=np.sqrt(m))

    def evaluate(self, imf, beta, amps):
        r = self.rows(imf, beta)
        v2 = r['m']*r['co'][0] + (r['star']*np.asarray(amps))@r['co'][1:]
        return dict(chi2=self.L.chi2_v(v2), einstein_residual=float((r['ae']@amps - r['ye'])/self.L.need), beta=float(beta))


# ------------------------------------------------------------------ the balanced joint objective (NL-1 and NK-1 part 2)
V_ACCEPT, CLUSTER_REF, KIN_LIMIT, N_LENSES = 21.9, 11.64, 128., 6


def balanced_weights(cluster_block, E, K):
    return [1./(cluster_block.N*2*CLUSTER_REF), 1./KIN_LIMIT, 1./N_LENSES] if E is not None else [1./(cluster_block.N*2*CLUSTER_REF)]


def galaxy_loss_balanced(block, gN_key='gN_rec'):
    loss = L2.SpeedLoss(block['A'], block[gN_key], block['R'], block['v'], block['sizes'])
    loss.w = loss.w/V_ACCEPT**2
    return loss
