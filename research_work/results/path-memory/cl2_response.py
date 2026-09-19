"""CL-2's response family: a universal spectrum of footprint widths, summed from CL-1's per-width operators
(protocol-cl2.md, section 2), and the nonnegative least-squares solve over blocks of observations.

    C(x) = sum_j Lambda_j int rho_b(y) exp[-|x-y|^2/(2 w_j^2)] d^3y,   Lambda_j >= 0,

with the unnormalised Gaussian of RUT-1 and CL-1, so Lambda_j w_j/G is the dimensionless writing strength of
member j. Two analytic members are known before any data: a log-flat weight gives a 1/r force between its
cut-offs, and a 1/w weight gives the Newtonian Green's function exactly, so a rescaled G is inside the family.

A Block holds an operator matrix A (observations x widths, per unit Lambda), the target y (observed minus the
Newtonian baseline), the errors s, and optional nuisance columns of its own (a cluster's boundary pressure).
The solve stacks blocks with weight 1/N_block each, so a large block cannot buy weight by its size.
"""
import numpy as np
from scipy.optimize import nnls

WIDTHS = np.geomspace(.1, 1e4, 26)                    # kpc, 5 per decade
WIDTH_GRIDS = dict(primary=WIDTHS, four_per_decade=np.geomspace(.1, 1e4, 21), six_per_decade=np.geomspace(.1, 1e4, 31),
                   wide=np.geomspace(.03, 3e4, 31))
CLUSTER_MIN_WIDTH = 10.                               # kpc: narrower footprints are not grid-converged in the clusters
GALAXY_MIN_WIDTH = .15                                # kpc: narrower footprints read the photometry's kinks (amendment 2)


class Block:
    def __init__(self, name, A, y, s, nuisance=None, meta=None):
        self.name = name
        self.A = np.asarray(A, float)
        self.y = np.asarray(y, float)
        self.s = np.asarray(s, float)
        self.nuisance = None if nuisance is None else np.asarray(nuisance, float)
        self.meta = meta or {}
        self.N = len(self.y)

    @property
    def n_nuisance(self):
        return 0 if self.nuisance is None else self.nuisance.shape[1]

    def predict(self, x, nu=None):
        out = self.A@x
        if self.nuisance is not None and nu is not None:
            out = out + self.nuisance@nu
        return out

    def chi2(self, x, nu=None):
        return float(np.sum(((self.predict(x, nu) - self.y)/self.s)**2))


def solve(blocks, weights=None, maxiter=100000):
    """Nonnegative least squares over the stacked blocks; returns amplitudes, per-block nuisances, chi2 per block
    and the Karush-Kuhn-Tucker residuals."""
    weights = [1./b.N for b in blocks] if weights is None else list(weights)
    nw = blocks[0].A.shape[1]
    n_nu = [b.n_nuisance for b in blocks]
    cols = nw + sum(n_nu)
    rows, rhs = [], []
    off = nw
    for b, wt in zip(blocks, weights):
        M = np.zeros((b.N, cols))
        M[:, :nw] = b.A
        if b.n_nuisance:
            M[:, off:off + b.n_nuisance] = b.nuisance
            off += b.n_nuisance
        rows.append(M*np.sqrt(wt)/b.s[:, None])
        rhs.append(b.y*np.sqrt(wt)/b.s)
    A = np.vstack(rows)
    y = np.concatenate(rhs)
    x, _ = nnls(A, y, maxiter=maxiter)
    grad = A.T@(A@x - y)
    scale = float(np.max(np.abs(A.T@y))) or 1.
    active = x > 0
    kkt = dict(max_abs_gradient_active=float(np.max(np.abs(grad[active]))/scale) if active.any() else 0.,
               min_gradient_inactive=float(np.min(grad[~active])/scale) if (~active).any() else 0., scale=scale)
    amp = x[:nw]
    nus, off = [], nw
    for k in n_nu:
        nus.append(x[off:off + k] if k else None)
        off += k
    chi2 = [b.chi2(amp, nu) for b, nu in zip(blocks, nus)]
    return dict(amplitudes=amp, nuisances=nus, chi2=chi2, kkt=kkt, widths_active=[int(i) for i in np.flatnonzero(active[:nw])])


def spectrum_summary(widths, amplitudes, G):
    return [dict(w_kpc=float(w), Lambda=float(a), w_over_ell=float(w*a/G)) for w, a in zip(widths, amplitudes) if a > 0]


# ---------------------------------------------------------------- analytic members (gate G1)
def log_flat_force(M, r, w_min, w_max):
    """g(r) for dLambda = d ln w between w_min and w_max, per unit A: (M/r)[e^{-r^2/2wmax^2} - e^{-r^2/2wmin^2}]."""
    r = np.asarray(r, float)
    return (M/r)*(np.exp(-r**2/(2*w_max**2)) - np.exp(-r**2/(2*w_min**2)))


def inverse_width_force(M, r, prefactor=np.sqrt(np.pi/2)):
    """g(r) for dLambda = w^{-1} d ln w over all widths, per unit A: sqrt(pi/2) M/r^2 -- the Newtonian member.
    prefactor=sqrt(pi) is gate G1's negative control."""
    r = np.asarray(r, float)
    return prefactor*M/r**2


def point_source_sum(source, r, widths, weights):
    """Sum of per-width forces of a SphericalSource at radii r with the given weights (already including d ln w)."""
    return sum(wt*source.g_mem(r, w) for w, wt in zip(widths, weights))
