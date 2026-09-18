"""Sample the declared *discrete* f(E,L) with the correct phase-space measure.

For L = R*v_theta, R*dR*dtheta*dv_r*dv_theta = dR*dtheta*dv_r*dL.
There is no extra R in the joint node mass. Positive-v_r nodes represent both
signs and theta is uniform, giving a factor 4*pi. Radial trapezoidal weights
match WarmAnnulus.summary(); rectangular, equal L weights deliberately match
its existing surface_density(). This does not repair the continuum model,
its energy cutoffs, equilibrium discretization, or any stability inference.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Any

import numpy as np

VERSION = 'annulus-phase-space-v2'


def _count(value: int, name: str, minimum: int = 1) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return int(value)


def _vector(value: Any, name: str, length: int | None = None) -> np.ndarray:
    a = np.asarray(value, dtype=float)
    if a.ndim != 1 or not np.all(np.isfinite(a)):
        raise ValueError(f'{name} must be a finite one-dimensional array')
    if length is not None and len(a) != length:
        raise ValueError(f'{name} must have length {length}')
    return a


def radial_weights(radii: Any) -> np.ndarray:
    """Positive trapezoidal nodal weights, including half-weight endpoints."""
    r = _vector(radii, 'radii')
    if len(r) < 2 or np.any(r <= 0) or np.any(np.diff(r) <= 0):
        raise ValueError('radii must be positive and strictly increasing, with at least two nodes')
    dx = np.diff(r)
    return np.r_[dx[0]/2, (dx[:-1] + dx[1:])/2, dx[-1]/2]


@dataclass(frozen=True)
class NodeDistribution:
    """Inspectable masses on (R,L,positive-v_r) nodes, with independent sign/angle draws."""

    radii: np.ndarray
    radius_index: np.ndarray
    angular_momentum: np.ndarray
    radial_speed: np.ndarray
    energy: np.ndarray
    node_mass: np.ndarray

    def __post_init__(self) -> None:
        radial_weights(self.radii)  # validate without changing the source
        size = len(self.node_mass)
        if size == 0:
            raise ValueError('no positive probability mass in the declared domain')
        idx = np.asarray(self.radius_index)
        if idx.ndim != 1 or len(idx) != size or idx.dtype.kind not in 'iu':
            raise ValueError('radius_index must contain one integer per node')
        if np.any(idx < 0) or np.any(idx >= len(self.radii)):
            raise ValueError('radius_index is outside the radial grid')
        for name in ('angular_momentum', 'radial_speed', 'energy', 'node_mass'):
            a = _vector(getattr(self, name), name, size)
            if name == 'node_mass' and np.any(a <= 0):
                raise ValueError('node masses must be strictly positive')
            if name == 'radial_speed' and np.any(a < 0):
                raise ValueError('radial_speed contains magnitudes, not signed velocities')
        total = float(np.sum(self.node_mass))
        if not np.isfinite(total) or total <= 0:
            raise ValueError('total node mass must be positive and finite')
        # Do not retain writable views into an equilibrium or caller-owned arrays.
        for name in self.__dataclass_fields__:
            a = np.array(getattr(self, name), copy=True)
            a.setflags(write=False)
            object.__setattr__(self, name, a)

    @property
    def mass(self) -> float:
        return float(np.sum(self.node_mass))

    @property
    def probabilities(self) -> np.ndarray:
        return self.node_mass/self.mass

    @property
    def node_radius(self) -> np.ndarray:
        return self.radii[self.radius_index]

    def radial_marginal(self) -> np.ndarray:
        return np.bincount(self.radius_index, weights=self.node_mass,
                           minlength=len(self.radii))/self.mass

    def observables(self) -> dict[str, np.ndarray]:
        return {'radius': self.node_radius, 'angular_momentum': self.angular_momentum,
                'energy': self.energy, 'radial_speed_squared': self.radial_speed**2}

    def moments(self) -> dict[str, dict[str, float]]:
        p = self.probabilities
        out = {}
        for name, a in self.observables().items():
            mean = float(np.sum(p*a))
            out[name] = {'mean': mean, 'variance': float(np.sum(p*(a - mean)**2))}
        return out

    def draw(self, n: int, seed: int | None = 0) -> tuple[np.ndarray, np.ndarray]:
        """Independent equal-mass draws; no global RNG, source mutation, or orbit resets."""
        n = _count(n, 'n', minimum=0)
        rng = np.random.default_rng(seed)
        pick = rng.choice(len(self.node_mass), size=n, p=self.probabilities)
        r = self.node_radius[pick]
        ell = self.angular_momentum[pick]
        vr = self.radial_speed[pick]*rng.choice([-1., 1.], size=n)
        theta = rng.uniform(0., 2*np.pi, n)
        ct, st = np.cos(theta), np.sin(theta)
        vt = ell/r
        x = np.column_stack((r*ct, r*st))
        v = np.column_stack((vr*ct - vt*st, vr*st + vt*ct))
        return x, v


def build_distribution(annulus: Any, n_L: int = 160, n_vr: int = 96,
                       reach: float = 5.) -> NodeDistribution:
    """Build masses matching annulus.surface_density on the same adapted quadrature.

    Accepts an existing WarmAnnulus object, solved or not. Verification below
    compares to its density for its current potential; it does not certify that
    the source is a converged equilibrium. No density or force parameter changes.
    """
    n_L, n_vr = _count(n_L, 'n_L', 2), _count(n_vr, 'n_vr', 2)
    if not np.isfinite(reach) or reach <= 0:
        raise ValueError('reach must be positive and finite')
    for name in ('dL', 'dE', 'w', 'tau_keep'):
        value = float(getattr(annulus, name))
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be positive and finite')
    if not np.isfinite(annulus.L0):
        raise ValueError('L0 must be finite')
    r = _vector(annulus.r, 'radii')
    wr = radial_weights(r)
    _vector(annulus.C, 'C', len(r))
    sigma = _vector(annulus.sigma, 'sigma', len(r))
    if np.any(sigma < 0):
        raise ValueError('sigma must be nonnegative')
    L, Ec, xg, wg = annulus._grids(annulus.C, n_L, n_vr, reach)
    L, Ec = _vector(L, 'L', n_L), _vector(Ec, 'Ec', n_L)
    xg, wg = _vector(xg, 'quadrature nodes', n_vr), _vector(wg, 'quadrature weights', n_vr)
    if np.any(np.abs(xg) > 1) or np.any(wg <= 0):
        raise ValueError('positive quadrature weights and nodes in [-1,1] are required')
    dl = np.diff(L)
    if np.any(dl <= 0) or not np.allclose(dl, dl[0], rtol=1e-12, atol=1e-14*abs(dl[0])):
        raise ValueError('the archived density rule requires a uniform increasing L grid')
    phi = _vector(annulus.phi_total(r, annulus.C), 'potential', len(r))
    chunks: dict[str, list[np.ndarray]] = {k: [] for k in
                                        ('radius_index', 'angular_momentum', 'radial_speed', 'energy', 'node_mass')}
    for i, R in enumerate(r):
        emin = phi[i] + .5*(L/R)**2
        hi = 2*(Ec + reach*annulus.dE - emin)
        lo = np.maximum(2*(Ec - reach*annulus.dE - emin), 0.)
        live = np.flatnonzero(hi > lo)
        if not len(live):
            continue
        a, b = np.sqrt(lo[live]), np.sqrt(hi[live])
        vr = .5*(b-a)[:, None]*xg + .5*(a+b)[:, None]
        wvr = .5*(b-a)[:, None]*wg
        ell = np.broadcast_to(L[live, None], vr.shape)
        en = emin[live, None] + .5*vr**2
        f = np.exp(-.5*((ell-annulus.L0)/annulus.dL)**2
                   - .5*((en-Ec[live, None])/annulus.dE)**2)
        mass = 4*np.pi*wr[i]*dl[0]*wvr*f  # R canceled against dv_theta=dL/R
        if not np.all(np.isfinite(mass)) or np.any(mass < 0):
            raise ValueError('nonfinite or negative node masses')
        keep = mass > 0  # exact underflow zeros have no probability
        if not np.any(keep):
            continue
        chunks['radius_index'].append(np.full(int(keep.sum()), i, dtype=np.int64))
        for name, arr in (('angular_momentum', ell), ('radial_speed', vr),
                          ('energy', en), ('node_mass', mass)):
            chunks[name].append(arr[keep])
    if not chunks['node_mass']:
        raise ValueError('no positive probability mass in the declared domain')
    return NodeDistribution(r, **{k: np.concatenate(v) for k, v in chunks.items()})


def sample(annulus: Any, n: int, seed: int | None = 0, n_L: int = 160,
           n_vr: int = 96, reach: float = 5.) -> tuple[np.ndarray, np.ndarray]:
    """Corrected drop-in function for an already constructed historical object."""
    _count(n, 'n', minimum=0)
    return build_distribution(annulus, n_L, n_vr, reach).draw(n, seed)


def verify_distribution(annulus: Any, nodes: NodeDistribution | None = None,
                        n_L: int = 160, n_vr: int = 96, reach: float = 5.,
                        tolerance: float = 1e-12) -> dict[str, Any]:
    """Numerical check against a separate density evaluation, not a stability gate."""
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError('tolerance must be positive and finite')
    nodes = build_distribution(annulus, n_L, n_vr, reach) if nodes is None else nodes
    r = _vector(annulus.r, 'radii')
    if not np.array_equal(nodes.radii, r):
        raise ValueError('node distribution belongs to a different radial grid')
    density = _vector(annulus.surface_density(annulus.C, n_L=n_L, n_vr=n_vr,
                                             reach=reach, sigma=annulus.sigma), 'density', len(r))
    if np.any(density < 0):
        raise ValueError('reference density must be nonnegative')
    radial_mass = 2*np.pi*r*density*radial_weights(r)
    total = float(np.sum(radial_mass))
    if not np.isfinite(total) or total <= 0:
        raise ValueError('reference mass must be positive and finite')
    p = radial_mass/total
    marginal_error = float(np.max(np.abs(nodes.radial_marginal() - p)))
    mass_error = abs(nodes.mass/total - 1)
    node_p = nodes.probabilities
    errors = {}
    for order in (1, 2):
        ref = float(np.sum(p*r**order))
        got = float(np.sum(node_p*nodes.node_radius**order))
        errors[f'R^{order}'] = abs(got/ref - 1) if ref else abs(got)
    worst = max(marginal_error, mass_error, *errors.values())
    return {'numerical_verification_passed': bool(worst <= tolerance),
            'sampler_version': VERSION, 'tolerance': tolerance,
            'radial_marginal_max_absolute_error': marginal_error,
            'mass_relative_error': float(mass_error), 'radius_moment_relative_errors': errors,
            'node_mass': nodes.mass, 'density_integral_mass': total,
            'node_count': len(nodes.node_mass), 'moments': nodes.moments(),
            'scientific_stability': 'not_evaluated',
            'scope': 'Agreement with the declared discrete density, not continuum convergence or finite-mass proof.'}
