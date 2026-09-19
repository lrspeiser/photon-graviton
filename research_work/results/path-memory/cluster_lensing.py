"""CL-F1: spherical cluster optics for the project's written field, not a halo.

Units: kpc, Msun, km/s; physical and reduced deflections in radians.
L1 is a conditional coupling, not a derived photon interaction. Sources end
at an explicit outer radius. Signed field densities are never clipped.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

import steady_field as SF


@dataclass(frozen=True)
class Geometry:
    """Explicit adopted distances; no redshift-to-distance cosmology is inferred."""

    lens_kpc: float
    source_kpc: float
    lens_source_kpc: float
    label: str

    def __post_init__(self):
        ds = (self.lens_kpc, self.source_kpc, self.lens_source_kpc)
        if not all(np.isfinite(x) and x > 0 for x in ds):
            raise ValueError('Distances must be finite and positive.')
        if self.source_kpc <= self.lens_kpc or self.lens_source_kpc > self.source_kpc:
            raise ValueError('Require a background source and D_ls <= D_s.')
        if not self.label:
            raise ValueError('Geometry needs a provenance/scenario label.')

    @property
    def ratio(self):
        return self.lens_source_kpc / self.source_kpc

    @property
    def sigma_crit(self):
        return SF.C_KMS**2 / (4 * np.pi * SF.G * self.lens_kpc * self.ratio)


class BaryonShells:
    """Mass-conserving uniform-volume shells from cumulative ordinary mass.

    r is strictly positive and increasing; the missing inner edge is zero.
    No matter is continued outside r[-1]. Kernels use the inherited midpoint
    quadrature; projection itself is exact for these finite shells.
    """

    def __init__(self, r, cumulative_mass):
        self.r, self.mass = np.asarray(r, float), np.asarray(cumulative_mass, float)
        if (self.r.ndim != 1 or len(self.r) == 0 or self.mass.shape != self.r.shape
                or not np.all(np.isfinite(self.r)) or not np.all(np.isfinite(self.mass))
                or self.r[0] <= 0 or np.any(np.diff(self.r) <= 0)
                or self.mass[0] < 0 or np.any(np.diff(self.mass) < 0)):
            raise ValueError('Need ordered positive radii and nondecreasing finite baryon mass.')
        self.inner = np.r_[0., self.r[:-1]]
        self.dm = np.diff(np.r_[0., self.mass])
        self.rho = self.dm / ((4 * np.pi / 3) * (self.r**3 - self.inner**3))

    def project(self, R):
        R = _radii(R)
        masses, sigmas = [], []
        # Chunking bounds memory for dense exported or ray-tracing grids.
        for part in np.array_split(R, max(1, int(np.ceil(len(R) / 256)))):
            b = part[:, None]
            outer_z = np.sqrt(np.maximum(self.r**2 - b**2, 0.))
            inner_z = np.sqrt(np.maximum(self.inner**2 - b**2, 0.))
            # Rationalized difference avoids subtracting nearly equal roots.
            dz = np.divide(self.r**2 - np.maximum(self.inner, b)**2,
                           outer_z + inner_z, out=np.zeros_like(outer_z),
                           where=(outer_z + inner_z) > 0)
            sigmas.extend(2 * np.maximum(dz, 0.) @ self.rho)

            def volume(a):
                safe_a = np.where(a > 0, a, 1.)
                x = np.minimum((b / safe_a)**2, 1.)
                with np.errstate(divide='ignore'):
                    return a**3 * (-np.expm1(1.5 * np.log1p(-x)))

            masses.extend((4 * np.pi / 3) * (volume(self.r) - volume(self.inner)) @ self.rho)
        return np.asarray(masses), np.asarray(sigmas)

    def enclosed(self, r):
        r = np.asarray(r, float)
        j = np.minimum(np.searchsorted(self.r, r), len(self.r) - 1)
        base = np.r_[0., self.mass[:-1]][j]
        vol = np.maximum(np.minimum(r, self.r[j])**3 - self.inner[j]**3, 0.)
        return base + (4 * np.pi / 3) * self.rho[j] * vol


def _radii(R):
    R = np.atleast_1d(np.asarray(R, float))
    if R.ndim != 1 or not len(R) or np.any(~np.isfinite(R)) or np.any(R <= 0):
        raise ValueError('Evaluation radii must be finite and positive.')
    return R


class WrittenClusterLens:
    def __init__(self, baryons, widths=(), amplitudes=(), n_project=4000):
        self.baryons = baryons
        self.widths, self.amplitudes = np.asarray(widths, float), np.asarray(amplitudes, float)
        if (self.widths.ndim != 1 or self.widths.shape != self.amplitudes.shape
                or np.any(~np.isfinite(self.widths)) or np.any(self.widths <= 0)
                or np.any(~np.isfinite(self.amplitudes)) or np.any(self.amplitudes < 0)):
            raise ValueError('Need matching positive widths and nonnegative finite amplitudes.')
        if not isinstance(n_project, (int, np.integer)) or n_project < 32:
            raise ValueError('n_project must be an integer >= 32.')
        grid = np.geomspace(min(baryons.r[0], baryons.r[-1] * 1e-6), baryons.r[-1], n_project)
        mass, _ = baryons.project(grid)
        self.projected = SF.ProjectedSource(grid, mass)

    def profile(self, R):
        R = _radii(R)
        mb, sb = self.baryons.project(R)
        mm, sm = np.zeros_like(R), np.zeros_like(R)
        for w, a in zip(self.widths, self.amplitudes):
            if a:
                for idx in np.array_split(np.arange(len(R)), max(1, int(np.ceil(len(R) / 256)))):
                    mm[idx] += a * self.projected.projected_equivalent_mass(R[idx], w)
                    sm[idx] += a * self.projected.equivalent_surface_density(R[idx], w)
        mass, sigma = mb + mm, sb + sm
        return dict(radius_kpc=R, baryon_projected_mass_Msun=mb,
                    field_equivalent_projected_mass_Msun=mm,
                    baryon_surface_density_Msun_kpc2=sb,
                    field_equivalent_surface_density_Msun_kpc2=sm,
                    sigma_Msun_kpc2=sigma, delta_sigma_Msun_kpc2=mass / (np.pi * R**2) - sigma,
                    alpha_rad=SF.deflection_from_projected_mass(mass, R))

    def observables(self, R, geometry):
        return observables(self.profile(R), geometry)


def observables(profile, geometry):
    R = profile['radius_kpc']
    alpha = profile['alpha_rad']
    kappa = profile['sigma_Msun_kpc2'] / geometry.sigma_crit
    gamma = profile['delta_sigma_Msun_kpc2'] / geometry.sigma_crit
    kbar = kappa + gamma
    lt, lr = 1 - kbar, 1 - kappa + gamma
    with np.errstate(divide='ignore', invalid='ignore'):
        reduced = gamma / (1 - kappa)
        magnification = 1 / (lt * lr)
    return dict(profile, theta_arcsec=R / geometry.lens_kpc * SF.ARCSEC,
                beta_arcsec=(R / geometry.lens_kpc - geometry.ratio * alpha) * SF.ARCSEC,
                kappa=kappa, gamma_t=gamma, reduced_shear_t=reduced,
                lambda_t=lt, lambda_r=lr, signed_magnification=magnification)


def bracketed_roots(fun, grid, values=None):
    """All sampled sign changes, not a certificate against tangent/missed roots."""
    grid = np.asarray(grid, float)
    if grid.ndim != 1 or len(grid) < 2 or np.any(~np.isfinite(grid)) or np.any(np.diff(grid) <= 0):
        raise ValueError('Root grid must be finite, increasing and have >= 2 points.')
    y = np.asarray(fun(grid) if values is None else values, float)
    if y.shape != grid.shape or np.any(~np.isfinite(y)):
        raise ValueError('Root scan values must be finite and match the grid.')
    roots = list(grid[y == 0])
    for i in np.flatnonzero(np.signbit(y[:-1]) != np.signbit(y[1:])):
        roots.append(brentq(lambda x: float(np.asarray(fun(np.array([x])))[0]),
                            grid[i], grid[i + 1], xtol=1e-11, rtol=1e-13))
    return sorted(set(float(x) for x in roots))


def solve_lens(lens, geometry, grid, beta_arcsec=5., sampled=None):
    grid = _radii(grid)
    if not np.isfinite(beta_arcsec):
        raise ValueError('Source offset must be finite.')
    p = lens.observables(grid, geometry) if sampled is None else sampled
    critical = {}
    for key in ('lambda_t', 'lambda_r'):
        roots = bracketed_roots(lambda r: lens.observables(r, geometry)[key], grid, p[key])
        critical[key] = [dict(radius_kpc=r, theta_arcsec=r / geometry.lens_kpc * SF.ARCSEC)
                         for r in roots]
    images = []
    if beta_arcsec != 0:
        def equation(x):
            return np.sign(x) * lens.observables(np.abs(x), geometry)['beta_arcsec'] - beta_arcsec
        for sign in (-1, 1):
            xs = sign * grid if sign > 0 else -grid[::-1]
            ys = sign * p['beta_arcsec'] - beta_arcsec
            if sign < 0:
                ys = ys[::-1]
            for x in bracketed_roots(equation, xs, ys):
                obs = lens.observables([abs(x)], geometry)
                det = float(obs['lambda_t'][0] * obs['lambda_r'][0])
                images.append(dict(theta_arcsec=x / geometry.lens_kpc * SF.ARCSEC,
                                   residual_arcsec=float(equation(np.array([x]))[0]),
                                   parity=int(np.sign(det)), signed_magnification=float(1 / det) if det else None))
    return dict(source_beta_arcsec=float(beta_arcsec), images=sorted(images, key=lambda x: x['theta_arcsec']),
                critical_curves=critical, aligned_source_rings_only=bool(beta_arcsec == 0),
                search_radius_kpc=[float(grid[0]), float(grid[-1])],
                completeness='Sign-bracketed roots only; tangent roots and radii outside this range are not certified.')
