"""Fixed-charge RK sampling and lattice-Coulomb fits for issue #6."""
from __future__ import annotations

import itertools
import math
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from finite_em_core import link_index
from finite_em_flux import FluxLattice

def lattice_green_function(L: int) -> NDArray[np.float64]:
    Gk = np.zeros((L, L, L), dtype=complex)
    for nx, ny, nz in itertools.product(range(L), repeat=3):
        if nx == ny == nz == 0:
            continue
        q2 = 4.0 * (
            math.sin(math.pi * nx / L) ** 2
            + math.sin(math.pi * ny / L) ** 2
            + math.sin(math.pi * nz / L) ** 2
        )
        Gk[nx, ny, nz] = 1.0 / q2
    return np.fft.ifftn(Gk).real


class ChargeWorm(FluxLattice):
    """Uniform finite-flux ensemble with a fixed +1 charge and mobile -1 charge."""

    def __init__(self, L: int, spin: int, rng: np.random.Generator):
        super().__init__(L=L, spin=spin, rng=rng)
        self.origin = (0, 0, 0)
        self.head = [0, 0, 0]
        self.head_attempts = 0
        self.head_accepts = 0

    def head_move(self) -> bool:
        mu = int(self.rng.integers(0, 3))
        sign = 1 if self.rng.random() < 0.5 else -1
        old = self.head.copy()
        link_site = old.copy()
        if sign > 0:
            link_site_tuple = tuple(link_site)
            delta = +1
            new = old.copy()
            new[mu] = (new[mu] + 1) % self.L
        else:
            link_site[mu] = (link_site[mu] - 1) % self.L
            link_site_tuple = tuple(link_site)
            delta = -1
            new = old.copy()
            new[mu] = (new[mu] - 1) % self.L
        self.head_attempts += 1
        value = int(self.E[link_site_tuple + (mu,)]) + delta
        if value < -self.spin or value > self.spin:
            return False
        self.E[link_site_tuple + (mu,)] = value
        self.head = new
        self.head_accepts += 1
        return True

    def expected_charge(self) -> NDArray[np.int16]:
        q = np.zeros((self.L, self.L, self.L), dtype=np.int16)
        q[self.origin] += 1
        q[tuple(self.head)] -= 1
        return q

    def worm_sweep(self, head_moves: int | None = None) -> None:
        for _ in range(3 * self.L**3):
            self.local_move()
        for _ in range(head_moves if head_moves is not None else self.L**3):
            self.head_move()
        for _ in range(max(1, self.L // 2)):
            self.winding_move()


def group_charge_histogram(
    histogram: NDArray[np.int64], green: NDArray[np.float64]
) -> list[dict[str, object]]:
    L = histogram.shape[0]
    groups: dict[float, dict[str, object]] = {}
    for x, y, z in itertools.product(range(L), repeat=3):
        # Cubic lattice Green values identify the periodic shells more accurately
        # than a naive Euclidean radius.
        key = round(float(green[x, y, z]), 12)
        group = groups.setdefault(
            key,
            {
                "green": float(green[x, y, z]),
                "count": 0,
                "degeneracy": 0,
                "minimum_image_r2": min(x, L - x) ** 2
                + min(y, L - y) ** 2
                + min(z, L - z) ** 2,
            },
        )
        group["count"] = int(group["count"]) + int(histogram[x, y, z])
        group["degeneracy"] = int(group["degeneracy"]) + 1
        group["minimum_image_r2"] = min(
            int(group["minimum_image_r2"]),
            min(x, L - x) ** 2 + min(y, L - y) ** 2 + min(z, L - z) ** 2,
        )
    return sorted(groups.values(), key=lambda row: float(row["green"]), reverse=True)


def fit_charge_potential(
    histogram: NDArray[np.int64], green: NDArray[np.float64]
) -> dict[str, object]:
    groups = group_charge_histogram(histogram, green)
    usable = [
        row
        for row in groups
        if int(row["minimum_image_r2"]) >= 2 and int(row["count"]) >= 30
    ]
    if len(usable) < 4:
        return {"fit_available": False, "shells": groups}
    total = float(np.sum(histogram))
    x = np.asarray([float(row["green"]) for row in usable])
    counts = np.asarray([float(row["count"]) for row in usable])
    degeneracy = np.asarray([float(row["degeneracy"]) for row in usable])
    per_site_probability = counts / degeneracy / total
    y = -np.log(per_site_probability)
    X = np.column_stack([np.ones(len(x)), -x])
    sw = np.sqrt(counts)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    prediction = X @ beta
    mean = np.average(y, weights=counts)
    ss_res = float(np.sum(counts * (y - prediction) ** 2))
    ss_tot = float(np.sum(counts * (y - mean) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    # Compare against a constant (screened/no-force at these distances) and a
    # linear minimum-image separation model (string tension/confinement).
    r = np.sqrt([float(row["minimum_image_r2"]) for row in usable])
    Xlin = np.column_stack([np.ones(len(r)), r])
    blin, *_ = np.linalg.lstsq(Xlin * sw[:, None], y * sw, rcond=None)
    predlin = Xlin @ blin
    rmse_coulomb = math.sqrt(ss_res / np.sum(counts))
    rmse_linear = math.sqrt(float(np.sum(counts * (y - predlin) ** 2)) / np.sum(counts))
    rmse_constant = math.sqrt(float(np.sum(counts * (y - mean) ** 2)) / np.sum(counts))

    return {
        "fit_available": True,
        "shells": groups,
        "fit_shell_count": len(usable),
        "potential_model": "V(r)=constant-A*G_lattice(r)",
        "coulomb_amplitude_A": float(beta[1]),
        "intercept": float(beta[0]),
        "weighted_r_squared": r2,
        "weighted_rmse_coulomb": rmse_coulomb,
        "weighted_rmse_linear_string": rmse_linear,
        "linear_string_intercept": float(blin[0]),
        "linear_string_slope": float(blin[1]),
        "weighted_rmse_constant": rmse_constant,
        "coulomb_preferred_to_linear_string": bool(rmse_coulomb < rmse_linear),
        "coulomb_preferred_to_constant": bool(rmse_coulomb < rmse_constant),
    }


def run_charge_ensemble(
    L: int,
    spin: int,
    thermal_sweeps: int,
    sample_sweeps: int,
    seed: int,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    lattice = ChargeWorm(L=L, spin=spin, rng=rng)
    for _ in range(thermal_sweeps):
        lattice.worm_sweep()
    histogram = np.zeros((L, L, L), dtype=np.int64)
    max_constraint = 0
    for sweep in range(sample_sweeps):
        lattice.worm_sweep()
        histogram[tuple(lattice.head)] += 1
        # Add inexpensive endpoint samples between full flux decorrelation sweeps.
        for _ in range(4 * L**3):
            lattice.head_move()
            histogram[tuple(lattice.head)] += 1
        if sweep % 20 == 0:
            max_constraint = max(
                max_constraint,
                int(np.max(np.abs(lattice.gauss() - lattice.expected_charge()))),
            )
    green = lattice_green_function(L)
    return {
        "lattice_size": L,
        "spin": spin,
        "seed": seed,
        "thermal_sweeps": thermal_sweeps,
        "sample_sweeps": sample_sweeps,
        "histogram_samples": int(np.sum(histogram)),
        "max_gauss_minus_charge_residual": max_constraint,
        "head_move_acceptance": lattice.head_accepts / max(1, lattice.head_attempts),
        "plaquette_move_acceptance": lattice.local_accepts / max(1, lattice.local_attempts),
        "potential_fit": fit_charge_potential(histogram, green),
    }
