"""RK flux, winding, transverse-correlation, and Wilson-shift sampling."""
from __future__ import annotations

import collections
import itertools
import math
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from finite_em_core import PAIRS, all_plaquettes, link_index

class FluxLattice:
    """Periodic finite electric-flux lattice with exact plaquette/winding moves."""

    def __init__(self, L: int, spin: int, rng: np.random.Generator):
        self.L = int(L)
        self.spin = int(spin)
        self.rng = rng
        self.E = np.zeros((L, L, L, 3), dtype=np.int8)
        self.local_attempts = 0
        self.local_accepts = 0
        self.loop_attempts = 0
        self.loop_accepts = 0

    def local_move(self) -> bool:
        L = self.L
        x, y, z = self.rng.integers(0, L, size=3)
        mu, nu = PAIRS[int(self.rng.integers(0, 3))]
        sign = 1 if self.rng.random() < 0.5 else -1
        base = [int(x), int(y), int(z)]
        pmu = base.copy()
        pmu[mu] = (pmu[mu] + 1) % L
        pnu = base.copy()
        pnu[nu] = (pnu[nu] + 1) % L
        refs = (
            (tuple(base), mu, +sign),
            (tuple(pmu), nu, +sign),
            (tuple(pnu), mu, -sign),
            (tuple(base), nu, -sign),
        )
        self.local_attempts += 1
        for site, axis, delta in refs:
            value = int(self.E[site + (axis,)]) + delta
            if value < -self.spin or value > self.spin:
                return False
        for site, axis, delta in refs:
            self.E[site + (axis,)] += delta
        self.local_accepts += 1
        return True

    def winding_move(self) -> bool:
        L = self.L
        mu = int(self.rng.integers(0, 3))
        sign = 1 if self.rng.random() < 0.5 else -1
        transverse = [int(self.rng.integers(0, L)) for _ in range(2)]
        coords: list[tuple[int, int, int]] = []
        for t in range(L):
            site = [0, 0, 0]
            site[mu] = t
            other = [a for a in range(3) if a != mu]
            site[other[0]] = transverse[0]
            site[other[1]] = transverse[1]
            coords.append(tuple(site))
        self.loop_attempts += 1
        if any(
            int(self.E[site + (mu,)]) + sign < -self.spin
            or int(self.E[site + (mu,)]) + sign > self.spin
            for site in coords
        ):
            return False
        for site in coords:
            self.E[site + (mu,)] += sign
        self.loop_accepts += 1
        return True

    def sweep(self, mix_winding: bool = True, loop_attempts: int | None = None) -> None:
        for _ in range(3 * self.L**3):
            self.local_move()
        if mix_winding:
            default_attempts = max(3, self.L * self.L // 2)
            for _ in range(loop_attempts if loop_attempts is not None else default_attempts):
                self.winding_move()

    def gauss(self) -> NDArray[np.int16]:
        div = np.zeros((self.L, self.L, self.L), dtype=np.int16)
        for mu in range(3):
            div += self.E[..., mu]
            div -= np.roll(self.E[..., mu], shift=1, axis=mu)
        return div

    def winding(self) -> tuple[int, int, int]:
        values = []
        for mu in range(3):
            selector = [slice(None), slice(None), slice(None), mu]
            selector[mu] = 0
            values.append(int(np.sum(self.E[tuple(selector)])))
        return tuple(values)

    def fourier_field(self, mode: tuple[int, int, int]) -> NDArray[np.complex128]:
        L = self.L
        n = np.asarray(mode, dtype=float)
        k = 2.0 * np.pi * n / L
        # np.fft.fftn uses exp(-2 pi i n.x/L) and returns link-origin fields.
        values = np.empty(3, dtype=complex)
        norm = math.sqrt(L**3)
        index = tuple(int(v) % L for v in mode)
        for mu in range(3):
            ft = np.fft.fftn(self.E[..., mu])
            values[mu] = ft[index] * np.exp(-0.5j * k[mu]) / norm
        return values

    def wilson_shift_allowed(self, a: int, b: int) -> float:
        """Average +/- closed-loop shift acceptance over a random rectangle."""
        L = self.L
        mu, nu = PAIRS[int(self.rng.integers(0, 3))]
        start = [int(v) for v in self.rng.integers(0, L, size=3)]
        changes: list[tuple[tuple[int, int, int], int, int]] = []
        pos = start.copy()
        for _ in range(a):
            changes.append((tuple(pos), mu, +1))
            pos[mu] = (pos[mu] + 1) % L
        for _ in range(b):
            changes.append((tuple(pos), nu, +1))
            pos[nu] = (pos[nu] + 1) % L
        for _ in range(a):
            pos[mu] = (pos[mu] - 1) % L
            changes.append((tuple(pos), mu, -1))
        for _ in range(b):
            pos[nu] = (pos[nu] - 1) % L
            changes.append((tuple(pos), nu, -1))

        allowed = 0
        for sign in (-1, +1):
            okay = True
            for site, axis, delta in changes:
                value = int(self.E[site + (axis,)]) + sign * delta
                if value < -self.spin or value > self.spin:
                    okay = False
                    break
            allowed += int(okay)
        return 0.5 * allowed


def transverse_projector(mode: tuple[int, int, int], L: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    k = 2.0 * np.pi * np.asarray(mode, dtype=float) / L
    q = 2.0 * np.sin(0.5 * k)
    qnorm = np.linalg.norm(q)
    if qnorm == 0:
        raise ValueError("zero mode has no transverse projector")
    qhat = q / qnorm
    return np.eye(3) - np.outer(qhat, qhat), qhat


def fit_flux_histogram(values: Sequence[int], L: int) -> dict[str, object]:
    counts = collections.Counter(int(v) for v in values)
    usable = sorted(w for w, c in counts.items() if c >= 20 and -4 <= w <= 4)
    if 0 not in usable or len(usable) < 3:
        return {
            "counts": {str(k): int(v) for k, v in sorted(counts.items())},
            "fit_available": False,
        }
    p0 = counts[0]
    x = []
    y = []
    weights = []
    for w in usable:
        x.append(float(w * w))
        y.append(-math.log(counts[w] / p0))
        weights.append(float(counts[w]))
    X = np.column_stack([np.ones(len(x)), np.asarray(x)])
    sw = np.sqrt(np.asarray(weights))
    beta, *_ = np.linalg.lstsq(X * sw[:, None], np.asarray(y) * sw, rcond=None)
    prediction = X @ beta
    mean = np.average(y, weights=weights)
    ss_res = float(np.sum(np.asarray(weights) * (np.asarray(y) - prediction) ** 2))
    ss_tot = float(np.sum(np.asarray(weights) * (np.asarray(y) - mean) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    slope = float(beta[1])
    return {
        "counts": {str(k): int(v) for k, v in sorted(counts.items())},
        "fit_available": True,
        "quadratic_slope": slope,
        "flux_stiffness_kappa": 2.0 * L * slope,
        "one_flux_free_energy": slope,
        "scaled_one_flux_free_energy_L_deltaF": L * slope,
        "weighted_r_squared": r2,
        "fit_windings": usable,
    }


def fit_wilson_loops(loop_data: dict[tuple[int, int], list[float]]) -> dict[str, object]:
    records = []
    for (a, b), samples in sorted(loop_data.items()):
        mean = float(np.mean(samples))
        # Jeffreys-sized floor prevents log(0) while retaining a rejection signal.
        floor = 0.5 / max(1, len(samples))
        records.append(
            {
                "a": a,
                "b": b,
                "area": a * b,
                "perimeter": 2 * (a + b),
                "overlap": max(mean, floor),
                "raw_overlap": mean,
                "samples": len(samples),
            }
        )
    y = -np.log([r["overlap"] for r in records])

    def fit(column: str) -> tuple[list[float], float]:
        x = np.asarray([r[column] for r in records], dtype=float)
        X = np.column_stack([np.ones(len(x)), x])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        return [float(v) for v in beta], float(np.sqrt(np.mean(resid**2)))

    perimeter_beta, perimeter_rmse = fit("perimeter")
    area_beta, area_rmse = fit("area")
    return {
        "loops": records,
        "perimeter_fit": {"intercept_slope": perimeter_beta, "rmse": perimeter_rmse},
        "area_fit": {"intercept_slope": area_beta, "rmse": area_rmse},
        "perimeter_preferred": bool(perimeter_rmse < area_rmse),
        "rmse_ratio_perimeter_over_area": perimeter_rmse / area_rmse if area_rmse else 0.0,
    }


def run_flux_ensemble(
    L: int,
    spin: int,
    thermal_sweeps: int,
    sample_sweeps: int,
    seed: int,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    lattice = FluxLattice(L=L, spin=spin, rng=rng)
    for _ in range(thermal_sweeps):
        lattice.sweep(mix_winding=True)

    modes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0)]
    correlations = {mode: np.zeros((3, 3), dtype=complex) for mode in modes}
    windings: list[int] = []
    loop_shapes = [(1, 1), (1, 2)]
    if L >= 6:
        loop_shapes += [(1, 3), (2, 2)]
    if L >= 8:
        loop_shapes += [(1, 4), (2, 3)]
    loop_data: dict[tuple[int, int], list[float]] = {shape: [] for shape in loop_shapes}
    max_gauss = 0

    for sweep in range(sample_sweeps):
        lattice.sweep(mix_winding=True)
        w = lattice.winding()
        windings.extend(w)
        for mode in modes:
            b = lattice.fourier_field(mode)
            correlations[mode] += np.outer(b, b.conj())
        for shape in loop_shapes:
            loop_data[shape].append(lattice.wilson_shift_allowed(*shape))
        if sweep % 25 == 0:
            max_gauss = max(max_gauss, int(np.max(np.abs(lattice.gauss()))))

    mode_results = []
    for mode in modes:
        C = correlations[mode] / sample_sweeps
        C = 0.5 * (C + C.conj().T)
        P, qhat = transverse_projector(mode, L)
        evals = np.linalg.eigvalsh(C.real)
        transverse_block = P @ C.real @ P
        transverse_evals = np.linalg.eigvalsh(transverse_block)
        transverse_nonzero = sorted(transverse_evals)[-2:]
        long_power = float(qhat @ C.real @ qhat)
        total_power = float(np.trace(C.real))
        mode_results.append(
            {
                "mode": list(mode),
                "lattice_q_squared": float(
                    np.sum((2.0 * np.sin(np.pi * np.asarray(mode) / L)) ** 2)
                ),
                "correlation_eigenvalues": [float(v) for v in evals],
                "longitudinal_power_fraction": long_power / total_power if total_power else 0.0,
                "two_transverse_eigenvalues": [float(v) for v in transverse_nonzero],
                "transverse_degeneracy_fractional_split": abs(
                    transverse_nonzero[1] - transverse_nonzero[0]
                )
                / max(1.0e-15, 0.5 * sum(transverse_nonzero)),
            }
        )

    flux_fit = fit_flux_histogram(windings, L)
    flux_fit["winding_variance"] = float(np.var(windings))
    flux_fit["variance_over_L"] = float(np.var(windings) / L)
    return {
        "lattice_size": L,
        "spin": spin,
        "seed": seed,
        "link_dimension": 2 * spin + 1,
        "thermal_sweeps": thermal_sweeps,
        "sample_sweeps": sample_sweeps,
        "max_gauss_residual": max_gauss,
        "local_move_acceptance": lattice.local_accepts / max(1, lattice.local_attempts),
        "winding_move_acceptance": lattice.loop_accepts / max(1, lattice.loop_attempts),
        "flux_sector_response": flux_fit,
        "transverse_correlations": mode_results,
        "wilson_shift_overlap": fit_wilson_loops(loop_data),
    }
