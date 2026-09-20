from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from functools import reduce
from math import gcd
from pathlib import Path
from typing import Sequence

import numpy as np


CHARGES: tuple[int, ...] = (-11, -5, -1, -1, 9, 9)
DEFAULT_LS = 12
R0 = 0.12
ETA = 0.04


def _pauli() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ident = np.eye(2, dtype=complex)
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sy = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    return ident, sx, sy, sz


I2, SX, SY, SZ = _pauli()
GAMMA1 = np.kron(SX, SX)
GAMMA2 = np.kron(SY, SX)
GAMMA3 = np.kron(SZ, SX)
GAMMA4 = np.kron(I2, SY)
GAMMA5 = np.kron(I2, SZ)
SPATIAL_GAMMAS = (GAMMA1, GAMMA2, GAMMA3)
GAMMAS = (*SPATIAL_GAMMAS, GAMMA4, GAMMA5)
WALL_CHIRALITY = 1.0j * GAMMA4 @ GAMMA5


@dataclass(frozen=True)
class SpeciesResult:
    charge: int
    multiplicity_index: int
    localization_ratio: float
    wilson_mass: float
    residual_gap: float
    analytic_gap: float
    analytic_relative_error: float
    fitted_velocity: float
    dispersion_max_relative_error: float
    left_wall_position: float
    right_wall_position: float
    left_chirality: float
    right_chirality: float
    nonzero_corner_minimum_gap: float
    spectral_pairing_residual: float
    finite_ls_slopes: dict[str, float]


def canonical_charge_tuple(charges: Sequence[int]) -> tuple[int, ...]:
    ordered = tuple(sorted(int(q) for q in charges))
    reversed_sign = tuple(sorted(-int(q) for q in charges))
    return min(ordered, reversed_sign)


def primitive(charges: Sequence[int]) -> bool:
    return reduce(gcd, (abs(int(q)) for q in charges)) == 1


def has_vectorlike_pair(charges: Sequence[int]) -> bool:
    values = set(int(q) for q in charges)
    return any(-q in values for q in values)


def anomaly_sums(charges: Sequence[int]) -> tuple[int, int]:
    return sum(int(q) for q in charges), sum(int(q) ** 3 for q in charges)


def search_odd_chiral_sets(max_species: int = 6, max_abs_charge: int = 11) -> dict[str, object]:
    """Exhaustively search primitive odd-charge U(1) chiral spectra.

    Conditions:
    * nonzero odd integer charges;
    * sum(q)=sum(q^3)=0;
    * no q,-q pair (no vectorlike bilinear);
    * gcd(|q|)=1;
    * spectra identified under permutation and overall charge reversal.
    """

    first_level: tuple[int, int] | None = None
    representatives: list[tuple[int, ...]] = []
    counts_by_species: dict[str, int] = {}

    for number in range(2, max_species + 1):
        values = [
            q
            for q in range(-max_abs_charge, max_abs_charge + 1)
            if q != 0 and abs(q) % 2 == 1
        ]
        found: set[tuple[int, ...]] = set()
        for charges in itertools.combinations_with_replacement(values, number):
            if anomaly_sums(charges) != (0, 0):
                continue
            if has_vectorlike_pair(charges):
                continue
            if not primitive(charges):
                continue
            found.add(canonical_charge_tuple(charges))
        counts_by_species[str(number)] = len(found)
        if found and first_level is None:
            smallest_max = min(max(abs(q) for q in spectrum) for spectrum in found)
            first_level = (number, smallest_max)
            representatives = sorted(
                spectrum
                for spectrum in found
                if max(abs(q) for q in spectrum) == smallest_max
            )

    # A second scan proves that the first solution needs |q|max=11, rather
    # than only showing that a solution exists inside the requested box.
    first_by_charge_bound: dict[str, int] = {}
    for charge_bound in range(1, max_abs_charge + 1, 2):
        values = [
            q
            for q in range(-charge_bound, charge_bound + 1)
            if q != 0 and abs(q) % 2 == 1
        ]
        found: set[tuple[int, ...]] = set()
        for charges in itertools.combinations_with_replacement(values, max_species):
            if anomaly_sums(charges) != (0, 0):
                continue
            if has_vectorlike_pair(charges) or not primitive(charges):
                continue
            found.add(canonical_charge_tuple(charges))
        first_by_charge_bound[str(charge_bound)] = len(found)

    return {
        "max_species": max_species,
        "max_abs_charge": max_abs_charge,
        "counts_by_species": counts_by_species,
        "counts_at_six_species_by_charge_bound": first_by_charge_bound,
        "first_solution_level": list(first_level) if first_level else None,
        "first_representatives": [list(item) for item in representatives],
    }


def clifford_residual() -> float:
    ident4 = np.eye(4, dtype=complex)
    residual = 0.0
    for i, gamma_i in enumerate(GAMMAS):
        for j, gamma_j in enumerate(GAMMAS):
            target = 2.0 * ident4 if i == j else np.zeros((4, 4), dtype=complex)
            residual = max(
                residual,
                float(np.linalg.norm(gamma_i @ gamma_j + gamma_j @ gamma_i - target)),
            )
    residual = max(
        residual,
        float(np.linalg.norm(WALL_CHIRALITY @ WALL_CHIRALITY - ident4)),
        float(np.linalg.norm(WALL_CHIRALITY - WALL_CHIRALITY.conj().T)),
    )
    return residual


def localization_ratio(charge: int) -> float:
    return R0 + ETA * (abs(int(charge)) - 1)


def wilson_mass(charge: int) -> float:
    # The open 4+1D Wilson slab is topological for -2 < m_5 < 0 in this
    # convention.  r_q=|1+m_5| is the wall-mode decay ratio at k=0.
    return -1.0 + localization_ratio(charge)


def slab_hamiltonian(
    momentum: Sequence[float],
    mass: float,
    slab_width: int,
    frame: np.ndarray | None = None,
) -> np.ndarray:
    momentum = np.asarray(momentum, dtype=float)
    if momentum.shape != (3,):
        raise ValueError("momentum must have exactly three components")
    if frame is None:
        frame = np.eye(3, dtype=float)
    frame = np.asarray(frame, dtype=float)
    if frame.shape != (3, 3):
        raise ValueError("frame must be a 3x3 matrix")

    kinetic = np.zeros((4, 4), dtype=complex)
    for local_axis, gamma in enumerate(SPATIAL_GAMMAS):
        coefficient = sum(
            frame[local_axis, coordinate_axis] * math.sin(momentum[coordinate_axis])
            for coordinate_axis in range(3)
        )
        kinetic += coefficient * gamma

    onsite_wilson = mass + 1.0 + sum(1.0 - math.cos(value) for value in momentum)
    onsite = kinetic + onsite_wilson * GAMMA5
    forward = (-GAMMA5 - 1.0j * GAMMA4) / 2.0

    dimension = 4 * slab_width
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for layer in range(slab_width):
        block = slice(4 * layer, 4 * (layer + 1))
        hamiltonian[block, block] = onsite
        if layer + 1 < slab_width:
            next_block = slice(4 * (layer + 1), 4 * (layer + 2))
            hamiltonian[block, next_block] = forward
            hamiltonian[next_block, block] = forward.conj().T
    return hamiltonian


def residual_gap(mass: float, slab_width: int) -> float:
    eigenvalues = np.linalg.eigvalsh(slab_hamiltonian((0.0, 0.0, 0.0), mass, slab_width))
    return float(np.min(np.abs(eigenvalues)))


def analytic_residual_gap(ratio: float, slab_width: int) -> float:
    # Exact asymptotic overlap formula for the isolated wall mode, with the
    # finite-normalization denominator retained.  It is an increasingly good
    # approximation as the two wall wavefunctions separate.
    numerator = (1.0 - ratio * ratio) * ratio**slab_width
    denominator = 1.0 - ratio ** (2 * slab_width)
    return float(numerator / denominator)


def low_energy_dispersion(mass: float, slab_width: int) -> tuple[float, float]:
    momenta = np.array([0.0, 0.006, 0.012, 0.018, 0.024, 0.030], dtype=float)
    energies: list[float] = []
    for momentum in momenta:
        eigenvalues = np.linalg.eigvalsh(
            slab_hamiltonian((float(momentum), 0.0, 0.0), mass, slab_width)
        )
        energies.append(float(np.min(np.abs(eigenvalues))))

    x = np.sin(momenta) ** 2
    y = np.asarray(energies) ** 2
    design = np.column_stack((np.ones_like(x), x))
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    intercept, slope = coefficients
    fitted = design @ coefficients
    scale = max(float(np.max(y)), 1.0e-30)
    max_relative_error = float(np.max(np.abs(fitted - y)) / scale)
    velocity = float(math.sqrt(max(slope, 0.0)))
    return velocity, max_relative_error


def wall_localized_chirality(mass: float, slab_width: int) -> tuple[float, float, float, float]:
    hamiltonian = slab_hamiltonian((0.0, 0.0, 0.0), mass, slab_width)
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    low_indices = np.argsort(np.abs(eigenvalues))[:4]
    low_vectors = eigenvectors[:, low_indices]

    wall_coordinate = np.repeat(np.linspace(-1.0, 1.0, slab_width), 4)
    wall_operator = np.diag(wall_coordinate)
    projected_wall = low_vectors.conj().T @ wall_operator @ low_vectors
    wall_values, rotation = np.linalg.eigh(projected_wall)
    localized_vectors = low_vectors @ rotation

    chirality_operator = np.kron(np.eye(slab_width, dtype=complex), WALL_CHIRALITY)
    chiralities = np.array(
        [
            float(
                np.real(
                    localized_vectors[:, index].conj()
                    @ chirality_operator
                    @ localized_vectors[:, index]
                )
            )
            for index in range(4)
        ]
    )

    left = wall_values < 0.0
    right = wall_values > 0.0
    return (
        float(np.mean(wall_values[left])),
        float(np.mean(wall_values[right])),
        float(np.mean(chiralities[left])),
        float(np.mean(chiralities[right])),
    )


def nonzero_corner_gap(mass: float, slab_width: int) -> tuple[int, float, dict[str, float]]:
    corner_gaps: dict[str, float] = {}
    low_corners = 0
    nonzero_minimum = math.inf
    for corner in itertools.product((0.0, math.pi), repeat=3):
        eigenvalues = np.linalg.eigvalsh(slab_hamiltonian(corner, mass, slab_width))
        gap = float(np.min(np.abs(eigenvalues)))
        label = "".join("1" if value else "0" for value in corner)
        corner_gaps[label] = gap
        if gap < 0.05:
            low_corners += 1
        if corner != (0.0, 0.0, 0.0):
            nonzero_minimum = min(nonzero_minimum, gap)
    return low_corners, float(nonzero_minimum), corner_gaps


def spectral_pairing_residual(mass: float, slab_width: int) -> float:
    eigenvalues = np.linalg.eigvalsh(slab_hamiltonian((0.11, 0.07, 0.03), mass, slab_width))
    return float(np.max(np.abs(eigenvalues + eigenvalues[::-1])))


def finite_width_scaling(mass: float, ratio: float) -> dict[str, float]:
    widths = np.array([6, 8, 10, 12, 14], dtype=float)
    gaps = np.array([residual_gap(mass, int(width)) for width in widths], dtype=float)
    slope, intercept = np.polyfit(widths, np.log(gaps), 1)
    return {
        "fit_log_gap_slope": float(slope),
        "predicted_log_ratio": float(math.log(ratio)),
        "slope_error": float(abs(slope - math.log(ratio))),
        "fit_intercept": float(intercept),
        "gap_L6": float(gaps[0]),
        "gap_L14": float(gaps[-1]),
    }
