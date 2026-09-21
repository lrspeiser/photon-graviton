#!/usr/bin/env python3
"""Stage 3C: derive a collective charging mode and test one shared cone.

The construction uses the already committed elementary cube bookkeeping:

* 8 cube vertices x 18 torsion-free connection channels;
* equivalently 24 positive-oriented links x 6 symmetric frame channels.

Both counts give 144 microscopic matter-geometry hinge channels.  A connected
relative-phase lock leaves one common phase.  Identical microscopic hinge
capacitances then add, giving U_common = Delta/144 without fitting either the
photon or gravity spectra.

The script next solves the common-cone condition using the actual Stage-6B
spin-1 transverse spectrum and reruns the finite compact gravity regulator at
the same derived charging coefficient.  The plaquette diagonal ratio v/t is
kept explicit because it has not yet been derived from the shared move set.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

from check_finite_dressed_frame_hamiltonian import (
    digits,
    encode,
    resolve,
    translate_perm,
)
from check_shared_junction_move import electromagnetic_coupling, gravity_coupling
from finite_em_dynamics_impl import axial_spectrum, fit_branch


FULL_V_VALUES = (0.2, 0.4, 0.6)
FULL_PRIMES = (11, 13, 17)
QUICK_V_VALUES = (0.4,)
QUICK_PRIMES = (11, 17)


def ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [ready(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def cube_laplacian() -> np.ndarray:
    adjacency = np.zeros((8, 8), dtype=float)
    for vertex in range(8):
        for bit in range(3):
            adjacency[vertex, vertex ^ (1 << bit)] = 1.0
    return np.diag(adjacency.sum(axis=1)) - adjacency


def complete_laplacian(size: int) -> np.ndarray:
    return size * np.eye(size) - np.ones((size, size))


def cycle_laplacian(size: int) -> np.ndarray:
    adjacency = np.zeros((size, size), dtype=float)
    for index in range(size):
        adjacency[index, (index + 1) % size] = 1.0
        adjacency[(index + 1) % size, index] = 1.0
    return np.diag(adjacency.sum(axis=1)) - adjacency


def multiplicities(values: np.ndarray, decimals: int = 10) -> dict[str, int]:
    unique, counts = np.unique(np.round(values, decimals), return_counts=True)
    return {f"{value:.10g}": int(count) for value, count in zip(unique, counts)}


def collective_bundle() -> dict[str, Any]:
    vertices = 8
    connection_channels = 18
    positive_oriented_links = 24
    symmetric_frame_channels = 6
    channel_count = vertices * connection_channels
    assert channel_count == positive_oriented_links * symmetric_frame_channels

    spatial = cube_laplacian()
    internal_complete = complete_laplacian(connection_channels)
    internal_cycle = cycle_laplacian(connection_channels)

    spatial_only = np.kron(spatial, np.eye(connection_channels))
    internal_only = np.kron(np.eye(vertices), internal_complete)
    bundle_complete = spatial_only + internal_only
    bundle_cycle = spatial_only + np.kron(np.eye(vertices), internal_cycle)

    complete_eigenvalues = np.linalg.eigvalsh(bundle_complete)
    cycle_eigenvalues = np.linalg.eigvalsh(bundle_cycle)
    common = np.ones(channel_count)

    # Each hinge has microscopic capacitance C0=1/Delta.  On the locked
    # subspace theta_i=Theta, the kinetic coefficient is N*C0, so the
    # Hamiltonian charging coefficient is Delta/N.
    capacitance_in_delta_inverse_units = np.eye(channel_count)
    effective_capacitance = float(
        common @ capacitance_in_delta_inverse_units @ common
    )
    charging_suppression = 1.0 / effective_capacitance

    return {
        "structural_counts": {
            "cube_vertices": vertices,
            "torsion_free_connection_channels_per_vertex": connection_channels,
            "vertex_connection_product": channel_count,
            "positive_oriented_links_on_periodic_2_cube": positive_oriented_links,
            "symmetric_frame_channels_per_link": symmetric_frame_channels,
            "link_frame_product": positive_oriented_links * symmetric_frame_channels,
        },
        "bundle_channels": channel_count,
        "complete_internal_lock": {
            "rank": int(np.linalg.matrix_rank(bundle_complete, tol=1.0e-9)),
            "nullity": int(
                channel_count
                - np.linalg.matrix_rank(bundle_complete, tol=1.0e-9)
            ),
            "relative_lock_gap": float(complete_eigenvalues[1]),
            "eigenvalue_multiplicities": multiplicities(complete_eigenvalues),
            "uniform_mode_residual": float(
                np.linalg.norm(bundle_complete @ common)
            ),
        },
        "cycle_internal_lock_control": {
            "rank": int(np.linalg.matrix_rank(bundle_cycle, tol=1.0e-9)),
            "nullity": int(
                channel_count
                - np.linalg.matrix_rank(bundle_cycle, tol=1.0e-9)
            ),
            "relative_lock_gap": float(cycle_eigenvalues[1]),
            "uniform_mode_residual": float(
                np.linalg.norm(bundle_cycle @ common)
            ),
            "interpretation": (
                "The 1/N charging result needs connectivity, not an all-to-all "
                "internal graph; the relative-mode gap changes but the unique "
                "uniform mode does not."
            ),
        },
        "negative_controls": {
            "spatial_lock_only_nullity": int(
                channel_count
                - np.linalg.matrix_rank(spatial_only, tol=1.0e-9)
            ),
            "internal_lock_only_nullity": int(
                channel_count
                - np.linalg.matrix_rank(internal_only, tol=1.0e-9)
            ),
            "interpretation": (
                "Without both spatial and internal locking there are 18 or 8 "
                "independent charging clocks, so one shared photon/gravity "
                "normalization is lost."
            ),
        },
        "effective_capacitance_in_C0_units": effective_capacitance,
        "derived_common_charging_factor_r": charging_suppression,
        "derived_common_charging_coefficient": "U_A=U_g=Delta/144",
    }


@functools.lru_cache(maxsize=None)
def photon_fit_at(x_rounded: float, v_rounded: float) -> dict[str, Any]:
    x = float(x_rounded)
    v_over_t = float(v_rounded)
    k_a = float(electromagnetic_coupling(x)[0])
    k_g = float(gravity_coupling(x)[0])
    r = 1.0 / 144.0
    u_over_t = r / k_a
    rows = [
        axial_spectrum(
            spin=1,
            lattice_size=lattice_size,
            u_over_t=u_over_t,
            v_over_t=v_over_t,
            eigenpairs=16,
        )
        for lattice_size in (6, 8, 10)
    ]
    fit = fit_branch(rows)
    physical_photon_speed = k_a * fit["linear_speed_c_gamma"]
    harmonic_gravity_speed = math.sqrt(r * k_g)
    return {
        "x_t_over_Delta": x,
        "v_over_t": v_over_t,
        "K_A_over_Delta": k_a,
        "K_g_over_Delta": k_g,
        "u_over_t_plaquette": u_over_t,
        "photon_fit": fit,
        "photon_speed_in_Delta_units": physical_photon_speed,
        "harmonic_gravity_speed_in_Delta_units": harmonic_gravity_speed,
        "speed_difference": physical_photon_speed - harmonic_gravity_speed,
        "speed_ratio": physical_photon_speed / harmonic_gravity_speed,
        "photon_rows": [
            {
                key: row[key]
                for key in (
                    "spin",
                    "lattice_size",
                    "hilbert_dimension",
                    "lattice_momentum",
                    "single_polarization_gap",
                    "signed_momentum_split",
                )
            }
            for row in rows
        ],
    }


def common_cone_root(v_over_t: float, iterations: int) -> dict[str, Any]:
    lower, upper = 0.10, 0.25
    lower_result = photon_fit_at(round(lower, 10), round(v_over_t, 10))
    upper_result = photon_fit_at(round(upper, 10), round(v_over_t, 10))
    if lower_result["speed_difference"] * upper_result["speed_difference"] >= 0.0:
        raise RuntimeError(f"common-cone root not bracketed for v/t={v_over_t}")

    for _ in range(iterations):
        midpoint = 0.5 * (lower + upper)
        result = photon_fit_at(round(midpoint, 10), round(v_over_t, 10))
        if lower_result["speed_difference"] * result["speed_difference"] <= 0.0:
            upper = midpoint
            upper_result = result
        else:
            lower = midpoint
            lower_result = result

    midpoint = 0.5 * (lower + upper)
    result = photon_fit_at(round(midpoint, 10), round(v_over_t, 10))
    result["root_bracket_width"] = upper - lower
    return result


def gravity_clock_hamiltonian(
    prime: int,
    lattice_size: int,
    k_physical: float,
    u_physical: float,
):
    """Compact clock with canonical angle theta=2*pi*q/p.

    The shift coefficient must be U_clock=U_phys*(p/2pi)^2 so the small
    momentum expansion gives U_phys*n^2/2 rather than introducing a separate
    prime-dependent clock normalization.
    """
    u_clock = u_physical * (prime / (2.0 * np.pi)) ** 2
    dimension = prime ** (lattice_size - 1)
    states = np.arange(dimension, dtype=np.int64)
    powers, reduced, coordinates = digits(states, prime, lattice_size)
    differences = (np.roll(coordinates, -1, axis=1) - coordinates) % prime

    rows = [states]
    cols = [states]
    data = [
        k_physical
        * np.sum(1.0 - np.cos(2.0 * np.pi * differences / prime), axis=1)
        + u_clock * lattice_size
    ]
    amplitude = -0.5 * u_clock
    for index in range(lattice_size - 1):
        plus = reduced.copy()
        plus[:, index] = (plus[:, index] + 1) % prime
        minus = reduced.copy()
        minus[:, index] = (minus[:, index] - 1) % prime
        rows.extend((states, states))
        cols.extend((encode(plus, prime, powers), encode(minus, prime, powers)))
        data.extend(
            (
                np.full(dimension, amplitude),
                np.full(dimension, amplitude),
            )
        )
    rows.extend((states, states))
    cols.extend(
        (
            encode(reduced - 1, prime, powers),
            encode(reduced + 1, prime, powers),
        )
    )
    data.extend(
        (
            np.full(dimension, amplitude),
            np.full(dimension, amplitude),
        )
    )

    hamiltonian = coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(dimension, dimension),
    ).tocsr()
    return hamiltonian, reduced, powers, u_clock


def fit_gap_rows(rows: list[dict[str, Any]]) -> dict[str, float]:
    momenta = np.asarray([row["lattice_momentum"] for row in rows])
    gaps = np.asarray([row["single_branch_gap"] for row in rows])
    power, logarithmic_prefactor = np.polyfit(
        np.log(momenta), np.log(gaps), 1
    )
    design = np.column_stack((momenta, momenta**3))
    speed, cubic = np.linalg.lstsq(design, gaps, rcond=None)[0]
    predicted = design @ np.asarray((speed, cubic))
    return {
        "gap_proportional_to_lattice_momentum_power": float(power),
        "log_fit_prefactor": float(np.exp(logarithmic_prefactor)),
        "linear_speed_c_g": float(speed),
        "cubic_cutoff_coefficient": float(cubic),
        "maximum_relative_linear_plus_cubic_fit_residual": float(
            np.max(np.abs(predicted - gaps) / gaps)
        ),
    }


def gravity_row(
    prime: int,
    lattice_size: int,
    k_physical: float,
    u_physical: float,
) -> dict[str, Any]:
    hamiltonian, reduced, powers, u_clock = gravity_clock_hamiltonian(
        prime, lattice_size, k_physical, u_physical
    )
    eigenpairs = min(24, hamiltonian.shape[0] - 2)
    energies, vectors = eigsh(
        hamiltonian,
        k=eigenpairs,
        which="SA",
        tol=1.0e-10,
        maxiter=60000,
    )
    order = np.argsort(energies)
    states = resolve(
        energies[order],
        vectors[:, order],
        translate_perm(reduced, prime, powers),
        lattice_size,
    )
    ground = min(states, key=lambda state: state[0])
    positive = min(
        (
            state
            for state in states
            if state[1] == 1 and state[0] - ground[0] > 1.0e-9
        ),
        key=lambda state: state[0],
    )
    negative = min(
        (
            state
            for state in states
            if state[1] == lattice_size - 1
            and state[0] - ground[0] > 1.0e-9
        ),
        key=lambda state: state[0],
    )
    positive_gap = positive[0] - ground[0]
    negative_gap = negative[0] - ground[0]
    lattice_momentum = 2.0 * np.sin(np.pi / lattice_size)
    return {
        "prime": prime,
        "lattice_size": lattice_size,
        "hilbert_dimension": hamiltonian.shape[0],
        "U_physical_over_Delta": u_physical,
        "K_physical_over_Delta": k_physical,
        "U_clock_over_Delta": u_clock,
        "lattice_momentum": lattice_momentum,
        "single_branch_gap": 0.5 * (positive_gap + negative_gap),
        "signed_momentum_split": abs(positive_gap - negative_gap),
    }


def finite_gravity_check(
    root: dict[str, Any], primes: tuple[int, ...]
) -> list[dict[str, Any]]:
    r = 1.0 / 144.0
    target_speed = root["harmonic_gravity_speed_in_Delta_units"]
    output = []
    for prime in primes:
        rows = [
            gravity_row(
                prime,
                lattice_size,
                root["K_g_over_Delta"],
                r,
            )
            for lattice_size in (3, 4, 5)
        ]
        fit = fit_gap_rows(rows)
        output.append(
            {
                "prime": prime,
                "rows": rows,
                "fit": fit,
                "harmonic_target_speed": target_speed,
                "relative_speed_error": (
                    fit["linear_speed_c_g"] / target_speed - 1.0
                ),
            }
        )
    return output


def build_report(quick: bool) -> dict[str, Any]:
    bundle = collective_bundle()
    v_values = QUICK_V_VALUES if quick else FULL_V_VALUES
    primes = QUICK_PRIMES if quick else FULL_PRIMES
    iterations = 10 if quick else 18

    roots = [common_cone_root(value, iterations) for value in v_values]
    gravity = {
        f"v_over_t={root['v_over_t']:.1f}": finite_gravity_check(root, primes)
        for root in roots
    }

    representative = min(
        roots, key=lambda root: abs(root["v_over_t"] - 0.4)
    )
    p5_rows = [
        gravity_row(
            5,
            lattice_size,
            representative["K_g_over_Delta"],
            1.0 / 144.0,
        )
        for lattice_size in (4, 5, 6)
    ]
    p5_fit = fit_gap_rows(p5_rows)

    root_speed_errors = [abs(root["speed_ratio"] - 1.0) for root in roots]
    photon_powers = [
        root["photon_fit"]["gap_proportional_to_lattice_momentum_power"]
        for root in roots
    ]
    photon_residuals = [
        root["photon_fit"][
            "maximum_relative_linear_plus_cubic_fit_residual"
        ]
        for root in roots
    ]
    finite_gravity_entries = [
        entry for entries in gravity.values() for entry in entries
    ]
    gravity_powers = [
        entry["fit"]["gap_proportional_to_lattice_momentum_power"]
        for entry in finite_gravity_entries
    ]
    gravity_residuals = [
        entry["fit"]["maximum_relative_linear_plus_cubic_fit_residual"]
        for entry in finite_gravity_entries
    ]

    monotonic_convergence = True
    for entries in gravity.values():
        speeds = [entry["fit"]["linear_speed_c_g"] for entry in entries]
        monotonic_convergence &= all(
            right > left for left, right in zip(speeds, speeds[1:])
        )

    checks = {
        "two_independent_structural_counts_equal_144": (
            bundle["structural_counts"]["vertex_connection_product"]
            == bundle["structural_counts"]["link_frame_product"]
            == 144
        ),
        "complete_lock_has_one_uniform_mode": (
            bundle["complete_internal_lock"]["nullity"] == 1
            and bundle["complete_internal_lock"]["uniform_mode_residual"]
            < 1.0e-12
        ),
        "topology_control_keeps_same_uniform_mode": (
            bundle["cycle_internal_lock_control"]["nullity"] == 1
            and bundle["cycle_internal_lock_control"][
                "uniform_mode_residual"
            ]
            < 1.0e-12
        ),
        "partial_lock_controls_have_multiple_clocks": (
            bundle["negative_controls"]["spatial_lock_only_nullity"] == 18
            and bundle["negative_controls"]["internal_lock_only_nullity"]
            == 8
        ),
        "charging_factor_is_one_over_144": abs(
            bundle["derived_common_charging_factor_r"] - 1.0 / 144.0
        )
        < 1.0e-15,
        "common_cone_root_found_for_each_tested_v": max(root_speed_errors)
        < (2.0e-3 if quick else 2.0e-5),
        "photon_roots_remain_z1_like": min(photon_powers) > 0.94
        and max(photon_powers) < 1.16,
        "photon_fit_residuals_below_0p002": max(photon_residuals) < 0.002,
        "finite_gravity_primes_are_linear": min(gravity_powers) > 0.94
        and max(gravity_powers) < 1.16,
        "finite_gravity_fit_residuals_below_0p003": max(gravity_residuals)
        < 0.003,
        "finite_gravity_speed_moves_toward_harmonic_target": monotonic_convergence,
        "coarse_p5_clock_is_negative_control": p5_fit[
            "gap_proportional_to_lattice_momentum_power"
        ]
        < 0.5,
        "diagonal_plaquette_ratio_is_not_silently_fixed": len(
            {round(root["x_t_over_Delta"], 5) for root in roots}
        )
        == len(roots),
    }

    report = {
        "module": (
            "phase_junction_network/microscopic/"
            "check_collective_capacitance.py"
        ),
        "status": (
            "Stage 3C PASS: a structurally counted 144-channel collective "
            "junction derives r=1/144 and admits shared photon/gravity cone "
            "roots; the plaquette diagonal ratio and full microscopic "
            "embedding remain open"
        ),
        "run_mode": "quick" if quick else "full",
        "collective_bundle": bundle,
        "common_cone_roots": roots,
        "finite_gravity_at_same_r": gravity,
        "coarse_prime_negative_control": {
            "prime": 5,
            "rows": p5_rows,
            "fit": p5_fit,
            "interpretation": (
                "The five-state clock is too coarse at the derived charging "
                "ratio. Larger finite local dimension is required before the "
                "compact frame regulator approaches the linear branch."
            ),
        },
        "parameter_count": {
            "charging_factor": (
                "r=1/144 is fixed by the declared elementary-cube bundle, not "
                "by photon or gravity spectral data."
            ),
            "common_cone": (
                "For each still-undetermined v/t, the common cone fixes the "
                "remaining move ratio x=t/Delta."
            ),
            "remaining_dimensionless_input": (
                "The diagonal plaquette ratio v/t has not yet been derived "
                "from the shared virtual paths."
            ),
            "overall_scale": "Delta remains the one dimensionful scale.",
        },
        "checks": checks,
        "stage_pass": all(checks.values()),
        "claim_boundary": {
            "established": [
                "an explicit 144 by 144 connected constraint matrix with one uniform mode",
                "the same channel count from cube-vertex/connection and link/frame bookkeeping",
                "a derived common charging suppression r=1/144",
                "common-cone roots using the actual finite photon spectrum",
                "linear finite gravity spectra at the same r for primes 11, 13, and 17",
                "failure of the coarse p=5 gravity clock at the same microscopic ratio",
            ],
            "not_established": [
                "that the 144 channels are the unique microscopic bundle",
                "the microscopic derivation of the internal locking graph strength",
                "the plaquette diagonal ratio v/t",
                "a full three-dimensional photon spectrum at the new roots",
                "exact equality of finite-p gravity and photon speeds before the p limit",
                "matter and companion coefficients from the same matrix",
                "physical alpha, G, nonlinear gravity, or formal novelty",
            ],
        },
        "next_gate": (
            "Derive the plaquette diagonal/flippability coefficient v from the "
            "same virtual move graph. Then choose no spectral input: the "
            "derived v and r must fix x, after which the full 3D photon, finite "
            "gravity, matter, and companion sectors must all be rerun without "
            "sector-specific clock rescaling."
        ),
    }

    if not report["stage_pass"]:
        raise AssertionError(json.dumps(ready(checks), indent=2, sort_keys=True))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "phase_junction_network/microscopic/"
            "collective_capacitance_results.json"
        ),
    )
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    report = build_report(bool(args.quick))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(ready(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "stage_pass": report["stage_pass"],
                "charging_factor": report["collective_bundle"][
                    "derived_common_charging_factor_r"
                ],
                "roots": [
                    {
                        "v_over_t": root["v_over_t"],
                        "x_t_over_Delta": root["x_t_over_Delta"],
                        "u_over_t": root["u_over_t_plaquette"],
                        "speed_ratio": root["speed_ratio"],
                    }
                    for root in report["common_cone_roots"]
                ],
                "failed_checks": [
                    key
                    for key, passed in report["checks"].items()
                    if not passed
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
