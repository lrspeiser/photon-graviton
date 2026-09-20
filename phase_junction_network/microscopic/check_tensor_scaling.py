#!/usr/bin/env python3
"""Finite-size and rotational-symmetry checks for the physical tensor branch."""
from __future__ import annotations
import argparse
import json
import numpy as np

from check_first_order_frame import fp_matrix, tt_basis


def lattice_momentum(mode, lattice_size):
    mode = np.asarray(mode, dtype=float)
    return 2.0 * np.sin(np.pi * mode / lattice_size)


def continuum_momentum(mode, lattice_size):
    mode = np.asarray(mode, dtype=float)
    return 2.0 * np.pi * mode / lattice_size


def tensor_frequency(mode, lattice_size):
    khat = lattice_momentum(mode, lattice_size)
    basis = tt_basis(khat)
    values = np.linalg.eigvalsh(basis.T @ fp_matrix(khat) @ basis)
    return float(np.sqrt(np.mean(values))), values


def power_fit(x, y):
    slope, intercept = np.polyfit(np.log(np.asarray(x)), np.log(np.asarray(y)), 1)
    return float(slope), float(np.exp(intercept))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sizes",
        default="12,16,24,32,48,64,96,128",
        help="comma-separated periodic lattice sizes",
    )
    parser.add_argument("--output")
    args = parser.parse_args()

    sizes = [int(value) for value in args.sizes.split(",")]
    axial_mode = (1, 0, 0)
    rows = []
    maximum_tt_split = 0.0
    for lattice_size in sizes:
        omega, values = tensor_frequency(axial_mode, lattice_size)
        continuum = float(np.linalg.norm(continuum_momentum(axial_mode, lattice_size)))
        exact_lattice = float(np.linalg.norm(lattice_momentum(axial_mode, lattice_size)))
        maximum_tt_split = max(
            maximum_tt_split,
            float(abs(values[1] - values[0]) / np.mean(values)),
        )
        rows.append(
            {
                "L": lattice_size,
                "omega": omega,
                "continuum_k": continuum,
                "lattice_k": exact_lattice,
                "omega_over_lattice_k": omega / exact_lattice,
                "omega_over_continuum_k": omega / continuum,
                "tt_eigenvalues": [float(x) for x in values],
            }
        )

    slope, prefactor = power_fit(
        [row["L"] for row in rows],
        [row["omega"] for row in rows],
    )

    anisotropy_rows = []
    for lattice_size in sizes:
        if lattice_size <= 8:
            continue
        omega_axis, _ = tensor_frequency((3, 0, 0), lattice_size)
        omega_mixed, _ = tensor_frequency((2, 2, 1), lattice_size)
        relative = abs(omega_axis - omega_mixed) / (
            0.5 * (omega_axis + omega_mixed)
        )
        anisotropy_rows.append(
            {
                "L": lattice_size,
                "axis_mode_omega": omega_axis,
                "mixed_mode_omega": omega_mixed,
                "relative_anisotropy": relative,
            }
        )
    anisotropy_slope, anisotropy_prefactor = power_fit(
        [row["L"] for row in anisotropy_rows],
        [row["relative_anisotropy"] for row in anisotropy_rows],
    )

    report = {
        "status": "finite-size tensor scaling verified",
        "axial_lowest_mode": rows,
        "gap_power_fit": {
            "omega_proportional_to_L_power": slope,
            "prefactor": prefactor,
            "target_power": -1.0,
        },
        "maximum_relative_tt_polarization_split": maximum_tt_split,
        "equal_continuum_norm_direction_test": {
            "modes": [[3, 0, 0], [2, 2, 1]],
            "rows": anisotropy_rows,
            "anisotropy_proportional_to_L_power": anisotropy_slope,
            "prefactor": anisotropy_prefactor,
            "target_leading_power": -2.0,
        },
        "interpretation": (
            "The physical tensor gap closes as 1/L, both polarizations remain "
            "degenerate, and cubic directional anisotropy vanishes as an "
            "O(L^-2) cutoff effect for fixed integer momentum modes."
        ),
    }

    assert abs(slope + 1.0) < 0.02
    assert maximum_tt_split < 1.0e-12
    assert abs(anisotropy_slope + 2.0) < 0.1
    assert max(abs(row["omega_over_lattice_k"] - 1.0) for row in rows) < 1.0e-12

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
