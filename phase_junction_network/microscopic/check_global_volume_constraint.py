#!/usr/bin/env python3
"""Check a fixed-global-volume constraint for the periodic conformal mode."""
from __future__ import annotations

import argparse
import json

import numpy as np
from scipy.linalg import null_space

from check_constrained_local_reduction import KINETIC_METRIC, TRACE_VECTOR


def inertia(values: np.ndarray, tolerance: float = 1.0e-10):
    return {
        "negative": int(np.count_nonzero(values < -tolerance)),
        "zero": int(np.count_nonzero(np.abs(values) <= tolerance)),
        "positive": int(np.count_nonzero(values > tolerance)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    args = parser.parse_args()

    poisson = np.block(
        [
            [np.zeros((6, 6)), np.eye(6)],
            [-np.eye(6), np.zeros((6, 6))],
        ]
    )
    symplectic = -poisson
    hamiltonian = np.block(
        [
            [np.zeros((6, 6)), np.zeros((6, 6))],
            [np.zeros((6, 6)), KINETIC_METRIC],
        ]
    )

    # Primary fixed-volume condition and its secondary trace-momentum condition.
    constraints = np.zeros((2, 12), dtype=float)
    constraints[0, :6] = TRACE_VECTOR
    constraints[1, 6:] = TRACE_VECTOR

    bracket = constraints @ poisson @ constraints.T
    constraint_surface = null_space(constraints)
    reduced_symplectic = (
        constraint_surface.T @ symplectic @ constraint_surface
    )
    reduced_hamiltonian = (
        constraint_surface.T @ hamiltonian @ constraint_surface
    )

    # The time derivative of V=t.h is t.T p = -1/2 t.p.
    preservation_map = TRACE_VECTOR @ KINETIC_METRIC
    expected_preservation = -0.5 * TRACE_VECTOR

    report = {
        "status": (
            "fixed-global-volume second-class pair removes the "
            "homogeneous conformal instability"
        ),
        "unconstrained_global_sector": {
            "phase_dimension": 12,
            "hamiltonian_inertia": inertia(
                np.linalg.eigvalsh(hamiltonian)
            ),
            "kinetic_eigenvalues": [
                float(value)
                for value in np.linalg.eigvalsh(KINETIC_METRIC)
            ],
        },
        "constraints": {
            "primary": "V = t_I h_I = 0",
            "secondary": "P_V = t_I p_I = 0",
            "poisson_bracket_matrix": bracket.tolist(),
            "bracket_rank": int(np.linalg.matrix_rank(bracket)),
            "preservation_residual": float(
                np.max(np.abs(preservation_map - expected_preservation))
            ),
        },
        "reduced_global_sector": {
            "phase_dimension": int(constraint_surface.shape[1]),
            "symplectic_rank": int(
                np.linalg.matrix_rank(reduced_symplectic)
            ),
            "hamiltonian_inertia": inertia(
                np.linalg.eigvalsh(reduced_hamiltonian)
            ),
            "minimum_hamiltonian_eigenvalue": float(
                np.min(np.linalg.eigvalsh(reduced_hamiltonian))
            ),
        },
        "degree_count_on_N_site_periodic_lattice": {
            "before_fixed_volume": "4N + 8 physical phase dimensions",
            "after_fixed_volume": "4N + 6 physical phase dimensions",
            "decomposition_after_fixed_volume": (
                "4(N-1) propagating phase dimensions plus "
                "10 homogeneous trace-free shear phase dimensions"
            ),
        },
        "interpretation": (
            "Fixing the total spatial volume is not first-class by itself. "
            "Hamiltonian preservation generates the homogeneous trace-momentum "
            "constraint, and the pair is second-class. It removes exactly the "
            "negative conformal canonical pair while leaving five positive "
            "global shear momenta and five zero-potential global shape moduli."
        ),
        "limitation": (
            "This controls the periodic homogeneous ghost at linear order. "
            "It does not yet derive a microscopic fixed-volume rule, remove "
            "the remaining global shear moduli, or solve the local vacuum-energy "
            "and nonlinear cosmological-constant problem."
        ),
    }

    assert report["constraints"]["bracket_rank"] == 2
    assert report["constraints"]["preservation_residual"] < 1.0e-12
    assert report["reduced_global_sector"]["phase_dimension"] == 10
    assert report["reduced_global_sector"]["symplectic_rank"] == 10
    assert report["reduced_global_sector"]["hamiltonian_inertia"] == {
        "negative": 0,
        "zero": 5,
        "positive": 5,
    }

    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
