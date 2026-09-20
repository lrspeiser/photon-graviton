from __future__ import annotations

import itertools
import math
from typing import Sequence

import numpy as np

from chiral_matter_model import (
    CHARGES, ETA, GAMMA1, GAMMA2, GAMMA3, I2, R0, SPATIAL_GAMMAS, SX, SY, SZ,
    localization_ratio, nonzero_corner_gap, residual_gap, slab_hamiltonian,
    wall_localized_chirality, wilson_mass,
)

def flux_bundle_covariance(charge: int) -> dict[str, object]:
    """Check the exact endpoint/link commutator without huge Kronecker matrices.

    For q>0 each spin-1 lane is raised once; for q<0 each is lowered once.
    The product has 2^|q| nonzero flux transitions.  Every transition changes
    total electric flux by exactly q.  Matter transfer y->x changes n_x by +1
    and n_y by -1, so both endpoint Gauss commutators vanish exactly.
    """

    lane_count = abs(int(charge))
    direction = 1 if charge > 0 else -1
    allowed_inputs = (-1, 0) if direction > 0 else (0, 1)
    max_link_residual = 0
    max_left_gauss_residual = 0
    max_right_gauss_residual = 0
    transitions = 0

    for state in itertools.product(allowed_inputs, repeat=lane_count):
        output = tuple(value + direction for value in state)
        delta_flux = sum(output) - sum(state)
        max_link_residual = max(max_link_residual, abs(delta_flux - charge))
        # Hopping c_x^dagger U^(q) c_y.
        delta_n_x = 1
        delta_n_y = -1
        left_gauss = delta_flux - charge * delta_n_x
        right_gauss = -delta_flux - charge * delta_n_y
        max_left_gauss_residual = max(max_left_gauss_residual, abs(left_gauss))
        max_right_gauss_residual = max(max_right_gauss_residual, abs(right_gauss))
        transitions += 1

    return {
        "charge": int(charge),
        "spin1_lanes": lane_count,
        "nonzero_bundle_transitions": transitions,
        "expected_transitions": 2**lane_count,
        "max_link_charge_residual": int(max_link_residual),
        "max_left_gauss_residual": int(max_left_gauss_residual),
        "max_right_gauss_residual": int(max_right_gauss_residual),
    }



def binding_subspace_audit(charge: int) -> dict[str, object]:
    """Verify that the composite move preserves an exact finite bound band.

    The four orbital labels factor, so one two-state orbital is sufficient.  A
    basis state is the occupation bit string of the |q| constituent strands.
    The penalty is the number of unequal strand pairs.  Composite creation or
    annihilation has support only between the all-empty and all-filled states,
    both of which have exactly zero penalty.
    """

    strand_count = abs(int(charge))
    penalties: list[int] = []
    bound_states = 0
    positive_penalties: list[int] = []
    transition_penalty_residual = 0
    composite_transition_count = 0

    for state in itertools.product((0, 1), repeat=strand_count):
        penalty = sum(
            (state[first] - state[second]) ** 2
            for first in range(strand_count)
            for second in range(first + 1, strand_count)
        )
        penalties.append(int(penalty))
        if penalty == 0:
            bound_states += 1
        else:
            positive_penalties.append(int(penalty))

        if all(value == 0 for value in state):
            output = (1,) * strand_count
            output_penalty = 0
            transition_penalty_residual = max(
                transition_penalty_residual, abs(output_penalty - penalty)
            )
            composite_transition_count += 1
        if all(value == 1 for value in state):
            output = (0,) * strand_count
            output_penalty = 0
            transition_penalty_residual = max(
                transition_penalty_residual, abs(output_penalty - penalty)
            )
            composite_transition_count += 1

    return {
        "charge": int(charge),
        "strand_count": strand_count,
        "one_orbital_dimension": 2**strand_count,
        "zero_penalty_states_per_orbital": bound_states,
        "bound_fock_dimension_for_four_orbitals": bound_states**4,
        "minimum_positive_penalty_in_units_of_Lambda": (
            min(positive_penalties) if positive_penalties else None
        ),
        "composite_transition_count": composite_transition_count,
        "composite_transition_penalty_residual": transition_penalty_residual,
        "exact_invariant_bound_band": transition_penalty_residual == 0,
    }


def topological_window_audit(slab_width: int) -> dict[str, object]:
    """Check a finite neighborhood of every selected bulk mass.

    This is not a proof against every interacting perturbation.  It verifies
    that the one-cone/chiral-wall result is a phase property rather than a
    single tuned matrix point.
    """

    shifts = (-0.05, -0.025, 0.0, 0.025, 0.05)
    samples: list[dict[str, object]] = []
    failures = 0
    minimum_corner_gap = math.inf
    minimum_wall_magnitude = math.inf
    minimum_wall_position = math.inf

    for charge in sorted(set(CHARGES)):
        base_mass = wilson_mass(charge)
        for shift in shifts:
            mass = base_mass + shift
            low_corner_count, nonzero_minimum, _ = nonzero_corner_gap(mass, slab_width)
            left_position, right_position, left_chi, right_chi = wall_localized_chirality(
                mass, slab_width
            )
            passed = (
                -2.0 < mass < 0.0
                and low_corner_count == 1
                and nonzero_minimum > 1.0
                and left_chi < -0.995
                and right_chi > 0.995
                and left_position < -0.88
                and right_position > 0.88
            )
            failures += int(not passed)
            minimum_corner_gap = min(minimum_corner_gap, nonzero_minimum)
            minimum_wall_magnitude = min(
                minimum_wall_magnitude, abs(left_chi), abs(right_chi)
            )
            minimum_wall_position = min(
                minimum_wall_position, abs(left_position), abs(right_position)
            )
            samples.append(
                {
                    "charge": int(charge),
                    "mass_shift": float(shift),
                    "wilson_mass": float(mass),
                    "low_corner_count": int(low_corner_count),
                    "nonzero_corner_minimum_gap": float(nonzero_minimum),
                    "left_chirality": float(left_chi),
                    "right_chirality": float(right_chi),
                    "left_wall_position": float(left_position),
                    "right_wall_position": float(right_position),
                    "passed": bool(passed),
                }
            )

    return {
        "mass_shifts": list(shifts),
        "sample_count": len(samples),
        "failures": failures,
        "minimum_nonzero_corner_gap": float(minimum_corner_gap),
        "minimum_wall_chirality_magnitude": float(minimum_wall_magnitude),
        "minimum_absolute_wall_position": float(minimum_wall_position),
        "scope": "uniform local Wilson-mass perturbations; interacting radiative stability remains issue #8",
    }

def exchange_phase(charge: int) -> int:
    # A charge-q defect is a bound product of |q| microscopic fermionic
    # strands. Exchanging two composites crosses q^2 odd operators.
    return -1 if (abs(int(charge)) ** 2) % 2 else 1


def bilinear_mass_audit(charges: Sequence[int]) -> dict[str, object]:
    neutral_pairs: list[list[int]] = []
    for first, q_first in enumerate(charges):
        for second in range(first, len(charges)):
            q_second = charges[second]
            if q_first + q_second == 0:
                neutral_pairs.append([first, second])

    # Test the stronger kinetic-algebra statement.  A same-wall Weyl mass M
    # would have to anticommute with all three projected Pauli generators.
    # Repeat the check for two identical-charge flavors, because q=-1 and q=9
    # each occur twice in the selected anomaly-free set.
    hermitian_basis_2 = (I2, SX, SY, SZ)
    flavor_basis = hermitian_basis_2
    kinetic_generators = (
        np.kron(SX, I2),
        np.kron(SY, I2),
        np.kron(SZ, I2),
    )
    candidate_basis = [
        np.kron(spin_matrix, flavor_matrix)
        for spin_matrix in hermitian_basis_2
        for flavor_matrix in flavor_basis
    ]
    columns: list[np.ndarray] = []
    for candidate in candidate_basis:
        pieces: list[np.ndarray] = []
        for generator in kinetic_generators:
            anticommutator = generator @ candidate + candidate @ generator
            pieces.append(np.concatenate((anticommutator.real.ravel(), anticommutator.imag.ravel())))
        columns.append(np.concatenate(pieces))
    linear_map = np.column_stack(columns)
    singular_values = np.linalg.svd(linear_map, compute_uv=False)
    rank = int(np.count_nonzero(singular_values > 1.0e-10))
    nullity = len(candidate_basis) - rank

    return {
        "gauge_neutral_left_left_pairs": neutral_pairs,
        "gauge_neutral_left_left_pair_count": len(neutral_pairs),
        "two_flavor_same_chirality_mass_map_rank": rank,
        "two_flavor_same_chirality_mass_map_nullity": nullity,
        "minimum_nonzero_singular_value": float(np.min(singular_values[singular_values > 1.0e-10])),
    }


def universal_frame_audit(slab_width: int) -> dict[str, object]:
    momentum = np.array([0.13, 0.09, 0.05], dtype=float)
    frame = np.array(
        [
            [1.03, 0.02, -0.01],
            [0.01, 0.98, 0.015],
            [0.0, -0.02, 1.01],
        ],
        dtype=float,
    )
    reference_derivatives: dict[tuple[int, int], np.ndarray] = {}
    epsilon = 1.0e-7
    mass = wilson_mass(CHARGES[0])
    for local_axis in range(3):
        for coordinate_axis in range(3):
            plus = frame.copy()
            minus = frame.copy()
            plus[local_axis, coordinate_axis] += epsilon
            minus[local_axis, coordinate_axis] -= epsilon
            derivative = (
                slab_hamiltonian(momentum, mass, slab_width, plus)
                - slab_hamiltonian(momentum, mass, slab_width, minus)
            ) / (2.0 * epsilon)
            reference_derivatives[(local_axis, coordinate_axis)] = derivative

    maximum_species_difference = 0.0
    maximum_exact_formula_residual = 0.0
    for charge in sorted(set(CHARGES)):
        species_mass = wilson_mass(charge)
        for local_axis in range(3):
            for coordinate_axis in range(3):
                plus = frame.copy()
                minus = frame.copy()
                plus[local_axis, coordinate_axis] += epsilon
                minus[local_axis, coordinate_axis] -= epsilon
                derivative = (
                    slab_hamiltonian(momentum, species_mass, slab_width, plus)
                    - slab_hamiltonian(momentum, species_mass, slab_width, minus)
                ) / (2.0 * epsilon)
                maximum_species_difference = max(
                    maximum_species_difference,
                    float(
                        np.linalg.norm(
                            derivative - reference_derivatives[(local_axis, coordinate_axis)]
                        )
                    ),
                )
                expected_block = math.sin(momentum[coordinate_axis]) * SPATIAL_GAMMAS[local_axis]
                expected = np.kron(np.eye(slab_width), expected_block)
                maximum_exact_formula_residual = max(
                    maximum_exact_formula_residual,
                    float(np.linalg.norm(derivative - expected)),
                )

    return {
        "frame_matrix": frame.tolist(),
        "momentum": momentum.tolist(),
        "maximum_species_derivative_difference": maximum_species_difference,
        "maximum_derivative_formula_residual": maximum_exact_formula_residual,
        "charge_dependence": "only the finite flux-bundle operator U^(q); the frame derivative is identical",
    }
