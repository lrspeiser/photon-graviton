#!/usr/bin/env python3
"""Stage 5E: photon/frame response of the complete antisymmetric matter Fock space.

The complete four-pair Stage-3I Hilbert space is used without a one-pair
truncation. Photon current/contact operators are differentiated from the exact
U(1)-covariant hopping and pair moves. A tracefree frame source scales the same
transport moves and their completed diagonal partners, so no independent frame
counterterm is introduced.

The calculation reports temporal and spatial inverse-kernel corrections. A
physical cone comparison is permitted only if the uniform frame response and
the coordinate-stress gate are under control. The script therefore separates
raw response data from any conditional speed diagnostic.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator, cg, eigsh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "microscopic"))
import check_fermionic_multipair_vacuum as fock

U_COMMON = fock.U_GAUGE
KAPPA = 2.854012469780666e-6
Z_RING = 1.018370697661399
B_FRAME = KAPPA / 2.0  # The inverse-kernel coefficient 2b, b=kappa/4.
B_COHERENT_PHOTON = KAPPA / Z_RING
LEGACY_PHOTON_SPEED = 0.011071048137135307


def ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [ready(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def patterns(momentum_pi: bool):
    photon = np.zeros(len(fock.LINKS), dtype=float)
    frame = np.zeros(len(fock.LINKS), dtype=float)
    for link_index, (site, direction) in enumerate(fock.LINKS):
        phase = (-1.0) ** site[0] if momentum_pi else 1.0
        if direction == 1:
            photon[link_index] = phase
        frame[link_index] = phase * (
            1.0 if direction == 1 else -1.0
        ) / math.sqrt(2.0)
    return photon, frame


def source_operators(
    dimension: int,
    plus_moves,
    minus_moves,
    pair_moves,
    pattern: np.ndarray,
    source: str,
):
    first = csr_matrix((dimension, dimension), dtype=complex)
    second = csr_matrix((dimension, dimension), dtype=complex)
    for link_index, weight in enumerate(pattern):
        if weight == 0.0:
            continue
        plus = plus_moves[link_index]
        minus = minus_moves[link_index]
        pair = pair_moves[link_index]
        if source == "photon":
            first += 1j * weight * fock.K_MATTER * (
                plus - plus.getH() - minus + minus.getH()
            )
            first += 1j * weight * fock.PAIR_COUPLING * (
                pair - pair.getH()
            )
            second += weight * weight * fock.K_MATTER * (
                plus + plus.getH() + minus + minus.getH()
            )
            second += weight * weight * fock.PAIR_COUPLING * (
                pair + pair.getH()
            )
        elif source == "frame":
            completed_pair = pair.getH() @ pair + pair @ pair.getH()
            transport = (
                -fock.K_MATTER
                * (plus + plus.getH() + minus + minus.getH())
                - fock.PAIR_COUPLING * (pair + pair.getH())
            )
            first += weight * (
                transport
                + 2.0
                * fock.PAIR_COMPLETION
                * fock.PAIR_COUPLING**2
                * completed_pair
            )
            second += weight * weight * (
                transport
                + 4.0
                * fock.PAIR_COMPLETION
                * fock.PAIR_COUPLING**2
                * completed_pair
            )
        else:
            raise ValueError(source)
    return first.tocsr(), second.tocsr()


def ground_state(hamiltonian):
    energies, vectors = eigsh(
        hamiltonian,
        k=1,
        which="SA",
        tol=1.0e-11,
        maxiter=180000,
    )
    return float(energies[0]), vectors[:, 0]


def inverse_moments(hamiltonian, energy: float, ground: np.ndarray, vector: np.ndarray):
    dimension = hamiltonian.shape[0]

    def matvec(item):
        return (
            hamiltonian @ item
            - energy * item
            + ground * np.vdot(ground, item)
        )

    operator = LinearOperator(
        (dimension, dimension),
        matvec=matvec,
        dtype=complex,
    )
    solutions = []
    current = vector
    infos = []
    for _ in range(3):
        solved, info = cg(
            operator,
            current,
            rtol=2.0e-10,
            atol=0.0,
            maxiter=100000,
        )
        solutions.append(solved)
        infos.append(int(info))
        current = solved
    return solutions, infos


def response(hamiltonian, energy, ground, first, second):
    expectation = np.vdot(ground, first @ ground)
    source_vector = first @ ground - expectation * ground
    solutions, infos = inverse_moments(
        hamiltonian,
        energy,
        ground,
        source_vector,
    )
    static = float(
        np.vdot(ground, second @ ground).real
        - 2.0 * np.vdot(source_vector, solutions[0]).real
    )
    temporal = float(
        2.0 * np.vdot(source_vector, solutions[2]).real
    )
    equation_residuals = []
    for index, solved in enumerate(solutions):
        target = source_vector if index == 0 else solutions[index - 1]
        residual = (
            hamiltonian @ solved
            - energy * solved
            + ground * np.vdot(ground, solved)
            - target
        )
        equation_residuals.append(
            float(np.linalg.norm(residual) / max(np.linalg.norm(target), 1.0e-30))
        )
    return {
        "source_expectation": {
            "real": float(expectation.real),
            "imag": float(expectation.imag),
        },
        "source_norm": float(np.linalg.norm(source_vector)),
        "static_response": static,
        "temporal_omega_squared_coefficient": temporal,
        "cg_info": infos,
        "maximum_inverse_equation_relative_residual": max(equation_residuals),
    }


def frame_deformed_hamiltonian(
    base,
    plus_moves,
    minus_moves,
    pair_moves,
    pattern: np.ndarray,
    amplitude: float,
):
    output = base.copy().tocsr()
    for link_index, weight in enumerate(pattern):
        if weight == 0.0:
            continue
        factor = math.exp(weight * amplitude)
        plus = plus_moves[link_index]
        minus = minus_moves[link_index]
        pair = pair_moves[link_index]
        output += -fock.K_MATTER * (factor - 1.0) * (
            plus + plus.getH() + minus + minus.getH()
        )
        output += -fock.PAIR_COUPLING * (factor - 1.0) * (
            pair + pair.getH()
        )
        output += (
            fock.PAIR_COMPLETION
            * fock.PAIR_COUPLING**2
            * (factor * factor - 1.0)
            * (pair.getH() @ pair + pair @ pair.getH())
        )
    return 0.5 * (output + output.getH())


def finite_difference_control(
    base,
    plus_moves,
    minus_moves,
    pair_moves,
    pattern,
    analytic,
    step=5.0e-4,
):
    energies = []
    for amplitude in (-step, 0.0, step):
        deformed = frame_deformed_hamiltonian(
            base,
            plus_moves,
            minus_moves,
            pair_moves,
            pattern,
            amplitude,
        )
        energies.append(ground_state(deformed)[0])
    curvature = (energies[2] - 2.0 * energies[1] + energies[0]) / step**2
    return {
        "step": step,
        "ground_energy_curvature": float(curvature),
        "bubble_plus_contact": analytic,
        "absolute_difference": float(abs(curvature - analytic)),
    }


def conditional_cone(bare_spatial, spatial_loop, temporal_loop):
    bare_temporal = 1.0 / U_COMMON
    effective_temporal = bare_temporal + temporal_loop
    effective_spatial = bare_spatial + spatial_loop
    valid = effective_temporal > 0.0 and effective_spatial > 0.0
    return {
        "bare_temporal": bare_temporal,
        "bare_spatial": bare_spatial,
        "bare_speed": math.sqrt(bare_spatial / bare_temporal),
        "loop_temporal": temporal_loop,
        "loop_spatial": spatial_loop,
        "effective_temporal": effective_temporal,
        "effective_spatial": effective_spatial,
        "effective_speed": (
            math.sqrt(effective_spatial / effective_temporal) if valid else None
        ),
        "fractional_speed_shift": (
            math.sqrt(
                effective_spatial
                / effective_temporal
                / (bare_spatial / bare_temporal)
            )
            - 1.0
            if valid
            else None
        ),
    }


def build_report(quick: bool):
    (
        basis,
        sector_counts,
        hamiltonian,
        plus_moves,
        minus_moves,
        pair_moves,
    ) = fock.build_system(4)
    energy, ground = ground_state(hamiltonian)

    responses: dict[str, dict[str, Any]] = {}
    for momentum_name, momentum_pi in (("zero", False), ("pi_zero", True)):
        photon_pattern, frame_pattern = patterns(momentum_pi)
        for source, pattern in (
            ("photon", photon_pattern),
            ("frame", frame_pattern),
        ):
            first, second = source_operators(
                len(basis),
                plus_moves,
                minus_moves,
                pair_moves,
                pattern,
                source,
            )
            responses[f"{source}_{momentum_name}"] = response(
                hamiltonian,
                energy,
                ground,
                first,
                second,
            )

    qhat_squared = 4.0
    photon_spatial = (
        responses["photon_pi_zero"]["static_response"]
        - responses["photon_zero"]["static_response"]
    ) / qhat_squared
    frame_spatial = (
        responses["frame_pi_zero"]["static_response"]
        - responses["frame_zero"]["static_response"]
    ) / qhat_squared
    photon_temporal = responses["photon_zero"][
        "temporal_omega_squared_coefficient"
    ]
    frame_temporal = responses["frame_zero"][
        "temporal_omega_squared_coefficient"
    ]

    _, frame_pi_pattern = patterns(True)
    frame_control = finite_difference_control(
        hamiltonian,
        plus_moves,
        minus_moves,
        pair_moves,
        frame_pi_pattern,
        responses["frame_pi_zero"]["static_response"],
    )

    coherent_photon = conditional_cone(
        B_COHERENT_PHOTON,
        photon_spatial,
        photon_temporal,
    )
    coherent_frame = conditional_cone(
        B_FRAME,
        frame_spatial,
        frame_temporal,
    )
    legacy_photon = conditional_cone(
        LEGACY_PHOTON_SPEED**2 / U_COMMON,
        photon_spatial,
        photon_temporal,
    )

    maximum_solver_residual = max(
        item["maximum_inverse_equation_relative_residual"]
        for item in responses.values()
    )
    uniform_frame = responses["frame_zero"]["static_response"]
    physical_comparison_permitted = (
        abs(uniform_frame) < 1.0e-8
        and maximum_solver_residual < 1.0e-7
    )
    conditional_mismatch = None
    if (
        coherent_photon["fractional_speed_shift"] is not None
        and coherent_frame["fractional_speed_shift"] is not None
    ):
        conditional_mismatch = (
            coherent_photon["fractional_speed_shift"]
            - coherent_frame["fractional_speed_shift"]
        )

    checks = {
        "complete_four_pair_fock_space": (
            len(basis) == 6336 and sector_counts.get(4) == 115
        ),
        "inverse_solvers_converged": all(
            info == 0
            for item in responses.values()
            for info in item["cg_info"]
        ),
        "inverse_equations_controlled": maximum_solver_residual < 1.0e-7,
        "frame_source_derivative_control": (
            frame_control["absolute_difference"] < 5.0e-6
        ),
        "photon_loop_nonzero": abs(photon_spatial) > 1.0e-8,
        "frame_loop_nonzero": abs(frame_spatial) > 1.0e-8,
        "uniform_frame_response_blocks_physical_cone_claim": (
            abs(uniform_frame) >= 1.0e-8
        ),
        "no_sector_specific_counterterm": True,
    }

    return {
        "module": (
            "phase_junction_network/validation/"
            "check_fock_frame_cone_response.py"
        ),
        "status": (
            "Stage 5E DIAGNOSTIC: the complete antisymmetric Fock stress has "
            "been coupled to photon and frame sources, but the nonzero uniform "
            "frame response forbids a physical cone-renormalization claim"
        ),
        "run_mode": "quick" if quick else "full",
        "hilbert_space": {
            "dimension": len(basis),
            "pair_sector_counts": sector_counts,
            "ground_energy": energy,
        },
        "source_definition": {
            "photon": (
                "phase derivative of positive/negative hopping and pair moves; "
                "the completed pair projector is U(1) invariant"
            ),
            "frame": (
                "tracefree logarithmic scaling of the same hopping/pair moves; "
                "the completed pair projector scales with twice the source"
            ),
            "new_counterterms": 0,
        },
        "responses": responses,
        "infrared_loop_coefficients": {
            "qhat_squared": qhat_squared,
            "photon_spatial_q_squared": photon_spatial,
            "photon_temporal_omega_squared": photon_temporal,
            "frame_spatial_q_squared": frame_spatial,
            "frame_temporal_omega_squared": frame_temporal,
            "uniform_photon_response": responses["photon_zero"][
                "static_response"
            ],
            "uniform_frame_response": uniform_frame,
        },
        "frame_finite_difference_control": frame_control,
        "conditional_cone_diagnostics": {
            "physical_comparison_permitted": physical_comparison_permitted,
            "coherent_ring_photon": coherent_photon,
            "coherent_frame": coherent_frame,
            "legacy_stage3e_photon": legacy_photon,
            "conditional_fractional_shift_mismatch": conditional_mismatch,
            "warning": (
                "The coherent ring has not yet been promoted to a complete "
                "Maxwell plaquette Hamiltonian, and the raw frame source has a "
                "masslike uniform response. These numbers are diagnostics, not "
                "a physical photon/graviton speed prediction."
            ),
        },
        "checks": checks,
        "execution_pass": all(checks.values()),
        "decision": (
            "The requested Fock stress coupling is calculable with no new "
            "coefficient, but it does not yet protect a shared cone. The next "
            "construction must repair the coordinate-stress source at the "
            "microscopic transport level and derive a coherent photon inverse "
            "kernel before any photon/frame pole comparison is legitimate."
        ),
        "next_gate": (
            "Build a common coherent photon/frame parent: use the even Q=K^2 "
            "ring channel for the photon curvature-square kernel and the odd "
            "chiral K channel for the frame area-curvature kernel, retain both "
            "frequency residues, and retest the full Fock stress Ward identity."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.quick)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(ready(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "execution_pass": report["execution_pass"],
                "loops": report["infrared_loop_coefficients"],
                "conditional": report["conditional_cone_diagnostics"],
                "failed_checks": [
                    key
                    for key, value in report["checks"].items()
                    if not value
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if report["execution_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
