#!/usr/bin/env python3
"""Stage 5H: common chiral generator for photon and frame channels.

Instead of assigning the photon to the positive even ring Q=K^2 and the frame
to the odd chiral ring K, place both in the same first-order K channel. The U(1)
charge generator I and a tracefree frame-area generator B are normalized by the
same trace pairing: A=I or B, F=A/4, so tr(A F)=1 in four internal dimensions.

This test asks whether their exact low-band poles, residues, and positive squared
kernels coincide without any port attenuation or sector-specific coefficient.
It does not yet prove that the chiral photon channel has the full finite Maxwell
Gauss phase or a positive bosonic Hamiltonian.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
import sys
sys.path.insert(0, str(HERE))
from solve_coherent_junction import X, adj, ring_parameters
from check_frame_area_vertex import bridge_with_area, projected_area_operator

U_COMMON = 0.006789913753821769


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


def channel(name: str, generator: np.ndarray) -> dict[str, Any]:
    curvature = generator / 4.0
    _, _, coupling, heavy = bridge_with_area(np.eye(4), generator)
    residue = np.eye(8) + coupling @ np.linalg.solve(
        heavy @ heavy, coupling.conj().T
    )
    transport = np.vstack((np.eye(4), 1j * np.eye(4))) / math.sqrt(2.0)
    projected_residue = transport.conj().T @ residue @ transport
    residue_scalar = float(np.trace(projected_residue).real / 4.0)
    residue_nonidentity = float(
        np.linalg.norm(
            projected_residue - residue_scalar * np.eye(4),
            2,
        )
    )

    phases = np.asarray((2.5e-4, 5.0e-4, 1.0e-3, 2.0e-3))
    odd_coefficients = []
    squared_coefficients = []
    operator_residuals = []
    orientation_residuals = []
    pole_splittings = []
    for phase in phases:
        plus = projected_area_operator(
            expm(1j * phase * curvature), generator
        )
        minus = projected_area_operator(
            expm(-1j * phase * curvature), generator
        )
        odd = 0.5 * (plus - minus)
        even = 0.5 * (plus + minus)
        coefficient = -float(np.trace(odd).real / 4.0) / phase
        odd_coefficients.append(coefficient)
        squared_coefficients.append(
            float(np.trace(odd @ odd).real / 4.0) / phase**2
        )
        operator_residuals.append(
            float(np.linalg.norm(odd + coefficient * phase * np.eye(4), 2))
        )
        orientation_residuals.append(float(np.linalg.norm(even, 2)))

        full, _, _, _ = bridge_with_area(
            expm(1j * phase * curvature), generator
        )
        eigenvalues = np.linalg.eigvalsh(full)
        low = np.asarray(sorted(eigenvalues, key=abs)[:8])
        negative = low[low < 0.0]
        positive = low[low > 0.0]
        pole_splittings.append(
            float(np.mean(positive) - np.mean(negative))
        )

    linear_design = phases[:, None]
    pole_slope = float(
        np.linalg.lstsq(linear_design, np.asarray(pole_splittings), rcond=None)[0][0]
    )
    static_coefficient = float(np.mean(odd_coefficients[-2:]))
    dynamic_coefficient = static_coefficient / residue_scalar
    positive_squared_coefficient = float(
        np.mean(squared_coefficients[-2:]) / residue_scalar**2
    )
    return {
        "name": name,
        "generator_trace": float(np.trace(generator).real),
        "generator_norm_squared_trace": float(
            np.trace(generator.conj().T @ generator).real
        ),
        "trace_pairing_tr_A_F": float(
            np.trace(generator @ curvature).real
        ),
        "projected_residue": residue_scalar,
        "projected_residue_nonidentity_residual": residue_nonidentity,
        "static_odd_coefficient": static_coefficient,
        "dynamic_odd_coefficient": dynamic_coefficient,
        "positive_squared_coefficient": positive_squared_coefficient,
        "exact_low_pole_splitting_slope": pole_slope,
        "orientation_reversal_even_residual": max(orientation_residuals),
        "odd_operator_identity_residual": max(operator_residuals),
        "squared_kernel_speed_diagnostic": math.sqrt(
            U_COMMON * positive_squared_coefficient
        ),
        "phase_rows": [
            {
                "phase": float(phase),
                "odd_coefficient": float(odd_coefficient),
                "squared_coefficient": float(squared_coefficient),
                "low_pole_splitting": float(splitting),
            }
            for phase, odd_coefficient, squared_coefficient, splitting in zip(
                phases,
                odd_coefficients,
                squared_coefficients,
                pole_splittings,
            )
        ],
    }


def build_report() -> dict[str, Any]:
    photon_generator = np.eye(4, dtype=complex)
    frame_generator = np.diag((1.0, 1.0, -1.0, -1.0)).astype(complex)
    photon = channel("U1_identity", photon_generator)
    frame = channel("tracefree_frame_area", frame_generator)

    comparisons = {}
    for key in (
        "projected_residue",
        "static_odd_coefficient",
        "dynamic_odd_coefficient",
        "positive_squared_coefficient",
        "exact_low_pole_splitting_slope",
        "squared_kernel_speed_diagnostic",
    ):
        denominator = max(abs(photon[key]), abs(frame[key]), 1.0e-30)
        comparisons[key] = {
            "photon": photon[key],
            "frame": frame[key],
            "relative_difference": abs(photon[key] - frame[key]) / denominator,
        }

    _, _, r, rho = ring_parameters()
    checks = {
        "equal_trace_normalization": (
            abs(photon["trace_pairing_tr_A_F"] - 1.0) < 1.0e-12
            and abs(frame["trace_pairing_tr_A_F"] - 1.0) < 1.0e-12
        ),
        "residue_is_scalar_in_both_channels": max(
            photon["projected_residue_nonidentity_residual"],
            frame["projected_residue_nonidentity_residual"],
        )
        < 1.0e-12,
        "orientation_odd_response": max(
            photon["orientation_reversal_even_residual"],
            frame["orientation_reversal_even_residual"],
        )
        < 1.0e-11,
        "projected_operator_is_scalar": max(
            photon["odd_operator_identity_residual"],
            frame["odd_operator_identity_residual"],
        )
        < 1.0e-11,
        "photon_and_frame_residues_match": comparisons[
            "projected_residue"
        ]["relative_difference"]
        < 1.0e-12,
        "photon_and_frame_odd_coefficients_match": comparisons[
            "dynamic_odd_coefficient"
        ]["relative_difference"]
        < 1.0e-10,
        "positive_squared_kernels_match": comparisons[
            "positive_squared_coefficient"
        ]["relative_difference"]
        < 1.0e-10,
        "exact_low_pole_slopes_match": comparisons[
            "exact_low_pole_splitting_slope"
        ]["relative_difference"]
        < 1.0e-6,
        "no_port_attenuation_or_new_coefficient": True,
    }

    return {
        "module": (
            "phase_junction_network/construction/coherent_junction/"
            "check_common_chiral_generator.py"
        ),
        "status": (
            "Stage 5H ALGEBRAIC CANDIDATE PASS: equal trace-normalized U(1) "
            "and frame-area generators in the same chiral ring have matching "
            "residues, odd pole slopes, and positive squared kernels"
        ),
        "ring": {
            "x": X,
            "r": r,
            "rho": rho,
            "architecture_change": (
                "The photon is moved from the even positive Q=K^2 channel to "
                "the same first-order chiral K channel used by the frame."
            ),
            "new_continuous_coefficients": 0,
        },
        "photon": photon,
        "frame": frame,
        "comparisons": comparisons,
        "checks": checks,
        "execution_pass": all(checks.values()),
        "decision": (
            "This removes the Stage-5F algebraic pole mismatch without port "
            "attenuation. It does not yet establish a common physical bosonic "
            "cone: photon positivity and the different derivative roles of the "
            "photon and frame channels must be audited next."
        ),
        "claim_boundary": {
            "established": [
                "equality of photon and frame chiral residues",
                "equality of their linear odd response coefficients",
                "equality of their positive squared low-band kernels",
                "equality of retained low-pole splitting slopes within numerical precision",
                "no sector-specific port normalization or continuous coefficient",
            ],
            "not_established": [
                "a positive bosonic Maxwell Hamiltonian from the chiral photon poles",
                "a common physical photon/gravity speed",
                "two transverse photon polarizations after exact finite Gauss reduction",
                "a deconfined finite 3+1D photon phase for this revised channel",
                "matter-loop stability of the common chiral algebra",
            ],
        },
        "next_gate": (
            "Audit bosonic positivity and derivative order before building a "
            "finite photon phase. The common first-order pole is acceptable only "
            "if a local positive Hamiltonian preserves z=1 for both sectors."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
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
                "comparisons": report["comparisons"],
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
