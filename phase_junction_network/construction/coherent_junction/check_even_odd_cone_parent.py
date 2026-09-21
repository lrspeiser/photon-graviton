#!/usr/bin/env python3
"""Stage 5F: compare the even photon and odd frame channels of one ring parent.

The positive completed ring Q=K^2 supplies an orientation-even curvature-square
channel. The chiral dilation K supplies the orientation-odd frame area-curvature
channel. Both use the inherited x and the same physical/reference architecture.
No photon or gravity coefficient is fitted in the baseline calculation.

The calculation also scans discrete low-port normalization denominators. That
scan is diagnostic only: a denominator is not accepted unless a microscopic
link/port construction derives it.
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
from solve_coherent_junction import X, adj, ring_dirac, ring_parameters
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


def positive_ring_bridge(holonomy: np.ndarray, denominator: float = 2.0):
    """Two low ports coupled to physical/reference positive completed rings."""
    size = len(holonomy)
    _, physical, _ = ring_dirac(holonomy)
    _, reference, _ = ring_dirac(np.eye(size))
    heavy_size = len(physical)
    heavy = np.block(
        [
            [physical, np.zeros_like(physical)],
            [np.zeros_like(reference), reference],
        ]
    )
    coupling = np.zeros((2 * size, 2 * heavy_size), dtype=complex)
    amplitude = X / math.sqrt(denominator)
    coupling[:size, :size] = amplitude * np.eye(size)
    coupling[:size, heavy_size : heavy_size + size] = amplitude * np.eye(size)
    coupling[size:, :size] = amplitude * np.eye(size)
    coupling[size:, heavy_size : heavy_size + size] = -amplitude * np.eye(size)
    full = np.block(
        [
            [np.zeros((2 * size, 2 * size)), coupling],
            [adj(coupling), heavy],
        ]
    )
    return full, coupling, heavy


def low_band_splitting(theta: float, denominator: float) -> float:
    holonomy = np.asarray([[np.exp(1j * theta)]])
    full, _, _ = positive_ring_bridge(holonomy, denominator)
    eigenvalues = np.linalg.eigvalsh(full)
    low = sorted(eigenvalues, key=abs)[:2]
    return float(abs(low[1] - low[0]))


def photon_coefficient(denominator: float, dynamic: bool = True) -> float:
    theta = np.asarray((5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3))
    values = []
    for phase in theta:
        holonomy = np.asarray([[np.exp(1j * phase)]])
        full, coupling, heavy = positive_ring_bridge(holonomy, denominator)
        if dynamic:
            eigenvalues = np.linalg.eigvalsh(full)
            low = sorted(eigenvalues, key=abs)[:2]
            values.append(abs(low[1] - low[0]))
        else:
            effective = -coupling @ np.linalg.solve(heavy, adj(coupling))
            eigenvalues = np.linalg.eigvalsh(effective)
            values.append(float(eigenvalues[-1] - eigenvalues[0]))
    design = np.column_stack((theta**2, theta**4))
    coefficients, *_ = np.linalg.lstsq(design, np.asarray(values), rcond=None)
    return float(coefficients[0])


def photon_residue(denominator: float) -> float:
    _, coupling, heavy = positive_ring_bridge(np.eye(1), denominator)
    residue = np.eye(2) + coupling @ np.linalg.solve(heavy @ heavy, adj(coupling))
    return float(np.trace(residue).real / 2.0)


def frame_channel() -> dict[str, float]:
    # B^2=I. F=B/4 is normalized so tr(BF)=1 in four spinor dimensions.
    area = np.diag((1.0, 1.0, -1.0, -1.0)).astype(complex)
    curvature = area / 4.0
    _, coupling, heavy = bridge_with_area(np.eye(4), area)
    residue = np.eye(8) + coupling @ np.linalg.solve(heavy @ heavy, adj(coupling))
    identity = np.eye(4)
    transport = np.vstack((identity, 1j * identity)) / math.sqrt(2.0)
    projected_residue = transport.conj().T @ residue @ transport
    residue_scalar = float(np.trace(projected_residue).real / 4.0)

    epsilon = 1.0e-5
    plus = projected_area_operator(expm(1j * epsilon * curvature), area)
    minus = projected_area_operator(expm(-1j * epsilon * curvature), area)
    derivative = (plus - minus) / (2.0 * epsilon)
    static_b = float(-np.trace(derivative).real / 4.0)
    dynamic_b = static_b / residue_scalar
    return {
        "static_b": static_b,
        "induced_residue": residue_scalar,
        "dynamic_b": dynamic_b,
        "bare_speed": math.sqrt(2.0 * U_COMMON * dynamic_b),
        "derivative_matrix_nonidentity_residual": float(
            np.linalg.norm(
                derivative
                + static_b * np.eye(4),
                2,
            )
        ),
    }


def continuous_match(frame_b: float) -> float:
    lower, upper = 2.0, 40.0
    for _ in range(80):
        midpoint = 0.5 * (lower + upper)
        if photon_coefficient(midpoint) > frame_b:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def build_report() -> dict[str, Any]:
    _, _, _, rho = ring_parameters()
    frame = frame_channel()
    baseline_dynamic = photon_coefficient(2.0, dynamic=True)
    baseline_static = photon_coefficient(2.0, dynamic=False)
    baseline_residue = photon_residue(2.0)

    scan = []
    for denominator in range(2, 33):
        coefficient = photon_coefficient(float(denominator), dynamic=True)
        ratio = math.sqrt(coefficient / frame["dynamic_b"])
        scan.append(
            {
                "denominator": denominator,
                "photon_dynamic_curvature": coefficient,
                "photon_to_frame_speed_ratio": ratio,
                "fractional_cone_mismatch": ratio - 1.0,
            }
        )
    best = min(scan, key=lambda row: abs(row["fractional_cone_mismatch"]))
    required = continuous_match(frame["dynamic_b"])

    baseline_ratio = math.sqrt(baseline_dynamic / frame["dynamic_b"])
    checks = {
        "frame_derivative_is_scalar_in_physical_spin_space": (
            frame["derivative_matrix_nonidentity_residual"] < 1.0e-12
        ),
        "positive_ring_residue_reproduces_dynamic_reduction": abs(
            baseline_dynamic - baseline_static / baseline_residue
        )
        / baseline_dynamic
        < 2.0e-4,
        "baseline_common_parent_cones_do_not_match": abs(baseline_ratio - 1.0) > 0.1,
        "no_integer_normalization_matches_exactly": abs(
            best["fractional_cone_mismatch"]
        )
        > 1.0e-4,
        "sixteen_denominator_is_near_but_not_derived": (
            best["denominator"] == 16
            and abs(best["fractional_cone_mismatch"]) < 0.01
        ),
        "required_denominator_is_not_integer": abs(required - round(required)) > 1.0e-3,
    }

    return {
        "module": (
            "phase_junction_network/construction/coherent_junction/"
            "check_even_odd_cone_parent.py"
        ),
        "status": (
            "Stage 5F REJECTION/SEARCH: the declared two-path even photon and "
            "odd frame channels of the coherent ring do not share a bare cone; "
            "a sixteen-denominator discrete candidate is close but unprotected"
        ),
        "ring": {
            "x": X,
            "rho": rho,
            "photon_channel": "orientation-even positive completed ring Q=K^2",
            "frame_channel": "orientation-odd chiral K area-curvature bridge",
            "new_continuous_coefficients": 0,
        },
        "baseline_two_path_photon": {
            "static_curvature_coefficient": baseline_static,
            "induced_residue": baseline_residue,
            "dynamic_curvature_coefficient": baseline_dynamic,
            "bare_speed": math.sqrt(2.0 * U_COMMON * baseline_dynamic),
        },
        "frame": frame,
        "baseline_photon_to_frame_speed_ratio": baseline_ratio,
        "integer_port_normalization_scan": scan,
        "best_integer_candidate": best,
        "continuous_denominator_required_for_exact_match": required,
        "checks": checks,
        "execution_pass": all(checks.values()),
        "decision": (
            "The current equal two-path coherent parent is not a common-cone "
            "theory. Weakening the photon low-port amplitude by an unexplained "
            "continuous factor would be a fit and is prohibited. Denominator 16 "
            "is a discrete near-match, but it must be derived from an explicit "
            "finite link/port Hilbert space before it can replace the rejected "
            "baseline."
        ),
        "claim_boundary": {
            "established": [
                "an explicit even photon parent built from Q=K^2",
                "an explicit odd frame parent built from K and the area operator",
                "both frequency residues retained in the low-band coefficients",
                "the declared two-path parent has a large bare cone mismatch",
                "integer denominator 16 gives the nearest discrete match but is not exact",
            ],
            "not_established": [
                "a microscopic derivation of denominator 16",
                "a complete finite-spin three-dimensional photon phase for the Q parent",
                "matter-loop stability of any revised port architecture",
                "a physical common photon/gravity cone",
            ],
        },
        "next_gate": (
            "Construct the finite photon link/port Hilbert space explicitly. "
            "Either derive the sixteen-channel normalization and rerun the "
            "Fock stress gate, or reject the coherent common-parent branch."
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
                "baseline_ratio": report[
                    "baseline_photon_to_frame_speed_ratio"
                ],
                "best_integer": report["best_integer_candidate"],
                "required_denominator": report[
                    "continuous_denominator_required_for_exact_match"
                ],
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
