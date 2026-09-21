#!/usr/bin/env python3
"""Stage 5G: test the equal-coupling coherent ring family for a common cone.

The four-stage result might be an accident of ring length. This check rebuilds
the same positive even photon channel and chiral odd frame channel for rings of
2 through 12 elementary transports. All local couplings retain the inherited x
and the same physical/reference 1/sqrt(2) port normalization.

The test is a bounded no-go for this architecture class, not for all possible
Phase Junction theories.
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


def ring_shift(holonomy: np.ndarray, stages: int) -> np.ndarray:
    size = len(holonomy)
    shift = np.zeros((stages * size, stages * size), dtype=complex)
    for stage in range(stages - 1):
        shift[
            (stage + 1) * size : (stage + 2) * size,
            stage * size : (stage + 1) * size,
        ] = np.eye(size)
    shift[:size, (stages - 1) * size :] = holonomy
    return shift


def ring_blocks(holonomy: np.ndarray, stages: int):
    _, alpha, r, _ = ring_parameters()
    shift = ring_shift(holonomy, stages)
    first_order = np.sqrt(alpha) * (np.eye(len(shift)) - r * shift)
    zero = np.zeros_like(first_order)
    chiral = np.block(
        [[zero, first_order], [first_order.conj().T, zero]]
    )
    positive = first_order.conj().T @ first_order
    return chiral, positive


def positive_bridge(holonomy: np.ndarray, stages: int):
    size = len(holonomy)
    _, physical = ring_blocks(holonomy, stages)
    _, reference = ring_blocks(np.eye(size), stages)
    heavy_size = len(physical)
    heavy = np.block(
        [[physical, np.zeros_like(physical)], [np.zeros_like(reference), reference]]
    )
    coupling = np.zeros((2 * size, 2 * heavy_size), dtype=complex)
    amplitude = X / math.sqrt(2.0)
    coupling[:size, :size] = amplitude * np.eye(size)
    coupling[:size, heavy_size : heavy_size + size] = amplitude * np.eye(size)
    coupling[size:, :size] = amplitude * np.eye(size)
    coupling[size:, heavy_size : heavy_size + size] = -amplitude * np.eye(size)
    full = np.block(
        [
            [np.zeros((2 * size, 2 * size)), coupling],
            [coupling.conj().T, heavy],
        ]
    )
    return full, coupling, heavy


def photon_curvature(stages: int) -> float:
    phases = np.asarray((5.0e-4, 1.0e-3, 2.0e-3, 4.0e-3))
    splittings = []
    for phase in phases:
        full, _, _ = positive_bridge(
            np.asarray([[np.exp(1j * phase)]]), stages
        )
        eigenvalues = np.linalg.eigvalsh(full)
        low = sorted(eigenvalues, key=abs)[:2]
        splittings.append(abs(low[1] - low[0]))
    design = np.column_stack((phases**2, phases**4))
    coefficients, *_ = np.linalg.lstsq(
        design, np.asarray(splittings), rcond=None
    )
    return float(coefficients[0])


def frame_bridge(holonomy: np.ndarray, area: np.ndarray, stages: int):
    size = len(holonomy)
    physical, _ = ring_blocks(holonomy, stages)
    reference, _ = ring_blocks(np.eye(size), stages)
    heavy_size = len(physical)
    heavy = np.block(
        [[physical, np.zeros_like(physical)], [np.zeros_like(reference), reference]]
    )
    coupling = np.zeros((2 * size, 2 * heavy_size), dtype=complex)
    amplitude = X / math.sqrt(2.0)
    coupling[:size, :size] = amplitude * area
    coupling[:size, heavy_size : heavy_size + size] = amplitude * area
    lower_start = stages * size
    coupling[
        size:, lower_start : lower_start + size
    ] = amplitude * np.eye(size)
    coupling[
        size:,
        heavy_size + lower_start : heavy_size + lower_start + size,
    ] = -amplitude * np.eye(size)
    full = np.block(
        [
            [np.zeros((2 * size, 2 * size)), coupling],
            [coupling.conj().T, heavy],
        ]
    )
    effective = -coupling @ np.linalg.solve(heavy, coupling.conj().T)
    return full, effective, coupling, heavy


def projected_area(holonomy: np.ndarray, area: np.ndarray, stages: int):
    _, effective, _, _ = frame_bridge(holonomy, area, stages)
    size = len(holonomy)
    cross = effective[:size, size:]
    return 0.5j * (cross - cross.conj().T)


def frame_curvature(stages: int) -> dict[str, float]:
    area = np.diag((1.0, 1.0, -1.0, -1.0)).astype(complex)
    curvature = area / 4.0
    _, _, coupling, heavy = frame_bridge(np.eye(4), area, stages)
    residue = np.eye(8) + coupling @ np.linalg.solve(
        heavy @ heavy, coupling.conj().T
    )
    transport = np.vstack((np.eye(4), 1j * np.eye(4))) / math.sqrt(2.0)
    projected_residue = transport.conj().T @ residue @ transport
    residue_scalar = float(np.trace(projected_residue).real / 4.0)
    epsilon = 1.0e-5
    derivative = (
        projected_area(expm(1j * epsilon * curvature), area, stages)
        - projected_area(expm(-1j * epsilon * curvature), area, stages)
    ) / (2.0 * epsilon)
    static = float(-np.trace(derivative).real / 4.0)
    dynamic = static / residue_scalar
    return {
        "static_curvature": static,
        "induced_residue": residue_scalar,
        "dynamic_curvature": dynamic,
        "derivative_nonidentity_residual": float(
            np.linalg.norm(derivative + static * np.eye(4), 2)
        ),
    }


def build_report() -> dict[str, Any]:
    rows = []
    for stages in range(2, 13):
        photon = photon_curvature(stages)
        frame = frame_curvature(stages)
        ratio = math.sqrt(photon / frame["dynamic_curvature"])
        rows.append(
            {
                "stages": stages,
                "winding_factor_r_to_N": ring_parameters()[2] ** stages,
                "photon_dynamic_curvature": photon,
                "frame_dynamic_curvature": frame["dynamic_curvature"],
                "frame_induced_residue": frame["induced_residue"],
                "photon_to_frame_speed_ratio": ratio,
                "fractional_cone_mismatch": ratio - 1.0,
                "frame_derivative_nonidentity_residual": frame[
                    "derivative_nonidentity_residual"
                ],
            }
        )
    best = min(rows, key=lambda row: abs(row["fractional_cone_mismatch"]))
    ratios = [row["photon_to_frame_speed_ratio"] for row in rows]
    checks = {
        "all_frame_vertices_are_scalar_in_physical_spin_space": max(
            row["frame_derivative_nonidentity_residual"] for row in rows
        )
        < 1.0e-11,
        "no_tested_ring_length_has_a_common_cone": min(
            abs(row["fractional_cone_mismatch"]) for row in rows
        )
        > 1.0,
        "speed_ratio_stays_above_two": min(ratios) > 2.0,
        "all_tested_mismatches_have_the_same_sign": all(
            row["fractional_cone_mismatch"] > 0.0 for row in rows
        ),
    }
    return {
        "module": (
            "phase_junction_network/construction/coherent_junction/"
            "check_ring_family_cone_no_go.py"
        ),
        "status": (
            "Stage 5G NO-GO: changing the number of equal-coupling coherent "
            "ring stages from 2 through 12 does not produce a common photon/"
            "frame cone"
        ),
        "architecture_class": {
            "ring_stages": [2, 12],
            "local_move_amplitude": X,
            "low_port_normalization": "X/sqrt(2) in both sectors",
            "photon": "orientation-even positive ring",
            "frame": "orientation-odd chiral area bridge",
            "new_coefficients": 0,
        },
        "rows": rows,
        "best_tested_ring": best,
        "checks": checks,
        "execution_pass": all(checks.values()),
        "decision": (
            "The Stage-5F mismatch is not a four-stage accident. Equal inherited "
            "low-port coupling across this coherent ring family keeps the photon "
            "more than twice as fast as the frame mode. Altering ring length is "
            "not a solution. A new representation-theoretic port projection or "
            "a different microscopic parent is required."
        ),
        "claim_boundary": {
            "established": [
                "a bounded no-go for equal-coupling rings of 2 through 12 stages",
                "the even and odd channels were both dynamically normalized",
                "the mismatch has the same sign at every tested ring length",
                "no continuous normalization was introduced",
            ],
            "not_established": [
                "a no-go for arbitrary graphs or port representations",
                "a no-go for tensor-product internal Hilbert spaces",
                "a no-go for symmetry-derived scalar projections",
                "a physical common photon/gravity cone",
            ],
        },
        "next_gate": (
            "Test explicit finite internal port representations, beginning with "
            "the common first-order chiral generator and the four-stage by "
            "four-spinor channel. Any projection must be derived from an exact "
            "symmetry and must pass the matter-stress Ward gate."
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
                "best": report["best_tested_ring"],
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
