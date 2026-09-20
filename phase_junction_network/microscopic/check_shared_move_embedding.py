#!/usr/bin/env python3
"""Stage 3B: test the Stage-3A shared move against the actual photon phase.

The naive Stage-3A common-Casimir hypothesis identifies the virtual defect
penalty Delta with the electromagnetic and frame charging coefficients. This
script maps that prediction into the Stage-6B finite photon Hamiltonian and
checks whether it lies in the demonstrated z~1 region.

A negative result is an intended scientific outcome: the executable passes
when it reproduces the incompatibility and quantifies the missing collective
charging mechanism.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from check_shared_junction_move import solve_common_cone
from finite_em_dynamics_impl import axial_spectrum, fit_branch, public_spectrum_row


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


def build_report() -> dict[str, Any]:
    common = solve_common_cone()
    k_over_delta = common["K_A_over_Delta"]
    naive_u_over_t = 1.0 / k_over_delta

    detuned_rows = []
    detuned_fits = {}
    for v_over_t in (0.2, 0.4, 0.6):
        rows = [
            axial_spectrum(
                spin=1,
                lattice_size=lattice_size,
                u_over_t=naive_u_over_t,
                v_over_t=v_over_t,
                eigenpairs=24,
            )
            for lattice_size in (6, 8, 10)
        ]
        key = f"u/t={naive_u_over_t:.9f}:v/t={v_over_t:.1f}"
        detuned_rows.extend(public_spectrum_row(row) for row in rows)
        detuned_fits[key] = fit_branch(rows)

    powers = [
        fit["gap_proportional_to_lattice_momentum_power"]
        for fit in detuned_fits.values()
    ]
    smallest_l10_gap = min(
        row["single_polarization_gap"]
        for row in detuned_rows
        if row["lattice_size"] == 10
    )

    demonstrated_u_over_t_max = 0.2
    representative_u_over_t = 0.1
    suppression_for_edge = demonstrated_u_over_t_max / naive_u_over_t
    suppression_for_representative = representative_u_over_t / naive_u_over_t

    gravity_self_dual_u_over_k = 1.0
    gravity_self_dual_suppression = (
        gravity_self_dual_u_over_k / naive_u_over_t
    )

    report = {
        "module": (
            "phase_junction_network/microscopic/"
            "check_shared_move_embedding.py"
        ),
        "status": (
            "Stage 3B diagnostic PASS: the naive one-hinge shared-Casimir "
            "embedding is rejected by the actual finite photon dynamics; "
            "a derived common collective charging suppression is required"
        ),
        "input_stage_3A": common,
        "mapping_to_stage_6B": {
            "identification": {
                "stage6B_plaquette_t": "K_A",
                "stage6B_electric_u": "U_A=Delta",
            },
            "implied_u_over_t": naive_u_over_t,
            "demonstrated_stage6B_u_over_t_range": [0.0, 0.2],
            "factor_above_demonstrated_upper_edge": (
                naive_u_over_t / demonstrated_u_over_t_max
            ),
        },
        "exact_spin1_axial_test_at_naive_ratio": {
            "v_over_t_values": [0.2, 0.4, 0.6],
            "lattice_sizes": [6, 8, 10],
            "rows": detuned_rows,
            "fits": detuned_fits,
            "minimum_power": min(powers),
            "maximum_power": max(powers),
            "minimum_L10_gap": smallest_l10_gap,
            "decision": (
                "The gaps remain order one in plaquette units and the fitted "
                "powers are near zero, not one. The naive common-Casimir point "
                "is not in the demonstrated photon phase."
            ),
        },
        "required_revision": {
            "common_collective_charging_factor": (
                "Introduce one derived factor r with U_A=U_g=r*Delta while "
                "retaining the shared t and Delta in the virtual processes."
            ),
            "maximum_r_to_enter_stage6B_edge": suppression_for_edge,
            "r_for_representative_u_over_t_0p1": (
                suppression_for_representative
            ),
            "optional_bundle_interpretation": {
                "if_r_equals_1_over_N_eff_minimum_N_for_edge": (
                    1.0 / suppression_for_edge
                ),
                "if_r_equals_1_over_N_eff_N_for_u_over_t_0p1": (
                    1.0 / suppression_for_representative
                ),
                "warning": (
                    "N_eff is only a target interpretation. It must be "
                    "derived from a finite capacitance or collective-mode "
                    "matrix; it is not obtained merely by declaring N hinges."
                ),
            },
            "gravity_self_dual_regulator": {
                "current_U_g_over_K_g": gravity_self_dual_u_over_k,
                "required_r_at_stage3A_K": gravity_self_dual_suppression,
                "mismatch_with_representative_photon_r": (
                    gravity_self_dual_suppression
                    / suppression_for_representative
                ),
                "decision": (
                    "The existing self-dual gravity clock normalization and a "
                    "representative Stage-6B photon point cannot both descend "
                    "from one r. The gravity regulator must be rerun at the "
                    "same derived charging ratio rather than independently "
                    "rescaled."
                ),
            },
        },
        "parameter_count": {
            "stage3A_naive": (
                "Delta and x=t/Delta; common cone fixes x and leaves Delta."
            ),
            "after_stage3B_obstruction": (
                "A common charging factor r is additionally required unless "
                "the full finite embedding generates it automatically."
            ),
            "forbidden_escape": (
                "Independent r_A and r_g would fit both sectors trivially but "
                "destroy the impedance prediction."
            ),
        },
        "checks": {
            "naive_ratio_outside_demonstrated_stage6B_region": (
                naive_u_over_t > demonstrated_u_over_t_max
            ),
            "naive_ratio_not_z1_in_exact_axial_blocks": max(powers) < 0.5,
            "naive_ratio_remains_gapped_at_L10": smallest_l10_gap > 1.0,
            "common_suppression_needed": suppression_for_edge < 0.02,
            "current_gravity_and_photon_normalizations_not_shared": (
                abs(
                    gravity_self_dual_suppression
                    / suppression_for_representative
                    - 1.0
                )
                > 1.0
            ),
        },
        "claim_boundary": {
            "established": [
                (
                    "the naive Stage-3A identification U_A=U_g=Delta maps "
                    "to u/t=14.5968888"
                ),
                (
                    "the exact spin-1 axial blocks at that ratio do not have "
                    "the Stage-6B z~1 scaling for v/t=0.2,0.4,0.6"
                ),
                (
                    "a common charging suppression of order 0.014 or smaller "
                    "is required to enter the demonstrated photon region"
                ),
                (
                    "the existing gravity self-dual normalization cannot be "
                    "kept independently if the same suppression is used"
                ),
            ],
            "not_established": [
                "a microscopic origin for the common suppression r",
                "that the full three-dimensional phase boundary ends at u/t=0.2",
                "that a collective hinge bundle necessarily gives r=1/N",
                "a revised shared-move model that passes both sectors",
                "physical alpha, G, or a common interacting cone",
            ],
        },
        "next_gate": (
            "Construct an explicit finite collective charging/capacitance "
            "matrix for the matter-geometry junction bundle, derive r rather "
            "than fit it, rerun the finite photon and gravity regulators at "
            "that same r, and reject the branch if separate clock "
            "normalizations are still required."
        ),
    }

    report["diagnostic_pass"] = all(report["checks"].values())
    if not report["diagnostic_pass"]:
        raise AssertionError(
            json.dumps(ready(report["checks"]), indent=2, sort_keys=True)
        )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "phase_junction_network/microscopic/"
            "shared_move_embedding_results.json"
        ),
    )
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
                "diagnostic_pass": report["diagnostic_pass"],
                "mapping": report["mapping_to_stage_6B"],
                "required_revision": report["required_revision"],
                "checks": report["checks"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
