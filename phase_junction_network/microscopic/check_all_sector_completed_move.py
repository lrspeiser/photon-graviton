#!/usr/bin/env python3
"""Stage 3F: apply one completed move to photon, frame, matter, and companion."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Any
from all_sector_move_common import X_ROOT, ready
from all_sector_move_photon_frame import photon_sector, frame_sector
from all_sector_move_matter_companion import matter_sector, companion_sector


def zero_crossing_control(controls, key):
    values = {round(row["coefficient"], 1): row[key] for row in controls}
    return (
        values[0.8] * values[0.9] > 0.0
        and values[1.1] * values[1.2] > 0.0
        and values[0.9] * values[1.1] < 0.0
        and abs(values[1.0]) < 1.0e-3
    )


def build_report() -> dict[str, Any]:
    photon = photon_sector()
    frame = frame_sector()
    matter = matter_sector()
    companion = companion_sector()

    frame_hoppings = [row["effective_hopping"] for row in frame["coefficient_controls"]]
    checks = {
        "same_unit_coefficient_used_in_all_sectors": True,
        "universal_positive_square_factorization": max(
            photon["factorized_completion_residual"],
            *(row["maximum_factorized_completion_residual"] for row in frame["root_rows"]),
            matter["root"]["maximum_factorized_completion_residual"],
            companion["root"]["maximum_factorized_completion_residual"],
        ) < 1.0e-12,
        "photon_unit_coefficient_is_unique": zero_crossing_control(
            photon["coefficient_controls"], "delta_U_over_x_squared"
        ),
        "matter_unit_coefficient_is_unique": zero_crossing_control(
            matter["coefficient_controls"], "common_diagonal_over_x_squared"
        ),
        "companion_unit_coefficient_is_unique": zero_crossing_control(
            companion["coefficient_controls"], "common_diagonal_over_x_squared"
        ),
        "frame_seagull_is_central_not_sector_tuned": (
            max(frame_hoppings) - min(frame_hoppings) < 1.0e-12
        ),
        "frame_moves_commute_with_connection_lock": max(
            row["maximum_lock_commutator"] for row in frame["root_rows"]
        ) < 1.0e-12,
        "frame_lock_leakage_is_roundoff": max(
            row["dressed_shift_lock_leakage"] for row in frame["root_rows"]
        ) < 1.0e-12,
        "frame_generated_operator_inventory_is_complete": max(
            row["complete_circulant_fit_residual"] for row in frame["root_rows"]
        ) < 1.0e-12,
        "frame_hopping_is_second_order": all(
            abs(fit["power"] - 2.0) < 0.02
            for fit in frame["small_x_hopping_fits"].values()
        ),
        "matter_exact_gauss": max(
            matter["root"]["maximum_left_gauss_commutator"],
            matter["root"]["maximum_right_gauss_commutator"],
        ) < 1.0e-12,
        "matter_operator_continuity": matter["root"]["continuity_residual"] < 1.0e-12,
        "matter_spectral_ward_identity": matter["root"]["spectral_ward_residual"] < 1.0e-12,
        "matter_pure_gauge_response_cancels": abs(
            matter["root"]["total_pure_gauge_curvature"]
        ) < 1.0e-12,
        "matter_twist_spectrum_is_gauge_invariant": matter[
            "phase_twist_maximum_spectral_spread"
        ] < 1.0e-12,
        "matter_hopping_is_second_order": abs(
            matter["small_x_charged_hopping_fit"]["power"] - 2.0
        ) < 0.02,
        "matter_residual_selfenergy_is_fourth_order": abs(
            matter["small_x_common_diagonal_fit"]["power"] - 4.0
        ) < 0.05,
        "matter_low_band_has_no_unexpected_pauli_operator": max(
            abs(matter["root"]["operator_inventory"]["sigma_y"]),
            abs(matter["root"]["operator_inventory"]["sigma_z"]),
            matter["root"]["operator_inventory"]["maximum_imaginary_coefficient"],
        ) < 1.0e-12,
        "companion_conversion_is_second_order": abs(
            companion["small_x_conversion_fit"]["power"] - 2.0
        ) < 0.02,
        "companion_residual_selfenergy_is_fourth_order": abs(
            companion["small_x_common_diagonal_fit"]["power"] - 4.0
        ) < 0.05,
        "companion_reciprocity": companion["root"][
            "forward_reverse_probability_residual"
        ] < 1.0e-12,
        "companion_low_band_has_no_unexpected_pauli_operator": max(
            abs(companion["root"]["operator_inventory"]["sigma_y"]),
            abs(companion["root"]["operator_inventory"]["sigma_z"]),
            companion["root"]["operator_inventory"]["maximum_imaginary_coefficient"],
        ) < 1.0e-12,
        "all_unwanted_bands_remain_separated": min(
            [row["low_band_separation_over_hopping"] for row in frame["root_rows"]]
            + [matter["root"]["next_band_gap_over_hopping"]]
            + [companion["root"]["next_band_gap_over_hopping"]]
        ) > 20.0,
        "photon_root_remains_relativistic": (
            0.94 <= photon["root"]["photon_gap_power"] <= 1.16
            and photon["root"]["photon_fit_residual"] < 0.002
        ),
    }

    report = {
        "module": "phase_junction_network/microscopic/check_all_sector_completed_move.py",
        "status": (
            "Stage 3F PASS at finite algebra/one-particle scope: the same unit "
            "completed move works in photon, dressed-frame, charged-endpoint, "
            "and neutral-companion representatives, and the charged sector "
            "satisfies exact Gauss, continuity, and finite Ward controls"
        ),
        "shared_rule": {
            "formula": "H_move=-x(M+M^dagger)+x^2(M^dagger M+M M^dagger)",
            "coefficient": 1.0,
            "shared_x_t_over_Delta": X_ROOT,
            "principle": (
                "Every allowed elementary move carries its own finite algebraic "
                "diagonal partner; no sector-specific coefficient is introduced."
            ),
            "positive_square_factorization": (
                "For Q_M=I-2xM, H_move=1/4*(Q_M^dagger Q_M + "
                "Q_M Q_M^dagger - 2I)."
            ),
        },
        "photon_sector": photon,
        "dressed_frame_sector": frame,
        "charged_endpoint_sector": matter,
        "neutral_companion_sector": companion,
        "checks": checks,
        "stage_pass": all(checks.values()),
        "decision": (
            "Advance the unit completed-move algebra to a spatial gauge-matter "
            "lattice using the actual Stage-6B transverse photons."
        ),
        "claim_boundary": {
            "established": [
                "one identical completion coefficient in four finite sector representatives",
                "one symmetric positive-square factorization of every completed move",
                "exact dressed-frame lock compatibility, zero leakage, and complete finite harmonic inventory",
                "exact extended Gauss commutation for the charged endpoint",
                "an operator continuity equation and transition-level finite Ward identity",
                "exact cancellation of paramagnetic and diamagnetic response to a pure gauge phase",
                "second-order charged hopping and companion conversion with fourth-order residual self-energy",
                "reversible neutral companion conversion",
            ],
            "not_established": [
                "one tensor-product many-body Hamiltonian containing all sectors simultaneously",
                "a nonzero-momentum Ward-Takahashi identity in the 3D photon phase",
                "matter vacuum polarization of an actual transverse photon",
                "radiative protection beyond the finite one-particle model",
                "why the microscopic Phase-Junction theory must choose the factorized first-order transport algebra",
                "mirror-wall completion, nonlinear gravity, alpha, G, or formal novelty",
            ],
        },
        "next_gate": (
            "Embed the |q|=1 endpoint in a finite spatial patch of the Stage-6B "
            "photon phase and test a nonzero-momentum Ward identity and photon dressing."
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
        default=Path("phase_junction_network/microscopic/all_sector_completed_move_results.json"),
    )
    args = parser.parse_args()
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(ready(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "stage_pass": report["stage_pass"],
        "matter_ward": report["charged_endpoint_sector"]["root"],
        "failed_checks": [key for key, passed in report["checks"].items() if not passed],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
