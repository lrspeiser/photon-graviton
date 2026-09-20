#!/usr/bin/env python3
"""Run finite charged chiral endpoint checks for Phase Junction issue #4."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from chiral_matter_model import (
    CHARGES, DEFAULT_LS, ETA, R0, SpeciesResult, analytic_residual_gap, anomaly_sums,
    clifford_residual, finite_width_scaling, has_vectorlike_pair,
    localization_ratio, low_energy_dispersion, nonzero_corner_gap, primitive,
    residual_gap, search_odd_chiral_sets, spectral_pairing_residual,
    wall_localized_chirality, wilson_mass,
)
from chiral_matter_audits import (
    bilinear_mass_audit, binding_subspace_audit, exchange_phase,
    flux_bundle_covariance, topological_window_audit, universal_frame_audit,
)

def build_results(slab_width: int) -> dict[str, object]:
    charge_search = search_odd_chiral_sets()
    linear_anomaly, cubic_anomaly = anomaly_sums(CHARGES)
    unique_counts: dict[int, int] = {}
    species_results: list[SpeciesResult] = []
    corner_details: dict[str, dict[str, float]] = {}

    for charge in CHARGES:
        unique_counts[charge] = unique_counts.get(charge, 0) + 1
        ratio = localization_ratio(charge)
        mass = wilson_mass(charge)
        gap = residual_gap(mass, slab_width)
        analytic = analytic_residual_gap(ratio, slab_width)
        relative_error = abs(gap - analytic) / max(gap, 1.0e-300)
        velocity, dispersion_error = low_energy_dispersion(mass, slab_width)
        left_position, right_position, left_chi, right_chi = wall_localized_chirality(
            mass, slab_width
        )
        low_corner_count, nonzero_minimum, corners = nonzero_corner_gap(mass, slab_width)
        corner_details[str(charge)] = corners
        if low_corner_count != 1:
            raise AssertionError(
                f"charge {charge}: expected exactly one low-energy physical corner, got {low_corner_count}"
            )
        species_results.append(
            SpeciesResult(
                charge=charge,
                multiplicity_index=unique_counts[charge],
                localization_ratio=ratio,
                wilson_mass=mass,
                residual_gap=gap,
                analytic_gap=analytic,
                analytic_relative_error=relative_error,
                fitted_velocity=velocity,
                dispersion_max_relative_error=dispersion_error,
                left_wall_position=left_position,
                right_wall_position=right_position,
                left_chirality=left_chi,
                right_chirality=right_chi,
                nonzero_corner_minimum_gap=nonzero_minimum,
                spectral_pairing_residual=spectral_pairing_residual(mass, slab_width),
                finite_ls_slopes=finite_width_scaling(mass, ratio),
            )
        )

    flux_checks = [flux_bundle_covariance(charge) for charge in sorted(set(CHARGES))]
    binding_checks = [binding_subspace_audit(charge) for charge in sorted(set(CHARGES))]
    mass_audit = bilinear_mass_audit(CHARGES)
    frame_audit = universal_frame_audit(slab_width)
    topological_window = topological_window_audit(slab_width)

    distinct_gaps = {
        str(charge): residual_gap(wilson_mass(charge), slab_width)
        for charge in sorted(set(CHARGES), key=lambda value: abs(value))
    }
    positive_distinct_gaps = [value for value in distinct_gaps.values() if value > 0.0]
    hierarchy_ratio = max(positive_distinct_gaps) / min(positive_distinct_gaps)
    duplicate_charge_gap_spreads = {
        str(charge): float(
            max(result.residual_gap for result in species_results if result.charge == charge)
            - min(result.residual_gap for result in species_results if result.charge == charge)
        )
        for charge in sorted(set(CHARGES))
        if sum(1 for value in CHARGES if value == charge) > 1
    }

    checks = {
        "clifford_algebra": clifford_residual() < 1.0e-12,
        "selected_charge_set": tuple(charge_search["first_representatives"][0]) == CHARGES,
        "linear_anomaly_cancels": linear_anomaly == 0,
        "cubic_anomaly_cancels": cubic_anomaly == 0,
        "no_vectorlike_pair": not has_vectorlike_pair(CHARGES),
        "all_endpoint_charges_odd": all(abs(charge) % 2 == 1 for charge in CHARGES),
        "all_exchange_phases_fermionic": all(exchange_phase(charge) == -1 for charge in CHARGES),
        "exact_bound_composite_subspace": all(
            check["zero_penalty_states_per_orbital"] == 2
            and check["bound_fock_dimension_for_four_orbitals"] == 16
            and check["exact_invariant_bound_band"]
            for check in binding_checks
        ),
        "exact_flux_bundle_gauss_covariance": all(
            check["max_link_charge_residual"] == 0
            and check["max_left_gauss_residual"] == 0
            and check["max_right_gauss_residual"] == 0
            for check in flux_checks
        ),
        "one_physical_weyl_corner_per_species": all(
            result.nonzero_corner_minimum_gap > 1.0 for result in species_results
        ),
        "opposite_wall_chiralities": all(
            result.left_chirality < -0.999999
            and result.right_chirality > 0.999999
            for result in species_results
        ),
        "wall_localization": all(
            result.left_wall_position < -0.90 and result.right_wall_position > 0.90
            for result in species_results
        ),
        "positive_linear_dispersion": all(
            0.995 < result.fitted_velocity < 1.005
            and result.dispersion_max_relative_error < 2.0e-3
            for result in species_results
        ),
        "particle_antiparticle_pairing": all(
            result.spectral_pairing_residual < 1.0e-12 for result in species_results
        ),
        "exponential_mass_protection": all(
            result.analytic_relative_error < 2.0e-5
            and result.finite_ls_slopes["slope_error"] < 3.0e-3
            for result in species_results
        ),
        "no_same_wall_bilinear_mass": mass_audit[
            "gauge_neutral_left_left_pair_count"
        ]
        == 0
        and mass_audit["two_flavor_same_chirality_mass_map_nullity"] == 0,
        "universal_frame_operator": frame_audit[
            "maximum_species_derivative_difference"
        ]
        < 1.0e-10
        and frame_audit["maximum_derivative_formula_residual"] < 2.0e-8,
        "topological_window_stability": topological_window["failures"] == 0,
        "nontrivial_charge_mass_hierarchy": hierarchy_ratio > 1.0e7,
        "duplicate_charge_mass_degeneracy": all(
            spread < 1.0e-15 for spread in duplicate_charge_gap_spreads.values()
        ),
    }

    failed = sorted(name for name, passed in checks.items() if not passed)
    result = {
        "protocol": "phase-junction finite chiral endpoint v1",
        "status": "PASS" if not failed else "FAIL",
        "interpretation": (
            "Finite anomaly-free chiral endpoint prototype passes its declared algebraic and spectral gates; "
            "it does not derive the Standard Model, observed masses, mirror-wall removal, or the interacting continuum."
        ),
        "finite_hilbert": {
            "spin1_dimension_per_electromagnetic_lane": 3,
            "maximum_parallel_lanes_per_geometric_link": max(abs(q) for q in CHARGES),
            "maximum_link_bundle_dimension": 3 ** max(abs(q) for q in CHARGES),
            "microscopic_fermionic_strands_for_charge_q": "|q|",
            "wilson_orbitals_per_strand_per_site_per_fifth_layer": 4,
            "unconstrained_fock_dimension_for_charge_q": "2^(4|q|)",
            "maximum_unconstrained_fock_dimension_per_species_site_layer": 2 ** (4 * max(abs(q) for q in CHARGES)),
            "binding_projector": "sum_{a<b,A} Lambda*(n_{a,A}-n_{b,A})^2",
            "bound_composite_fock_dimension_per_species_site_layer": 2**4,
            "fifth_layers": slab_width,
            "physical_sector": "Gauss law plus exact equal-occupation strand binding; all selected |q| are odd",
        },
        "charges": list(CHARGES),
        "charge_search": charge_search,
        "anomalies": {
            "sum_q": linear_anomaly,
            "sum_q_cubed": cubic_anomaly,
            "primitive": primitive(CHARGES),
            "vectorlike_pair_present": has_vectorlike_pair(CHARGES),
        },
        "charge_statistics_relation": {
            "composite_operator": "D_q,A^dagger = product_{a=1}^{|q|} f_q,a,A^dagger",
            "binding_constraint": "n_q,a,A = n_q,b,A for every strand pair a,b",
            "exchange_phases": {str(q): exchange_phase(q) for q in sorted(set(CHARGES))},
            "consequence": "exchanging two composites gives (-1)^(q^2)=-1 for every selected odd charge",
        },
        "binding_subspace_checks": binding_checks,
        "flux_bundle_checks": flux_checks,
        "frame_coupling": frame_audit,
        "topological_window_stability": topological_window,
        "mass_protection": {
            "mechanism": "opposite-chirality wall overlap only; no gauge-invariant same-wall bilinear",
            "shared_relation": "r_q = r0 + eta*(|q|-1), m5(q) = -1 + r_q",
            "r0": R0,
            "eta": ETA,
            "slab_width": slab_width,
            "gap_hierarchy_ratio": hierarchy_ratio,
            "distinct_gaps": distinct_gaps,
            "duplicate_charge_gap_spreads": duplicate_charge_gap_spreads,
            "bilinear_audit": mass_audit,
        },
        "species": [result.__dict__ for result in species_results],
        "brillouin_corner_gaps": corner_details,
        "checks": checks,
        "failed_checks": failed,
        "claim_boundary": {
            "established": [
                "finite local endpoint/link Hilbert space",
                "exact quantum-link Gauss covariance of charge-bundle hopping",
                "fermionic endpoint exchange under the parity constraint",
                "anomaly-free same-chirality charge spectrum in the searched finite class",
                "one wall-localized Weyl cone per species and no physical-corner doublers",
                "stability throughout a finite local Wilson-mass neighborhood",
                "exponentially protected finite-slab gaps and a charge-linked hierarchy",
                "universal bare frame operator for every species",
                "particle/antiparticle spectral pairing",
            ],
            "not_established": [
                "observed Standard Model charge assignments or generations",
                "observed fermion masses or Yukawa couplings",
                "a purely 3+1D regulator with the remote mirror wall removed",
                "nonperturbative mirror-wall symmetric mass generation",
                "an interacting QED continuum or radiative stability",
                "common microscopic derivation of r0, eta, gauge stiffness, and gravity stiffness",
            ],
        },
    }
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("phase_junction_network/microscopic/chiral_matter_results.json"),
        help="JSON output path",
    )
    parser.add_argument(
        "--slab-width",
        type=int,
        default=DEFAULT_LS,
        help="finite fifth-direction width",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.slab_width < 6:
        raise ValueError("slab width must be at least 6 for the frozen checks")
    results = build_results(args.slab_width)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": results["status"],
        "failed_checks": results["failed_checks"],
        "charges": results["charges"],
        "gap_hierarchy_ratio": results["mass_protection"]["gap_hierarchy_ratio"],
    }, indent=2, sort_keys=True))
    return 0 if results["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
