#!/usr/bin/env python3
"""Finite spin-link Coulomb diagnostics for the Phase Junction network.

The executable deliberately separates exact finite-link/static Coulomb claims
from the still-open relativistic photon and interacting-QED gates.
"""
from __future__ import annotations

import argparse
import collections
import json
import math
from pathlib import Path

import numpy as np

from finite_em_core import exact_periodic_cube_check, local_rk_check
from finite_em_flux import run_flux_ensemble
from finite_em_charge import run_charge_ensemble

def finite_size_summary(flux_runs: list[dict[str, object]]) -> dict[str, object]:
    by_spin: dict[int, list[dict[str, object]]] = collections.defaultdict(list)
    for run in flux_runs:
        by_spin[int(run["spin"])].append(run)
    result: dict[str, object] = {}
    for spin, runs in sorted(by_spin.items()):
        runs = sorted(runs, key=lambda row: int(row["lattice_size"]))
        Ls = np.asarray([float(row["lattice_size"]) for row in runs])
        variance_scaled = np.asarray(
            [float(row["flux_sector_response"]["variance_over_L"]) for row in runs]
        )
        scaled_gap = np.asarray(
            [
                float(row["flux_sector_response"].get("scaled_one_flux_free_energy_L_deltaF", np.nan))
                for row in runs
            ]
        )
        max_long = max(
            float(mode["longitudinal_power_fraction"])
            for row in runs
            for mode in row["transverse_correlations"]
        )
        axis_splits_by_run = [
            np.mean([
                float(mode["transverse_degeneracy_fractional_split"])
                for mode in row["transverse_correlations"]
                if mode["mode"] in ([1, 0, 0], [0, 1, 0], [0, 0, 1])
            ])
            for row in runs
        ]
        max_split = max(
            float(mode["transverse_degeneracy_fractional_split"])
            for row in runs
            for mode in row["transverse_correlations"]
            if mode["mode"] in ([1, 0, 0], [0, 1, 0], [0, 0, 1])
        )
        # A constant variance/L is the expected Coulomb winding law.  Record the
        # fractional range instead of forcing a noisy exponent from three sizes.
        frac_range = float((np.max(variance_scaled) - np.min(variance_scaled)) / np.mean(variance_scaled))
        valid_gaps = scaled_gap[np.isfinite(scaled_gap)]
        gap_frac_range = (
            float((np.max(valid_gaps) - np.min(valid_gaps)) / np.mean(valid_gaps))
            if len(valid_gaps) >= 2 and np.mean(valid_gaps) != 0
            else None
        )
        result[str(spin)] = {
            "lattice_sizes": [int(v) for v in Ls],
            "winding_variance_over_L": [float(v) for v in variance_scaled],
            "winding_variance_over_L_fractional_range": frac_range,
            "L_times_one_flux_free_energy": [float(v) for v in scaled_gap],
            "L_times_one_flux_free_energy_fractional_range": gap_frac_range,
            "maximum_longitudinal_power_fraction": max_long,
            "maximum_axis_transverse_fractional_split": max_split,
            "maximum_axis_mean_transverse_fractional_split": float(max(axis_splits_by_run)),
        }
    return result


def ensemble_run_parameters(quick: bool) -> tuple[list[int], int, int]:
    if quick:
        return [4, 6], 80, 280
    return [4, 6, 8], 120, 520


def run_all_ensembles(quick: bool, seed: int) -> dict[str, object]:
    sizes, thermal, samples = ensemble_run_parameters(quick)

    # The charge worm is the most allocation-sensitive stage, so run its largest
    # lattice first.  This ordering is deterministic and avoids allocator
    # fragmentation after thousands of small FFT temporaries.
    charge_specs = ([(sizes[-1], 60, 180)] if quick else [(8, 220, 1600), (6, 220, 1200)])
    charge_runs = [
        run_charge_ensemble(
            L=L,
            spin=1,
            thermal_sweeps=charge_thermal,
            sample_sweeps=charge_samples,
            seed=seed + 9000 + L,
        )
        for L, charge_thermal, charge_samples in charge_specs
    ]

    flux_runs: list[dict[str, object]] = []
    for spin in (1, 2):
        run_sizes = sizes if spin == 1 else sizes[-2:]
        for L in run_sizes:
            flux_runs.append(
                run_flux_ensemble(
                    L=L,
                    spin=spin,
                    thermal_sweeps=thermal,
                    sample_sweeps=samples,
                    seed=seed + 1000 * spin + L,
                )
            )
    return {"flux_runs": flux_runs, "charge_runs": charge_runs}


def charge_finite_size_summary(charge_runs: list[dict[str, object]]) -> dict[str, object]:
    runs = sorted(
        [run for run in charge_runs if int(run["spin"]) == 1],
        key=lambda row: int(row["lattice_size"]),
    )
    slopes = [float(run["potential_fit"].get("linear_string_slope", math.nan)) for run in runs]
    amplitudes = [float(run["potential_fit"].get("coulomb_amplitude_A", math.nan)) for run in runs]
    r2 = [float(run["potential_fit"].get("weighted_r_squared", math.nan)) for run in runs]
    decreasing = (
        len(slopes) >= 2
        and all(math.isfinite(value) for value in slopes)
        and slopes[-1] < 0.80 * slopes[0]
    )
    return {
        "lattice_sizes": [int(run["lattice_size"]) for run in runs],
        "coulomb_amplitudes": amplitudes,
        "coulomb_weighted_r_squared": r2,
        "effective_linear_string_slopes": slopes,
        "effective_string_slope_decreases_with_size": decreasing,
    }


def acceptance_report(
    local: list[dict[str, object]],
    flux_runs: list[dict[str, object]],
    charge_runs: list[dict[str, object]],
    quick: bool,
) -> dict[str, object]:
    finite = finite_size_summary(flux_runs)
    spin1 = finite["1"]
    local_pass = all(
        float(row["link_algebra_error"]) < 1e-12
        and float(row["gauss_commutator_max_abs"]) < 1e-12
        and float(row["rk_minimum_eigenvalue"]) > -1e-10
        and float(row["rk_row_sum_max_abs"]) < 1e-12
        for row in local
    )
    transverse_pass = (
        float(spin1["maximum_longitudinal_power_fraction"]) < 1e-12
        and float(spin1["maximum_axis_mean_transverse_fractional_split"]) < (0.30 if quick else 0.25)
    )
    winding_pass = float(spin1["winding_variance_over_L_fractional_range"]) < (0.80 if quick else 0.45)
    asymptotic_charge_runs = [
        run for run in charge_runs
        if int(run["spin"]) == 1 and int(run["lattice_size"]) >= 6
    ]
    charge_summary = charge_finite_size_summary(asymptotic_charge_runs)
    minimum_r2 = 0.03 if quick else 0.35
    charge_pass = bool(asymptotic_charge_runs) and all(
        bool(run["potential_fit"].get("fit_available"))
        and float(run["potential_fit"]["coulomb_amplitude_A"]) > 0
        and float(run["potential_fit"]["weighted_r_squared"]) > minimum_r2
        and bool(run["potential_fit"]["coulomb_preferred_to_constant"])
        for run in asymptotic_charge_runs
    )
    if not quick:
        charge_pass = charge_pass and bool(
            charge_summary["effective_string_slope_decreases_with_size"]
        )
    representation_pass = "2" in finite and float(finite["2"]["maximum_longitudinal_power_fraction"]) < 1e-12
    return {
        "exact_finite_link_algebra_and_positive_rk_parent": local_pass,
        "two_transverse_static_components_no_longitudinal_component": transverse_pass,
        "coulomb_winding_scaling": winding_pass,
        "opposite_charge_lattice_coulomb_potential": charge_pass,
        "charge_finite_size_summary": charge_summary,
        "larger_representation_static_universality_control": representation_pass,
        "relativistic_z1_photon_dynamics": False,
        "dynamical_charged_matter_and_qed_ward_identities": False,
        "issue_6_closed": False,
        "stage_pass": bool(
            local_pass
            and transverse_pass
            and winding_pass
            and charge_pass
            and representation_pass
        ),
    }


def build_report(args: argparse.Namespace) -> dict[str, object]:
    quick = bool(args.quick)
    # Run stochastic lattices before dense local LAPACK.  The exact 146,327-state
    # cube is intentionally a separate process/invocation (`--exact-only`) so its
    # large Python tuple arena cannot perturb the scaling run.
    ensembles = run_all_ensembles(quick=quick, seed=args.seed)
    local = [local_rk_check(spin) for spin in (1, 2)]
    flux_runs = ensembles["flux_runs"]
    charge_runs = ensembles["charge_runs"]
    acceptance = acceptance_report(local, flux_runs, charge_runs, quick=quick)
    report = {
        "status": (
            "finite spin-1 RK Coulomb ground-state stage; static deconfinement passes, "
            "relativistic quantum dynamics and interacting QED remain open"
        ),
        "model": {
            "links": "integer spin-S finite U(1) quantum links, E=S_z and U=S_+/max|S_+|",
            "gauss_law": "G_x=sum_out E-sum_in E-Q_x",
            "hamiltonian": (
                "H(u,t,v)=u/2 sum_l E_l^2 - t sum_p(W_p+W_p^dagger) "
                "+ v sum_p D_p, D_p=sqrt(W_p^dagger W_p)+sqrt(W_p W_p^dagger)"
            ),
            "rk_point": "u=0 and v=t; each local term is a weighted graph Laplacian",
            "rishon_number_role": "fixed by the selected spin-S irreducible link representation",
            "determinant_term": "not required for U(1); no extra U(N) determinant sector is present",
        },
        "run_mode": "quick" if quick else "full",
        "seed": args.seed,
        "local_operator_checks": local,
        "companion_exact_result_file": "finite_em_exact_cube_results.json",
        "flux_ensembles": flux_runs,
        "finite_size_summary": finite_size_summary(flux_runs),
        "charge_ensembles": charge_runs,
        "charge_finite_size_summary": charge_finite_size_summary(charge_runs),
        "acceptance": acceptance,
        "claim_boundary": {
            "established": [
                "exact finite U(1) Gauss symmetry",
                "positive finite spin-1 frustration-free 3+1D parent Hamiltonian",
                "Coulomb winding-sector scaling in the finite-link RK ground state",
                "rank-two transverse equal-time correlation tensor",
                "opposite-charge interaction fitted by the periodic lattice Green function",
                "static behavior persists in a larger finite representation",
            ],
            "not_established": [
                "a z=1 quantum photon away from the RK point",
                "1/L dynamical photon-gap scaling",
                "dynamical charged fermions",
                "vacuum polarization or charge renormalization",
                "interacting QED Ward identities",
                "QED precision bounds on Lorentz-violating cutoff operators",
            ],
        },
    }
    if not acceptance["stage_pass"]:
        raise AssertionError(json.dumps(acceptance, indent=2, sort_keys=True))
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--quick", action="store_true", help="smaller deterministic CI run")
    parser.add_argument("--exact-gap", action="store_true", help="Lanczos gap on the 146327-state cube")
    parser.add_argument("--exact-only", action="store_true", help="run only the exact 2x2x2 gauge-reduced cube")
    parser.add_argument("--flux-only", action="store_true", help="run the frozen full flux/transverse scaling suite")
    parser.add_argument("--charge-only", type=int, choices=(6, 8), help="run one frozen full charged-sector lattice")
    parser.add_argument("--seed", type=int, default=26092026)
    args = parser.parse_args()
    if args.exact_only:
        exact = exact_periodic_cube_check(L=2, spin=1, compute_gap=args.exact_gap)
        if int(exact["gauge_reduced_component_dimension"]) != 146327:
            raise AssertionError("unexpected exact gauge-sector dimension")
        if int(exact["sampled_gauss_residual_max_abs"]) != 0:
            raise AssertionError("Gauss law failed in the exact cube")
        text = json.dumps(exact, indent=2, sort_keys=True)
        print(text)
        if args.output:
            args.output.write_text(text + "\n", encoding="utf-8")
        return
    if args.flux_only:
        runs = []
        for spin in (1, 2):
            sizes = (4, 6, 8) if spin == 1 else (6, 8)
            for L in sizes:
                runs.append(
                    run_flux_ensemble(
                        L=L,
                        spin=spin,
                        thermal_sweeps=120,
                        sample_sweeps=520,
                        seed=args.seed + 1000 * spin + L,
                    )
                )
        summary = finite_size_summary(runs)
        local_checks = [local_rk_check(spin) for spin in (1, 2)]
        spin1 = summary["1"]
        spin2 = summary["2"]
        wilson_ok = all(
            bool(run["wilson_shift_overlap"]["perimeter_preferred"])
            for run in runs
        )
        flux_pass = (
            all(float(row["gauss_commutator_max_abs"]) < 1e-12 for row in local_checks)
            and float(spin1["winding_variance_over_L_fractional_range"]) < 0.15
            and float(spin1["L_times_one_flux_free_energy_fractional_range"]) < 0.30
            and float(spin1["maximum_longitudinal_power_fraction"]) < 1e-12
            and float(spin1["maximum_axis_mean_transverse_fractional_split"]) < 0.25
            and float(spin2["winding_variance_over_L_fractional_range"]) < 0.15
            and float(spin2["maximum_longitudinal_power_fraction"]) < 1e-12
            and wilson_ok
        )
        report = {
            "status": "full finite-link RK flux and transverse scaling run",
            "seed": args.seed,
            "local_operator_checks": local_checks,
            "flux_ensembles": runs,
            "finite_size_summary": summary,
            "acceptance": {"stage_pass": flux_pass, "wilson_perimeter_all_runs": wilson_ok},
        }
        if not flux_pass:
            raise AssertionError(json.dumps(report["acceptance"], indent=2, sort_keys=True))
        text = json.dumps(report, indent=2, sort_keys=True)
        print(text)
        if args.output:
            args.output.write_text(text + "\n", encoding="utf-8")
        return
    if args.charge_only is not None:
        L = int(args.charge_only)
        thermal, samples = ((220, 1200) if L == 6 else (220, 1600))
        report = run_charge_ensemble(
            L=L,
            spin=1,
            thermal_sweeps=thermal,
            sample_sweeps=samples,
            seed=args.seed + 9000 + L,
        )
        fit = report["potential_fit"]
        charge_pass = (
            int(report["max_gauss_minus_charge_residual"]) == 0
            and bool(fit.get("fit_available"))
            and float(fit["coulomb_amplitude_A"]) > 0
            and float(fit["weighted_r_squared"]) > 0.35
            and bool(fit["coulomb_preferred_to_constant"])
        )
        report["acceptance"] = {"stage_pass": charge_pass}
        if not charge_pass:
            raise AssertionError(json.dumps(report["acceptance"], indent=2, sort_keys=True))
        text = json.dumps(report, indent=2, sort_keys=True)
        print(text)
        if args.output:
            args.output.write_text(text + "\n", encoding="utf-8")
        return
    # The integrated invocation is the deterministic CI regression.  Full frozen
    # runs are split with --flux-only, --charge-only, and --exact-only so each
    # large stage starts in a fresh process.
    args.quick = True
    report = build_report(args)
    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
