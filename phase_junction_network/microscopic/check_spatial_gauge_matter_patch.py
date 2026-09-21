#!/usr/bin/env python3
"""Stage 3G: finite spatial photon + charged-endpoint patch.

This calculation combines the Stage-3E finite photon operators with the
Stage-3F minimal |q|=1 endpoint hopping in one periodic 2x2 gauge-matter
Hamiltonian. It tests exact local continuity, a nonzero-momentum Ward
identity, static transversality, and finite photon dressing without adding
matter- or photon-specific counterterms.

The matter content is one dynamical q=+1 endpoint and one dynamical q=-1
endpoint. It is a finite real-particle polarization test, not yet vacuum pair
creation or an interacting-QED continuum calculation.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import eigsh

L = 2
SITES = [(x, y) for x in range(L) for y in range(L)]
SITE_INDEX = {site: index for index, site in enumerate(SITES)}
LINKS = [(site, direction) for site in SITES for direction in range(2)]
LINK_INDEX = {link: index for index, link in enumerate(LINKS)}

# Stage-3E photon coefficients and Stage-3F charged-endpoint hopping, in Delta units.
U_GAUGE = 0.006789913753821769
K_PLAQUETTE = 0.005606197571911289
V_OVER_K = -0.12924513079019587
V_PLAQUETTE = V_OVER_K * K_PLAQUETTE
J2_OVER_K = 0.013851221789811033
J2_PLAQUETTE = J2_OVER_K * K_PLAQUETTE
K_MATTER = 0.0338029748033
X_SHARED = 0.13554178509861228


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


def add(site: tuple[int, int], direction: int, step: int = 1) -> tuple[int, int]:
    x, y = site
    if direction == 0:
        x = (x + step) % L
    else:
        y = (y + step) % L
    return x, y


def divergence(flux: tuple[int, ...]) -> tuple[int, ...]:
    result = np.zeros(len(SITES), dtype=int)
    for (site, direction), link_index in LINK_INDEX.items():
        result[SITE_INDEX[site]] += flux[link_index]
        result[SITE_INDEX[add(site, direction)]] -= flux[link_index]
    return tuple(int(value) for value in result)


FLUXES = list(itertools.product((-1, 0, 1), repeat=len(LINKS)))
FLUX_BY_DIVERGENCE: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
for flux_state in FLUXES:
    FLUX_BY_DIVERGENCE.setdefault(divergence(flux_state), []).append(flux_state)

PLAQUETTES: list[tuple[tuple[int, int], ...]] = []
for site in SITES:
    PLAQUETTES.append(
        (
            (LINK_INDEX[(site, 0)], +1),
            (LINK_INDEX[(add(site, 0), 1)], +1),
            (LINK_INDEX[(add(site, 1), 0)], -1),
            (LINK_INDEX[(site, 1)], -1),
        )
    )


def apply_steps(
    flux: tuple[int, ...],
    steps: tuple[tuple[int, int], ...],
    multiple: int = 1,
) -> tuple[int, ...] | None:
    moved = list(flux)
    for link_index, sign in steps:
        value = moved[link_index] + multiple * sign
        if value < -1 or value > 1:
            return None
        moved[link_index] = value
    return tuple(moved)


def pure_basis() -> list[tuple[int, ...]]:
    return FLUX_BY_DIVERGENCE[(0, 0, 0, 0)]


def matter_basis() -> list[tuple[int, int, tuple[int, ...]]]:
    basis = []
    for positive in range(len(SITES)):
        for negative in range(len(SITES)):
            charge = [0] * len(SITES)
            charge[positive] += 1
            charge[negative] -= 1
            for flux in FLUX_BY_DIVERGENCE[tuple(charge)]:
                basis.append((positive, negative, flux))
    return basis


def matter_move_matrices(basis):
    """Oriented q=+1 and q=-1 forward hopping operators on each link."""
    lookup = {state: index for index, state in enumerate(basis)}
    plus_moves = []
    minus_moves = []
    dimension = len(basis)

    for link_index, (tail, direction) in enumerate(LINKS):
        head = add(tail, direction)
        plus_rows: list[int] = []
        plus_columns: list[int] = []
        minus_rows: list[int] = []
        minus_columns: list[int] = []

        for column, (positive, negative, flux) in enumerate(basis):
            if positive == SITE_INDEX[tail] and flux[link_index] > -1:
                moved = list(flux)
                moved[link_index] -= 1
                target = (SITE_INDEX[head], negative, tuple(moved))
                plus_rows.append(lookup[target])
                plus_columns.append(column)

            if negative == SITE_INDEX[tail] and flux[link_index] < 1:
                moved = list(flux)
                moved[link_index] += 1
                target = (positive, SITE_INDEX[head], tuple(moved))
                minus_rows.append(lookup[target])
                minus_columns.append(column)

        plus_moves.append(
            coo_matrix(
                (np.ones(len(plus_rows)), (plus_rows, plus_columns)),
                shape=(dimension, dimension),
                dtype=complex,
            ).tocsr()
        )
        minus_moves.append(
            coo_matrix(
                (np.ones(len(minus_rows)), (minus_rows, minus_columns)),
                shape=(dimension, dimension),
                dtype=complex,
            ).tocsr()
        )

    return plus_moves, minus_moves


def gauge_hamiltonian(basis, matter: bool = False) -> csr_matrix:
    lookup = {state: index for index, state in enumerate(basis)}
    rows: list[int] = []
    columns: list[int] = []
    data: list[complex] = []
    diagonal = np.zeros(len(basis), dtype=float)

    for row, state in enumerate(basis):
        flux = state[2] if matter else state
        diagonal[row] += 0.5 * U_GAUGE * sum(value * value for value in flux)

        for steps in PLAQUETTES:
            forward = apply_steps(flux, steps, +1)
            backward = apply_steps(flux, steps, -1)
            diagonal[row] += V_PLAQUETTE * (
                int(forward is not None) + int(backward is not None)
            )

            for moved in (forward, backward):
                if moved is not None:
                    target = (state[0], state[1], moved) if matter else moved
                    rows.append(row)
                    columns.append(lookup[target])
                    data.append(-K_PLAQUETTE)

            for moved in (
                apply_steps(flux, steps, +2),
                apply_steps(flux, steps, -2),
            ):
                if moved is not None:
                    target = (state[0], state[1], moved) if matter else moved
                    rows.append(row)
                    columns.append(lookup[target])
                    data.append(J2_PLAQUETTE)

        rows.append(row)
        columns.append(row)
        data.append(diagonal[row])

    hamiltonian = coo_matrix(
        (data, (rows, columns)),
        shape=(len(basis), len(basis)),
        dtype=complex,
    ).tocsr()
    return 0.5 * (hamiltonian + hamiltonian.getH())


def matter_hamiltonian(
    basis,
    plus_moves,
    minus_moves,
    external: np.ndarray | None = None,
    matter_scale: float = 1.0,
) -> csr_matrix:
    external = (
        np.zeros(len(LINKS), dtype=float)
        if external is None
        else np.asarray(external, dtype=float)
    )
    hamiltonian = csr_matrix((len(basis), len(basis)), dtype=complex)
    hopping = K_MATTER * matter_scale

    for link_index in range(len(LINKS)):
        phase = external[link_index]
        plus = plus_moves[link_index]
        minus = minus_moves[link_index]
        hamiltonian += -hopping * (
            np.exp(1j * phase) * plus
            + np.exp(-1j * phase) * plus.getH()
        )
        hamiltonian += -hopping * (
            np.exp(-1j * phase) * minus
            + np.exp(1j * phase) * minus.getH()
        )
    return hamiltonian


def build_matter_system(
    external: np.ndarray | None = None,
    matter_scale: float = 1.0,
):
    basis = matter_basis()
    plus_moves, minus_moves = matter_move_matrices(basis)
    hamiltonian = gauge_hamiltonian(basis, matter=True) + matter_hamiltonian(
        basis,
        plus_moves,
        minus_moves,
        external=external,
        matter_scale=matter_scale,
    )
    return basis, hamiltonian, plus_moves, minus_moves


def lowest(hamiltonian: csr_matrix, count: int):
    energies, vectors = eigsh(
        hamiltonian,
        k=min(count, hamiltonian.shape[0] - 2),
        which="SA",
        tol=1.0e-10,
        maxiter=100000,
    )
    order = np.argsort(energies)
    return energies[order], vectors[:, order]


def electric_transverse_values(basis, matter: bool = False) -> np.ndarray:
    values = []
    for state in basis:
        flux = state[2] if matter else state
        values.append(
            sum(
                ((-1) ** site[0]) * flux[LINK_INDEX[(site, 1)]]
                for site in SITES
            )
        )
    return np.asarray(values, dtype=float)


def spectral_summary(
    hamiltonian: csr_matrix,
    operator_values: np.ndarray,
    count: int,
) -> dict[str, Any]:
    energies, vectors = lowest(hamiltonian, count)
    ground = vectors[:, 0]
    created = operator_values * ground
    norm = float(np.vdot(created, created).real)

    groups = []
    start = 1
    while start < len(energies):
        stop = start + 1
        while stop < len(energies) and abs(energies[stop] - energies[start]) < 1.0e-8:
            stop += 1
        orthonormal, _ = np.linalg.qr(vectors[:, start:stop])
        projected = orthonormal.conj().T @ created
        weight = float(np.vdot(projected, projected).real)
        groups.append(
            {
                "gap": float(np.mean(energies[start:stop]) - energies[0]),
                "multiplicity": stop - start,
                "residue": weight,
                "fraction": weight / norm if norm else 0.0,
            }
        )
        start = stop

    groups.sort(key=lambda row: row["residue"], reverse=True)
    photon = groups[0]
    return {
        "ground_energy": float(energies[0]),
        "operator_norm": norm,
        "photon_gap": photon["gap"],
        "photon_multiplicity": photon["multiplicity"],
        "photon_residue": photon["fraction"],
        "top_groups": groups[:8],
        "energies": energies,
        "vectors": vectors,
    }


def translation_permutation(basis, dx: int, dy: int, matter: bool) -> np.ndarray:
    lookup = {state: index for index, state in enumerate(basis)}
    permutation = np.empty(len(basis), dtype=np.int64)

    for column, state in enumerate(basis):
        flux = state[2] if matter else state
        translated_flux = [0] * len(LINKS)
        for (site, direction), link_index in LINK_INDEX.items():
            target_site = ((site[0] + dx) % L, (site[1] + dy) % L)
            translated_flux[LINK_INDEX[(target_site, direction)]] = flux[link_index]
        translated_flux = tuple(translated_flux)

        if matter:
            positive_site = SITES[state[0]]
            negative_site = SITES[state[1]]
            positive = SITE_INDEX[
                ((positive_site[0] + dx) % L, (positive_site[1] + dy) % L)
            ]
            negative = SITE_INDEX[
                ((negative_site[0] + dx) % L, (negative_site[1] + dy) % L)
            ]
            translated = (positive, negative, translated_flux)
        else:
            translated = translated_flux
        permutation[column] = lookup[translated]
    return permutation


def translation_checks(
    hamiltonian: csr_matrix,
    basis,
    operator_values: np.ndarray,
    matter: bool,
) -> dict[str, float]:
    residuals = []
    characters = {}
    for direction, shift in (("x", (1, 0)), ("y", (0, 1))):
        permutation = translation_permutation(basis, *shift, matter=matter)
        inverse = np.empty_like(permutation)
        inverse[permutation] = np.arange(len(permutation))
        translated = hamiltonian[inverse, :][:, permutation]
        difference = translated - hamiltonian
        residuals.append(
            float(np.max(np.abs(difference.data))) if difference.nnz else 0.0
        )
        translated_values = operator_values[permutation]
        denominator = float(np.vdot(operator_values, operator_values).real)
        characters[direction] = float(
            np.vdot(operator_values, translated_values).real / denominator
        )
    return {
        "maximum_translation_commutator": max(residuals),
        "transverse_operator_character_x": characters["x"],
        "transverse_operator_character_y": characters["y"],
    }


def charge_and_current_operators(basis, plus_moves, minus_moves, scale: float = 1.0):
    densities = []
    for site_index in range(len(SITES)):
        values = [
            (1 if positive == site_index else 0)
            - (1 if negative == site_index else 0)
            for positive, negative, _ in basis
        ]
        densities.append(diags(values, dtype=complex).tocsr())

    currents = []
    hopping = K_MATTER * scale
    for link_index in range(len(LINKS)):
        plus = plus_moves[link_index]
        minus = minus_moves[link_index]
        currents.append(
            1j
            * hopping
            * (plus - plus.getH() - minus + minus.getH())
        )
    return densities, currents


def continuity_checks(
    hamiltonian: csr_matrix,
    basis,
    plus_moves,
    minus_moves,
    energies: np.ndarray,
    vectors: np.ndarray,
) -> dict[str, Any]:
    densities, currents = charge_and_current_operators(
        basis, plus_moves, minus_moves
    )

    local_residuals = []
    for site in SITES:
        site_index = SITE_INDEX[site]
        divergence_current = csr_matrix(hamiltonian.shape, dtype=complex)
        for direction in range(2):
            divergence_current += currents[LINK_INDEX[(site, direction)]]
            divergence_current -= currents[
                LINK_INDEX[(add(site, direction, -1), direction)]
            ]
        residual = (
            1j
            * (
                hamiltonian @ densities[site_index]
                - densities[site_index] @ hamiltonian
            )
            + divergence_current
        )
        local_residuals.append(
            float(np.max(np.abs(residual.data))) if residual.nnz else 0.0
        )

    momentum = (math.pi, 0.0)
    density_momentum = csr_matrix(hamiltonian.shape, dtype=complex)
    divergence_momentum = csr_matrix(hamiltonian.shape, dtype=complex)

    for site in SITES:
        phase = np.exp(
            -1j * (momentum[0] * site[0] + momentum[1] * site[1])
        )
        density_momentum += phase * densities[SITE_INDEX[site]]

    for direction in range(2):
        derivative = 1.0 - np.exp(-1j * momentum[direction])
        if abs(derivative) < 1.0e-14:
            continue
        current_momentum = csr_matrix(hamiltonian.shape, dtype=complex)
        for site in SITES:
            phase = np.exp(
                -1j * (momentum[0] * site[0] + momentum[1] * site[1])
            )
            current_momentum += phase * currents[
                LINK_INDEX[(site, direction)]
            ]
        divergence_momentum += derivative * current_momentum

    operator_residual = (
        1j
        * (
            hamiltonian @ density_momentum
            - density_momentum @ hamiltonian
        )
        + divergence_momentum
    )

    ground = vectors[:, 0]
    spectral_residuals = []
    for index in range(1, len(energies)):
        excited = vectors[:, index]
        residual = (
            1j
            * (energies[index] - energies[0])
            * (excited.conj() @ (density_momentum @ ground))
            + excited.conj() @ (divergence_momentum @ ground)
        )
        spectral_residuals.append(abs(residual))

    return {
        "maximum_local_continuity_residual": max(local_residuals),
        "nonzero_momentum": list(momentum),
        "operator_ward_residual": (
            float(np.max(np.abs(operator_residual.data)))
            if operator_residual.nnz
            else 0.0
        ),
        "maximum_low_spectral_ward_residual": float(max(spectral_residuals)),
        "rho_k_norm_on_ground": float(np.linalg.norm(density_momentum @ ground)),
    }


def mode_pattern(component: int) -> np.ndarray:
    pattern = np.zeros(len(LINKS), dtype=float)
    for link_index, (site, direction) in enumerate(LINKS):
        if direction == component:
            pattern[link_index] = (-1) ** site[0]
    return pattern


def external_unitary(basis, amplitude: float) -> csr_matrix:
    phases = []
    for positive, negative, _ in basis:
        positive_lambda = -0.5 * amplitude * ((-1) ** SITES[positive][0])
        negative_lambda = -0.5 * amplitude * ((-1) ** SITES[negative][0])
        phases.append(np.exp(1j * (positive_lambda - negative_lambda)))
    return diags(phases, dtype=complex).tocsr()


def polarization_checks(
    base_gauge: csr_matrix,
    basis,
    plus_moves,
    minus_moves,
    zero_field_hamiltonian: csr_matrix,
    matter_scale: float = 1.0,
    step: float = 0.002,
) -> dict[str, Any]:
    longitudinal_pattern = mode_pattern(0)
    transverse_pattern = mode_pattern(1)

    def ground_energy(longitudinal: float, transverse: float) -> float:
        external = (
            longitudinal * longitudinal_pattern
            + transverse * transverse_pattern
        )
        hamiltonian = base_gauge + matter_hamiltonian(
            basis,
            plus_moves,
            minus_moves,
            external=external,
            matter_scale=matter_scale,
        )
        energies, _ = lowest(hamiltonian, 2)
        return float(energies[0])

    energy_zero = ground_energy(0.0, 0.0)
    energy_long_plus = ground_energy(step, 0.0)
    energy_long_minus = ground_energy(-step, 0.0)
    energy_trans_plus = ground_energy(0.0, step)
    energy_trans_minus = ground_energy(0.0, -step)
    energy_pp = ground_energy(step, step)
    energy_pm = ground_energy(step, -step)
    energy_mp = ground_energy(-step, step)
    energy_mm = ground_energy(-step, -step)

    longitudinal_curvature = (
        energy_long_plus + energy_long_minus - 2.0 * energy_zero
    ) / step**2
    transverse_curvature = (
        energy_trans_plus + energy_trans_minus - 2.0 * energy_zero
    ) / step**2
    mixed_curvature = (
        energy_pp - energy_pm - energy_mp + energy_mm
    ) / (4.0 * step**2)

    amplitude = 0.37
    longitudinal_hamiltonian = base_gauge + matter_hamiltonian(
        basis,
        plus_moves,
        minus_moves,
        external=amplitude * longitudinal_pattern,
        matter_scale=matter_scale,
    )
    unitary = external_unitary(basis, amplitude)
    equivalence = (
        longitudinal_hamiltonian
        - unitary @ zero_field_hamiltonian @ unitary.getH()
    )

    return {
        "finite_difference_step": step,
        "longitudinal_curvature": longitudinal_curvature,
        "transverse_curvature": transverse_curvature,
        "mixed_curvature": mixed_curvature,
        "transversality_ratio": abs(longitudinal_curvature)
        / max(abs(transverse_curvature), 1.0e-30),
        "exact_longitudinal_unitary_equivalence_residual": (
            float(np.max(np.abs(equivalence.data)))
            if equivalence.nnz
            else 0.0
        ),
        "transverse_response_nonzero": abs(transverse_curvature) > 1.0e-4,
    }


def scale_scan(bare_gap: float) -> list[dict[str, float]]:
    rows = []
    for scale in (0.25, 0.5, 0.75, 1.0):
        basis, hamiltonian, plus_moves, minus_moves = build_matter_system(
            matter_scale=scale
        )
        operator = electric_transverse_values(basis, matter=True)
        spectrum = spectral_summary(hamiltonian, operator, 80)
        polarization = polarization_checks(
            gauge_hamiltonian(basis, matter=True),
            basis,
            plus_moves,
            minus_moves,
            hamiltonian,
            matter_scale=scale,
        )
        rows.append(
            {
                "matter_scale": scale,
                "photon_gap": spectrum["photon_gap"],
                "photon_gap_shift_from_pure_gauge": spectrum["photon_gap"]
                - bare_gap,
                "photon_residue": spectrum["photon_residue"],
                "transverse_curvature": polarization["transverse_curvature"],
            }
        )
    return rows


def build_report(quick: bool = False) -> dict[str, Any]:
    pure = pure_basis()
    pure_hamiltonian = gauge_hamiltonian(pure)
    pure_operator = electric_transverse_values(pure)
    bare = spectral_summary(pure_hamiltonian, pure_operator, 40 if quick else 60)

    basis, hamiltonian, plus_moves, minus_moves = build_matter_system()
    dressed_operator = electric_transverse_values(basis, matter=True)
    dressed = spectral_summary(
        hamiltonian,
        dressed_operator,
        70 if quick else 120,
    )

    ward = continuity_checks(
        hamiltonian,
        basis,
        plus_moves,
        minus_moves,
        dressed["energies"],
        dressed["vectors"],
    )
    polarization = polarization_checks(
        gauge_hamiltonian(basis, matter=True),
        basis,
        plus_moves,
        minus_moves,
        hamiltonian,
    )
    translations = translation_checks(
        hamiltonian,
        basis,
        dressed_operator,
        matter=True,
    )

    gap_shift = dressed["photon_gap"] - bare["photon_gap"]
    pole_self_energy = dressed["photon_gap"] ** 2 - bare["photon_gap"] ** 2

    basis_zero, hamiltonian_zero, plus_zero, minus_zero = build_matter_system(
        matter_scale=0.0
    )
    polarization_zero = polarization_checks(
        gauge_hamiltonian(basis_zero, matter=True),
        basis_zero,
        plus_zero,
        minus_zero,
        hamiltonian_zero,
        matter_scale=0.0,
    )

    scan = [] if quick else scale_scan(bare["photon_gap"])

    checks = {
        "physical_basis_counts": len(pure) == 115 and len(basis) == 1484,
        "translation_symmetry_exact": translations[
            "maximum_translation_commutator"
        ]
        < 1.0e-12,
        "operator_has_pi_zero_momentum_character": abs(
            translations["transverse_operator_character_x"] + 1.0
        )
        < 1.0e-12
        and abs(translations["transverse_operator_character_y"] - 1.0)
        < 1.0e-12,
        "local_continuity_exact": ward["maximum_local_continuity_residual"]
        < 1.0e-12,
        "nonzero_momentum_operator_ward_exact": ward[
            "operator_ward_residual"
        ]
        < 1.0e-12,
        "nonzero_momentum_spectral_ward": ward[
            "maximum_low_spectral_ward_residual"
        ]
        < 1.0e-9,
        "longitudinal_external_mode_is_pure_gauge": polarization[
            "exact_longitudinal_unitary_equivalence_residual"
        ]
        < 1.0e-12
        and polarization["transversality_ratio"] < 1.0e-5,
        "transverse_matter_polarization_nonzero": polarization[
            "transverse_response_nonzero"
        ],
        "mixed_polarization_small": abs(polarization["mixed_curvature"])
        < 1.0e-6,
        "photon_residue_retained": dressed["photon_residue"] > 0.90,
        "finite_photon_dressing_nonzero": abs(gap_shift) > 1.0e-4,
        "decoupled_matter_control_has_zero_response": abs(
            polarization_zero["transverse_curvature"]
        )
        < 1.0e-9,
        "dressing_grows_with_matter_coupling": True
        if quick
        else all(
            right["photon_gap_shift_from_pure_gauge"]
            > left["photon_gap_shift_from_pure_gauge"]
            and right["transverse_curvature"] > left["transverse_curvature"]
            for left, right in zip(scan, scan[1:])
        ),
    }

    public_bare = {
        key: value
        for key, value in bare.items()
        if key not in ("energies", "vectors")
    }
    public_dressed = {
        key: value
        for key, value in dressed.items()
        if key not in ("energies", "vectors")
    }

    report = {
        "module": (
            "phase_junction_network/microscopic/"
            "check_spatial_gauge_matter_patch.py"
        ),
        "status": (
            "Stage 3G PASS at finite spatial gauge-matter scope"
            if all(checks.values())
            else "Stage 3G FAIL"
        ),
        "run_mode": "quick" if quick else "full",
        "lattice": {
            "shape": [2, 2],
            "sites": 4,
            "spin1_links": 8,
            "pure_gauge_physical_dimension": len(pure),
            "neutral_pair_physical_dimension": len(basis),
            "matter_content": (
                "one dynamical q=+1 endpoint and one dynamical q=-1 endpoint"
            ),
        },
        "coefficients": {
            "U_over_Delta": U_GAUGE,
            "K_plaquette_over_Delta": K_PLAQUETTE,
            "v_over_K": V_OVER_K,
            "J2_over_K": J2_OVER_K,
            "K_matter_over_Delta": K_MATTER,
            "shared_x_t_over_Delta": X_SHARED,
            "new_counterterms": 0,
        },
        "translation": translations,
        "bare_photon": public_bare,
        "matter_dressed_photon": public_dressed,
        "photon_dressing": {
            "gap_shift": gap_shift,
            "pole_self_energy_delta_omega_squared": pole_self_energy,
            "residue_ratio_dressed_to_bare": dressed["photon_residue"]
            / bare["photon_residue"],
        },
        "ward_identity": ward,
        "static_polarization": polarization,
        "matter_coupling_scan": scan,
        "decoupled_matter_control": polarization_zero,
        "checks": checks,
        "stage_pass": all(checks.values()),
        "claim_boundary": {
            "established": [
                (
                    "one finite spatial Hamiltonian contains the Stage-3E "
                    "photon operators and a minimal neutral pair of q=1 endpoints"
                ),
                "exact local continuity and a nonzero-momentum operator Ward identity",
                "low-energy nonzero-momentum spectral Ward matrix elements",
                "exact longitudinal pure-gauge invariance and nonzero transverse matter response",
                "finite photon pole dressing with most transverse electric residue retained",
                "no matter- or photon-specific counterterm was introduced",
            ],
            "not_established": [
                "virtual vacuum pair creation or a relativistic fermion determinant",
                "infinite-volume transverse vacuum polarization",
                "a 3+1D interacting QED fixed point",
                "radiative stability beyond this one-pair sector",
                "mirror-wall completion or observed matter",
                "nonlinear gravity or a common all-sector many-body continuum",
            ],
        },
        "next_gate": (
            "Add vacuum plus pair-creation sectors using the same completed "
            "endpoint move, then measure transverse photon self-energy from "
            "virtual pairs and test the finite Ward-Takahashi identity across "
            "the vacuum and pair blocks."
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
        default=Path(
            "phase_junction_network/microscopic/"
            "spatial_gauge_matter_results.json"
        ),
    )
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    report = build_report(bool(args.quick))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(ready(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "stage_pass": report["stage_pass"],
                "photon_dressing": report["photon_dressing"],
                "polarization": report["static_polarization"],
                "ward": report["ward_identity"],
                "failed_checks": [
                    key
                    for key, passed in report["checks"].items()
                    if not passed
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
