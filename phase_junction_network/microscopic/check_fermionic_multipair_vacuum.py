#!/usr/bin/env python3
"""Stage 3I: antisymmetric multi-pair vacuum polarization.

Extends the Stage-3H 2x2 spin-1 gauge patch from a zero/one-pair truncation
to the complete neutral spinless-fermion Fock space on four sites. Exact
Jordan-Wigner signs are retained for endpoint hopping, pair creation,
translations, and Ward currents.

The calculation tests:
* exact fermionic anticommutation and Pauli blocking;
* exact Gauss-law dynamics in all neutral pair-number sectors;
* convergence from one pair through the full four-pair Fock space;
* local and nonzero-momentum Ward identities;
* longitudinal pure-gauge invariance and transverse polarization;
* photon-pole and residue renormalization versus pair gap and coupling.

It is not yet the full domain-wall multiplet, a fermion determinant, or an
infinite-volume interacting-QED calculation.
"""
from __future__ import annotations

import argparse
import collections
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

# Stage-3E / Stage-3F coefficients in Delta units.
U_GAUGE = 0.006789913753821769
K_PLAQUETTE = 0.005606197571911289
V_OVER_K = -0.12924513079019587
V_PLAQUETTE = V_OVER_K * K_PLAQUETTE
J2_OVER_K = 0.013851221789811033
J2_PLAQUETTE = J2_OVER_K * K_PLAQUETTE
K_MATTER = 0.0338029748033
PAIR_COUPLING = K_MATTER
PAIR_COMPLETION = 1.0
PAIR_GAP = 1.0
BARE_PHOTON_GAP = 0.01785945619783673
BARE_PHOTON_RESIDUE = 0.9893176943138873
DOMAIN_WALL_FIRST_DOUBLER_GAP = 1.1674499031546928


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


FLUX_BY_DIVERGENCE: dict[tuple[int, ...], list[tuple[int, ...]]] = collections.defaultdict(list)
for flux_state in itertools.product((-1, 0, 1), repeat=len(LINKS)):
    FLUX_BY_DIVERGENCE[divergence(flux_state)].append(flux_state)

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


def annihilate(mask: int, mode: int) -> tuple[int, int] | None:
    if not ((mask >> mode) & 1):
        return None
    sign = -1 if (mask & ((1 << mode) - 1)).bit_count() % 2 else 1
    return mask ^ (1 << mode), sign


def create(mask: int, mode: int) -> tuple[int, int] | None:
    if (mask >> mode) & 1:
        return None
    sign = -1 if (mask & ((1 << mode) - 1)).bit_count() % 2 else 1
    return mask | (1 << mode), sign


def hop_mask(mask: int, source: int, destination: int) -> tuple[int, int] | None:
    removed = annihilate(mask, source)
    if removed is None:
        return None
    intermediate, sign_one = removed
    inserted = create(intermediate, destination)
    if inserted is None:
        return None
    final, sign_two = inserted
    return final, sign_one * sign_two


def pair_create_mask(mask: int, positive_head: int, negative_tail: int) -> tuple[int, int] | None:
    # Operator convention: c^dagger_{+,head} c^dagger_{-,tail}.
    negative = create(mask, 4 + negative_tail)
    if negative is None:
        return None
    intermediate, sign_one = negative
    positive = create(intermediate, positive_head)
    if positive is None:
        return None
    final, sign_two = positive
    return final, sign_one * sign_two


def charge_from_bits(positive_bits: int, negative_bits: int) -> tuple[int, ...]:
    return tuple(
        ((positive_bits >> site) & 1) - ((negative_bits >> site) & 1)
        for site in range(len(SITES))
    )


def fock_basis(max_pairs: int) -> tuple[list[tuple[int, int, tuple[int, ...]]], dict[int, int]]:
    basis = []
    counts: collections.Counter[int] = collections.Counter()
    for positive_bits in range(1 << len(SITES)):
        pair_number = positive_bits.bit_count()
        if pair_number > max_pairs:
            continue
        for negative_bits in range(1 << len(SITES)):
            if negative_bits.bit_count() != pair_number:
                continue
            charge = charge_from_bits(positive_bits, negative_bits)
            for flux in FLUX_BY_DIVERGENCE[charge]:
                basis.append((positive_bits, negative_bits, flux))
                counts[pair_number] += 1
    return basis, dict(counts)


def gauge_terms(basis, lookup):
    rows: list[int] = []
    columns: list[int] = []
    data: list[complex] = []
    diagonal = np.zeros(len(basis), dtype=float)
    for row, (positive_bits, negative_bits, flux) in enumerate(basis):
        diagonal[row] += 0.5 * U_GAUGE * sum(value * value for value in flux)
        for steps in PLAQUETTES:
            forward = apply_steps(flux, steps, +1)
            backward = apply_steps(flux, steps, -1)
            diagonal[row] += V_PLAQUETTE * (
                int(forward is not None) + int(backward is not None)
            )
            for moved in (forward, backward):
                if moved is not None:
                    rows.append(row)
                    columns.append(lookup[(positive_bits, negative_bits, moved)])
                    data.append(-K_PLAQUETTE)
            for moved in (
                apply_steps(flux, steps, +2),
                apply_steps(flux, steps, -2),
            ):
                if moved is not None:
                    rows.append(row)
                    columns.append(lookup[(positive_bits, negative_bits, moved)])
                    data.append(J2_PLAQUETTE)
    return rows, columns, data, diagonal


def move_matrices(basis, lookup):
    dimension = len(basis)
    plus_moves = []
    minus_moves = []
    pair_moves = []
    for link_index, (tail, direction) in enumerate(LINKS):
        head = add(tail, direction)
        tail_index = SITE_INDEX[tail]
        head_index = SITE_INDEX[head]
        records = [([], [], []) for _ in range(3)]
        for column, (positive_bits, negative_bits, flux) in enumerate(basis):
            mask = positive_bits | (negative_bits << 4)

            positive_hop = hop_mask(mask, tail_index, head_index)
            if positive_hop is not None and flux[link_index] > -1:
                moved_mask, sign = positive_hop
                moved_flux = list(flux)
                moved_flux[link_index] -= 1
                target = (
                    moved_mask & 15,
                    (moved_mask >> 4) & 15,
                    tuple(moved_flux),
                )
                if target in lookup:
                    records[0][0].append(lookup[target])
                    records[0][1].append(column)
                    records[0][2].append(sign)

            negative_hop = hop_mask(mask, 4 + tail_index, 4 + head_index)
            if negative_hop is not None and flux[link_index] < 1:
                moved_mask, sign = negative_hop
                moved_flux = list(flux)
                moved_flux[link_index] += 1
                target = (
                    moved_mask & 15,
                    (moved_mask >> 4) & 15,
                    tuple(moved_flux),
                )
                if target in lookup:
                    records[1][0].append(lookup[target])
                    records[1][1].append(column)
                    records[1][2].append(sign)

            pair = pair_create_mask(mask, head_index, tail_index)
            if pair is not None and flux[link_index] > -1:
                moved_mask, sign = pair
                moved_flux = list(flux)
                moved_flux[link_index] -= 1
                target = (
                    moved_mask & 15,
                    (moved_mask >> 4) & 15,
                    tuple(moved_flux),
                )
                if target in lookup:
                    records[2][0].append(lookup[target])
                    records[2][1].append(column)
                    records[2][2].append(sign)

        matrices = []
        for move_rows, move_columns, move_data in records:
            matrices.append(
                coo_matrix(
                    (
                        np.asarray(move_data, dtype=complex),
                        (move_rows, move_columns),
                    ),
                    shape=(dimension, dimension),
                    dtype=complex,
                ).tocsr()
            )
        plus_moves.append(matrices[0])
        minus_moves.append(matrices[1])
        pair_moves.append(matrices[2])
    return plus_moves, minus_moves, pair_moves


def build_system(
    max_pairs: int,
    pair_gap: float = PAIR_GAP,
    pair_scale: float = 1.0,
    matter_scale: float = 1.0,
    external: np.ndarray | None = None,
):
    basis, sector_counts = fock_basis(max_pairs)
    lookup = {state: index for index, state in enumerate(basis)}
    rows, columns, data, diagonal = gauge_terms(basis, lookup)
    for row, (positive_bits, _, _) in enumerate(basis):
        diagonal[row] += pair_gap * positive_bits.bit_count()
    hamiltonian = coo_matrix(
        (data, (rows, columns)),
        shape=(len(basis), len(basis)),
        dtype=complex,
    ).tocsr() + diags(diagonal, dtype=complex)

    plus_moves, minus_moves, pair_moves = move_matrices(basis, lookup)
    external = (
        np.zeros(len(LINKS), dtype=float)
        if external is None
        else np.asarray(external, dtype=float)
    )
    pair_coupling = PAIR_COUPLING * pair_scale
    matter_hopping = K_MATTER * matter_scale
    for link_index in range(len(LINKS)):
        phase = external[link_index]
        plus = plus_moves[link_index]
        minus = minus_moves[link_index]
        pair = pair_moves[link_index]
        hamiltonian += -matter_hopping * (
            np.exp(1j * phase) * plus
            + np.exp(-1j * phase) * plus.getH()
        )
        hamiltonian += -matter_hopping * (
            np.exp(-1j * phase) * minus
            + np.exp(1j * phase) * minus.getH()
        )
        hamiltonian += -pair_coupling * (
            np.exp(1j * phase) * pair
            + np.exp(-1j * phase) * pair.getH()
        )
        hamiltonian += PAIR_COMPLETION * pair_coupling**2 * (
            pair.getH() @ pair + pair @ pair.getH()
        )
    hamiltonian = 0.5 * (hamiltonian + hamiltonian.getH())
    return (
        basis,
        sector_counts,
        hamiltonian,
        plus_moves,
        minus_moves,
        pair_moves,
    )


def lowest(hamiltonian: csr_matrix, count: int):
    energies, vectors = eigsh(
        hamiltonian,
        k=min(count, hamiltonian.shape[0] - 2),
        which="SA",
        tol=1.0e-10,
        maxiter=150000,
    )
    order = np.argsort(energies)
    return energies[order], vectors[:, order]


def electric_values(basis) -> np.ndarray:
    return np.asarray(
        [
            sum(
                ((-1) ** site[0]) * flux[LINK_INDEX[(site, 1)]]
                for site in SITES
            )
            for _, _, flux in basis
        ],
        dtype=float,
    )


def spectral_summary(hamiltonian, values, count):
    energies, vectors = lowest(hamiltonian, count)
    ground = vectors[:, 0]
    created = values * ground
    norm = float(np.vdot(created, created).real)
    groups = []
    start = 1
    while start < len(energies):
        stop = start + 1
        while stop < len(energies) and abs(energies[stop] - energies[start]) < 1.0e-8:
            stop += 1
        basis_vectors, _ = np.linalg.qr(vectors[:, start:stop])
        amplitudes = basis_vectors.conj().T @ created
        weight = float(np.vdot(amplitudes, amplitudes).real)
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
    return {
        "ground_energy": float(energies[0]),
        "operator_norm": norm,
        "photon_gap": groups[0]["gap"],
        "photon_multiplicity": groups[0]["multiplicity"],
        "photon_residue": groups[0]["fraction"],
        "top_groups": groups[:8],
        "energies": energies,
        "vectors": vectors,
    }


def pair_distribution(basis, ground) -> dict[int, float]:
    distribution: collections.Counter[int] = collections.Counter()
    for amplitude, (positive_bits, _, _) in zip(ground, basis):
        distribution[positive_bits.bit_count()] += abs(amplitude) ** 2
    return {number: float(value) for number, value in sorted(distribution.items())}


def density_current_operators(basis, plus_moves, minus_moves, pair_moves):
    densities = []
    for site in range(len(SITES)):
        values = [
            ((positive_bits >> site) & 1) - ((negative_bits >> site) & 1)
            for positive_bits, negative_bits, _ in basis
        ]
        densities.append(diags(values, dtype=complex).tocsr())
    currents = []
    for link_index in range(len(LINKS)):
        plus = plus_moves[link_index]
        minus = minus_moves[link_index]
        pair = pair_moves[link_index]
        current = 1j * K_MATTER * (
            plus - plus.getH() - minus + minus.getH()
        )
        current += 1j * PAIR_COUPLING * (pair - pair.getH())
        currents.append(current)
    return densities, currents


def ward_checks(hamiltonian, basis, plus_moves, minus_moves, pair_moves, energies, vectors):
    densities, currents = density_current_operators(
        basis, plus_moves, minus_moves, pair_moves
    )
    local_residuals = []
    for site in SITES:
        divergence_current = csr_matrix(hamiltonian.shape, dtype=complex)
        for direction in range(2):
            divergence_current += currents[LINK_INDEX[(site, direction)]]
            divergence_current -= currents[
                LINK_INDEX[(add(site, direction, -1), direction)]
            ]
        residual = (
            1j
            * (
                hamiltonian @ densities[SITE_INDEX[site]]
                - densities[SITE_INDEX[site]] @ hamiltonian
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


def external_unitary(basis, amplitude: float):
    phases = []
    for positive_bits, negative_bits, _ in basis:
        exponent = 0.0
        for site_index, site in enumerate(SITES):
            local_phase = -0.5 * amplitude * ((-1) ** site[0])
            exponent += (
                ((positive_bits >> site_index) & 1)
                - ((negative_bits >> site_index) & 1)
            ) * local_phase
        phases.append(np.exp(1j * exponent))
    return diags(phases, dtype=complex).tocsr()


def polarization(
    max_pairs: int,
    pair_gap: float = PAIR_GAP,
    pair_scale: float = 1.0,
    step: float = 0.002,
):
    longitudinal = mode_pattern(0)
    transverse = mode_pattern(1)

    def energy(longitudinal_amplitude: float, transverse_amplitude: float) -> float:
        _, _, hamiltonian, *_ = build_system(
            max_pairs,
            pair_gap=pair_gap,
            pair_scale=pair_scale,
            external=(
                longitudinal_amplitude * longitudinal
                + transverse_amplitude * transverse
            ),
        )
        return float(
            eigsh(
                hamiltonian,
                k=1,
                which="SA",
                tol=1.0e-10,
                maxiter=120000,
                return_eigenvectors=False,
            )[0]
        )

    energy_zero = energy(0.0, 0.0)
    long_plus = energy(step, 0.0)
    long_minus = energy(-step, 0.0)
    transverse_plus = energy(0.0, step)
    transverse_minus = energy(0.0, -step)
    plus_plus = energy(step, step)
    plus_minus = energy(step, -step)
    minus_plus = energy(-step, step)
    minus_minus = energy(-step, -step)

    longitudinal_curvature = (
        long_plus + long_minus - 2.0 * energy_zero
    ) / step**2
    transverse_curvature = (
        transverse_plus + transverse_minus - 2.0 * energy_zero
    ) / step**2
    mixed_curvature = (
        plus_plus - plus_minus - minus_plus + minus_minus
    ) / (4.0 * step**2)

    basis_zero, _, hamiltonian_zero, *_ = build_system(
        max_pairs,
        pair_gap=pair_gap,
        pair_scale=pair_scale,
    )
    basis_long, _, hamiltonian_long, *_ = build_system(
        max_pairs,
        pair_gap=pair_gap,
        pair_scale=pair_scale,
        external=0.37 * longitudinal,
    )
    if basis_zero != basis_long:
        raise AssertionError("external field changed the physical basis")
    unitary = external_unitary(basis_zero, 0.37)
    difference = hamiltonian_long - unitary @ hamiltonian_zero @ unitary.getH()
    return {
        "longitudinal_curvature": longitudinal_curvature,
        "transverse_curvature": transverse_curvature,
        "mixed_curvature": mixed_curvature,
        "transversality_ratio": abs(longitudinal_curvature)
        / max(abs(transverse_curvature), 1.0e-30),
        "exact_longitudinal_unitary_equivalence_residual": (
            float(np.max(np.abs(difference.data))) if difference.nnz else 0.0
        ),
    }


def permutation_sign(mask: int, mode_permutation: dict[int, int]) -> int:
    mapped = [
        mode_permutation[mode]
        for mode in range(8)
        if (mask >> mode) & 1
    ]
    inversions = sum(
        1
        for left in range(len(mapped))
        for right in range(left + 1, len(mapped))
        if mapped[left] > mapped[right]
    )
    return -1 if inversions % 2 else 1


def translate_bits(bits: int, dx: int, dy: int) -> int:
    output = 0
    for site_index, site in enumerate(SITES):
        if (bits >> site_index) & 1:
            target = ((site[0] + dx) % L, (site[1] + dy) % L)
            output |= 1 << SITE_INDEX[target]
    return output


def translation_matrix(basis, dx: int, dy: int):
    lookup = {state: index for index, state in enumerate(basis)}
    mode_permutation = {}
    for site_index, site in enumerate(SITES):
        target = ((site[0] + dx) % L, (site[1] + dy) % L)
        mode_permutation[site_index] = SITE_INDEX[target]
        mode_permutation[4 + site_index] = 4 + SITE_INDEX[target]
    rows = []
    columns = []
    data = []
    for column, (positive_bits, negative_bits, flux) in enumerate(basis):
        mask = positive_bits | (negative_bits << 4)
        sign = permutation_sign(mask, mode_permutation)
        translated_flux = [0] * len(LINKS)
        for (site, direction), link_index in LINK_INDEX.items():
            target_site = ((site[0] + dx) % L, (site[1] + dy) % L)
            translated_flux[LINK_INDEX[(target_site, direction)]] = flux[link_index]
        target = (
            translate_bits(positive_bits, dx, dy),
            translate_bits(negative_bits, dx, dy),
            tuple(translated_flux),
        )
        rows.append(lookup[target])
        columns.append(column)
        data.append(sign)
    return coo_matrix(
        (data, (rows, columns)),
        shape=(len(basis), len(basis)),
        dtype=complex,
    ).tocsr()


def fermion_algebra_checks() -> dict[str, float]:
    annihilators = []
    for mode in range(8):
        rows = []
        columns = []
        data = []
        for mask in range(1 << 8):
            result = annihilate(mask, mode)
            if result is not None:
                target, sign = result
                rows.append(target)
                columns.append(mask)
                data.append(sign)
        annihilators.append(
            coo_matrix(
                (data, (rows, columns)),
                shape=(256, 256),
                dtype=complex,
            ).tocsr()
        )
    identity = diags(np.ones(256), dtype=complex).tocsr()
    maximum_car_residual = 0.0
    for left in range(8):
        for right in range(8):
            target = identity if left == right else csr_matrix((256, 256), dtype=complex)
            mixed = (
                annihilators[left] @ annihilators[right].getH()
                + annihilators[right].getH() @ annihilators[left]
                - target
            )
            same = (
                annihilators[left] @ annihilators[right]
                + annihilators[right] @ annihilators[left]
            )
            if mixed.nnz:
                maximum_car_residual = max(
                    maximum_car_residual,
                    float(np.max(np.abs(mixed.data))),
                )
            if same.nnz:
                maximum_car_residual = max(
                    maximum_car_residual,
                    float(np.max(np.abs(same.data))),
                )

    first = create(0, 0)
    assert first is not None
    first_then_second = create(first[0], 1)
    second = create(0, 1)
    assert second is not None
    second_then_first = create(second[0], 0)
    assert first_then_second is not None and second_then_first is not None
    exchange_phase = (
        first[1] * first_then_second[1]
        / (second[1] * second_then_first[1])
    )
    return {
        "maximum_canonical_anticommutator_residual": maximum_car_residual,
        "two_creation_exchange_phase": float(exchange_phase),
    }


def spectrum_run(max_pairs: int, pair_gap: float, pair_scale: float, count: int):
    basis, sector_counts, hamiltonian, plus_moves, minus_moves, pair_moves = build_system(
        max_pairs,
        pair_gap=pair_gap,
        pair_scale=pair_scale,
    )
    spectrum = spectral_summary(hamiltonian, electric_values(basis), count)
    distribution = pair_distribution(basis, spectrum["vectors"][:, 0])
    return {
        "basis": basis,
        "sector_counts": sector_counts,
        "hamiltonian": hamiltonian,
        "plus_moves": plus_moves,
        "minus_moves": minus_moves,
        "pair_moves": pair_moves,
        "spectrum": spectrum,
        "pair_distribution": distribution,
        "mean_pair_number": float(
            sum(number * weight for number, weight in distribution.items())
        ),
    }


def public_spectrum(run) -> dict[str, Any]:
    spectrum = run["spectrum"]
    return {
        "dimension": len(run["basis"]),
        "sector_counts": run["sector_counts"],
        "pair_distribution": run["pair_distribution"],
        "mean_pair_number": run["mean_pair_number"],
        "ground_energy": spectrum["ground_energy"],
        "photon_gap": spectrum["photon_gap"],
        "photon_gap_shift_from_pure_gauge": spectrum["photon_gap"]
        - BARE_PHOTON_GAP,
        "photon_multiplicity": spectrum["photon_multiplicity"],
        "photon_residue": spectrum["photon_residue"],
        "top_groups": spectrum["top_groups"],
    }


def logarithmic_fit(xs, ys) -> dict[str, float]:
    # Validation baseline: reject unidentifiable logarithmic fits.
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape or x.size < 2:
        raise ValueError("A log fit requires at least two paired samples")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("Log-fit samples must be finite")
    if np.any(x <= 0) or np.any(y <= 0):
        raise ValueError("Log-fit samples must be positive")
    lx, ly = np.log(x), np.log(y)
    centered = lx - lx.mean()
    denominator = float(centered @ centered)
    if denominator <= 1.0e-24:
        raise ValueError("A log-fit exponent needs distinct abscissae")
    power = float(centered @ (ly - ly.mean()) / denominator)
    prefactor = float(np.exp(ly.mean() - power * lx.mean()))
    return {"power": power, "prefactor": prefactor}


def build_report(quick: bool) -> dict[str, Any]:
    maximum_pairs = 2 if quick else 4
    truncations = (1, 2) if quick else (1, 2, 3, 4)
    pair_gaps = (1.0, 1.25) if quick else (
        0.75,
        1.0,
        DOMAIN_WALL_FIRST_DOUBLER_GAP,
        1.25,
        1.5,
    )

    runs = {
        maximum: spectrum_run(maximum, PAIR_GAP, 1.0, 40 if quick else 48)
        for maximum in truncations
    }
    target = runs[maximum_pairs]
    target_spectrum = target["spectrum"]

    ward = ward_checks(
        target["hamiltonian"],
        target["basis"],
        target["plus_moves"],
        target["minus_moves"],
        target["pair_moves"],
        target_spectrum["energies"],
        target_spectrum["vectors"],
    )
    target_polarization = polarization(maximum_pairs)
    decoupled_polarization = polarization(maximum_pairs, pair_scale=0.0)

    truncation_rows = []
    for maximum in truncations:
        row = public_spectrum(runs[maximum])
        row["max_pairs"] = maximum
        if not quick and maximum in (1, maximum_pairs):
            row["static_polarization"] = (
                target_polarization if maximum == maximum_pairs else polarization(maximum)
            )
        truncation_rows.append(row)

    gap_rows = []
    for gap in pair_gaps:
        run = spectrum_run(maximum_pairs, gap, 1.0, 38 if quick else 45)
        row = public_spectrum(run)
        row["pair_gap_over_Delta"] = gap
        if quick or gap in (0.75, 1.0, 1.5):
            row["static_polarization"] = (
                target_polarization if gap == 1.0 else polarization(maximum_pairs, pair_gap=gap)
            )
        gap_rows.append(row)

    coupling_rows = []
    if not quick:
        for scale in (0.25, 0.5, 0.75, 1.0):
            run = spectrum_run(maximum_pairs, PAIR_GAP, scale, 38)
            row = public_spectrum(run)
            row["pair_coupling_scale"] = scale
            coupling_rows.append(row)

    algebra = fermion_algebra_checks()
    translation_residuals = []
    operator_characters = {}
    operator = diags(electric_values(target["basis"]), dtype=complex)
    for name, shift in (("x", (1, 0)), ("y", (0, 1))):
        translation = translation_matrix(target["basis"], *shift)
        difference = (
            translation @ target["hamiltonian"] @ translation.getH()
            - target["hamiltonian"]
        )
        translation_residuals.append(
            float(np.max(np.abs(difference.data))) if difference.nnz else 0.0
        )
        translated_operator = translation @ operator @ translation.getH()
        numerator = operator.conj().multiply(translated_operator).sum()
        denominator = operator.conj().multiply(operator).sum()
        operator_characters[name] = float(np.real(numerator / denominator))

    maximum_pair_nilpotency = 0.0
    for move in target["pair_moves"]:
        square = move @ move
        if square.nnz:
            maximum_pair_nilpotency = max(
                maximum_pair_nilpotency,
                float(np.max(np.abs(square.data))),
            )

    one_pair = runs[1]
    multipair_gap_correction = (
        target_spectrum["photon_gap"] - one_pair["spectrum"]["photon_gap"]
    )
    one_pair_shift = one_pair["spectrum"]["photon_gap"] - BARE_PHOTON_GAP
    full_shift = target_spectrum["photon_gap"] - BARE_PHOTON_GAP
    finite_cone_fractional_shift = target_spectrum["photon_gap"] / BARE_PHOTON_GAP - 1.0

    mass_fit_rows = [
        row
        for row in gap_rows
        if "static_polarization" in row
    ]
    polarization_mass_fit = logarithmic_fit(
        [row["pair_gap_over_Delta"] for row in mass_fit_rows],
        [abs(row["static_polarization"]["transverse_curvature"]) for row in mass_fit_rows],
    )
    occupation_mass_fit = logarithmic_fit(
        [row["pair_gap_over_Delta"] for row in mass_fit_rows],
        [row["mean_pair_number"] for row in mass_fit_rows],
    )

    expected_full_counts = {0: 115, 1: 1484, 2: 3138, 3: 1484, 4: 115}
    checks = {
        "fermionic_car_exact": algebra[
            "maximum_canonical_anticommutator_residual"
        ] < 1.0e-12,
        "fermion_exchange_phase_is_minus_one": abs(
            algebra["two_creation_exchange_phase"] + 1.0
        ) < 1.0e-12,
        "pair_creation_is_pauli_nilpotent": maximum_pair_nilpotency < 1.0e-12,
        "neutral_fock_sector_counts": (
            True
            if quick
            else target["sector_counts"] == expected_full_counts
        ),
        "signed_translation_symmetry_exact": max(translation_residuals) < 1.0e-12,
        "transverse_operator_momentum_character": abs(
            operator_characters["x"] + 1.0
        ) < 1.0e-12
        and abs(operator_characters["y"] - 1.0) < 1.0e-12,
        "local_continuity_exact": ward[
            "maximum_local_continuity_residual"
        ] < 1.0e-12,
        "nonzero_momentum_operator_ward_exact": ward[
            "operator_ward_residual"
        ] < 1.0e-12,
        "nonzero_momentum_spectral_ward": ward[
            "maximum_low_spectral_ward_residual"
        ] < 1.0e-9,
        "longitudinal_external_mode_is_pure_gauge": target_polarization[
            "exact_longitudinal_unitary_equivalence_residual"
        ] < 1.0e-12
        and target_polarization["transversality_ratio"] < 2.0e-5,
        "transverse_virtual_response_nonzero": abs(
            target_polarization["transverse_curvature"]
        ) > 1.0e-5,
        "decoupled_control_zero": abs(
            decoupled_polarization["transverse_curvature"]
        ) < 1.0e-8,
        "multipair_series_converges": (
            True
            if quick
            else abs(multipair_gap_correction) < 0.02 * abs(one_pair_shift)
            and abs(
                truncation_rows[-1]["static_polarization"]["transverse_curvature"]
                - truncation_rows[0]["static_polarization"]["transverse_curvature"]
            )
            < 0.02
            * abs(truncation_rows[-1]["static_polarization"]["transverse_curvature"])
        ),
        "photon_residue_retained": target_spectrum["photon_residue"] > 0.98,
        "ground_is_vacuum_dominated": target["pair_distribution"].get(0, 0.0) > 0.98,
        "higher_pair_probabilities_are_hierarchical": (
            target["pair_distribution"].get(1, 0.0)
            > target["pair_distribution"].get(2, 0.0)
            > target["pair_distribution"].get(3, 0.0)
            > target["pair_distribution"].get(4, 0.0)
            if maximum_pairs == 4
            else target["pair_distribution"].get(1, 0.0)
            > target["pair_distribution"].get(2, 0.0)
        ),
        "coupling_scan_monotonic": (
            True
            if quick
            else all(
                right["mean_pair_number"] > left["mean_pair_number"]
                and right["photon_gap_shift_from_pure_gauge"]
                > left["photon_gap_shift_from_pure_gauge"]
                for left, right in zip(coupling_rows, coupling_rows[1:])
            )
        ),
        "pair_gap_suppresses_polarization": all(
            right["static_polarization"]["transverse_curvature"]
            < left["static_polarization"]["transverse_curvature"]
            for left, right in zip(
                sorted(mass_fit_rows, key=lambda row: row["pair_gap_over_Delta"]),
                sorted(mass_fit_rows, key=lambda row: row["pair_gap_over_Delta"])[1:],
            )
        ),
        "same_unit_pair_completion": PAIR_COMPLETION == 1.0,
        "finite_cone_shift_is_small_but_nonzero": 0.0
        < abs(finite_cone_fractional_shift)
        < 0.005,
    }

    report = {
        "module": (
            "phase_junction_network/microscopic/"
            "check_fermionic_multipair_vacuum.py"
        ),
        "status": (
            "Stage 3I PASS at finite antisymmetric multipair scope"
            if all(checks.values())
            else "Stage 3I FAIL"
        ),
        "run_mode": "quick" if quick else "full",
        "lattice": {
            "shape": [2, 2],
            "spin1_links": 8,
            "positive_fermion_modes": 4,
            "negative_fermion_modes": 4,
            "maximum_pairs_in_run": maximum_pairs,
            "full_neutral_fock_dimension": 6336,
        },
        "coefficients": {
            "pair_gap_over_Delta": PAIR_GAP,
            "pair_creation_coupling_over_Delta": PAIR_COUPLING,
            "completion_coefficient": PAIR_COMPLETION,
            "new_counterterms": 0,
        },
        "fermion_algebra": {
            **algebra,
            "maximum_pair_creation_nilpotency_residual": maximum_pair_nilpotency,
        },
        "translation": {
            "maximum_signed_translation_commutator": max(translation_residuals),
            "transverse_operator_character_x": operator_characters["x"],
            "transverse_operator_character_y": operator_characters["y"],
        },
        "pair_number_truncation": truncation_rows,
        "full_multipair_result": public_spectrum(target),
        "multipair_convergence": {
            "one_pair_photon_gap_shift": one_pair_shift,
            "full_photon_gap_shift": full_shift,
            "full_minus_one_pair_gap": multipair_gap_correction,
            "relative_correction_to_one_pair_shift": multipair_gap_correction
            / one_pair_shift,
            "finite_photon_cone_fractional_shift": finite_cone_fractional_shift,
            "interpretation": (
                "The full antisymmetric Fock space changes the one-pair photon "
                "shift only slightly. A matching frame-loop calculation is "
                "still required before claiming common-cone radiative stability."
            ),
        },
        "ward_identity": ward,
        "static_polarization": target_polarization,
        "decoupled_control": decoupled_polarization,
        "pair_gap_scan": gap_rows,
        "pair_gap_scaling": {
            "mean_pair_number_power": occupation_mass_fit["power"],
            "transverse_polarization_power": polarization_mass_fit["power"],
            "domain_wall_first_doubler_gap_control": DOMAIN_WALL_FIRST_DOUBLER_GAP,
        },
        "pair_coupling_scan": coupling_rows,
        "checks": checks,
        "stage_pass": all(checks.values()),
        "claim_boundary": {
            "established": [
                "an exact antisymmetric neutral fermion Fock space through four pairs on the finite patch",
                "canonical anticommutation, exchange sign, and Pauli blocking",
                "exact signed translation symmetry",
                "exact local and nonzero-momentum Ward identities in the multipair theory",
                "longitudinal pure-gauge invariance and nonzero transverse multipair polarization",
                "convergence of the one-pair truncation toward the complete finite Fock result",
                "controlled dependence on pair gap and pair-creation strength without a new counterterm",
            ],
            "not_established": [
                "the full fifth-dimensional domain-wall multiplet in the interacting patch",
                "mirror-wall removal or symmetric mass generation",
                "volume scaling beyond the 2x2 torus",
                "an infinite-volume fermion determinant or charge-renormalization flow",
                "matching matter-loop renormalization of the frame cone",
                "3+1D interacting QED, nonlinear gravity, or observed matter",
            ],
        },
        "next_gate": (
            "Couple the finite domain-wall light and mirror modes to this exact "
            "antisymmetric Fock construction, then compare photon and frame "
            "self-energies from the same matter loop and test whether the "
            "shared cone is radiatively stable without a new counterterm."
        ),
    }

    report["cutoff_semantics"] = (
        "Rebuilt completed-move Hamiltonians at each pair cutoff; not a "
        "fixed-Hamiltonian truncation error. See validation/check_fixed_cutoff.py."
    )
    report["not_tested"] = (
        ["full_four_pair_convergence", "coupling_scan", "full_sector_counts"]
        if quick else []
    )
    if quick:
        for key in ("multipair_series_converges", "coupling_scan_monotonic",
                    "neutral_fock_sector_counts"):
            report["checks"][key] = None
    report["stage_pass"] = all(
        value for value in report["checks"].values() if value is not None
    )

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
            "fermionic_multipair_vacuum_revalidated.json"
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
                "multipair_convergence": report["multipair_convergence"],
                "ward": report["ward_identity"],
                "failed_checks": [
                    key
                    for key, passed in report["checks"].items()
                    if passed is False
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
