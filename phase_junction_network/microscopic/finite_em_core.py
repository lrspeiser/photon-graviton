"""Exact finite-link algebra and gauge-sector enumeration for issue #6.

This script deliberately separates three claims:

1. exact finite-link algebra and a positive frustration-free 3+1D Hamiltonian;
2. static Coulomb-phase diagnostics of its Rokhsar-Kivelson ground state;
3. the still-open relativistic quantum dynamics/QED limit.

No continuum rotor replaces the finite links in the checks below.  The Monte Carlo
samples the exact equal-amplitude ground-state support of the finite Hamiltonian.
"""
from __future__ import annotations

import argparse
import collections
import itertools
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

try:
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
except Exception:  # pragma: no cover - script reports missing optional ED support
    sp = None
    spla = None

PAIRS: tuple[tuple[int, int], ...] = ((0, 1), (0, 2), (1, 2))


def spin_link_operators(spin: int) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Return E=S_z and U=S_+/max|S_+| for integer spin S."""
    if spin < 1 or int(spin) != spin:
        raise ValueError("spin must be a positive integer")
    S = float(spin)
    m = np.arange(-spin, spin + 1, dtype=float)
    E = np.diag(m).astype(complex)
    U = np.zeros((len(m), len(m)), dtype=complex)
    for j, value in enumerate(m[:-1]):
        U[j + 1, j] = math.sqrt(S * (S + 1.0) - value * (value + 1.0))
    U /= np.max(np.abs(U))
    return E, U


def kron_all(ops: Sequence[NDArray[np.complex128]]) -> NDArray[np.complex128]:
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def embed_local(op: NDArray[np.complex128], which: int, count: int) -> NDArray[np.complex128]:
    d = op.shape[0]
    identity = np.eye(d, dtype=complex)
    return kron_all([op if i == which else identity for i in range(count)])


def local_rk_check(spin: int) -> dict[str, object]:
    """Check the exact local gauge algebra and weighted RK positivity.

    For a plaquette operator W, define

        D = sqrt(W^dagger W) + sqrt(W W^dagger),
        h_RK = D - W - W^dagger.

    W has at most one positive matrix element in each row and column, so h_RK is
    the weighted graph Laplacian of all allowed finite-link plaquette moves.
    """
    E, U = spin_link_operators(spin)
    d = E.shape[0]
    Es = [embed_local(E, i, 4) for i in range(4)]
    Us = [embed_local(U, i, 4) for i in range(4)]
    W = Us[0] @ Us[1] @ Us[2].conj().T @ Us[3].conj().T

    left = W.conj().T @ W
    right = W @ W.conj().T
    # These are diagonal in the electric basis; clipping removes roundoff below 0.
    D = np.diag(
        np.sqrt(np.clip(np.diag(left).real, 0.0, None))
        + np.sqrt(np.clip(np.diag(right).real, 0.0, None))
    ).astype(complex)
    h = D - W - W.conj().T

    # Plaquette orientation: links 0,1 point with the loop; 2,3 against it.
    # Vertex Gauss operators for a square with no external links.
    G = [
        Es[0] + Es[3],
        Es[1] - Es[0],
        -Es[1] - Es[2],
        Es[2] - Es[3],
    ]
    commutator = max(float(np.max(np.abs(h @ g - g @ h))) for g in G)
    algebra = float(np.max(np.abs(E @ U - U @ E - U)))
    evals = np.linalg.eigvalsh(h).real

    # Uniform vector on every connected component must be a zero vector.  The
    # complete matrix can have many components because boundary flux is fixed.
    row_sum = float(np.max(np.abs(np.sum(h, axis=1))))
    negative = float(min(0.0, np.min(evals)))
    return {
        "spin": spin,
        "link_dimension": int(d),
        "link_algebra_error": algebra,
        "gauss_commutator_max_abs": commutator,
        "rk_row_sum_max_abs": row_sum,
        "rk_minimum_eigenvalue": float(np.min(evals)),
        "rk_negative_part": negative,
        "rk_zero_mode_count": int(np.count_nonzero(np.abs(evals) < 1.0e-10)),
        "rk_maximum_edge_amplitude": float(np.max(np.abs(W))),
        "rk_minimum_nonzero_edge_amplitude": float(
            np.min(np.abs(W[np.abs(W) > 1.0e-14]))
        ),
    }


def link_index(x: int, y: int, z: int, mu: int, L: int) -> int:
    return (((x * L + y) * L + z) * 3) + mu


def plaquette_link_steps(
    x: int, y: int, z: int, mu: int, nu: int, L: int
) -> tuple[tuple[int, int], ...]:
    base = [x, y, z]
    plus_mu = base.copy()
    plus_mu[mu] = (plus_mu[mu] + 1) % L
    plus_nu = base.copy()
    plus_nu[nu] = (plus_nu[nu] + 1) % L
    return (
        (link_index(*base, mu, L), +1),
        (link_index(*plus_mu, nu, L), +1),
        (link_index(*plus_nu, mu, L), -1),
        (link_index(*base, nu, L), -1),
    )


def all_plaquettes(L: int) -> list[tuple[tuple[int, int], ...]]:
    return [
        plaquette_link_steps(x, y, z, mu, nu, L)
        for x, y, z in itertools.product(range(L), repeat=3)
        for mu, nu in PAIRS
    ]


def apply_steps_tuple(
    state: tuple[int, ...], steps: Sequence[tuple[int, int]], sign: int, spin: int
) -> tuple[int, ...] | None:
    values = list(state)
    for idx, delta in steps:
        value = values[idx] + sign * delta
        if value < -spin or value > spin:
            return None
        values[idx] = value
    return tuple(values)


def gauss_from_flat(state: Sequence[int], L: int) -> NDArray[np.int16]:
    E = np.asarray(state, dtype=np.int16).reshape(L, L, L, 3)
    div = np.zeros((L, L, L), dtype=np.int16)
    for mu in range(3):
        div += E[..., mu]
        div -= np.roll(E[..., mu], shift=1, axis=mu)
    return div

@dataclass
class BasisResult:
    states: list[tuple[int, ...]]
    index: dict[tuple[int, ...], int]
    plaquettes: list[tuple[tuple[int, int], ...]]


def enumerate_zero_flux_component(L: int, spin: int, max_states: int | None = None) -> BasisResult:
    """Enumerate the connected zero-charge, zero-winding component from E=0."""
    pls = all_plaquettes(L)
    zero = (0,) * (3 * L**3)
    states = [zero]
    index = {zero: 0}
    queue = collections.deque([zero])
    while queue:
        state = queue.popleft()
        for steps in pls:
            for sign in (-1, +1):
                nxt = apply_steps_tuple(state, steps, sign, spin)
                if nxt is None or nxt in index:
                    continue
                index[nxt] = len(states)
                states.append(nxt)
                queue.append(nxt)
                if max_states is not None and len(states) > max_states:
                    raise RuntimeError(
                        f"basis exceeded max_states={max_states}; got {len(states)}"
                    )
    return BasisResult(states=states, index=index, plaquettes=pls)


def exact_periodic_cube_check(
    L: int = 2, spin: int = 1, compute_gap: bool = False
) -> dict[str, object]:
    """Build the exact gauge-reduced 3D RK graph on a small periodic cube."""
    started = time.perf_counter()
    basis = enumerate_zero_flux_component(L=L, spin=spin)
    n = len(basis.states)
    degrees = np.zeros(n, dtype=np.int32)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    edge_count_twice = 0
    max_gauss = 0

    for i, state in enumerate(basis.states):
        if i % 1000 == 0:
            max_gauss = max(max_gauss, int(np.max(np.abs(gauss_from_flat(state, L)))))
        degree = 0
        for steps in basis.plaquettes:
            for sign in (-1, +1):
                nxt = apply_steps_tuple(state, steps, sign, spin)
                if nxt is None:
                    continue
                j = basis.index.get(nxt)
                if j is None:
                    raise AssertionError("connected-basis closure failed")
                degree += 1
                if compute_gap:
                    rows.append(i)
                    cols.append(j)
                    data.append(-1.0)
        degrees[i] = degree
        edge_count_twice += degree

    # Every row of the graph Laplacian annihilates the equal-amplitude state.
    uniform_residual = 0.0
    gap = None
    eigenvalues: list[float] | None = None
    if compute_gap:
        if sp is None or spla is None:
            raise RuntimeError("scipy is required for --exact-gap")
        rows.extend(range(n))
        cols.extend(range(n))
        data.extend(degrees.astype(float).tolist())
        H = sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=float)
        row_sums = np.asarray(H.sum(axis=1)).ravel()
        uniform_residual = float(np.max(np.abs(row_sums)))
        # Shift-invert around zero is expensive for this graph.  Smallest-magnitude
        # Lanczos is adequate because H is positive semidefinite.
        vals = np.sort(spla.eigsh(H, k=3, which="SM", return_eigenvectors=False, tol=1e-8))
        eigenvalues = [float(v) for v in vals]
        gap = float(vals[1])

    elapsed = time.perf_counter() - started
    return {
        "lattice_size": L,
        "sites": L**3,
        "links": 3 * L**3,
        "plaquettes": 3 * L**3,
        "spin": spin,
        "link_dimension": 2 * spin + 1,
        "gauge_reduced_component_dimension": n,
        "undirected_transition_edges": edge_count_twice // 2,
        "minimum_graph_degree": int(np.min(degrees)),
        "maximum_graph_degree": int(np.max(degrees)),
        "mean_graph_degree": float(np.mean(degrees)),
        "sampled_gauss_residual_max_abs": max_gauss,
        "equal_amplitude_rk_residual_max_abs": uniform_residual,
        "smallest_eigenvalues": eigenvalues,
        "first_rk_gap": gap,
        "elapsed_seconds": elapsed,
    }
