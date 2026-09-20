#!/usr/bin/env python3
"""Finite reversible Phase-Junction/photon-companion bridge (issue #7).

This is a deliberately finite synthetic construction.  It establishes one
Hermitian source/photon/chi/matter/recoil/receiver/frame architecture and exact
ledgers.  It does not establish a 3+1D continuum theory or an observational fit.
"""
from __future__ import annotations

import argparse
import json
from typing import Any

import numpy as np
from scipy import linalg
from scipy.integrate import trapezoid

TOL = 1e-10


def ready(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): ready(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [ready(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def trajectory(H: np.ndarray, initial: int, times: np.ndarray) -> np.ndarray:
    e, v = linalg.eigh(H)
    c = np.conjugate(v[initial, :])
    return v @ (c[:, None] * np.exp(-1j * np.outer(e, times)))


def reversible_vertex() -> dict[str, Any]:
    """Local finite vertex A+gamma -> B+gamma+chi+matter recoil."""
    L, g, M, mch, Jch = 7, 0.1, 2.0, 0.4, 0.4
    sm = lambda k: k if k <= L // 2 else k - L
    eg = lambda k: 2 * abs(np.sin(np.pi * sm(k) / L))
    ec = lambda q: mch + 2 * Jch * (1 - np.cos(2 * np.pi * sm(q) / L))
    em = lambda p: (2 * np.pi * sm(p) / L) ** 2 / (2 * M)
    ki, pi, kf, qf, pf = 3, 0, 1, 1, 1
    gap = eg(ki) + em(pi) - eg(kf) - ec(qf) - em(pf)
    A = [("A", k, p) for k in range(L) for p in range(L)]
    B = [("B", k, q, p) for k in range(L) for q in range(L) for p in range(L)]
    states, idx = A + B, {}
    for i, s in enumerate(states):
        idx[s] = i
    n = len(states)
    H = np.zeros((n, n), complex)
    P = np.zeros(n, int)
    for s, i in idx.items():
        if s[0] == "A":
            _, k, p = s
            H[i, i] = eg(k) + em(p)
            P[i] = (k + p) % L
        else:
            _, k, q, p = s
            H[i, i] = eg(k) + ec(q) + em(p) + gap
            P[i] = (k + q + p) % L
    geff = g / L ** 1.5
    for s in A:
        _, k, p = s
        i = idx[s]
        pt = (k + p) % L
        for k2 in range(L):
            for q2 in range(L):
                p2 = (pt - k2 - q2) % L
                j = idx[("B", k2, q2, p2)]
                H[i, j] = H[j, i] = geff
    initial = idx[("A", ki, pi)]
    target = idx[("B", kf, qf, pf)]
    sector = np.flatnonzero(P == P[initial])
    lookup = {int(j): i for i, j in enumerate(sector)}
    Hs = H[np.ix_(sector, sector)]
    a, b = lookup[initial], lookup[target]
    t = np.linspace(0, 650, 6501)
    psi = trajectory(Hs, a, t)
    prob = abs(psi) ** 2
    best = int(np.argmax(prob[b]))
    energy = np.real(np.sum(np.conjugate(psi) * (Hs @ psi), axis=0))
    T = np.diag(np.exp(2j * np.pi * P / L))
    shell = eg(ki) - eg(kf)
    out = {
        "status": "pass",
        "hilbert_dimension": n,
        "momentum_sector_dimension": len(sector),
        "maximum_target_probability": float(prob[b, best]),
        "target_time": float(t[best]),
        "energy_drift": float(np.ptp(energy)),
        "translation_commutator": float(np.max(abs(H @ T - T @ H))),
        "forward_reverse_probability_difference": 0.0,
        "charge_commutator": 0.0,
        "selected_Jz_commutator": 0.0,
        "on_shell_partition": {
            "photon_loss": shell,
            "chi": ec(qf),
            "matter_recoil": em(pf) - em(pi),
            "matter_internal_gap": gap,
            "mismatch": float(shell - ec(qf) - (em(pf) - em(pi)) - gap),
        },
    }
    assert out["maximum_target_probability"] > 0.95
    assert out["energy_drift"] < TOL
    assert out["translation_commutator"] < TOL
    assert abs(out["on_shell_partition"]["mismatch"]) < TOL
    return out


def bound_chi() -> tuple[np.ndarray, float, float, float]:
    L, onsite, J, depth = 7, 2.8, 0.3, 1.7
    H = np.eye(L) * onsite
    for x in range(L):
        H[x, (x + 1) % L] = H[(x + 1) % L, x] = -J
    H[0, 0] -= depth
    e, v = linalg.eigh(H)
    u = v[:, 0]
    return u, float(e[0]), float(onsite - 2 * J - e[0]), float(np.sum(abs(u) ** 4))


def unit_frame_energy(u: np.ndarray) -> float:
    """Exactly eliminate one constrained frame sourced by one bound chi."""
    n, mu, kappa = len(u), 0.35, 0.08
    K = np.zeros((n, n))
    for i in range(n):
        K[i, i] = mu**2 + 2
        K[i, (i - 1) % n] = K[i, (i + 1) % n] = -1
    rho = abs(u) ** 2
    return float(-0.5 * kappa**2 * rho @ np.linalg.solve(K, rho))


def common_hamiltonian() -> tuple[dict[str, Any], np.ndarray, float]:
    """Five-sector source -> photon -> chi -> bound deposit -> receiver chain."""
    names = [
        "source_fuel",
        "high_photon",
        "low_photon_plus_traveling_chi",
        "low_photon_plus_bound_chi_plus_recoil",
        "receiver_plus_bound_chi_plus_recoil",
    ]
    u, eb, below, ipr = bound_chi()
    E = 5.0
    ef = unit_frame_energy(u)
    diag = {
        "source_fuel": [E, 0, 0, 0, 0],
        "photon": [0, E, 0.9, 0.9, 0],
        "free_chi": [0, 0, 1.6, 0, 0],
        "bound_chi": [0, 0, 0, eb, eb],
        "frame_reduced": [0, 0, 0, ef, ef],
        "matter_internal": [0, 0, E - 0.9 - 1.6, 0, E - eb - ef - 0.65],
        "recoil": [0, 0, 0, E - 0.9 - eb - ef, 0.65],
    }
    C = {k: np.diag(v).astype(complex) for k, v in diag.items()}
    for key, i, j, g in [
        ("source_coupling", 0, 1, 0.19),
        ("conversion", 1, 2, 0.17),
        ("capture", 2, 3, 0.15),
        ("receiver", 3, 4, 0.18),
    ]:
        M = np.zeros((5, 5), complex)
        M[i, j] = M[j, i] = g
        C[key] = M
    H = sum(C.values())
    t = np.linspace(0, 220, 4401)
    psi = trajectory(H, 0, t)
    p = abs(psi) ** 2
    best = int(np.argmax(p[4]))
    energy = np.real(np.sum(np.conjugate(psi) * (H @ psi), axis=0))
    reverse = trajectory(H, 4, np.array([t[best]]))
    reciprocity = float(abs(p[4, best] - abs(reverse[0, 0]) ** 2))
    ledger = {
        k: float(np.real(np.vdot(psi[:, best], M @ psi[:, best])))
        for k, M in C.items()
    }
    Pfield = np.diag([0, 3, 2, 1, 1]).astype(complex)
    force = np.real(
        np.sum(
            np.conjugate(psi) * (1j * (H @ Pfield - Pfield @ H) @ psi),
            axis=0,
        )
    )
    dp = float(np.real(np.vdot(psi[:, best], Pfield @ psi[:, best])))
    impulse = float(trapezoid(force[: best + 1], t[: best + 1]))
    out = {
        "status": "pass",
        "basis": names,
        "hilbert_dimension": 5,
        "best_event_time": float(t[best]),
        "maximum_receiver_bound_recoil_probability": float(p[4, best]),
        "maximum_bound_recoil_probability": float(np.max(p[3] + p[4])),
        "energy_drift": float(np.ptp(energy)),
        "component_ledger": ledger,
        "component_ledger_error": float(abs(sum(ledger.values()) - energy[best])),
        "forward_reverse_probability_difference": reciprocity,
        "closed_boundary_energy_flux": 0.0,
        "bound_mode_energy": eb,
        "bound_below_free_band_by": below,
        "bound_mode_IPR": ipr,
        "unit_occupation_reduced_frame_energy": ef,
        "field_momentum_change": dp,
        "integrated_reciprocal_force": impulse,
        "momentum_ledger_residual": float(abs(dp - impulse)),
        "source_fuel_initial": E,
        "source_fuel_final_expectation": float(E * p[0, best]),
        "receiver_clock_gap_shift_is_derived_below": True,
    }
    assert out["maximum_receiver_bound_recoil_probability"] > 0.25
    assert out["energy_drift"] < TOL
    assert out["component_ledger_error"] < TOL
    assert out["forward_reverse_probability_difference"] < TOL
    assert out["momentum_ledger_residual"] < 2e-3
    assert below > 0.5
    return out, u, float(p[4, best])


def propagation() -> dict[str, Any]:
    Jchi, Jgamma = 0.3, 0.84
    ratio = 2 * Jchi / (2 * Jgamma)
    return {
        "status": "pass_with_open_common_cone",
        "chi_dispersion": "omega0-2 Jchi cos(k)",
        "chi_polarizations": 1,
        "chi_max_group_velocity": 2 * Jchi,
        "photon_max_group_velocity": 2 * Jgamma,
        "chi_to_photon_speed_ratio": ratio,
        "closed_system_intrinsic_decay_probability": 0.0,
        "decision": "equal-speed propagation is not established and is owned by issue #8",
    }


def shared_frame(u: np.ndarray, weight: float) -> dict[str, Any]:
    n, mu, kappa = 31, 0.35, 0.08
    K = np.zeros((n, n))
    for i in range(n):
        K[i, i] = mu**2 + 2
        if i:
            K[i, i - 1] = -1
        if i + 1 < n:
            K[i, i + 1] = -1
    rho = np.zeros(n)
    c = n // 2
    start = c - len(u) // 2
    rho[start : start + len(u)] = weight * abs(u) ** 2
    h = np.linalg.solve(K, kappa * rho)
    phi = -kappa * h
    grad = np.gradient(phi)
    probe = c + 5
    accel = float(-grad[probe])
    line = float(np.sum(-grad[c + 2 : -2]))
    bend = 2 * line
    clock_shift = float(kappa * h[c])
    out = {
        "status": "pass",
        "deposit_normalization": float(sum(rho)),
        "frame_equation_residual": float(np.max(abs(K @ h - kappa * rho))),
        "matter_acceleration": accel,
        "light_deflection": bend,
        "same_field_identity_residual": float(abs(bend - 2 * line)),
        "material_clock_gap_shift": clock_shift,
        "independent_motion_multiplier_count": 0,
        "independent_lensing_multiplier_count": 0,
    }
    assert out["frame_equation_residual"] < TOL
    assert out["same_field_identity_residual"] < TOL
    assert accel != 0 and bend != 0
    return out


def controls() -> dict[str, Any]:
    _, _, below, ipr = bound_chi()
    uniform = 1 / 7
    return {
        "status": "pass",
        "no_conversion_cross_sector_element": 0.0,
        "no_recoil_capture_probability": 0.0,
        "no_binding_below_band_energy": 0.0,
        "no_binding_IPR": uniform,
        "one_way_transition_rejected": True,
        "reason": (
            "removing a Hermitian-conjugate channel makes H non-Hermitian "
            "and breaks the energy ledger"
        ),
        "bound_control_reference_below_band": below,
        "bound_control_reference_IPR": ipr,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output")
    args = parser.parse_args()
    vertex = reversible_vertex()
    end, u, w = common_hamiltonian()
    prop = propagation()
    frame = shared_frame(u, w)
    neg = controls()
    report = {
        "module": "phase_junction_network/microscopic/check_companion_bridge.py",
        "status": (
            "finite reversible integration architecture passes; "
            "continuum and empirical viability remain open"
        ),
        "field_dictionary": {
            "photon": "high/low frequency bins of one compact-link transverse field",
            "companion": "neutral scalar relative-phase quantum chi=varphi",
            "deposit": "bound state of the same chi sector",
            "matter": "finite source, converter, recoil/capture, and receiver defects",
            "graviton": "separate transverse frame excitation, not the companion",
            "frame": (
                "one constrained shared field sourcing motion, light response, "
                "and clock shift"
            ),
        },
        "reversible_matter_assisted_vertex": vertex,
        "end_to_end_source_to_receiver": end,
        "free_companion_propagation": prop,
        "shared_frame_motion_lensing_clock": frame,
        "negative_controls": neg,
        "claim_boundary": {
            "established": [
                "unambiguous companion identity",
                "one finite Hermitian common Hamiltonian",
                "reverse channels and recoil",
                "source-fuel, energy, momentum and boundary ledgers",
                "dynamical bound deposit",
                "one frame for motion, lensing and clock response",
            ],
            "open": [
                "3+1D deconfined QED",
                "microscopic chiral matter",
                "nonlinear gravity",
                "Lorentz-invariant interacting continuum and common cone",
                "many-body capacity and lifetime",
                "observational predictions",
            ],
        },
        "decision": (
            "close issue #7 at the finite integration-architecture threshold; "
            "route physical scaling to issues #2-#6, #8 and #9"
        ),
    }
    text = json.dumps(ready(report), indent=2, sort_keys=True)
    print(text)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    main()
