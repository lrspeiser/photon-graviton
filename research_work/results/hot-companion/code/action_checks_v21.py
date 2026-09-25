"""Round 21, step 4: checks on the relativistic action proposed for the law (supplied 25 September 2026).

The proposal introduces a response field chi, constrained along the companion's net energy-flow direction s^mu,
    S_chi = (1/8 pi G) int sqrt(-g) [ -1/2 P^{mu nu} d_mu chi d_nu chi + lambda (s^mu d_mu chi + f sqrt(a I)) ],
with I -> |g_N| + S in the weak field (the companion's occupation), f the release factor, P^{mu nu} acting across the
flow, and a metric in which matter and light both see Phi_N + chi (no slip). Four weak-field checks, with numbers:

  1. the sign: with s pointing outward (the flow's direction), s.grad chi = -f sqrt(a I) makes chi fall outward, so
     -grad chi points outward: the extra pull repels. The constraint needs s.grad chi = +f sqrt(a I).
  2. the constraint field: varying chi, with matter coupled to Phi_N + chi, gives div(lambda s) = -8 pi G rho across a
     spherical source, so lambda = -2 G M(r) / r^2 = -2 |g_N|, whatever the sign in 1. Its stress-energy is of order
     |lambda| |grad chi| / 8 pi G: the mass-equivalent it adds, relative to the source, is (v_f / c)^2 ln(r2 / r1).
  3. the reaction on warm matter: if I contains S = G int k rho / d^2 computed from the matter, varying the matter's
     positions adds a force on warm matter with density (k / 8 pi) grad int (lambda df/dS ...) / d^2, the same form as
     round 2's action (README section 8.2), whose reaction round 3 excluded (it pushed cluster stars outward by 32-73%
     of gravity, section 10.6). The ratio of the two reaction strengths is computed here.
  4. gravitational waves: if light and matter felt chi but gravitational waves did not (chi only in the matter's
     metric), GW170817's waves and light would be delayed differently by the extra potential along the way (the
     Shapiro delay). The Milky Way's share alone is computed here, with the law's constants and the companion's reach.

    python code/action_checks_v21.py --output run-action-checks-v21/action_checks_v21.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.integrate import quad

G_SI = 6.674e-11
C_SI = 2.99792458e8
MSUN = 1.989e30
KPC = 3.0857e19
A_LAW = 6.297889390049439e-11          # round 12
GD_LAW = 2.027e-10
U_KMS = 169.4


def rotation(M_msun, r_kpc, sign):
    """Circular speed (km/s) of a point mass with the constraint's sign (+1 attracts, -1 as written)."""
    gN = G_SI * M_msun * MSUN / (r_kpc * KPC) ** 2
    extra = np.exp(-gN / GD_LAW) * np.sqrt(A_LAW * gN)
    g = gN + sign * extra
    return float(np.sign(g) * np.sqrt(abs(g) * r_kpc * KPC) / 1e3), float(gN), float(extra)


def lambda_mass_fraction(M_msun, r1_kpc, r2_kpc):
    """Mass-equivalent of the constraint field's stress-energy between r1 and r2, over M: with lambda = 2 G M / r^2 and
    |grad chi| = f sqrt(a G M) / r, energy density |lambda| |grad chi| / 8 pi G = M f sqrt(a G M) / (4 pi r^3)."""
    GM = G_SI * M_msun * MSUN
    f = lambda r: np.exp(-(GM / (r * KPC) ** 2) / GD_LAW)
    integral = quad(lambda r: f(r) / r, r1_kpc, r2_kpc, limit=200)[0]
    return float(np.sqrt(A_LAW * GM) / C_SI ** 2 * integral), float(np.sqrt(np.sqrt(A_LAW * GM)) / 1e3)


def reaction_ratio(gN, S):
    """lambda-term reaction over round 2's: 2 gN f sqrt(a) / (2 sqrt(gN + S))  over  int_0^gN e^(-q/gd) sqrt(a / (q + S)) dq."""
    ours = gN * np.exp(-gN / GD_LAW) * np.sqrt(A_LAW / (gN + S))
    round2 = quad(lambda q: np.exp(-q / GD_LAW) * np.sqrt(A_LAW / (q + S)), 0.0, gN)[0]
    return float(ours / round2)


def shapiro_mw(M_msun, reach_kpc, r_sun_kpc=8.2):
    """Extra one-way delay of light over gravitational waves from the Milky Way's chi alone, if waves ignored chi.
    Outside the matter chi = sqrt(a G M) ln(r / reach) (zero beyond the companion's reach); the delay is
    (2 / c^3) int |chi| dl along a line leaving the Galaxy from the Sun (taken radially, a lower estimate)."""
    GM = G_SI * M_msun * MSUN
    v2 = np.sqrt(A_LAW * GM)
    R = reach_kpc * KPC
    integral = quad(lambda r: np.log(R / r), r_sun_kpc * KPC, R, limit=200)[0]
    return float(2 * v2 / C_SI ** 3 * integral)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    out = dict(experiment='round 21: weak-field checks on the proposed relativistic action', constants=dict(a=A_LAW, g_d=GD_LAW, u_kms=U_KMS))
    # 1. the sign, for a Milky-Way-like mass
    M = 6.0e10
    rows = []
    for r in (5.0, 10.0, 20.0, 40.0):
        vp, gN, ex = rotation(M, r, +1); vm, _, _ = rotation(M, r, -1)
        vn = float(np.sqrt(gN * r * KPC) / 1e3)
        rows.append(dict(r_kpc=r, v_newton=vn, v_constraint_plus=vp, v_constraint_as_written=vm))
        print(f"1. r = {r:4.0f} kpc: Newton {vn:6.1f} km/s; constraint s.grad chi = +f sqrt(aI): {vp:6.1f}; as written (-): {vm:6.1f}")
    out['sign'] = dict(M_msun=M, rows=rows)
    # 2. the constraint field's mass-equivalent
    frac, vf = lambda_mass_fraction(M, 1.0, 1000.0)
    out['lambda_mass_fraction'] = dict(M_msun=M, r_kpc=[1.0, 1000.0], fraction=frac, v_flat_kms=vf)
    print(f"2. lambda's stress-energy between 1 kpc and 1 Mpc: {frac:.2e} of the source mass (v_f = {vf:.0f} km/s)")
    # 3. the reaction on warm matter, relative to round 2's
    gN = 1.0e-10
    rr = []
    for s_over in (0.01, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0):
        q = reaction_ratio(gN, s_over * gN)
        rr.append(dict(S_over_gN=s_over, ratio=q))
    out['reaction_ratio_vs_round2'] = dict(g_N=gN, rows=rr)
    print('3. reaction strength over round 2\'s, at S/g_N = ' + ', '.join(f"{d['S_over_gN']:g}: {d['ratio']:.2f}" for d in rr))
    # 4. the Shapiro delay from the Milky Way's chi
    reach = U_KMS * 1.0227121650537077 * 13.0
    dt = shapiro_mw(M, reach)
    out['shapiro_mw'] = dict(M_msun=M, reach_kpc=reach, delay_s=dt, delay_years=dt / 3.156e7, gw170817_observed_s=1.7)
    print(f"4. Milky Way alone, reach {reach:.0f} kpc: light would lag gravitational waves by {dt:.2e} s ({dt / 3.156e7:.1f} years); "
          f"GW170817 observed 1.7 s")
    args.output.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
