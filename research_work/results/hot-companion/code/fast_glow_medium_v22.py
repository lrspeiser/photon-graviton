"""Round 22, step 2b: what a fast, even (isotropic) glow from collisions would require of the companion's medium.

Round 10 (README section 20.1, condition 3) found that the companion's crests must move much more slowly than its energy:
v_phase <= 0.05 u, or the coherent power that pays for the pull would spoil v^4 = G M a. A star crossing another
system's flow at a relative speed w therefore meets that flow's crests at the Doppler-shifted frequency
    omega' = omega_0 |1 - w_r / v_phase|,
hundreds of times the companion's own frequency omega_0. That is the frequency at which a collision drives each star, and
so the natural home of a distinct "hot" glow. If that glow is to spread evenly at v_h = 600 km/s (round 21), the
medium's group velocity must rise from u at omega_0 to about 3.5 u at omega'. This script computes:
  1. omega'/omega_0 for the Bullet's speeds and several crest speeds;
  2. the average logarithmic slope of v_g(omega) that this requires;
  3. why a single power law omega ~ k^n cannot give both (it has v_g/v_phase = n everywhere, so a slow crest at omega_0,
     n >= 20, forces v_g to rise by a factor 460^0.95 ~ 340 at omega', far more than 3.5);
  4. one two-regime dispersion that does give both, omega = u q + beta q^2 above a wavenumber k_0 where the frequency
     crosses zero (q = k - k_0): the curvature beta it needs.
It is a statement of requirements, not a derivation: nothing in the project yet fixes the medium's dispersion.

    python code/fast_glow_medium_v22.py --output run-crossing-frame-v22/fast_glow_medium_v22.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

U = 169.4


def two_regime(vph_over_u, ratio, vh_over_u):
    """omega = u q + beta q^2 (q = k - k0 >= 0), units u = k0 = 1. The carrier q1 has omega1/(1 + q1) = vph; find beta such
    that at omega' = ratio * omega1 the group velocity u + 2 beta q' equals vh."""
    def carrier(beta):
        return brentq(lambda q: (q + beta * q * q) / (1 + q) - vph_over_u, 1e-9, 10.0)

    def miss(beta):
        q1 = carrier(beta); w1 = q1 + beta * q1 * q1
        wp = ratio * w1
        qp = (-1 + np.sqrt(1 + 4 * beta * wp)) / (2 * beta)
        return (1 + 2 * beta * qp) - vh_over_u
    beta = brentq(miss, 1e-6, 100.0)
    q1 = carrier(beta)
    wp = ratio * (q1 + beta * q1 * q1)
    qp = (-1 + np.sqrt(1 + 4 * beta * wp)) / (2 * beta)
    return dict(beta_u_over_k0=float(beta), carrier_q_over_k0=float(q1), crest_crossing_q_over_k0=float(qp),
                group_velocity_at_carrier_over_u=float(1 + 2 * beta * q1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    vh_over_u = 600.0 / U
    rows = []
    for w in (3000.0, 3900.0):
        for f in (0.05, 0.02, 0.01):
            vph = f * U
            ratio = abs(1 - w / vph)
            alpha = np.log(vh_over_u) / np.log(ratio)
            n_needed = 1 / (1 - alpha)
            powerlaw_rise = ratio ** (1 - 1 / 20.0)          # v_g ~ omega^(1 - 1/n) with n = 20 (v_g/v_phase = 20)
            tr = two_regime(f, ratio, vh_over_u)
            rows.append(dict(w_kms=w, vphase_over_u=f, crest_crossing_ratio=float(ratio), slope_needed=float(alpha),
                             powerlaw_n_needed=float(n_needed), powerlaw_n20_group_velocity_rise=float(powerlaw_rise), **tr))
            print(f"w {w:5.0f} km/s, v_phase {f:.2f} u: omega'/omega0 = {ratio:7.0f}; needed d ln v_g / d ln omega = {alpha:.3f} "
                  f"(power law n = {n_needed:.2f}, i.e. v_g/v_phase = {n_needed:.2f}, not >= 20; with n = 20 v_g would rise x{powerlaw_rise:.0f}); "
                  f"two-regime curvature beta = {tr['beta_u_over_k0']:.3f} u/k0, carrier at q = {tr['carrier_q_over_k0']:.4f} k0")
    out = dict(experiment='round 22: what an even fast glow requires of the medium', u_kms=U, v_h_kms=600.0, rows=rows)
    args.output.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
