"""Round 17, step B: the independent review's protected-rhythm candidate, reproduced for one piece, and why it does not
carry over to pieces that share one wave.

The candidate: a quiet oscillation D mixed by the motion (strength g) into two radiating families B+ and B- whose own
frequencies sit at +Delta and -Delta, each radiating at gamma through its own channel:
    dD/dt  = -gamma0 D - i g (B+ + B-)
    dB+/dt = -(gamma + i Delta) B+ - i g D,      dB-/dt = -(gamma - i Delta) B- - i g D.
The quiet mode's complex rate is the eigenvalue of this 3 x 3 matrix that tends to -gamma0 as g -> 0. Its real part is
the added decay (the released glow), its imaginary part the shift of the rhythm. The review's numbers:
  * the added decay grows as g^2 (x 4.00008, 4.00030, 4.00120 per doubling of g from 0.0025 to 0.02),
  * the shift stays below 1.1e-18 (the two families' pulls cancel exactly, by symmetry),
  * opposite 1% deviations of the two couplings shift the rhythm by about 2% of the added decay,
  * collisions at rate nu = gamma (each family's coherence lost at gamma + nu), with Delta = gamma, keep 80% of the
    extra glow (a single resonant family keeps 50%).
Also computed: the first-order sensitivity of the quiet mode's rate to a common complex shift eps of both families'
frequencies, d(rate)/d(eps), which vanishes at Delta = gamma (sum over +- of 1/(gamma +- i Delta)^2 = 0), and the same
sensitivity when the two families radiate into ONE shared channel, so that the surroundings act on the channel (on the
sum B+ + B-) instead of on each family's frequency. In the shared case the sensitivity is proportional to the square of
the channel's own response, chi^2 with chi = sum over +- of 1/(gamma +- i Delta) = 2 gamma/(gamma^2 + Delta^2), which is
real and never zero: no choice of Delta cancels it.

    python code/paired_channel_v17.py --output run-rhythm-budget-v17/paired_channel_v17.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np


def matrix(g, gamma0=1e-3, gamma=1.0, Delta=1.0, gp=None, gm=None, nu=0.0, eps_p=0.0, eps_m=0.0):
    """The 3 x 3 rate matrix of (D, B+, B-). gp, gm: the two couplings (default g); nu: collision rate added to both
    families' decay; eps_p, eps_m: complex shifts of the families' frequencies by their surroundings."""
    gp = g if gp is None else gp; gm = g if gm is None else gm
    return np.array([[-gamma0, -1j * gp, -1j * gm],
                     [-1j * gp, -(gamma + nu + 1j * (Delta + eps_p)), 0],
                     [-1j * gm, 0, -(gamma + nu + 1j * (-Delta + eps_m))]], complex)


def quiet_rate(A, gamma0):
    ev = np.linalg.eigvals(A)
    return ev[np.argmin(np.abs(ev + gamma0))]


def shared_matrix(g, gamma0=1e-3, gamma=1.0, Delta=1.0, f=0.2, eps=0.0):
    """The same two families, each losing gamma (1 - f) inside, but radiating their share gamma f into ONE shared
    channel: each family feels the radiation reaction of the channel's total, gamma f (B+ + B-). The surroundings act on
    the channel: gamma f -> gamma f + i eps (eps real: a reactive shift; eps = -i e: extra damping)."""
    gi, Gr = gamma * (1 - f), gamma * f + 1j * eps
    return np.array([[-gamma0, -1j * g, -1j * g],
                     [-1j * g, -(gi + 1j * Delta + Gr), -Gr],
                     [-1j * g, -Gr, -(gi - 1j * Delta + Gr)]], complex)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    g0, gam, Dl = 1e-3, 1.0, 1.0
    out = {}
    # 1. growth with g and the shift
    gs = [0.0025, 0.005, 0.01, 0.02]
    rates = [quiet_rate(matrix(g, g0, gam, Dl), g0) for g in gs]
    added = [-(r.real) - g0 for r in rates]
    out['added_decay'] = dict(g=gs, added=added, per_doubling=[added[i + 1] / added[i] for i in range(3)],
                              shift=[float(r.imag) for r in rates], expected_small_g=[2 * g * g * gam / (gam ** 2 + Dl ** 2) for g in gs])
    print('added decay per doubling of g:', ', '.join(f'{x:.5f}' for x in out['added_decay']['per_doubling']))
    print('largest |shift|:', max(abs(s) for s in out['added_decay']['shift']))
    # 2. opposite 1% coupling deviations
    g = 0.01
    r = quiet_rate(matrix(g, g0, gam, Dl, gp=g * 1.01, gm=g * 0.99), g0)
    out['asymmetry_1pc'] = dict(shift=float(r.imag), added=float(-r.real - g0), ratio=float(abs(r.imag) / (-r.real - g0)),
                                expected_ratio=2 * 0.01 * Dl / gam)
    print(f"1% opposite coupling deviations: shift/added decay = {out['asymmetry_1pc']['ratio']:.4f} (2 eps Delta/gamma = {2 * 0.01 * Dl / gam:.4f})")
    # 3. collisions nu = gamma
    free = -quiet_rate(matrix(g, g0, gam, Dl), g0).real - g0
    coll = -quiet_rate(matrix(g, g0, gam, Dl, nu=gam), g0).real - g0
    free1 = -quiet_rate(matrix(g, g0, gam, 0.0, gm=0.0), g0).real - g0     # one resonant family, same coupling
    coll1 = -quiet_rate(matrix(g, g0, gam, 0.0, gm=0.0, nu=gam), g0).real - g0
    out['collisions'] = dict(kept_pair=float(coll / free), kept_single=float(coll1 / free1))
    print(f"collisions at nu = gamma: the pair keeps {coll / free:.3f} of its extra glow, one resonant family {coll1 / free1:.3f}")
    # 4. first-order sensitivity to the surroundings
    h = 1e-7
    sens = {}
    for Dl_ in (0.0, 0.5, 1.0, 2.0):
        base = quiet_rate(matrix(g, g0, gam, Dl_), g0)
        sep_re = (quiet_rate(matrix(g, g0, gam, Dl_, eps_p=h, eps_m=h), g0) - base) / h      # reactive shift of both
        sep_im = (quiet_rate(matrix(g, g0, gam, Dl_, eps_p=-1j * h, eps_m=-1j * h), g0) - base) / h   # extra damping of both
        b2 = quiet_rate(shared_matrix(g, g0, gam, Dl_), g0)
        sh_re = (quiet_rate(shared_matrix(g, g0, gam, Dl_, eps=h), g0) - b2) / h
        sh_im = (quiet_rate(shared_matrix(g, g0, gam, Dl_, eps=-1j * h), g0) - b2) / h
        # the separate-channel sensitivity counts the surroundings acting on each family's whole rate; for a like-for-like
        # comparison with the shared channel, where they act on the radiated share f only, it is scaled by f below
        added_sep = -base.real - g0; added_sh = -b2.real - g0
        sens[str(Dl_)] = dict(separate_reactive=[sep_re.real / added_sep, sep_re.imag / added_sep],
                              separate_dissipative=[sep_im.real / added_sep, sep_im.imag / added_sep],
                              shared_added_decay_over_separate=float(added_sh / added_sep),
                              shared_reactive=[sh_re.real / added_sh, sh_re.imag / added_sh],
                              shared_dissipative=[sh_im.real / added_sh, sh_im.imag / added_sh])
        print(f"Delta = {Dl_}: d(rate)/d(eps) per unit added decay [decay, shift]: separate channels, reactive "
              f"{sep_re.real / added_sep:+.3f} {sep_re.imag / added_sep:+.3f}, dissipative {sep_im.real / added_sep:+.3f} {sep_im.imag / added_sep:+.3f} | "
              f"one shared channel, reactive {sh_re.real / added_sh:+.3f} {sh_re.imag / added_sh:+.3f}, dissipative {sh_im.real / added_sh:+.3f} {sh_im.imag / added_sh:+.3f}")
    out['sensitivity'] = sens
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
