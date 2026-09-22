#!/usr/bin/env python3
"""JR-10 stage 4: what the locked formula does and does not derive.

Stage 3 locked a law. This stage separates the parts that are exact linearised GR
from the part that is an empirical fit, bounds what stress engineering can ever do
for the energy problem, and checks whether the driving variable is even identified.

    python derive_stress.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np

import model as M

MSUN, C = 1.98892e33, 2.99792458e10
YR = 3.155693e7


def tolman(w_r, w_t):
    """tau = 1 + w_r + 2 w_t, the factor multiplying rho in the source of Phi."""
    return 1.0 + w_r + 2.0 * w_t


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--stage1', type=Path, default=Path('../run-v1/required-stress.json'))
    ap.add_argument('--stage3', type=Path, default=Path('../run-v3-bounded/bounded-law.json'))
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    s1 = json.loads(args.stage1.read_text())
    s3 = json.loads(args.stage3.read_text())
    rows = s1['rows']

    # ---- 1. What is exact, and what is fitted --------------------------------
    exact = dict(
        statement='These follow from linearised GR with no fitting whatsoever.',
        items=[
            'Phi is sourced by (rho + T^k_k) and Psi by rho, so gamma_chi = Psi/Phi = 1/tau '
            'with tau = 1 + w_r + 2 w_t.',
            'Light deflection carries the companion at weight (1 + gamma_chi)/2 relative to '
            'the static-dust value, which is exactly the R11 form.',
            'A free-streaming radial null flux has w_r = 1, w_t = 0, hence tau = 2 exactly.',
            'A static radial field gradient has w_r = 1, w_t = -1, hence tau = 0 exactly: '
            'zero active gravitational mass, so no force on slow stars, while its energy '
            'density still deflects light.',
            'Writing w_t = -f with w_r = 1 interpolates the two limits as tau = 2(1-f), and '
            'f in [0,1] keeps every fitted value inside the energy conditions by construction.'])
    fitted = dict(
        statement='This is an empirical law found by searching a declared grid. It is not derived.',
        law='f(z) = sigmoid(a + b ln z),  tau = 2(1-f),  gamma_chi = 1/tau',
        parameters=s3['best']['parameters'],
        leave_one_out_angle_RMS=s3['best']['loo_angle_fractional_RMS'],
        static_dust_baseline_RMS=s3['baseline_static_dust']['angle_fractional_RMS'],
        look_elsewhere_p=s3['permutation_test']['look_elsewhere_corrected_p'])

    # ---- 2. Can stress engineering ever solve the energy problem? ------------
    # The dominant energy condition bounds |p| <= rho, so each w lies in [-1, 1].
    grid = np.linspace(-1, 1, 401)
    WR, WT = np.meshgrid(grid, grid)
    T = tolman(WR, WT)
    tau_max = float(T.max())
    tau_min_positive = float(T[T > 0].min())
    # rho needed for a fixed stellar force scales as 1/tau, so the best possible
    # energy saving relative to static dust is a factor tau_max.
    shortfall = 4.89e4          # median starlight shortfall measured in JR-9 / consistent with RC-1
    energy_bound = dict(
        dominant_energy_condition='|p| <= rho, so w_r and w_t each lie in [-1, 1]',
        max_tolman_factor=tau_max,
        best_possible_energy_saving_factor=tau_max,
        measured_starlight_shortfall=shortfall,
        remaining_shortfall_after_best_stress=float(shortfall / tau_max),
        verdict='Stress engineering cannot solve the energy problem. The largest factor any '
                'physically allowed stress state can save is tau_max, four. The measured '
                'shortfall is about 5e4. Anisotropic stress fixes the DIRECTION of the '
                'response, not its cost. The energy arrow remains open and is not addressed '
                'by anything in JR-10.')

    # ---- 3. Is the driving variable identified? ------------------------------
    z = np.array([r['zFG'] for r in rows])
    bE = np.array([r['einstein_radius_kpc'] for r in rows])
    Dl = np.array([M.make_lens(r['name']).Dl for r in rows])
    pairs = {}
    for a, b in (('zFG', 'einstein_radius_kpc'), ('zFG', 'Dl_kpc'), ('einstein_radius_kpc', 'Dl_kpc')):
        va = {'zFG': z, 'einstein_radius_kpc': bE, 'Dl_kpc': Dl}[a]
        vb = {'zFG': z, 'einstein_radius_kpc': bE, 'Dl_kpc': Dl}[b]
        rk = lambda v: np.argsort(np.argsort(v))
        pairs[f'{a} vs {b}'] = dict(pearson=float(np.corrcoef(va, vb)[0, 1]),
                                    rank=float(np.corrcoef(rk(va), rk(vb))[0, 1]))
    identification = dict(
        problem='Across these six systems, lens redshift, lens distance and physical Einstein '
                'radius are collinear, so the law predicts well without identifying which '
                'variable drives it.',
        collinearity=pairs,
        redshift_span=[float(z.min()), float(z.max())],
        what_would_break_it='Lenses at matched redshift with different Einstein radii, or '
                            'matched Einstein radii at different redshifts. Neither exists in '
                            'this six-system set.',
        consequence='Deriving a redshift mechanism now would be deriving a variable the data '
                    'has not identified. The next acquisition, not the next derivation, is '
                    'what advances this.')

    # ---- 4. The local prediction the law makes, which is testable now --------
    local = dict(
        statement='The law sends tau to 2 exactly as z goes to zero, which is the '
                  'free-streaming radial null flux state, giving gamma_chi = 1/2.',
        consequence='Around nearby galaxies the companion should deflect light at '
                    '(1 + 1/2)/2 = 3/4 of the rate a static-dust companion of the same '
                    'stellar force would produce.',
        why_it_is_independent='This is a statement about local weak lensing at fixed measured '
                              'rotation. No lens in the fitted set sits near z = 0, and the '
                              'six SLACS angles cannot produce it.',
        caution='The z -> 0 limit follows from the log-linear form, which the six systems '
                'cannot distinguish from the linear form inside the observed range. The '
                'linear form does not give tau = 2 at z = 0. This prediction is therefore '
                'conditional on a form choice made after seeing the fits.')

    # ---- 5. Kinematics and SPARC are untouched, and why that cuts both ways --
    scope = dict(
        kinematics='No V_rms prediction changes for any gamma_chi, because R10 g_chi was '
                   'fitted to stellar kinematics and already is the Phi gradient.',
        sparc='All 149 SPARC galaxies have no lensing measurement, so not one rotation curve '
              'moves and none of them can constrain this law.',
        discipline='Because the law touches only the observable it was fitted to, the '
                   'discipline has to come from elsewhere: the energy-condition band, which '
                   'every required value satisfies; the leave-one-out score; and the '
                   'permutation test over the whole search grid. It does not come from a '
                   'second observable, and that is a genuine weakness.')

    payload = dict(experiment='JR-10 stage 4', exact_general_relativity=exact,
                   empirical_law=fitted, energy_bound=energy_bound,
                   identification=identification, local_prediction=local, scope=scope,
                   code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out / 'derivation-and-bounds.json').write_text(
        json.dumps(payload, indent=2, allow_nan=False) + '\n')

    print('EXACT (no fitting):')
    for i in exact['items']:
        print('  -', i)
    print(f"\nENERGY BOUND")
    print(f"  max Tolman factor allowed by the dominant energy condition: {tau_max:.1f}")
    print(f"  measured starlight shortfall: {shortfall:.2e}")
    print(f"  shortfall remaining after the best possible stress state: "
          f"{energy_bound['remaining_shortfall_after_best_stress']:.2e}")
    print(f"  => {energy_bound['verdict'][:96]}...")
    print(f"\nIDENTIFICATION (n=6 collinearity):")
    for k, v in pairs.items():
        print(f"  {k:<44} pearson {v['pearson']:+.3f}  rank {v['rank']:+.3f}")
    print('\nWROTE', out / 'derivation-and-bounds.json')


if __name__ == '__main__':
    main()
