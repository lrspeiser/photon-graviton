#!/usr/bin/env python3
"""JR-9 B: run the declared screening statistic on the real 46-element joint vector.

Builds the observed / baseline / candidate vectors in the units the telescope
actually reports: 40 seeing-convolved annular V_rms values in km/s across six
SLACS lenses, then six Einstein angles in arcsec. The covariance carries each
lens's released kinematic covariance; the lensing block uses a declared scan of
assumed fractional angle uncertainty, because the catalog SIE summaries ship no
error and JR-1's 5 percent was an optimization scale.

The candidate is the within-budget spherical redistribution from JR-9D: for each
lens, the shell that closes its Einstein-angle gap at the lowest kinematic cost
while staying inside the frozen companion's own mass budget.

CIRCULARITY, STATED PLAINLY: that shell was sized *by requiring the angle to
match*, so the lensing part of the template reproduces the lensing residual by
construction and its fitted amplitude near one carries no information. The
informative output is the motion-block amplitude and the amplitude difference:
does the resolved stellar kinematics independently support the same correction
that repairs the lens?

Nuisance marginalization uses the numerical Jacobian of all 46 predictions with
respect to every lens's (stellar offset, anisotropy) pair, with the published
stellar-mass error and the declared anisotropy width as the prior covariance.

    python joint_screen.py --output-dir <fresh dir>
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from scipy.linalg import block_diag

import model as M
import response_kernels as rk
from joint_template_test import evaluate

ANGLE_SIGMA_SCAN = (0.01, 0.02, 0.03, 0.05, 0.10)
KERNELS = M.ROOT / 'research_work/results/cross-prediction-response/run-v2-kernels/response-kernels.json'


def budget(lens, dm):
    _, A, _, rt, _ = rk.base_terms(lens, dm)
    return A * rt / rk.G


def pick_shell(lens, entry):
    """Cheapest angle-closing shell that stays inside the companion's own mass budget."""
    cap = budget(lens, entry['kinematic_best']['dm_dex'])
    inside = [c for c in entry['shell_costs'] if abs(c['shell_mass_Msun']) <= cap]
    if not inside:
        return None, cap
    return min(inside, key=lambda c: c['delta_penalty']), cap


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    if out.exists():
        raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True)

    kern = {e['name']: e for e in json.loads(KERNELS.read_text())['lenses']}
    lenses = [M.make_lens(n) for n in M.LENSES]

    observed, baseline, candidate, labels, covs = [], [], [], [], []
    jac_cols, prior_var, chosen = [], [], []
    ang_obs, ang_base, ang_cand = [], [], []

    for lens in lenses:
        nu = M.JR1_NUISANCE[lens.name]
        dm0, b0 = float(nu['dm']), float(nu['beta'])
        v_base = M.vrms(lens, dm0, b0)
        a_base = M.einstein_angle(lens, dm0)
        shell, cap = pick_shell(lens, kern[lens.name])
        if shell is None:
            raise RuntimeError('No within-budget repair for ' + lens.name)
        m, r_s = shell['shell_mass_Msun'], shell['r_kpc']
        v_cand = rk.vrms_with_shell(lens, dm0, b0, m, r_s)
        a_cand = rk.angle_with_shell(lens, dm0, m, r_s)
        observed.append(lens.y)
        baseline.append(v_base)
        candidate.append(v_cand)
        labels.append(np.full(len(lens.y), 'motion'))
        covs.append(np.asarray(lens.cov))
        ang_obs.append(lens.theta)
        ang_base.append(a_base)
        ang_cand.append(a_cand)
        chosen.append(dict(name=lens.name, shell_r_kpc=r_s, shell_r_over_Re=shell['r_over_Re'],
                           shell_mass_Msun=m, shell_over_stellar_mass=shell['shell_over_stellar_mass'],
                           shell_over_companion_budget=float(m / cap),
                           delta_penalty=shell['delta_penalty']))

    n_kin = sum(len(l.y) for l in lenses)
    n_tot = n_kin + len(lenses)
    observed = np.r_[np.concatenate(observed), np.array(ang_obs)]
    baseline = np.r_[np.concatenate(baseline), np.array(ang_base)]
    candidate = np.r_[np.concatenate(candidate), np.array(ang_cand)]
    labels = np.r_[np.concatenate(labels), np.full(len(lenses), 'lensing')]

    # Numerical nuisance Jacobian: 46 predictions vs each lens's (dm, beta).
    offset = 0
    for i, lens in enumerate(lenses):
        nu = M.JR1_NUISANCE[lens.name]
        dm0, b0 = float(nu['dm']), float(nu['beta'])
        nb = len(lens.y)
        for which, step, sigma in (('dm', 1e-3, lens.mass_log_error),
                                   ('beta', 1e-3, M.BETA_PRIOR_SIGMA)):
            col = np.zeros(n_tot)
            if which == 'dm':
                vp, vm = M.vrms(lens, dm0 + step, b0), M.vrms(lens, dm0 - step, b0)
                ap_, am = M.einstein_angle(lens, dm0 + step), M.einstein_angle(lens, dm0 - step)
                col[n_kin + i] = (ap_ - am) / (2 * step)
            else:
                vp, vm = M.vrms(lens, dm0, b0 + step), M.vrms(lens, dm0, b0 - step)
            col[offset:offset + nb] = (vp - vm) / (2 * step)
            jac_cols.append(col)
            prior_var.append(sigma ** 2)
        offset += nb

    B = np.column_stack(jac_cols)
    S = np.diag(prior_var)

    payload = dict(
        experiment='JR-9B', protocol='../PROTOCOL.md',
        scope='Screening diagnostic on the real joint vector. Not a full image or '
              'dynamical likelihood, and not a p-value.',
        vector=dict(n_total=int(n_tot), n_motion=int(n_kin), n_lensing=len(lenses),
                    motion_units='km/s seeing-convolved annular V_rms',
                    lensing_units='arcsec catalog SIE Einstein angle'),
        baseline='Published JR-1 R10 prediction at its own nuisance convention',
        candidate='JR-9D within-budget angle-closing spherical redistribution',
        candidate_shells=chosen,
        circularity_note='The candidate shell mass was solved from the observed angle, so the '
                         'lensing-block amplitude is fitted by construction and is not evidence. '
                         'Read the motion-block amplitude and the amplitude difference.',
        runs=[], code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        vendored_tool_sha256=hashlib.sha256(
            (Path(__file__).parent / 'joint_template_test.py').read_bytes()).hexdigest())

    for frac in ANGLE_SIGMA_SCAN:
        ang_cov = np.diag((frac * np.array(ang_obs)) ** 2)
        C = block_diag(*covs, ang_cov)
        plain = evaluate(observed, baseline, candidate, C, labels)
        marg = evaluate(observed, baseline, candidate, C, labels, B, S)
        payload['runs'].append(dict(assumed_fractional_angle_sigma=frac,
                                    without_nuisance_marginalization=plain,
                                    with_nuisance_marginalization=marg))
        print(f"sigma_theta={frac:.0%}  "
              f"sqrt(tWt)={np.sqrt(marg['expected_squared_SNR_for_candidate_vs_baseline']):7.2f}  "
              f"A={marg.get('common_diagnostic_amplitude', float('nan')):+.3f}"
              f"+-{marg.get('common_amplitude_standard_error', float('nan')):.3f}  "
              f"A_motion={marg.get('diagnostic_amplitude_motion', float('nan')):+.3f}  "
              f"A_lens={marg.get('diagnostic_amplitude_lensing', float('nan')):+.3f}  "
              f"diff={marg.get('amplitude_difference_in_standard_errors', float('nan')):+.2f} sigma",
              flush=True)

    (out / 'joint-screen.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    np.savez(out / 'joint-vector.npz', observed=observed, baseline=baseline, candidate=candidate,
             block=labels, covariance=block_diag(*covs, np.diag((0.03 * np.array(ang_obs)) ** 2)),
             nuisance_jacobian=B, nuisance_covariance=S)
    print('WROTE', out / 'joint-screen.json')


if __name__ == '__main__':
    main()
