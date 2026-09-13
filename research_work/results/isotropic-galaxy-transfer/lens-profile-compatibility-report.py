"""Summarize compatibility without turning a forced fit into a prediction."""
import json
from pathlib import Path
P=Path(__file__).resolve().parent;d=json.loads((P/'lens-profile-compatibility-results.json').read_text());rows=d['rows']
lines=['# Can redistribution reconcile stellar motion and lensing?','',
       'Under the retained spherical geometry, optical distances, constant stellar mass-to-light shape, and orbital-anisotropy bounds, the tested positive mixture family does not reconcile three of six lens galaxies. Other systems have less tension; one admits a low-residual joint fit at the orbital bound. This is a conditional inverse calculation, not a prediction or a proof that all redistribution is impossible.','',
       '## Fixed assumptions and added flexibility','',
       'The exact-one-third retention and its three capture constants are unchanged. The existing deposited inventory is mixed among seven conservative scale factors 0.1, 0.25, 0.5, 1, 2, 4 and 10. Weights are nonnegative and sum to one. This is much more flexible than the earlier shared migration law; each target fits its own weights, so a successful case would only establish a possible target profile. All stellar bins now enter the full covariance fit.','',
       'The catalog lens angle is imposed exactly by eliminating stellar mass from the lens equation at that angle. This is an algebraic target constraint, not a lens-angle uncertainty model or a new predicted angle. The implied stellar mass must satisfy the previous mass bounds. A zero lens residual in this calculation is therefore not evidence of success on its own.','',
       '## Stellar-motion cost of enforcing the lens angle','',
       'First population proxy (Chabrier); smaller covariance-weighted residual sums indicate closer stellar-motion agreement. No calibrated significance is claimed.','',
       '| Galaxy | Original profile: motions only | Original profile: exact lens | Flexible profile: exact lens | Largest motion residual / plotted error, flexible |','|---|---:|---:|---:|---:|']
names=list(dict.fromkeys(r['Name'] for r in rows))
for name in names:
    rr={r['fit']:r for r in rows if r['Name']==name and r['population']=='Chabrier'}
    a,b,c=[rr[k] for k in ['original_free','original_exact_lens','mixture_exact_lens']]
    lines.append(f"| {name} | {a['stellar_chi2']:.2f} | {b['stellar_chi2']:.2f} | {c['stellar_chi2']:.2f} | {c['max_marginal_sigma']:.2f} |")
lines += ['', 'The Salpeter cases give the same qualitative separation, with flexible exact-lens scores 177.13, 29.92, 131.32, 251.83, 3.84 and 7.84 in the listed order. Both sets of every fit and velocity prediction are archived in the JSON. Largest marginal residuals are descriptive bin errors, not independent significances; covariance is used in the summed score.','',
          'J0037, J1204 and J1402 retain severe stellar-motion mismatches even with target-fitted redistribution. J1112 already fits its lens nearly unchanged but retains a stellar-profile mismatch. J1630 also has little lens/motion conflict at this level; it is not a perfect velocity fit. J1621 permits a close fit within the tested family, but requires beta at the upper bound 0.45. None of these statements verifies a nonnegative, stable phase-space distribution.','',
          '## Where the inverse fit tries to move the deposits','',
          '| Galaxy | Mixture weights above 1% (weight at radius scale) | Orbital bound? | M_new / M_old at one effective radius |','|---|---|---|---:|']
for r in rows:
    if r['fit']!='mixture_exact_lens' or r['population']!='Chabrier':continue
    description=', '.join(f'{100*w:.1f}% at {s:g}×' for w,s in zip(r['weights'],d['scales']) if w>.01)
    change=r['cumulative_mass_changes'][1]
    lines.append(f"| {r['Name']} | {description} | {r['beta_boundary']} | {change['new_mass']/change['baseline_mass']:.4f} |")
lines += ['', 'A scale of 10 means moving each radius to ten times its former value, not putting all mass in a shell at ten effective radii. Several weights sit at the largest available scale, so these are restricted-family optima, not measured preferred radii. J1402 mixes compact and extended material while hitting the orbital bound. The different preferred mixtures do not yet define a universal law. Full cumulative corrections are recorded at 0.5, 1, 2 and 5 effective radii.','',
          '## Alternative diagnosis: geometry or lensing response','',
          'Holding each original motions-only fit fixed, the multiplier required on D_ls/D_s (or equivalently on the entire bending amplitude) at the catalog angle is:', '', '| Galaxy | Required multiplier |','|---|---:|']
for r in rows:
    if r['fit']=='original_free' and r['population']=='Chabrier':lines.append(f"| {r['Name']} | {1/(1+r['lens_equation_fractional_residual']):.4f} |")
lines += ['', 'This is a required-value diagnostic, not an adopted formula. It affects all bending, not only the companion term. Three problematic galaxies need roughly 17–21% less bending at the observed angle, while two others need almost no change. A common multiplier is therefore not an exact repair. Changes to optical geometry must also remain consistent with the redshift/brightness branch, and changes to gravitational response require an independently specified field model. Do not interpret these numbers as photon-energy conversion fractions.','',
          '## What follows','',
          'Preserve the one-third law while next testing the common optical geometry and non-spherical stellar/orbital assumptions. In these systems reducing the companion concentration alone cannot compensate for the ordinary stellar contribution required by the motions within the tested family. This is a reason to audit the light-to-mass/orbit/lensing mapping, not yet evidence for a modified light-bending law. Do not add per-galaxy bending multipliers to claim success.','',
          'The fixed photometric proxy and geometry uncertainties were not marginalized; arbitrary density functions, radial anisotropy and nonspherical dynamics were not exhausted. Thus the result does not establish the stronger claim that no reasonable ordinary-matter model can fit. Nor does it test clusters. All six goals remain open.','',
          '## Verification','',
          'The extended cumulative-mass setup reproduces the prior fixed-parameter stellar velocities to maximum '+f"{max(v['baseline_max_velocity_drift_kms'] for v in d['verification']):.3g}"+' km/s. Mass inventory is preserved by the exact dilation/CDF identity. Optimizations use multiple starts, then direct Jeans coefficients for final evaluation/refinement, not only beta interpolation. Exact-lens residuals are checked below 1e-9 fractional; nonnegative weights and unit sum are checked. Hashes track the frozen capture, optical and photometric inputs. No profile was claimed physically supported merely because its Jeans moments could be computed.','',
          'Run lens-profile-compatibility.py then lens-profile-compatibility-report.py. The pre-run specification is lens-profile-compatibility-protocol.md.']
(P/'lens-profile-compatibility-report.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[10:23]))
