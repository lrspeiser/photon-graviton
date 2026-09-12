"""Build the control report from all frozen cases, including failures."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
p=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
assert len(r['runs'])==6
assert {(v['truth_b'],v['seed']) for v in r['runs']}=={(b,s) for b in [0,1] for s in [701,702,703]}
lines=['# Can the timing estimator distinguish stretched from unstretched events?', '',
    'This is a small, frozen synthetic control experiment toward the joint redshift, timing and brightness test. It uses 98 previously exposed DES cadence/error patterns, with newly seeded artificial fluxes. No measured supernova flux is fitted. This tests the measuring method, not the physical companion hypothesis.', '',
    '## Why this test is needed', '',
    'A method can give stable numbers and still favor the wrong answer. We therefore create events whose timing rule is known: b=0 gives no duration stretch; b=1 multiplies duration by 1+z. The estimator is not told which answer was injected. It fits observer-time light curves without dividing times by 1+z.', '',
    'The law W proportional to (1+z)^b is a known phenomenological statistical parameterization, not a new physical derivation. An intrinsic source-duration trend can mimic b; these controls deliberately have no such trend.', '',
    '## Frozen design', '',
    'Three seeds (701, 702, 703) each supply paired b=0 and b=1 cases. The paired cases share intrinsic-width draws, peak offsets and noise draws, isolating the change of propagation exponent. The source has an asymmetric Gaussian shape plus a shoulder, outside the estimator family. Intrinsic log-width scatter is 0.1 and nominal peak-to-median-error ratio is 20. No cadence slots are removed.', '',
    'Each case uses 2048 Sobol shape samples and 321 duration-grid nodes. Before running, we specified refinement of the case with the largest absolute error from its injected exponent to 8192 samples with an independent scramble. This diagnostic case selection is not a validation split; it cannot guarantee all cases have converged.', '',
    'The small-control thresholds are absolute median bias at most 0.15 for each truth, median b separation at least 0.7, successful interior fits with negligible width-boundary probability, and refinement changes at most 0.05 in b and 0.1 in the mean centered event log-likelihood. These are engineering acceptance criteria, not a statistical significance claim or thresholds for observational agreement.', '',
    '| Injected b | Seed | True finite-sample slope | Recovered b | Difference from injected b |', '| ---: | ---: | ---: | ---: | ---: |']
for v in r['runs']:
    lines.append(f"| {v['truth_b']} | {v['seed']} | {v['true_finite_sample_slope']:.5f} | {v['fit']['b']:.5f} | {v['bias_relative_to_injected_b']:+.5f} |")
lines+=['', 'The true finite-sample slope differs from the injected population exponent because randomly drawn source widths happen to correlate with redshift. It is diagnostic information available only for artificial data, not a correction applied to real observations.', '',
    'Median recovered b for no stretch: '+str(r['median_b_by_truth']['0'])+'. For injected stretch: '+str(r['median_b_by_truth']['1'])+'. Median separation: '+str(r['median_separation'])+'.', '',
    '## Numerical check and outcome', '',
    'Refined case: '+r['refinement']['selected_case']+'. Change in b: '+str(r['refinement']['b_change'])+'. Mean centered event log-likelihood change: '+str(r['refinement']['mean_event_log_change'])+'. Both estimates are retained in the result file.', '',
    '**Small-control gate: '+('PASS' if r['control_gate_pass'] else 'FAIL')+'.** Median-bias gate: '+str(r['median_bias_gate_pass'])+'; refinement gate: '+str(r['refinement']['gate_pass'])+'; optimization/boundary gate: '+str(r['optimizer_gate_pass'])+'.', '',
    '## What this does not establish', '',
    'Three realizations per truth cannot establish calibrated uncertainty coverage or broad robustness. The samples use one shoulder shape, one signal level, no intrinsic evolution, and simplified band matching. Survey selection, actual filter transmission, source population variation, amplitude-prior sensitivity and additional noise levels remain necessary. Artificial amplitude normalization removes the real brightness-selection problem.', '',
    'Because each event amplitude is marginalized freely, this test supplies no evidence for the physical brightness law. Joint redshift/timing/brightness still needs calibrated flux predictions, source luminosity/distance information and a common transport mechanism. None of the six scientific demonstrations is marked complete.', '',
    '## Reproduce', '', '```powershell',
    'python research_work/results/timing-null-stretch-controls/run.py',
    'python research_work/results/timing-null-stretch-controls/summarize.py', '```', '']
lines += ['## Residual estimator shift', '']
for truth in [0,1]:
    differences=[v['difference_from_true_sample_slope'] for v in r['runs'] if v['truth_b']==truth]
    mean=sum(differences)/len(differences)
    lines.append(f'Mean recovered-minus-true-finite-sample slope for b={truth}: {mean:+.5f}.')
lines += ['', 'Every base-fit recovered slope exceeds its corresponding true finite-sample slope in these six controls. The small-control median criterion can pass despite a residual positive shift. Three paired samples cannot establish its general size or origin; subsequent source-shape, prior and selection sensitivity tests must address it rather than treating this pass as an unbiasedness guarantee.', '']
(HERE/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
print(json.dumps({'control_gate_pass':r['control_gate_pass'],'medians':r['median_b_by_truth'],'separation':r['median_separation']}))
