"""Report shape-quadrature convergence without changing the original gate."""
from pathlib import Path
import json
import sys
import numpy as np
from scipy.special import logsumexp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'timing-population'))
from estimator import log_width_weights

def main():
    r=json.loads((HERE/'timing-population-results.json').read_text())
    details=[]
    for comparison in r['comparisons']:
        left=comparison['left'];right=comparison['right']
        with np.load(HERE/(left+'-event-likelihoods.npz')) as a, np.load(HERE/(right+'-event-likelihoods.npz')) as b:
            assert np.array_equal(a['ids'],b['ids'])
            assert np.array_equal(a['log_width'],b['log_width'])
            lc=a['event_log_likelihoods'];rc=b['event_log_likelihoods']
            lc=lc-lc.max(axis=1)[:,None];rc=rc-rc.max(axis=1)[:,None]
            fit=r['runs'][right]
            post=rc+log_width_weights(b['log_width'],b['redshift'],[fit['a'],fit['b'],np.log(fit['sigma'])])
            post=np.exp(post-logsumexp(post,axis=1)[:,None])
            errors=np.sum(post*abs(lc-rc),axis=1)
            assert np.isclose(errors.mean(),comparison['mean_posterior_weighted_centered_event_log_change'])
            order=np.argsort(errors)[::-1]
            details.append({'left':left,'right':right,'median':float(np.median(errors)),
                'p90':float(np.quantile(errors,.9)),'maximum':float(errors.max()),
                'events_above_0_1':int(np.sum(errors>.1)),
                'largest_five':[{'id':str(b['ids'][i]),'weighted_log_change':float(errors[i])} for i in order[:5]]})
    (HERE/'event-diagnostics.json').write_text(json.dumps(details,indent=2)+'\n',encoding='utf-8',newline='\n')
    lines=['# Supernova timing: shape-integration refinement','',
        'This calculation repairs a numerical prerequisite for the joint redshift, timing and brightness demonstration. It uses artificial light curves at the same 98 previously exposed DES cadence patterns. It does not measure real supernova timing, test brightness, or validate the companion mechanism.','',
        '## What changed and why','',
        'The previous pilot failed its numerical probability-accuracy threshold. We kept its synthetic seed, source shape, noise, priors and thresholds fixed, and increased shape/peak integration from 512 samples to 2048 and 8192. Both comparisons now use the same 321-node observer-time width grid, isolating shape integration. An independent 8192-point scramble checks sensitivity to the quadrature realization. The protocol was written before this run. The original failed result is retained.','',
        'The likelihood formulas are known statistical marginalization and numerical quadrature, not new physical equations. Times are not divided by redshift. All events remain included.','',
        '| Calculation | Recovered timing exponent b | Intrinsic log-width scatter |','| --- | ---: | ---: |']
    for name,fit in r['runs'].items():lines.append(f"| {name} | {fit['b']:.6f} | {fit['sigma']:.6f} |")
    lines+=['','Injected exponent: 1. The finite-sample true-width slope is '+str(r['slope_of_true_injected_widths_in_finite_sample'])+'. Agreement in one realization cannot establish population bias or uncertainty coverage.','',
        '| Comparison | Change in b (limit 0.05) | Mean event log-likelihood change (limit 0.1) |','| --- | ---: | ---: |']
    for c in r['comparisons']:lines.append(f"| {c['left']} to {c['right']} | {c['absolute_b_change']:.6f} | {c['mean_posterior_weighted_centered_event_log_change']:.6f} |")
    lines+=['', '**Original numerical pilot gate: '+('PASS' if r['numerical_pilot_gate_pass'] else 'FAIL')+'.** No threshold was relaxed. The gate also requires successful interior optimization and negligible population probability outside the width bounds.','',
        'Event-level diagnostics are saved separately. The gate concerns a mean over events; passing does not assert every individual event meets 0.1.','',
        '## What remains before interpreting observations','',
        'Separate width-grid and prior-boundary sensitivity, repeated injections of both b=0 and b=1 with different source shapes and signal levels, calibrated uncertainty coverage, actual filter transmission, source evolution and survey selection remain unverified. The estimator integrates an arbitrary brightness amplitude for each event, so this timing calculation cannot itself test the modelÃ¢â‚¬â„¢s brightness prediction. A joint test must retain physical flux calibration and independently constrain source luminosities and distances.','',
        'If this numerical gate fails, local adaptive integration or further convergence work is needed before observational inference. If it passes, proceed to the remaining convergence and injection checks; do not describe the model as matching real supernovae.','',
        '## Reproduce','', 'From the repository root:', '', '```powershell',
        'python research_work/results/timing-population/run.py --protocol research_work/results/timing-quadrature-refinement/protocol.json --output research_work/results/timing-quadrature-refinement',
        'python research_work/results/timing-quadrature-refinement/check_width_grid.py', 'python research_work/results/timing-quadrature-refinement/summarize.py','```','']
    lines += ['## Individual-event limits', '', '| Comparison | Events above 0.1 | Largest change |', '| --- | ---: | ---: |']
    for d in details: lines.append(f"| {d['left']} to {d['right']} | {d['events_above_0_1']}/98 | {d['maximum']:.6f} |")
    lines += ['', 'The mean-based gate passes even though some individual curves remain less stable. These events are retained and should receive attention in subsequent sensitivity checks; no individual-event accuracy guarantee is claimed.', '']
    width=json.loads((HERE/'width-grid-results.json').read_text())
    lines += ['## Separate duration-grid check', '', f"The nested 161-to-321-node check uses identical refined shape samples. Change in b: {width['absolute_b_change']:.8g}; mean event log-likelihood change: {width['mean_event_log_change']:.6f}. Original limits: 0.05 and 0.1. Result: {'PASS' if width['nested_grid_check_pass'] else 'FAIL'}. This does not check a still finer 641-node grid.", '', 'Reproduce with `python research_work/results/timing-quadrature-refinement/check_width_grid.py` before running this summary script.', '']
    (HERE/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
    print(json.dumps({'numerical_gate_pass':r['numerical_pilot_gate_pass'],'event_diagnostics':details}))

if __name__=='__main__':main()
