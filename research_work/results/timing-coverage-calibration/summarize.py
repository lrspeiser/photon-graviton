"""Report complete calibration without turning a screening pass into validation."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
p=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
assert len(r['cases'])==160 and len(r['refinements'])==8
lines=['# Expanded synthetic timing calibration results','','All 160 frozen base cases and eight selected numerical refinements completed. These are artificial fluxes on exposed DES cadence patterns, not measurements of astronomical time dilation.','','| Shape | SNR | True b | Mean bias | Monte Carlo SE | Nominal 95% coverage | Exact 95% interval for coverage | Invalid fits |','| --- | ---: | ---: | ---: | ---: | --- | --- | ---: |']
criteria=p['calibration']['diagnostic_thresholds'];screens=[]
for s in r['summary']:
    v=s['coverage']['0.95'];lo,hi=v['exact_binomial_95_interval']
    lines.append(f"| {s['family']} | {s['snr']} | {s['truth_b']} | {s['mean_bias']:+.4f} | {s['bias_monte_carlo_se']:.4f} | {v['covered']}/{v['total']} | [{lo:.3f}, {hi:.3f}] | {s['invalid_fits']} |")
    screens.append(abs(s['mean_bias'])<=criteria['maximum_absolute_mean_bias'] and v['fraction']>=criteria['minimum_nominal_95_coverage'] and s['invalid_fits']<=criteria['maximum_optimizer_or_boundary_failures'])
lines+=['','Screening thresholds: absolute mean bias <=0.1, at least 85% nominal-95 coverage, and zero invalid fits in each cell. These were declared before execution. They are not precision validation of 95% coverage. Per-cell screening results: '+str(screens)+'.','','## Selected numerical refinements','','| Refined case | Change in b | Mean event log-likelihood change | Pass |','| --- | ---: | ---: | --- |']
for label,v in r['refinements'].items():lines.append(f"| {label} | {v['b_change']:.5f} | {v['mean_event_log_change']:.5f} | {v['numerical_pass']} |")
lines+=['','Selected refinements retain both estimates and are not substituted selectively into the base coverage counts. Refining the largest estimation error does not guarantee that every other case has converged numerically.','','## Interpretation','','The exact binomial intervals show how uncertain these coverage fractions remain with 20 seeds per cell. Seeds are paired across cells; aggregate counts cannot be treated as independent. Failed fits count as uncovered and are explicitly listed. The bias Monte Carlo standard error describes variation over these synthetic draws, not total astrophysical uncertainty.','','Source evolution, real filter transmission, spectroscopic selection, amplitude-prior sensitivity and physically calibrated brightness remain outside this batch. A numerical or screening pass does not establish that photon-companion conversion explains observed redshift, timing or brightness. All six scientific goals remain open.','']
(HERE/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
print(json.dumps({'cell_screening':screens,'numerical_refinements_pass':all(v['numerical_pass'] for v in r['refinements'].values())}))
