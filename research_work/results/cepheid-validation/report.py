"""Publish reserved-validation results with conditional scope and exposure status."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
s=json.loads((HERE/'results.json').read_text());p=json.loads((HERE/'protocol.json').read_text())
bins=json.loads((HERE/'validation-bins.json').read_text());v=json.loads((HERE/'verification.json').read_text())
names={'original_ordinary':'Original ordinary matter','original_completion':'Primary: original added field',
       'balanced_ordinary':'Adjusted ordinary matter','balanced_completion':'Secondary: added field with frozen mass adjustment'}
table='\n'.join(f"| {names[model]} | {score['rms_kms']:.2f} | {score['bias_kms']:.2f} | {score['underpredicted_bins']} / {score['bins']} |" for model,score in s['scores'].items())
bin_lines=[]
for r in bins:
    if r.get('moment_valid'):
        bin_lines.append(f"| {r['R_low_kpc']}–{r['R_high_kpc']} | {r['n']} | {r['jeans_proxy_kms']:.2f} | {r['original_completion_vc_kms']:.2f} | {r['balanced_completion_vc_kms']:.2f} |")
    else:bin_lines.append(f"| {r['R_low_kpc']}–{r['R_high_kpc']} | {r['n']} | Not scored | — | — |")
bin_table='\n'.join(bin_lines)
matched_table='\n'.join(f"| {names[model]} | {score['rms_kms']:.2f} | {s['scores'][model]['rms_kms']:.2f} |" for model,score in v['training_metrics_on_same_bin_indices'].items())
valid=[r for r in bins if r.get('moment_valid')]
R=np.array([r['R_mean_kpc'] for r in valid])
fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
ax.scatter(R,[r['jeans_proxy_kms'] for r in valid],color='#203a52',s=40,label='Reserved validation: conditional Jeans estimate',zorder=5)
for model,color,ls in [('original_ordinary','#999999',':'),('balanced_ordinary','#b19b81',':'),('original_completion','#b95727','-'),('balanced_completion','#286c91','-')]:
    ax.plot(R,[r[model+'_vc_kms'] for r in valid],color=color,ls=ls,label=names[model])
for r in valid:ax.annotate(str(r['n']),(r['R_mean_kpc'],r['jeans_proxy_kms']),xytext=(0,7),textcoords='offset points',ha='center',fontsize=8)
ax.set(xlabel='Distance from Galactic axis (kpc)',ylabel='Circular-speed estimate or prediction (km/s)',title='Reserved Cepheid validation: parameters fixed before evaluation')
ax.grid(alpha=.2);ax.legend(fontsize=8,loc='lower left')
fig.savefig(HERE/'comparison.png',dpi=170);plt.close(fig)
primary=s['scores']['original_completion'];secondary=s['scores']['balanced_completion']
improvement=100*(1-secondary['rms_kms']/primary['rms_kms'])
report=f'''# Frozen Cepheid validation results

**The reserved-star test is complete for the declared approximate pipeline.** From 433 measurement/mode candidates, the fixed spatial and velocity cuts retain **{s['selected_validation_stars']} validation stars**, with **{len(valid)} scored radial bins**. No parameter was fitted to these outcomes. The final test role remains untouched.

| Frozen model | Validation bin RMS (km/s) | Mean predicted minus inferred speed (km/s) | Bins underpredicted |
|---|---:|---:|---:|
{table}

The secondary mass choice changes RMS relative to the primary by a {improvement:.1f}% reduction (negative would mean worsening). This quantifies transfer to held-out stars under the same assumptions, not physical validation of the mass choices or photon origin. The primary predicts lower speeds in {primary['underpredicted_bins']} of {primary['bins']} scored bins; residual disagreement remains.

![Reserved validation comparison](comparison.png)

## What was frozen before opening outcomes

Commit **e5cf484** contains the protocol, complete evaluation program and component adapter before the first validation-star distance/velocity calculation. Hashes cover those programs, the calibration and frame, field coefficients, bar/nuclear caches, original catalog, role assignments, and the previously fitted component scales. The adapter's component class definitions were checked to match the prior calculation exactly.

The primary is the original conservative additional field. The secondary retains the previously fitted component multipliers: stellar disks **{p['balanced_component_scales'][0]:.6f}**, gas disks **{p['balanced_component_scales'][1]:.6f}**, central stars **{p['balanced_component_scales'][2]:.6f}**. The black hole and spatial component profiles remain fixed. No potential-stretch parameter is used. Both ordinary-matter-only counterparts are included as controls.

Distances use the same published period–Wesenheit calibration and the same solar frame as training. The same moment-error scenario includes an independent 7-percent distance uncertainty. Selection remains `6 <= R <= 18 kpc`, `|phi| <= 30 degrees`, `|z| <= 0.5 kpc`, and `|vz| <= 100 km/s`. Each one-kpc bin needs at least five stars and physically valid error-corrected moments. Every bin is listed below; there is no residual-based exclusion or merging.

## Per-bin outcome

All numerical speed columns are km/s. The inferred speed is a conditional population Jeans estimate, not a directly observed acceleration.

| Radius interval (kpc) | Validation stars | Inferred from motions | Primary prediction | Secondary prediction |
|---|---:|---:|---:|---:|
{bin_table}

For context, the following uses training metrics restricted to the same bin indices. The stars and exact radii still differ, so these are descriptive comparisons rather than a paired measurement test.

| Model | Training RMS on matched bins (km/s) | Validation RMS (km/s) |
|---|---:|---:|
{matched_table}

## Formula status and interpretation

The period–Wesenheit calibration, coordinate transforms, moment-error propagation, and simplified Jeans equation are **known methods**, not unique formulas of this project. The gravitational calculation uses the **previously archived empirical response in a known conservative field structure**. Its identification with photon-generated companions remains a hypothesis. The secondary ordinary-matter normalizations were fitted on training/exposed data and have not been established by independent mass constraints.

This opening is narrower than full theory validation: the larger stellar likelihood is unfinished, but the frozen pipeline can still be tested honestly on reserved stars. The test answers whether its fixed conditional predictions carry over, without claiming that shared calibration, equilibrium, radial-profile, selection or mass assumptions are correct. No acceptable-fit threshold based on complete observational uncertainty was available, so these RMS values are not turned into sigma-level acceptance or rejection.

The reserved stars share Gaia calibration and the same Galactic population with training. Published aggregate Cepheid curves had already been examined. Consequently this is **a held-out star subset, not a new independent catalog or a test of all systematic errors**. It says nothing new about photon production, supernova event timing, lossless transport, capture, long-lived storage, lensing, or total energy supply. Those remain separate requirements; the total photon-supply budget remains deferred.

## Verification and exposure record

All frozen hashes remain unchanged. Selected identifiers are unique, belong to validation, and are disjoint from both training and final test. Spatial cuts, bin counts, minimum-count rules and stored metrics are verified. Coarse/refined prediction differences are at most **{max(s['max_coarse_refined_change_kms'].values()):.4f} km/s**, below the declared 0.1 km/s numerical tolerance. Original empirical parameters and the final test assignment are unchanged.

**These validation outcomes are now exposed.** They must not be described as fresh validation after any future tuning. The final test should only be opened after the remaining method and model decisions are fixed; this result does not authorize repeated tuning against that sample.

Reproduce with `evaluate.py`, `verify.py`, then `report.py`, using the recorded local inputs. `freeze.py` refuses to overwrite an existing protocol. Raw and derived validation-star rows remain in the ignored data cache. The tracked files retain summaries, code, per-bin predictions and provenance.
'''
(HERE/'report.md').write_text(report,encoding='utf-8',newline='\n')
print(json.dumps(dict(secondary_RMS_reduction_percent=improvement),indent=2))
