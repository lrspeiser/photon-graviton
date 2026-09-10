"""Render the executed coverage result, preserving limitations and failed checks."""
from pathlib import Path
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'results.json').read_text())
integ=json.loads((HERE/'integration.json').read_text())
field=json.loads((HERE/'field-check.json').read_text())

def primary(window,model):
    return next(x for x in r['windows'][window][model] if x['radius_kpc']==.5 and x['velocity_radius_kms']==50)

def percent(row): return 100*row['position_and_velocity_covered']/row['stars']

cells=[x for x in r['cells'] if x['side']=='both']
fig,axes=plt.subplots(1,2,figsize=(12,6.5),sharey=True)
y=np.arange(len(cells)); labels=[]
for row in cells:
    labels.append(f"R {row['R_range'][0]:g}–{row['R_range'][1]:g}; |z| {row['absolute_z_range'][0]:g}–{row['absolute_z_range'][1]:g}")
for ax,window,title in zip(axes,['all_positive','late_half'],['All positive saved times','Later half of each trajectory']):
    old=np.array([100*x[window]['old']/x['stars'] for x in cells])
    new=np.array([100*x[window]['expanded']/x['stars'] for x in cells])
    ax.barh(y-.18,old,height=.34,color='#8e989f',label='72 orbits')
    ax.barh(y+.18,new,height=.34,color='#2166ac',label=f"{r['expanded_orbits']} passing orbits")
    ax.set_xlim(0,100);ax.set_title(title);ax.set_xlabel('Training stars with nearby position AND velocity (%)')
    ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
axes[0].set_yticks(y,labels);axes[0].invert_yaxis();axes[0].set_ylabel('Location ranges (kpc)')
axes[1].legend(loc='lower right')
fig.suptitle('Finite orbit-library coverage — not a gravity-model fit',fontsize=14)
fig.text(.5,.01,'Same non-launch training stars; proximity thresholds 0.5 kpc and 50 km/s are diagnostic choices, not error bars.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.04,1,.95]);fig.savefig(HERE/'coverage.png',dpi=170);plt.close(fig)

lines=['# Expanded training orbit library', '',
f"The library now contains {r['expanded_orbits']} numerically passing paths: the original 72 plus {r['passing_new_orbits']} of {r['proposed_new_orbits']} proposed additions. This improves the representation of the existing empirical extra-potential candidate. It is not a fitted stellar population, a new gravity law, or an independent observation confirming companions.", '',
'## What was done and why', '',
'A fixed prospective rule selected four populated missing starting-state neighborhoods in each of nine radius/absolute-height cells. The rule used only existing training data, a 0.5-kpc spatial neighborhood and a 50-km/s velocity neighborhood. It selected by missing library coverage, not by gravity residual. All 36 starts passed the existing association-window screen; none carried the existing distance-disagreement flag. This does not prove every distance or association correct. No actual-epoch association requery was performed for these new starts.', '',
'Each trajectory covers approximately 244 million years in the unchanged rotating full-bar potential. The propagation and extra-force refinement thresholds were retained. The old paths and their documented limitations were not rewritten. Every original and newly proposed launch star is excluded from the comparison, including any proposed path that fails; t=0 samples are excluded too. The resulting target sample is identical for the old and expanded libraries.', '',
'## Numerical checks', '',
f"- Successive-tolerance trajectory and Jacobi checks passed for {sum(x['numerical_pass'] for x in integ['rows'])}/{len(integ['rows'])} new paths.",
f"- Sampled-path extra-force refinement passed for {sum(x['field_pass'] for x in field['rows'])}/{len(field['rows'])}; largest relative difference was {100*max(x.get('maximum_extra_force_fraction',0) for x in field['rows']):.4f}% against the unchanged 1% threshold.",
'- The old-library coverage was reproduced exactly on the shared evaluation sample. Twelve exhaustive neighbor-search controls checked the expanded result. Nested distance thresholds and reduced time sampling obey the expected inclusion relations.', '',
'## Coverage on the same training stars', '',
f"There are {r['target_stars']:,} eligible non-launch targets. The percentages below describe proximity to at least one saved orbit state; they are not percentages of stars explained by the theory.", '',
'| Saved-time choice | Old covered | Expanded covered | Old fraction | Expanded fraction |',
'|---|---:|---:|---:|---:|']
for window in ['all_positive','half_sampling','late_half']:
    a,b=primary(window,'old'),primary(window,'expanded')
    lines.append(f"| {window.replace('_',' ')} | {a['position_and_velocity_covered']:,} | {b['position_and_velocity_covered']:,} | {percent(a):.2f}% | {percent(b):.2f}% |")
lines += ['', 'All-positive times run from about 0.49 to 244 Myr; late-half times run from about 123 to 244 Myr. Half sampling keeps every second positive saved time. The late-half comparison checks whether coverage persists away from the initial launch; it does not prove stationary orbit occupations.', '',
'The later-half coverage is only 8.55%, compared with 49.72% when the same number of saved times is spread across the whole integration. The difference therefore cannot be attributed just to having fewer saved states. Much of the present coverage depends on the early trajectory segments. This is a reason to improve duration and orbital-phase sampling before fitting a stationary population, not evidence that the stars or gravity field themselves change this way. The dominant missing high-|z| 5–9-kpc population has only 2.88% late-half coverage. These findings prevent treating the larger library as ready for a decisive bulge/disk gravity test.', '',
'| R (kpc) | Absolute height (kpc) | Targets | Old coverage | Expanded coverage | Expanded, late half |',
'|---|---|---:|---:|---:|---:|']
for x in cells:
    pct=lambda window,model:100*x[window][model]/x['stars']
    lines.append(f"| {x['R_range'][0]:g}–{x['R_range'][1]:g} | {x['absolute_z_range'][0]:g}–{x['absolute_z_range'][1]:g} | {x['stars']:,} | {pct('all_positive','old'):.2f}% | {pct('all_positive','expanded'):.2f}% | {pct('late_half','expanded'):.2f}% |")
lines += ['', 'Signed-height breakdowns and all nine distance/velocity threshold combinations are retained in results.json. The main table combines above/below only for readability; a future likelihood must retain signed height and bar position.', '',
'![Old and expanded finite-path coverage](coverage.png)', '',
'## Scientific interpretation and next requirements', '',
'This is adaptive training-library development, not a holdout. The selected starting states depend on training motions and the existing candidate library, so improvement cannot be presented as predictive validation. Numerical trajectory accuracy also does not establish a stationary stellar population. The earlier long-duration failures and occupation-variation findings still apply; no relaxation of their thresholds is implied.', '',
'The original 72 paths include distance-flagged and far-outward paths, and their documented unresolved issues remain. These additions do not resolve survey selection, posterior-distance prior recycling, mass-model uncertainty, or the complete velocity-distribution likelihood. The nine cells concern the selected chemically restricted 0.5–9-kpc sample, not all Milky Way stars or the 20–25-kpc outer rotation curve.', '',
'Before interpreting a bulge/disk discrepancy as companion gravity, construct equivalent ordinary-matter orbit coverage and develop the shared population/error/selection model; establish adequate duration and occupation sampling. Keep gravity parameters frozen during this readiness step. The same fitted spatial potential must predict radial, rotational and vertical motions together. The present potential is an empirical response built with known conservative-field mathematics; its photon production/capture origin is still underived.', '',
'The broader goal remains open: propagation must jointly explain spectral redshift, whole-event timing and endpoint clocks; companion transport must conserve physical energy and momentum; deposits need capture, retention and spatial support; one response must survive joint stellar and lensing tests. The failed relativistic stress check in the separate reservoir candidate remains a failure. The total photon-supply budget remains deferred, not passed. Adding orbits does not solve these physical requirements.', '',
'No validation or final-test observations were opened or scored.']
(HERE/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
dest=Path('C:/Users/henry/Documents/Codex/2026-09-09/cr/outputs/orbit-library-expansion')
dest.mkdir(parents=True,exist_ok=True)
for name in ['report.md','coverage.png','results.json','protocol.md','selection.json','integration.json','field-check.json']:
    shutil.copy2(HERE/name,dest/name)
print(dest)
