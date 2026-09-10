"""Create a source-linked account and plots from completed paired diagnostics."""
from pathlib import Path
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'results.json').read_text())
checks=json.loads((HERE/'checks.json').read_text())
datasets={kind:json.loads((HERE/f'{kind}-integration.json').read_text()) for kind in ['ordinary','full']}

def item(kind,window):return next(x for x in r['tables'] if x['kind']==kind and x['window']==window)

fig,axes=plt.subplots(1,2,figsize=(12,5.5))
groups=[(['early_equal250','middle_equal250','late_equal250'],['Early\n0–122 Myr','Middle\n428–550 Myr','Late\n856–978 Myr'],'Equal 250 saved times per orbit'),
        (['first_quarter500','full_history_equal500','full_history2000'],['First 244 Myr\n500 times','Full 978 Myr\n500 times','Full 978 Myr\n2,000 times'],'Accumulated coverage and sampling control')]
for ax,(windows,labels,title) in zip(axes,groups):
    xx=np.arange(3)
    for offset,kind,color,label in [(-.18,'ordinary','#8e989f','Ordinary matter'),(.18,'full','#2166ac','Empirical extra potential')]:
        vals=[100*item(kind,w)['fraction'] for w in windows]
        bars=ax.bar(xx+offset,vals,width=.34,color=color,label=label)
        ax.bar_label(bars,labels=[f'{v:.1f}%' for v in vals],fontsize=9,padding=3)
    ax.set_xticks(xx,labels);ax.set_ylim(0,100);ax.set_title(title,fontsize=11)
    ax.set_ylabel('Training targets with a nearby position and velocity (%)')
    ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
axes[1].legend(loc='upper right',fontsize=9)
fig.suptitle(f"{len(r['common_passing_pairs'])} common passing orbit pairs — a coverage diagnostic, not a gravity ranking",fontsize=13)
fig.text(.5,.015,'Same 27,606 non-launch training targets; fixed 0.5-kpc / 50-km/s proximity thresholds; no fitted populations.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.04,1,.93]);fig.savefig(HERE/'coverage.png',dpi=170);plt.close(fig)

lines=['# Paired long-duration orbit diagnostics','',
'Longer integration helps cumulative coverage but does not solve the late-window representation problem. For the eight common passing pairs, the candidate covers 21.31% of targets over the first 244 Myr and 33.03% over 978 Myr at the same 500 saved times per orbit. Yet coverage during the last equal-size window is only 1.02% for the candidate and 0.49% for ordinary matter. This is evidence that duration alone has not produced a reliable population representation; it is not a comparison of gravitational likelihoods.', '',
'This calculation tests whether longer integrations improve the finite orbit representation needed for the bulge/disk comparison. Nine fixed training-derived starts were followed for approximately 978 million years under both ordinary-matter gravity and the existing empirical extra-potential candidate. The initial stars, ordinary-matter components, rotating frame, physical parameters and numerical gates are shared. No gravity parameter or stellar likelihood was fitted.','',
'## Selection and numerical checks','',
'The starts are expansion indices 0, 4, 8, 12, 16, 20, 24, 28 and 32: the first populated-gap launch in each of nine radius/absolute-height cells. Selection used the candidate library, so even a matched comparison of these starts is not neutral model selection. All passed the prior association-window screen and none carried the existing distance-disagreement flag; exact-epoch associations and distance likelihoods are not newly established.','',
'| Model | Complete saved paths | Trajectory/Jacobi passes |','|---|---:|---:|']
for kind,data in datasets.items():
    lines.append(f"| {kind} | {sum(x['complete_trajectory'] for x in data['rows'])}/9 | {sum(x['numerical_pass'] for x in data['rows'])}/9 |")
lines += ['',f"The candidate extra-force refinement passed for {sum(x['field_pass'] for x in checks['rows'])}/9 probes; its largest sampled relative difference was {100*max(x.get('maximum_extra_force_difference',0) for x in checks['rows']):.4f}% against the unchanged 1% threshold. Its first 244-Myr segment agreed with the archived short path within the existing position/velocity gates for {sum(x['short_path_pass'] for x in checks['rows'])}/9 probes. Ordinary components are the same previously audited fields; this is not a new certification of their mass assumptions or discretization.", '',
f"The common passing indices are {r['common_passing_pairs']}. Excluded paired indices are {r['excluded_pairs']}. Every attempt is retained. The common passing subset is selected on numerical reliability and cannot be treated as an unbiased stellar sample or used to rank gravity. Leaving an interpolation domain would be a computational limitation, not evidence of physical escape.", '',
'## Does coverage persist later in the journey?','',
'The same 27,606 eligible training stars are used throughout, excluding all 108 original/expansion launch stars. A star is covered only if some saved state lies within both 0.5 kpc in position and 50 km/s in velocity. These thresholds are diagnostic choices, not measurement error bars. Coverage is not the fraction of stars explained by a gravitational theory. All initial-time samples are excluded.', '',
'| Time selection | Samples/orbit | Ordinary covered | Candidate covered | Ordinary fraction | Candidate fraction |',
'|---|---:|---:|---:|---:|---:|']
for w in ['early_equal250','middle_equal250','late_equal250','first_quarter500','full_history_equal500','full_history2000']:
    a,b=item('ordinary',w),item('full',w)
    lines.append(f"| {w.replace('_',' ')} | {a['saved_times_per_orbit']} | {a['covered']:,} | {b['covered']:,} | {100*a['fraction']:.2f}% | {100*b['fraction']:.2f}% |")
lines += ['', 'The first three rows have equal numbers of saved states: roughly 0.49–122, 428–550 and 856–978 Myr. The last three distinguish adding physical duration from merely adding sample count. Full-history equal500 uses every fourth positive saved time; full-history2000 uses every positive saved time. A cumulative union can only gain coverage when points are added, so its increase alone is not evidence for stable orbital populations.', '',
'![Coverage by time window and duration](coverage.png)','',
'## Orbit occupation, with failures retained','',
'Known statistical definition, not a new physical formula: TV = one half of the sum over bins of the absolute difference between the early and late time fractions. Zero means identical binned fractions; one means no overlap. The radius/signed-height/bar-angle bins are exactly those in the earlier duration audit. TV depends on those bins and is not a probability of model failure. Small values do not prove stationarity.', '',
'| Start | Model | Numerical pass | TV at 244 Myr | TV at 489 Myr | TV at 978 Myr | R/height-only TV at 978 Myr | Bar-angle range (turns) |',
'|---:|---|---|---:|---:|---:|---:|---:|']
for index in r['selected_indices']:
    for kind,data in datasets.items():
        row=next(x for x in data['rows'] if x['seed_index']==index)
        if not row['complete_trajectory']:
            lines.append(f"| {index} | {kind} | False | unavailable | unavailable | unavailable | unavailable | unavailable |")
        else:
            t=row['occupation_windows']
            lines.append(f"| {index} | {kind} | {row['numerical_pass']} | {t[0]['spatial_TV']:.3f} | {t[1]['spatial_TV']:.3f} | {t[2]['spatial_TV']:.3f} | {t[2]['R_z_TV']:.3f} | {row['bar_angle_range_turns']:.3f} |")
lines += ['', 'Occupation values on rows marked False are descriptive outputs of unverified long paths; they cannot support precise dynamical conclusions. Bar-angle ranges use sampled unwrapped positions. Large ranges alone do not establish mixing, and restricted ranges can correspond to physical libration. Half-time occupation sensitivities, full integration attempts, signed-height coverage cells and source hashes are retained in the JSON results.', '',
'## Why a billion years is not automatically sufficient', '',
'Supplementary kinematic diagnosis added after the integrations began, without changing their selection: in the existing rotating frame, the known identity is d(phi_bar)/dt = (x v_y - y v_x)/R^2 - Omega_bar. If this instantaneous relative rate remained constant, a full relative turn would take 2 pi divided by its absolute value. This is a standard rotating-frame calculation, not a new formula from our theory or a measured return period. Real rates change along noncircular orbits.', '',
'| Start | Initial radius (kpc) | Initial relative angular rate (km/s/kpc) | Constant-rate full-turn timescale (Myr) |',
'|---:|---:|---:|---:|']
for x in r['initial_phase_scales']:
    lines.append(f"| {x['seed_index']} | {x['radius_kpc']:.3f} | {x['instantaneous_relative_angular_rate_kms_per_kpc']:.3f} | {x['constant_rate_phase_scale_Myr']:.0f} |")
lines += ['', 'Two initial rates imply relative-phase timescales longer than the integration. A star and the bar can move at nearly the same angular rate, so their relative position changes slowly even while the star moves rapidly around the Galaxy. This explains why a duration expressed only as years is not a universal phase-coverage criterion. These initial timescales do not establish an actual resonance, orbital period, capture mechanism or the duration required for statistical convergence.', '',
'## Scope and scientific consequences','',
'These nine probes test duration and phase representation. They do not replace the 108-orbit library with a fitted population or provide equally complete libraries for both gravity models. Neither equal orbit weights nor a few observed launches specify the actual Galactic distribution of stars. The targets are the existing chemically restricted 0.5–9-kpc training sample; no claim extends this test to every Galactic population, the 20–25-kpc outer disk, or external galaxies.', '',
'A meaningful gravity comparison still needs a sufficiently sampled orbital population, the same population flexibility across models, justified ordinary-matter uncertainty, survey selection and a non-duplicated treatment of Gaia/StarHorse distance information. It must predict radial, rotational and vertical velocity distributions jointly before scoring reserved observations. The candidate field is an empirical conservative response, not a derived photon-deposition law.', '',
'The next methodological step is an orbit-population construction with a documented duration/phase-convergence criterion: determine whether its predicted spatial and velocity distributions remain stable when the integration window is lengthened or shifted. The present results do not justify selecting one common duration for all paths. A few favorable returning paths or additional cumulative samples must not stand in for that check.', '',
'The larger objective remains active. An acceptable nonexpanding photon/companion mechanism must also account for spectral redshift and event timing together, actual source/detector clocks, physical energy and momentum, capture and supported deposits, and lensing. The failed stress-energy gate in the separate reservoir construction remains unresolved; the cosmic photon supply remains deferred, not passed. No validation or final-test outcome was opened or scored.']
(HERE/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
dest=Path('C:/Users/henry/Documents/Codex/2026-09-09/cr/outputs/paired-orbit-duration');dest.mkdir(parents=True,exist_ok=True)
for name in ['report.md','coverage.png','protocol.md','ordinary-integration.json','full-integration.json','checks.json','results.json']:
    shutil.copy2(HERE/name,dest/name)
print(dest)
