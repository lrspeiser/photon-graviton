from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent
load=lambda s:json.loads((H/s).read_text())
base=load('predictions.json');ref=load('predictions-refined.json');higher=load('predictions-L64.json')
key=lambda r:(r['observable'],r['model'],r['R_kpc'],r['observed'])
assert [key(r) for r in base]==[key(r) for r in ref]==[key(r) for r in higher]
verification={}
for obs in ['rotation','vertical']:
    ix=[i for i,x in enumerate(base) if x['observable']==obs]
    verification[obs]={'max_sampling_refinement_prediction_change':max(abs(base[i]['predicted']-ref[i]['predicted']) for i in ix),
        'max_bar_L40_to_L64_prediction_change':max(abs(base[i]['predicted']-higher[i]['predicted']) for i in ix)}
    assert verification[obs]['max_sampling_refinement_prediction_change']<.01
    assert verification[obs]['max_bar_L40_to_L64_prediction_change']<1
results=load('results-refined.json')
verification['max_conservative_closed_work_kms2']=max(abs(l['work_kms2'][k]) for l in results['loops'] for k in ['baryons','spherical','potential_composition'])
assert verification['max_conservative_closed_work_kms2']<.01
assert min(abs(l['work_kms2']['direct_multiplier']) for l in results['loops'])>100
rot={m:np.array([r['predicted'] for r in ref if r['model']==m and r['observable']=='rotation']) for m in ['spherical','potential_composition']}
verification['max_two_conservative_rotation_difference_kms']=float(np.max(abs(rot['spherical']-rot['potential_composition'])))
assert verification['max_two_conservative_rotation_difference_kms']<.001
verification['interpretation']='Finite numerical comparisons only; not continuum error bounds, mass-model uncertainty or observational validation.'
(H/'verification.json').write_text(json.dumps(verification,indent=2)+'\n',encoding='utf-8')
fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
colors={'baryons':'#757575','spherical':'#0072b2','potential_composition':'#d55e00'}
labels={'baryons':'Ordinary matter','spherical':'Spherical extra well','potential_composition':'Deepen existing well shape'}
inp=load('../milky-way-capture/inputs.json')
for ax,obs in zip(axes,['rotation','vertical']):
    for model in colors:
        subset=sorted([r for r in ref if r['model']==model and r['observable']==obs],key=lambda r:r['R_kpc'])
        ax.plot([r['R_kpc'] for r in subset],[r['predicted'] for r in subset],label=labels[model],color=colors[model],ls='--' if model=='spherical' else '-',lw=2.3,zorder=3 if model=='spherical' else 2)
    subset=[r for r in ref if r['model']=='baryons' and r['observable']==obs]
    er=[x['Kz_error'] for x in inp['bovy']['rows']] if obs=='vertical' else [[x['err_minus_kms'] for x in inp['eilers']['rows']],[x['err_plus_kms'] for x in inp['eilers']['rows']]]
    ax.errorbar([r['R_kpc'] for r in subset],[r['observed'] for r in subset],yerr=er,fmt='.',color='black',alpha=.55,label='Published inference',capsize=1)
    ax.set_xlabel('Distance from Galactic center (kpc)');ax.grid(alpha=.15)
    ax.set_title('Same rotation prediction' if obs=='rotation' else 'Different pull at 1.1 kpc above disk')
    ax.set_ylabel('Circular speed (km/s)' if obs=='rotation' else r'$|K_z|/(2\pi G)$  ($M_\odot$/pc$^2$)')
axes[0].legend(fontsize=8)
fig.suptitle('Frozen empirical gravity fit: geometry remains an independent assumption',fontsize=12)
fig.savefig(H/'comparison.png',dpi=170)
print(json.dumps(verification,indent=2))
