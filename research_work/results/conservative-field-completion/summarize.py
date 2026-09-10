"""Compare numerical sensitivities and plot the exposed observational checks."""
from pathlib import Path
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
def load(name):return json.loads((HERE/name).read_text())
coarse=load('predictions.json');fine=load('predictions-refined.json');outer=load('predictions-refined-outer400.json')
checks={}
for obs in ['rotation','vertical']:
    select=lambda rows:np.array([r['predicted'] for r in rows if r['observable']==obs and r['model']=='conservative_completion'])
    c,f,o=map(select,[coarse,fine,outer])
    checks[obs]=dict(max_coarse_fine_prediction_change=float(np.max(abs(c-f))),
                     max_200_400_kpc_prediction_change=float(np.max(abs(o-f))))
    seen=np.array([r['observed'] for r in fine if r['observable']==obs and r['model']=='conservative_completion'])
    checks[obs]['median_prediction_observation_ratio']=float(np.median(f/seen))
assert checks['rotation']['max_coarse_fine_prediction_change']<1
assert checks['vertical']['max_coarse_fine_prediction_change']<2
assert checks['rotation']['max_200_400_kpc_prediction_change']<.01
assert checks['vertical']['max_200_400_kpc_prediction_change']<.01
analytic=load('analytic-verification.json')
mass=load('results-refined.json')['disk_mass_Msun']
checks['refined_disk_mass_relative_error']=abs(mass/analytic['disk_mass_independent_Msun']-1)
assert checks['refined_disk_mass_relative_error']<1e-3
old=json.loads((HERE.parent/'gravity-geometry-audit/predictions-refined.json').read_text())
for obs in ['rotation','vertical']:
    current=np.array([r['predicted'] for r in fine if r['model']=='ordinary' and r['observable']==obs])
    before=np.array([r['predicted'] for r in old if r['model']=='baryons' and r['observable']==obs])
    checks[obs]['max_new_old_ordinary_baseline_difference']=float(np.max(abs(current-before)))
checks['analytic']=analytic
checks['development_note']='The razor-thin Kuzmin limit initially failed the 1% force check at ell=128 (1.44%). Exact kernel integration reduced it to 1.19%; ell=512 with 2048 angular nodes reduced it to 0.247%. The MW finite-thickness comparison uses separately refined ell=64/128 response grids.'
paths=[HERE/'run.py',HERE/'verify.py',HERE.parent/'milky-way-capture/inputs.json',HERE.parent/'joint-galaxy-audit/results.json',HERE.parents[2]/'research_work/data-cache/bar-field/bar-L64.npz',HERE.parents[2]/'research_work/data-cache/bar-field/nuclei-L16.npz']
checks['input_sha256']={str(p.relative_to(HERE.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(HERE/'verification.json').write_text(json.dumps(checks,indent=2)+'\n',newline='\n')
fig,axes=plt.subplots(1,2,figsize=(11,4.5))
for ax,obs,ylabel in zip(axes,['rotation','vertical'],['Circular speed (km/s)','Vertical pull / (2 pi G), solar masses/pc²']):
    for label,color in [('ordinary','#78828b'),('conservative_completion','#037e88')]:
        rows=[r for r in fine if r['observable']==obs and r['model']==label]
        rr=np.array([r['R_kpc'] for r in rows]);ii=np.argsort(rr)
        ax.plot(rr[ii],np.array([r['predicted'] for r in rows])[ii],'.-',color=color,label=label.replace('_',' '))
    ax.scatter(rr,[r['observed'] for r in rows],s=17,color='#222222',label='Published inferred values',zorder=3)
    ax.set(xlabel='Radius from Galactic axis (kpc)',ylabel=ylabel)
    ax.legend(fontsize=8)
axes[0].set_title('Rotation remains much closer to the data')
axes[1].set_title('Pull at height 1.1 kpc is overpredicted')
fig.suptitle('Frozen empirical coefficients in a known conservative field equation')
fig.text(.5,.015,'Previously exposed, model-dependent force summaries; no new stellar holdouts, error-covariance fit, or photon-origin claim.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.04,1,.94));fig.savefig(HERE/'comparison.png',dpi=170);plt.close(fig)
print(json.dumps(checks,indent=2))
