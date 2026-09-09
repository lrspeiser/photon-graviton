"""Publish saved exploratory results without fitting additional parameters."""
from pathlib import Path
import csv,json,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
source=ROOT/'research_work/results/redshift-priority'
out=HERE/'analysis'
for src,dst in [('predictions.csv','sky-tile-predictions.csv'),('maser-comparison.csv','maser-comparison.csv'),('coarse-sky-predictions.csv','coarse-sky-predictions.csv'),('observed-time-stretch-targets.csv','observed-time-stretch-targets.csv'),('observation-factor-diagnostics.csv','observation-factor-diagnostics.csv')]:
    shutil.copyfile(source/src,out/dst)
rows=list(csv.DictReader((source/'maser-comparison.csv').open()))
d=np.array([float(r['distance_mpc']) for r in rows])
z=np.array([float(r['observed_cmb_z']) for r in rows])
p=np.array([float(r['exponential_z']) for r in rows])
assert len(rows)==6
assert np.allclose(p,np.expm1(.0002488993286382367*d),atol=1e-14,rtol=0)
fig,axes=plt.subplots(2,1,figsize=(7.2,6.5),layout='constrained')
x=np.linspace(0,140,250)
axes[0].plot(x,np.expm1(.0002488993286382367*x),color='#226b9b',label='Exponential: rate fixed from SBF sample')
axes[0].plot(x,75.17651383428588*x/299792.458,'--',color='#a04266',label='Linear: rate fixed from SBF sample')
axes[0].scatter(d,z,color='#222222',s=27,label='Maser observations')
axes[0].set(xlabel='Adopted geometric distance (Mpc)',ylabel='CMB-frame redshift z')
axes[0].legend(fontsize=8)
axes[1].barh([r['name'] for r in rows],299792.458*(p-z),color='#226b9b')
axes[1].axvline(0,color='#555555',lw=1)
axes[1].set(xlabel='c (predicted z - observed z), km/s')
fig.savefig(out/'maser-comparison.png',dpi=200)
plt.close(fig)
cv=json.loads((source/'results.json').read_text())
check={name:float(np.sqrt(np.mean([(299792.458*(float(r[name+'_oof_z'])-float(r['observed_z'])))**2 for r in csv.DictReader((source/'predictions.csv').open())]))) for name in ['linear','constant','smooth']}
assert all(abs(check[n]-cv['out_of_fold'][n]['rms'])<1e-9 for n in check)
(out/'extended-checks.json').write_text(json.dumps(dict(maser_rows=6,sky_tile_rows=164,cv_rms_recomputed=check,no_new_fitting=True),indent=2)+'\n',newline='\n')

coarse=json.loads((source/'coarse-sky-results.json').read_text())
coarse_rows=list(csv.DictReader((source/'coarse-sky-predictions.csv').open()))
assert len(coarse_rows)==164
coarse_check={name:float(np.sqrt(np.mean([(299792.458*(float(r[name+'_oof_z'])-float(r['observed_z'])))**2 for r in coarse_rows]))) for name in ['linear','constant','smooth']}
assert all(abs(coarse_check[n]-coarse['out_of_fold'][n]['rms'])<1e-9 for n in coarse_check)
fig,ax=plt.subplots(figsize=(7.2,3.5),layout='constrained')
x=np.arange(3);names=['linear','constant','smooth']
for offset,values,color,label in [(-.18,check,'#226b9b','65 sky tiles'),(.18,coarse_check,'#a04266','8 larger sky regions')]:
    bars=ax.bar(x+offset,[values[n] for n in names],width=.34,color=color,label=label)
    ax.bar_label(bars,fmt='%.1f',padding=3,fontsize=8)
ax.set(xticks=x,xticklabels=['Linear control','Constant rate','Smooth rate'],ylabel='Out-of-fold RMS, c times redshift residual (km/s)',ylim=(0,540))
ax.legend(loc='upper left',fontsize=8)
fig.savefig(out/'grouped-validation.png',dpi=200);plt.close(fig)
clock=json.loads((ROOT/'research_work/results/universal-clock-coupling/results.json').read_text())
assert len(clock['records'])==10
assert max(abs(r['wavelength_stretch']-r['instantaneous_event_stretch']) for r in clock['records'])<2e-9
(out/'revision-checks.json').write_text(json.dumps(dict(coarse_rows=164,coarse_rms_recomputed=coarse_check,clock_examples=10,clock_timing_identity_checked=True,new_candidate_outcomes_scored=False,no_new_fitting=True),indent=2)+'\n',newline='\n')
