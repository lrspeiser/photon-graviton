"""Standalone scientific summary of scoped measurements, not a validation seal."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
result=json.loads((HERE/'results.json').read_text())
obs=json.loads((HERE/'observations.json').read_text())
with (HERE/'redshift-predictions.csv').open() as f:
    redshift=list(csv.DictReader(f))
with (HERE/'spectral-aging-predictions.csv').open() as f:
    timing=list(csv.DictReader(f))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'electromagnetic-audit'})
fig,axs=plt.subplots(2,2,figsize=(12.8,9.2),layout='constrained')
blue,orange,green='#245EA8','#C35D26','#24765C'
ax=axs[0,0]
for name,label,color,marker in [('SBF_group','164 exposed groups',blue,'o'),('maser_galaxy','6 exposed maser galaxies',orange,'D')]:
    rows=[r for r in redshift if r['sample']==name]
    ax.scatter([float(r['D_mpc']) for r in rows],[float(r['observed_z']) for r in rows],s=16,c=color,marker=marker,label=label,alpha=.75)
d=np.linspace(0,140,200)
ax.plot(d,np.expm1(result['alpha_per_mpc']*d),color=green,label='Fixed fractional-loss prediction')
ax.set(xlabel='Adopted distance (Mpc)',ylabel='Observed redshift',title='A. Distance trend: empirical rate, residuals remain')
ax.legend(fontsize=8,loc='upper left')

ax=axs[0,1]
z=np.array([float(r['z']) for r in timing]); age=np.array([float(r['observed_aging_rate']) for r in timing]); sigma=np.array([float(r['sigma']) for r in timing])
ax.errorbar(z,age,yerr=sigma,fmt='o',ms=3,color=blue,alpha=.75,label='35 published spectral-aging estimates')
zz=np.linspace(0,.65,100)
ax.plot(zz,1/(1+zz),color=green,label='Matched wavelength/event stretch')
ax.axhline(1,color=orange,ls='--',label='Stationary loss, unchanged source evolution')
ax.set(xlabel='Measured redshift',ylabel='Apparent aging rate',title='B. Timing: source-template assumptions retained')
ax.legend(fontsize=8,loc='lower left')

ax=axs[1,0]
lo,hi=obs['radio_methanol']['independent_centroid_groups']
nu=np.linspace(lo['frequency_GHz'],hi['frequency_GHz'],100)
A=np.log1p(lo['z'])
for p,label,color in [(-1,'Fixed absolute loss',orange),(0,'Fixed fractional loss',green),(1,'Fractional rate proportional to energy',blue)]:
    Ap=A if p==0 else np.expm1(p*A)/p
    logS=np.full_like(nu,A) if p==0 else np.log1p(p*Ap*(nu/nu[0])**p)/p
    ax.plot(nu,np.expm1(logS),label=label,color=color)
ax.errorbar([lo['frequency_GHz'],hi['frequency_GHz']],[lo['z'],hi['z']],yerr=[lo['z_error'],hi['z_error']],fmt='ko',ms=5,label='Two centroid groups; errors smaller than markers')
ax.set(xlabel='Laboratory line frequency (GHz)',ylabel='Redshift',title='C. Radio: calibrate low line, predict high line')
ax.legend(fontsize=7.5,loc='lower left')

ax=axs[1,1]
profiles=result['thermal']['profiles']
for model,label,color,marker in [('number_preserved','Photon number retained',orange,'o'),('extra_photon_removal','Extra removal restores normalization',green,'s')]:
    r=[v for v in profiles if v['model']==model]
    ax.plot(range(len(r)),[v['chi2_diagonal'] for v in r],marker=marker,color=color,label=label)
ax.set_yscale('log')
ax.set_xticks(range(5),['1','0.9924','0.99','0.5','0.4545'])
ax.set(xlabel='Assumed photon energy survival',ylabel='Diagonal residual chi-square (43 channels)',title='D. FIRAS: specified initially thermal, fixed-volume bath')
ax.legend(fontsize=8)
fig.suptitle('Photon-companion model: available-data audit',fontsize=17,weight='bold')
fig.supxlabel('Exploratory comparisons, not a completed theory. Independent sample validation and full covariance remain outstanding.',fontsize=9)
fig.savefig(HERE/'summary.png',dpi=180)
fig.savefig(HERE/'summary.svg',metadata={'Date':None})
svg=HERE/'summary.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')
print(str(HERE/'summary.png'))
