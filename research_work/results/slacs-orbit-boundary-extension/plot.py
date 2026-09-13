from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE=Path(__file__).resolve().parent
result=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
data={r['Name']:r for r in json.loads((HERE.parent/'slacs-resolved-input-audit/results.json').read_text(encoding='utf-8'))['systems']}
names=list(dict.fromkeys(r['Name'] for r in result['rows']))
fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
colors={'baryons':'#2563a8','empirical_extra':'#d36126'}
for ax,name in zip(axes.flat,names):
 obs=data[name];left=np.array(obs['inner_arcsec']);right=np.array(obs['outer_arcsec']);mid=(left+right)/2
 ax.errorbar(mid,obs['vrms_kms'],yerr=np.sqrt(np.diag(obs['covariance_kms_squared'])),fmt='o',color='#222222',ms=4,capsize=2,zorder=3)
 for label in colors:
  row=next(r for r in result['rows'] if r['Name']==name and r['model']==label)
  ax.hlines(row['predicted_vrms_kms'],left,right,color=colors[label],lw=2)
 ax.set_title(name,fontsize=11);ax.set_xlabel('Projected radius (arcsec)');ax.set_ylabel('Annular Vrms (km/s)');ax.grid(alpha=.15)
fig.suptitle('Measured stellar motions and fitted predictions — six training lenses',fontsize=15)
fig.legend(handles=[Line2D([],[],color='#222222',marker='o',ls='',label='Measured; diagonal 1-sigma bars'),Line2D([],[],color=colors['baryons'],lw=2,label='Mass-follows-light baseline'),Line2D([],[],color=colors['empirical_extra'],lw=2,label='Fixed empirical extra force')],loc='outside lower center',ncols=3,frameon=False)
fig.savefig(HERE/'training-radial-profiles.png',dpi=160)
plt.close(fig)
