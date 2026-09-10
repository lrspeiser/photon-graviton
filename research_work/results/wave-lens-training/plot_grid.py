from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent/'parameter-grid'
summary=json.loads((H/'summary.json').read_text());assert summary['grid_complete']
best=summary['best_successful_cell']
rows=json.loads((H/f'cell-{best}'/'predictions.json').read_text())
fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
scores=np.array([r['selection_score'] for r in summary['cells']]).reshape(3,3)
axes[0].imshow(scores,cmap='viridis_r',aspect='auto')
axes[0].set(xticks=range(3),xticklabels=['0.3','0.6','1.0'],yticks=range(3),yticklabels=['0.7','1.0','1.4'],
            xlabel='Source / stellar mass',ylabel='Particle mass (10^-24 eV/c²)',title='Training log-error score (lower is better)')
for i in range(3):
    for j in range(3):axes[0].text(j,i,f'{scores[i,j]:.4f}',ha='center',va='center',color='white' if scores[i,j]>.06 else 'black')
for ax,pred,obs,label in [(axes[1],'sigma_pred_km_s','sigma_observed_km_s','Dispersion (km/s)'),
                          (axes[2],'theta_pred_arcsec','theta_SIE_arcsec','Einstein angle (arcsec)')]:
    for model,color,name in [('baryons','#7a8794','Ordinary matter'),('stationary_wave','#c95728','Selected wave setting')]:
        sub=[r for r in rows if r['model']==model]
        ax.scatter([r[obs] for r in sub],[r[pred] for r in sub],s=20,alpha=.8,label=name,color=color)
    lo=min(min(r[obs],r[pred]) for r in rows)*.85;hi=max(max(r[obs],r[pred]) for r in rows)*1.05
    ax.plot([lo,hi],[lo,hi],'--',color='black',lw=1)
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='Measured '+label.lower(),ylabel='Predicted '+label.lower())
axes[1].legend(fontsize=8)
fig.suptitle('32 training lenses — fixed photometric masses; selected on these same data',fontsize=12)
fig.savefig(H/'training-comparison.png',dpi=170)
plt.close(fig)
