"""Plot frozen extrapolations; markers are calculated radii, lines guide the eye."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
results=json.loads((P/'outer-companion-predictions-results.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
for row in results['rows']:
    v=row['outer_predictions'];x=[p['radius_kpc']/1000 for p in v]
    axes[0].plot(x,[p['circular_speed_kms'] for p in v],marker='o',label=row['Name'],linewidth=2.5 if row['Name']=='J1621+3931' else 1)
    if row['Name']=='J1621+3931':
        axes[1].plot(x,[p['companion_kappa'] for p in v],'o-',label='Companion convergence')
        axes[1].plot(x,[p['reduced_shear'] for p in v],'s-',label='Total reduced tangential shear')
axes[0].set(xscale='log',yscale='log',xlabel='Physical radius (Mpc)',ylabel='Implied circular speed (km/s)',title='Frozen outer gravity profiles')
axes[0].legend(fontsize=8,ncol=2)
axes[1].set(xscale='log',xlabel='Projected radius (Mpc)',ylabel='Dimensionless lensing strength',title='J1621+3931: extended-reservoir signature')
axes[1].axhline(0,color='black',linewidth=.7);axes[1].legend(fontsize=8)
for ax in axes:ax.grid(alpha=.2)
fig.suptitle('Model extrapolations, not outer observations; no refitting',fontsize=12)
fig.savefig(P/'outer-companion-predictions.png',dpi=170)
