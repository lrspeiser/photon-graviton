"""Plot frozen J1621 branch contrasts; no observed outer data."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
data=json.loads((p/'reservoir-branch-outer-results.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
for row in data['rows']:
 scale=row['capture_scale_Re']
 if scale not in [.1,.3,10,30,100]:continue
 v=row['outer_predictions'];x=[r['radius_kpc'] for r in v]
 label='Stellar-only control' if scale==.1 else f'ac/Re = {scale:g}'
 axes[0].plot(x,[r['circular_speed_kms'] for r in v],'.-',label=label)
 axes[1].plot(x,[r['reduced_shear'] for r in v],'.-',label=label)
for ax in axes:
 ax.set_xscale('log');ax.set_xlabel('Projected / spherical radius (kpc)');ax.grid(alpha=.2)
axes[0].set_yscale('log');axes[0].set_ylabel('Implied circular speed (km/s)');axes[0].legend(fontsize=8)
axes[1].set_ylabel('Reduced tangential shear');axes[1].set_ylim(-.025,.05);axes[1].axhline(0,color='grey',lw=.7)
axes[1].set_title('Outer differences (inner values exceed plotted range)',fontsize=9)
fig.suptitle('J1621+3931: frozen branches, not observed outer data\nMarkers are calculated radii; lines only connect them',fontsize=11)
fig.savefig(p/'reservoir-branch-outer.png',dpi=180)
