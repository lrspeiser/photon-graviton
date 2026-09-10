"""Plot the mass-normalization tradeoff; not a likelihood or confidence band."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
profile=json.loads((HERE/'profile.json').read_text())
result=json.loads((HERE/'results.json').read_text())
lam=np.array([r['lambda_mass'] for r in profile])
fits={r['calibration_objective']:r['lambda_mass'] for r in result['cases']}
fig,axes=plt.subplots(1,2,figsize=(10.5,4.4))
for ax,obs,label in zip(axes,['rotation','vertical'],['Rotation RMS error (km/s)','Vertical-force RMS (surface-equivalent units)']):
    ax.plot(lam,[r['metrics'][obs]['rmse'] for r in profile],color='#037e88',lw=2)
    for name,color in [('rotation','#8564ad'),('vertical','#c76d37')]:
        ax.axvline(fits[name],color=color,ls='--',label=f'{name} preferred mass')
    ax.axvline(1,color='#777777',ls=':',label='Original mass')
    ax.set(xlim=(.55,1.3),ylim=(0,80),xlabel='Common ordinary-matter mass multiplier',ylabel=label)
    ax.legend(fontsize=8)
fig.suptitle('One mass adjustment cannot optimize both Milky Way comparisons')
fig.text(.5,.015,'Fixed component shapes and empirical coefficients. Exposed, model-dependent force summaries; no confidence intervals or fresh holdouts.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.04,1,.94));fig.savefig(HERE/'mass-tradeoff.png',dpi=170);plt.close(fig)
