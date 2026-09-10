"""Show the preregistered predictions and newly evaluated Cepheid table."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
obs=json.loads((HERE/'observations.json').read_text())['rows']
pred=json.loads((HERE/'predictions.json').read_text())
fig,(ax,res)=plt.subplots(2,1,figsize=(8.5,6.5),sharex=True,gridspec_kw={'height_ratios':[2,1]})
for name,label,color in [('ordinary','Ordinary matter','#7b848c'),('conservative_completion','Primary frozen field','#037e88'),('rotation_mass_adjusted_completion','Earlier rotation-calibrated mass','#9464a1')]:
    rows=[r for r in pred if r['model']==name]
    R=[r['R_kpc'] for r in rows]
    ax.plot(R,[r['predicted_kms'] for r in rows],'.-',color=color,label=label)
    res.plot(R,[r['residual_kms'] for r in rows],'.-',color=color)
ax.errorbar([r['R_kpc'] for r in obs],[r['vc_kms'] for r in obs],yerr=[r['error_kms'] for r in obs],fmt='o',color='#222222',ms=4,capsize=3,label='Cepheid table; bootstrap errors')
ax.set(ylabel='Circular speed (km/s)',title='Frozen models checked against the 2026 Cepheid rotation table')
ax.legend(fontsize=9);res.axhline(0,color='black',lw=.8)
res.set(xlabel='Radius from Galactic axis (kpc)',ylabel='Prediction − table (km/s)')
fig.text(.5,.015,'Protocol committed before table acquisition. Shared Gaia systematics; error bars are not the complete uncertainty.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.04,1,1));fig.savefig(HERE/'comparison.png',dpi=170);plt.close(fig)
