"""Reproduce the frozen constant-rate benchmark; never alter adopted data."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'analysis'
OUT.mkdir(exist_ok=True)
C = 299792.458
source = ROOT / 'redshift_paper/all_164_groups.csv'
digest = hashlib.sha256(source.read_bytes()).hexdigest()
assert digest == '8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0'
rows = list(csv.DictReader(source.open(newline='', encoding='utf-8-sig')))
d = np.array([float(r['catalog_distance_mpc']) for r in rows])
z = np.array([float(r['observed_cmb_z']) for r in rows])
split = np.array([r['split'] for r in rows])
train = split == 'train'
objective = lambda k: np.sum((C * (np.expm1(k*d[train]/C)-z[train]))**2)
fit = minimize_scalar(objective, bounds=(0,150), method='bounded', options={'xatol':1e-9})
k = min([0., fit.x, 150.], key=objective)
alpha = k/C
pred = np.expm1(alpha*d)
linear_alpha = np.dot(d[train], z[train])/np.dot(d[train], d[train])
linear = linear_alpha*d
saved = json.loads((ROOT/'research_work/results/conversion-first/results.json').read_text())
assert abs(alpha-saved['fitted']['alpha_per_mpc']) < 1e-12
assert len(rows) == len({r['group_pgc'] for r in rows}) == 164
old = {r['pgc']:r for r in csv.DictReader((ROOT/'research_work/results/conversion-first/predictions.csv').open())}
assert max(abs(pred[i]-float(old[r['pgc']]['predicted_conversion_z'])) for i,r in enumerate(rows)) < 1e-10
metrics = {}
for name in ['train','validation','test']:
    m = split == name
    residual = C*(pred[m]-z[m])
    metrics[name] = dict(n=int(m.sum()), rmse_km_s=float(np.sqrt(np.mean(residual**2))),
        mae_km_s=float(np.mean(abs(residual))), bias_km_s=float(residual.mean()),
        linear_rmse_km_s=float(np.sqrt(np.mean((C*(linear[m]-z[m]))**2))))
result = dict(source_sha256=digest, alpha_per_mpc=alpha, c_alpha_km_s_per_mpc=k,
    distance_range_mpc=[float(d.min()),float(d.max())], metrics=metrics,
    status='Exploratory calibration on previously exposed catalog; not a blind physical prediction',
    saved_predictions_reproduced=True)
(OUT/'metrics.json').write_text(json.dumps(result, indent=2)+'\n')
with (OUT/'observed-predicted-redshift.csv').open('w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['pgc','group_pgc','split','distance_mpc','observed_cmb_z','predicted_z','residual_pred_minus_obs_km_s','linear_control_z','transferred_energy_fraction'])
    for i,r in enumerate(rows):
        w.writerow([r['pgc'],r['group_pgc'],r['split'],d[i],z[i],pred[i],C*(pred[i]-z[i]),linear[i],pred[i]/(1+pred[i])])
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
colors={'train':'#226b9b','validation':'#bf7421','test':'#38866a'}
fig,axs=plt.subplots(2,1,figsize=(7.2,7.4),layout='constrained')
grid=np.linspace(0,100,400)
for name,color in colors.items():
    m=split==name
    axs[0].scatter(d[m],z[m],s=20,color=color,alpha=.8,label=f'{name} ({m.sum()})')
    axs[1].scatter(z[m],pred[m],s=20,color=color,alpha=.8)
axs[0].plot(grid,np.expm1(alpha*grid),color='#242424',label='Fitted exponential')
axs[0].plot(grid,linear_alpha*grid,'--',color='#a04266',label='Fitted linear control')
axs[0].set(xlabel='Adopted catalog distance (Mpc)',ylabel='Observed redshift z',title='A. Distance relation: 164 previously exposed groups')
axs[0].legend(fontsize=8,ncol=2)
axs[1].plot([0,.028],[0,.028],color='#555555',linestyle='--')
axs[1].set(xlabel='Observed CMB-frame redshift z',ylabel='Predicted conversion redshift z',title='B. Observed versus predicted (dashed line: equality)')
fig.savefig(OUT/'redshift-comparison.png',dpi=200)
plt.close(fig)
fig,ax=plt.subplots(figsize=(7.2,3.5),layout='constrained')
for name,color in colors.items():
    m=split==name
    ax.scatter(d[m],C*(pred[m]-z[m]),s=23,color=color,alpha=.85,label=name)
ax.axhline(0,color='#555555',lw=1)
ax.set(xlabel='Adopted catalog distance (Mpc)',ylabel='c (predicted z - observed z), km/s',title='Residuals: all 164 groups retained')
ax.legend(fontsize=8,ncol=3)
fig.savefig(OUT/'redshift-residuals.png',dpi=200)
plt.close(fig)
print(json.dumps(result,indent=2))
