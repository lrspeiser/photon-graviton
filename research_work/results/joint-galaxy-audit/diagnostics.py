"""Post-result robustness checks. No candidate is selected on these scores."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
H=Path(__file__).resolve().parent
load=lambda f:json.loads((H/f).read_text())
save=lambda f,v:(H/f).write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n')
rows=load('milky-way-predictions.json'); out={}; G=4.30091727003628e-6
for variant in ['I','II']:
    rr=[x for x in rows if x['baryons']==variant and x['model']=='baryons' and x['observable']=='vc_kms']
    R=np.array([x['R_kpc'] for x in rr]); y=np.array([x['observed'] for x in rr]); vb=np.array([x['predicted'] for x in rr]); cal=R<=15
    fit=minimize_scalar(lambda k:np.mean((np.sqrt(vb[cal]**2+k*R[cal]**2)-y[cal])**2),bounds=(0,10000),method='bounded')
    yp=np.sqrt(vb*vb+fit.x*R*R)
    out[variant]=dict(k_kms2_per_kpc2=float(fit.x),M100_Msun=float(fit.x*100**3/G),inner_RMSE_kms=float(np.sqrt(np.mean((yp[cal]-y[cal])**2))),outer_RMSE_kms=float(np.sqrt(np.mean((yp[~cal]-y[~cal])**2))))
save('uniform-normalization-sensitivity.json',out)

# Paired bootstrap of reserved GALAXIES (not points); parameters remain fixed.
sr=load('galaxy-rotation-predictions.json'); names=sorted({x['galaxy'] for x in sr if x['split']=='test'})
err={}
for model in ['baryons','power']:
    err[model]=np.array([np.mean([(x['predicted_kms']-x['observed_kms'])**2 for x in sr if x['split']=='test' and x['galaxy']==n and x['model']==model]) for n in names])
rng=np.random.default_rng(9092026); indices=rng.integers(0,len(names),size=(5000,len(names)))
delta=np.sqrt(err['power'][indices].mean(axis=1))-np.sqrt(err['baryons'][indices].mean(axis=1))
save('robustness.json',{'status':'Descriptive paired galaxy bootstrap, fixed fitted parameters; not a full nuisance/systematic uncertainty or blinded validation',
    'n_galaxies':len(names),'draws':5000,'power_minus_baryons_RMSE_95_kms':np.quantile(delta,[.025,.975]).tolist(),
    'power_better_galaxies':int(np.sum(err['power']<err['baryons']))})

profiles=load('capture-profiles.json'); metrics=[]
for profile in profiles:
    R=np.array(profile['R_kpc']); rho=np.array(profile['rho_Msun_kpc3'])
    metrics.append(dict(baryons=profile['baryons'],kind=profile['kind'],rho_25_over_5=float(np.interp(25,R,rho)/np.interp(5,R,rho)),
                        shell_mass_per_kpc_25_over_5=float(25**2*np.interp(25,R,rho)/(5**2*np.interp(5,R,rho)))))
save('capture-shape-diagnostics.json',metrics)
print(json.dumps({'uniform_sensitivity':out,'capture_shapes':metrics},indent=2))
