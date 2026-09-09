from pathlib import Path
import json,numpy as np
from scipy.optimize import minimize
# Execute immutable loader/calculation definitions but omit expensive first-run fit loop.
src=(Path(__file__).parent/'run.py').read_text();head=src.split('results={};predrows=[]')[0];exec(compile(head,'run.py','exec'))
fold={n:int(hashlib.sha256(n.encode()).hexdigest()[:8],16)%5 for n in matched}
# Controls distinguish source information from spatial template. Boundary-driven diffusion is normalized to total energy budget.
def ext_mass_frac(r,rh,a):
 x=np.minimum(r,rh)/rh
 if a>=10:return x**3 # uniform limit, fractional error <0.1%
 y=x/a;Y=1/a
 return (y*np.cosh(y)-np.sinh(y))/(Y*np.cosh(Y)-np.sinh(Y))
def pred_control(m,p,d,uc0):
 l,e=templates(d);rh=10*d['rd'];r=d['r']
 if 'constant' in m:e=e*uc0/d['uc']
 if 'boundary' in m:
  a=float(m.split('_')[-1]);x=np.minimum(r/rh,1);e=e*ext_mass_frac(r,rh,a)/(x**3)
 return np.sqrt(d['vb2']+(10**p[0]*l+10**p[1]*e if m.startswith('mix') else 10**p[0]*e))
controls={};rows=[]
for m in ['external_constant','mix_constant','external_boundary_0.1','external_boundary_0.3','external_boundary_1.0']:
 cv={};params=[]
 for k in range(5):
  tr=[n for n in matched if fold[n]!=k];te=[n for n in matched if fold[n]==k];uc0=np.median([data[n]['uc'] for n in tr]);st=[9.,9.] if m.startswith('mix') else [9.]
  o=minimize(lambda p:np.mean([np.mean(np.log10(pred_control(m,p,data[n],uc0)/data[n]['v'])**2) for n in tr]),st,bounds=[(0,15)]*len(st),method='L-BFGS-B',options={'ftol':1e-13,'gtol':1e-9,'maxiter':1000});params.append(o.x.tolist())
  for n in te:cv[n]=pred_control(m,o.x,data[n],uc0)
 controls[m]=dict(cv_scores=scores(lambda n:cv[n],matched),parameters_by_fold=params)
 for n in matched:
  for r,v,pv in zip(data[n]['r'],data[n]['v'],cv[n]):rows.append(dict(model=m,galaxy=n,r_kpc=r,observed_kms=v,predicted_kms=pv))
# Consistent pointwise depth proxy at each measured radius, for disks and clusters.
a0=3085.677581491367
cr=list(csv.DictReader((P.parent/'companion_wave_test/cluster_comparison.csv').open()))
def pred_a(p,d,kind):
 gb=d['vb2']/d['r'];dep=d['vdepth'] if kind=='global' else np.sqrt(2*d['vb2']);gx=10**p[0]*a0*(gb/a0)**p[1]*(dep/200)**p[2];return np.sqrt(d['vb2']+d['r']*gx)
def fit_a(ns,kind):
 o=minimize(lambda p:np.mean([np.mean(np.log10(pred_a(p,data[n],kind)/data[n]['v'])**2) for n in ns]),[-.22,.35,.3],bounds=[(-3,3),(0,1),(-3,3)],method='L-BFGS-B',options={'ftol':1e-12,'gtol':1e-8,'maxiter':600});return o.x
def cluster_ratios(p,fs=.02):
 out=[]
 for c in cr:
  r=float(c['R500_Mpc'])*1000;mt=float(c['M500_1e14Msun'])*1e14;mb=mt*(float(c['fgas500'])+fs);gb=G*mb/r**2;dep=np.sqrt(2*G*mb/r);gx=10**p[0]*a0*(gb/a0)**p[1]*(dep/200)**p[2];out.append((mt-mb)/(gx*r*r/G))
 return out
response={};rng=np.random.default_rng(2047)
for kind in ['global','pointwise']:
 p=fit_a(split['train'],kind);ratios=cluster_ratios(p);boots=[]
 for i in range(120):
  ns=rng.choice(split['train'],len(split['train']),replace=True).tolist();pp=fit_a(ns,kind);boots.append([*pp,np.median(cluster_ratios(pp))])
 response[kind]=dict(parameters=p.tolist(),scores={s:scores(lambda n:pred_a(p,data[n],kind),ns) for s,ns in split.items()},cluster_ratios=ratios,cluster_median=float(np.median(ratios)),bootstrap_parameter_and_cluster_median_95=np.quantile(boots,[.025,.975],axis=0).tolist(),cluster_fstar_sensitivity={str(fs):float(np.median(cluster_ratios(p,fs))) for fs in [.01,.02,.03]})
 print(kind,response[kind],flush=True)
# SED systematic energy sensitivity and distance uncertainty diagnostic, not full covariance.
# Recompute illumination with all line-of-sight distances +/-20% independently via MC. Same positions held.
keys=list(dust);xyz=np.array([dust[k]['xyz'] for k in keys]);ll=np.array([dust[k]['lum'] for k in keys])*LSUN;idx={k:i for i,k in enumerate(keys)};irr=[]
for rep in range(200):
 xx=xyz*np.exp(rng.normal(0,.2,len(keys)))[:,None];vals=[]
 for n in matched:
  j=idx[canon(n)];sep=np.linalg.norm(xx-xx[j],axis=1);sep[j]=np.inf;sep=np.maximum(sep,.03);xxk=KAPPA*sep
  with np.errstate(invalid='ignore'):flux=ll*xxk*np.exp(-2*xxk)/(4*np.pi*(sep*1000*KPC)**2)
  flux[j]=0.;vals.append(np.sum(flux)/C/data[n]['uc'])
 irr.append(vals)
irr=np.array(irr);sens=dict(assumption='independent 20% log-distance perturbations; illustrative, not published distance errors',median_across_galaxies_of_each_individual_95_interval=np.median(np.quantile(irr,[.025,.975],axis=0),axis=1).tolist())
res=dict(controls=controls,response=response,illumination_distance_sensitivity=sens)
(P/'followup_results.json').write_text(json.dumps(res,indent=2))
with (P/'control_predictions.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(res,indent=2))
