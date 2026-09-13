"""Single shared deposited-profile homology; train once, transfer without refitting."""
from pathlib import Path
import json,hashlib
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar

P=Path(__file__).resolve().parent;ROOT=P.parents[2];OLD=P.parent/'isotropic-galaxy-transfer'
paths=[OLD/'third-radiation-retention-results.json',OLD/'model-comparison-predictions.json',OLD/'model-comparison-results.json',ROOT/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt',P/'homology-retention-protocol.md']
fit=json.loads(paths[0].read_text())['models']['attenuated']; saved=json.loads(paths[1].read_text()); controls=json.loads(paths[2].read_text())
cat={}
for line in paths[3].read_text().splitlines():
 f=line.split()
 if len(f)==19:
  try:cat[f[0]]=(float(f[11]),float(f[7])*1e9)
  except ValueError:pass
data=[]
for x in saved:
 if x['model']!='companion_third':continue
 b=next(v for v in saved if v['model']=='baryons' and v['galaxy']==x['galaxy'])
 assert b['R_kpc']==x['R_kpc'] and b['observed_kms']==x['observed_kms']
 data.append(dict(name=x['galaxy'],split=x['split'],r=np.array(x['R_kpc']),y=np.array(x['observed_kms']),vb2=np.array(b['predicted_kms'])**2,saved=np.array(x['predicted_kms']),rd=cat[x['galaxy']][0],L=cat[x['galaxy']][1]))
G=4.30091727003628e-6
def profiles(n,order):
 mu,w=leggauss(order);answer=[]
 for d in data:
  a=d['rd']*fit['scale_to_disk'];r=np.r_[0.,np.geomspace(1e-7*a,max(2.01*d['r'].max(),30*a),n)];x=r/a
  t=x[:,None]*mu;B2=1+x[:,None]**2*(1-mu*mu);B=np.sqrt(B2)
  tau=fit['k0_per_kpc']*a*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3))
  J=np.exp(-np.maximum(tau,0))@w/2;X=(d['L']/1e9)/d['rd']**2;eta=X**(1/3)/(1+X**(1/3))
  rho=2*fit['C_Msun_kpc3']*eta*J/(1+x*x)**2
  mass=4*np.pi*cumulative_trapezoid(r*r*rho,r,initial=0)
  assert np.all(np.diff(mass)>=0)
  answer.append(PchipInterpolator(r,mass,extrapolate=False))
 return answer
def velocities(s,masses):return [np.sqrt(d['vb2']+G*s*m(d['r']/s)/d['r']) for d,m in zip(data,masses)]
def loss(s,masses):
 v=velocities(s,masses)
 return float(np.mean([np.mean(np.log10(p/d['y'])**2) for d,p in zip(data,v) if d['split']=='train']))
def score(v):
 return {split:dict(n=sum(d['split']==split for d in data),RMSE_kms=float(np.sqrt(np.mean([np.mean((p-d['y'])**2) for d,p in zip(data,v) if d['split']==split]))),log_RMS=float(np.sqrt(np.mean([np.mean(np.log10(p/d['y'])**2) for d,p in zip(data,v) if d['split']==split])))) for split in ['train','validation','test']}
coarse=profiles(2048,96);grid=np.linspace(.5,1,101);ys=[loss(s,coarse) for s in grid];candidates=[(.5,ys[0]),(1.,ys[-1])]
for lo,hi in [(grid[0],grid[1]),(grid[-2],grid[-1])]:
 f=minimize_scalar(lambda s:loss(s,coarse),bounds=(lo,hi),method='bounded',options={'xatol':1e-10});assert f.success;candidates.append((float(f.x),float(f.fun)))
for i in range(1,len(grid)-1):
 if ys[i]<=ys[i-1] and ys[i]<=ys[i+1]:
  f=minimize_scalar(lambda s:loss(s,coarse),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-9});assert f.success;candidates.append((float(f.x),float(f.fun)))
s,objective=min(candidates,key=lambda v:v[1]);fine=profiles(4096,192)
vr=velocities(1,fine);vs=velocities(s,fine);vc=velocities(s,coarse)
delta=max(float(np.max(abs(p-d['saved']))) for d,p in zip(data,vr));assert delta<.05,delta
resolution=max(float(np.max(abs(a-b))) for a,b in zip(vs,vc));assert resolution<.05,resolution
rows=[];radial=[];improved={}
for split in ['train','validation','test']:
 dd=[(d,a,b) for d,a,b in zip(data,vr,vs) if d['split']==split]
 improved[split]=int(sum(np.mean(np.log10(b/d['y'])**2)<np.mean(np.log10(a/d['y'])**2) for d,a,b in dd))
for d,a,b in zip(data,vr,vs):rows.append(dict(galaxy=d['name'],split=d['split'],R_kpc=d['r'].tolist(),observed_kms=d['y'].tolist(),reference_kms=a.tolist(),modified_kms=b.tolist()))
for name,values in [('reference',vr),('modified',vs)]:
 for split in ['train','validation','test','all_exposed']:
  for k in range(3):
   errors=[]
   for d,v in zip(data,values):
    mask=np.digitize(d['r']/d['rd'],[1,3])==k
    if mask.any() and (split=='all_exposed' or d['split']==split):errors.append(v[mask]-d['y'][mask])
   radial.append(dict(model=name,split=split,bin=k,galaxies=len(errors),bias_kms=float(np.mean([e.mean() for e in errors])),RMSE_kms=float(np.sqrt(np.mean([np.mean(e*e) for e in errors])))))
out=dict(scope='One extra shared parameter trained on exposed training galaxies; unchanged one-third exponent; no lensing refit or physical formation claim.',s=s,stored_fraction=s,unretained_fraction=1-s,boundary=bool(s in [.5,1.]),fit_objective=objective,candidates=candidates,scan=[dict(s=float(a),loss=b) for a,b in zip(grid,ys)],scores=dict(reference=score(vr),modified=score(vs),coarse_modified=score(vc)),improved_galaxies_log_loss=improved,radial=radial,reference_max_difference_from_archived_kms=delta,refinement_max_prediction_change_kms=resolution,controls=controls['summary'],source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},rows=rows)
out['refined_scan']=[dict(s=float(a),loss=loss(a,fine)) for a in grid]
out['refined_grid_best_s']=min(out['refined_scan'],key=lambda x:x['loss'])['s']
out['fixed_diagnostics']=[dict(s=a,scores=score(velocities(a,fine))) for a in [.95,.9,.8,.5]]
for diagnostic in out['fixed_diagnostics']:
 values=velocities(diagnostic['s'],fine);bs=[];up=0;nn=0
 for d,v,b in zip(data,values,vr):
  mask=d['r']/d['rd']>=3
  if mask.any():bs.append(float(np.mean(v[mask]-d['y'][mask])));up+=int(np.sum(v[mask]>b[mask]));nn+=int(mask.sum())
 diagnostic['outer_bias_kms']=float(np.mean(bs));diagnostic['outer_points_with_increased_speed']=up;diagnostic['outer_points']=nn
ff=minimize_scalar(lambda a:loss(a,fine),bounds=(.995,1.),method='bounded',options={'xatol':1e-10});assert ff.success
out['fine_endpoint_polish_s']=float(ff.x)
out['training_log_loss_improvement_fine_at_frozen_s']=loss(1,fine)-loss(s,fine)
if s in [.5,1.]:assert out['refined_grid_best_s']==s
(P/'homology-retention-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print('s',s,'reference reproduction',delta,'refinement',resolution);print(json.dumps(out['scores'],indent=2));print('improved',improved)
print([x for x in radial if x['split']=='all_exposed'])
