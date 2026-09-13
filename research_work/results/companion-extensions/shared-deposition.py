"""Bounded shared inverse deposition law and frozen observable replay."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parent;L=P.parent/'isotropic-galaxy-transfer'
halos=json.loads((P/'free-nfw-results.json').read_text())['rows']
inputs=json.loads((L/'capacity-reference-optics-results.json').read_text())['rows']
bins=json.loads((P.parent/'slacs-resolved-input-audit/results.json').read_text())['systems']
profiles=json.loads((P.parent/'slacs-light-profile-audit/results.json').read_text())['rows']
data=[]
def F(x):return np.log1p(x)-x/(1+x)
for h in halos:
 a=next(r for r in inputs if r['Name']==h['Name'] and r['retention_mapping']['population']=='Chabrier')
 v=next(r for r in bins if r['Name']==h['Name']);pr=next(r for r in profiles if r['Name']==h['Name'])
 Re=h['nfw_rs_kpc']/h['nfw_rs_Re'];M=a['retention_mapping']['population_mass_Msun']
 radius=(np.array(v['inner_arcsec'])+v['outer_arcsec'])/2/pr['computed_equal_area_half_light_arcsec']*Re
 data.append(dict(name=h['Name'],M=M,eta=a['retention_mapping']['eta'],Re=Re,w=4.30091727003628e-6*M/Re/300**2,radius=radius,target=h['nfw_mass_amplitude_Msun']*F(radius/h['nfw_rs_kpc'])/M))
def params(q,d):
 k,b=np.exp(q[:2]);p=q[2] if len(q)>2 else 0.;s=q[3] if len(q)>3 else 0.
 return k*d['eta']*d['M']*d['w']**p,b*d['Re']*d['w']**s
def fit(n,train):
 def residual(q):
  return np.concatenate([(np.arcsinh(params(q,d)[0]*F(d['radius']/params(q,d)[1])/d['M'])-np.arcsinh(d['target']))/np.sqrt(len(d['radius'])) for d in train])
 lo=[np.log(1e-4),np.log(.01)]+[-4]*(n-2);hi=[np.log(1e5),np.log(1000)]+[4]*(n-2)
 runs=[least_squares(residual,[np.log(k),np.log(b)]+[0.]*(n-2),bounds=(lo,hi),max_nfev=2000,ftol=1e-10,xtol=1e-10,gtol=1e-10) for k in [1,100,10000] for b in [.1,3,100]]
 good=[r for r in runs if r.success];assert good
 best=min(good,key=lambda r:sum(r.fun**2))
 return dict(q=best.x.tolist(),loss=float(np.mean(best.fun**2)),successful_starts=len(good),total_starts=len(runs),on_bound=[bool(min(abs(v-l),abs(v-h))<1e-4) for v,l,h in zip(best.x,lo,hi)])
families=[]
for n in [2,3,4]:
 full=fit(n,data);omitted=[]
 for i,d in enumerate(data):
  run=fit(n,[v for j,v in enumerate(data) if i!=j]);A,rs=params(run['q'],d)
  omitted.append(dict(name=d['name'],fit=run,A=A,rs=rs))
 families.append(dict(n=n,full=full,omitted=omitted))
 print('Fitted shared family',n,flush=True)
out=dict(families=families,observables=[])
# Reuse production setup only; no optimizer, recalibration or mutation of saved fits.
prefix=(P/'free-nfw.py').read_text().split('    seeds=[];reproduction=[]')[0]
replay=r'''
    frozen=next(v for v in halo_rows if v['Name']==name)
    d=next(v for v in shared_data if v['name']==name)
    h=frozen['h'];b0=frozen['beta0'];bi=frozen['beta_infinity'];ra=Re*frozen['orbit_transition_radius_Re'];mass=frozen['stellar_mass_Msun']
    weight=h*meanH/(1+h*meanH)
    beta=b0+(bi-b0)*r*r/(r*r+ra*ra)
    factor=np.exp(2*b0*np.log(r/model.a)+(bi-b0)*np.log((r*r+ra*ra)/(model.a*model.a+ra*ra)))
    starforce=mass/1e11*((1-weight)*model.forces[0]+weight*model.forces[1])
    starbend=mass/1e11*((1-weight)*s0+weight*sh)
    for fam in shared_families:
      for mode in ['full','omitted']:
        if mode=='full':A,rs=shared_params(fam['full']['q'],d)
        else:
          pred=next(v for v in fam['omitted'] if v['name']==name);A,rs=pred['A'],pred['rs']
        force=starforce+G*A*nfw_fraction(r/rs)/r**2
        pressure=-cumulative_trapezoid((nu*force*factor)[::-1],r[::-1],initial=0)[::-1]/factor
        speed=np.sqrt(np.trapezoid(r*r*pressure*(model.W-beta*model.T),r,axis=1)/model.den)
        err=y-speed;chi=float(err@cho_solve(full_fac,err))
        bending=starbend+A*unit_deflection(np.log(rs/Re))
        shared_output['observables'].append(dict(name=name,n=fam['n'],mode=mode,A=float(A),rs_kpc=float(rs),chi2=chi,lens_fractional_error=float(bending/required-1),predicted_kms=speed.tolist()))
    print('Replayed',name,flush=True)
'''
env=dict(__file__=str(P/'free-nfw.py'),__name__='shared_replay',halo_rows=halos,shared_data=data,shared_families=families,shared_params=params,shared_output=out)
exec(compile(prefix+replay,'frozen-shared-replay','exec'),env)
out['summary']=[]
for fam in families:
 for mode in ['full','omitted']:
  rows=[v for v in out['observables'] if v['n']==fam['n'] and v['mode']==mode]
  out['summary'].append(dict(n=fam['n'],mode=mode,motion_chi2=sum(v['chi2'] for v in rows),lens_fractional_RMS=float(np.sqrt(np.mean([v['lens_fractional_error']**2 for v in rows])))))
out['baseline_motion_chi2']=sum(v['all_motion_chi2'] for v in halos)
paths=[P/'shared-deposition-protocol.md',Path(__file__),P/'free-nfw.py',P/'free-nfw-results.json',L/'capacity-reference-optics-results.json',P.parent/'slacs-resolved-input-audit/results.json',P.parent/'slacs-light-profile-audit/results.json']
out['source_sha256']={str(p.relative_to(P.parents[2])).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(P/'shared-deposition-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(out['summary'])
