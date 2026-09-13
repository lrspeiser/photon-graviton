"""Self-consistent deposited-field feedback in the local exchange closure."""
from pathlib import Path
import json,hashlib
import numpy as np
from numba import njit
P=Path(__file__).resolve().parent;G=4.30091e-6
files=['stream-capture-results.json','halo-deposition-map-results.json','third-retention-optics-results.json','reservoir-exchange-results.json']
streams=json.loads((P/files[0]).read_text())['rows'];halos=json.loads((P/files[1]).read_text())['rows'];stars=json.loads((P/files[2]).read_text())['rows'];prior=json.loads((P/files[3]).read_text())['results']
x=np.linspace(0,1,401);names=['binding','released_binding','occupancy_penalty']
params=np.array([(s,v) for s in np.geomspace(.03,1,25) for v in np.geomspace(10,3000,25)])

def setup(row,branch,n):
 h=next(v for v in halos if all(v[k]==row[k] for k in ['Name','geometry','population']))
 st=next(v for v in stars if v['Name']==row['Name'] and v['model']=='attenuated_'+row['population'])
 M=next(v['original_deposit_enclosed_Msun'] for v in h['samples'] if v['aperture']=='5Re')
 R=row['R_kpc'];edge=np.linspace(0,1,n+1);r=(edge[:-1]+edge[1:])*R/2
 dm=np.diff(np.interp(edge,row['r_over_R'],row['case_cdfs'][branch]));a=R/5/1.8153;Ms=st['mass_Msun']
 vol=4*np.pi/3*np.diff((edge*R)**3);rho_star=Ms*a/(2*np.pi*r*(r+a)**3)
 travel=dm*(1+(r/row['a_kpc'])**2)**2;travel/=travel.sum()
 return dm,r,M,Ms,a,vol*rho_star,travel

@njit
def potential_at(q,r,m):
 pref=np.cumsum(m);inv=np.cumsum(m/r);total=inv[-1];ans=np.empty(len(q))
 for i in range(len(q)):
  j=np.searchsorted(r,q[i],side='right')-1
  ans[i]=total if j<0 else pref[j]/q[i]+total-inv[j]
 return ans

@njit
def profile(dm,f,s):
 n=len(dm);edge=np.linspace(0,1,n+1)
 C=np.empty(n+1);C[0]=0;C[1:]=np.cumsum(dm*f)
 E=np.empty(n+1);E[0]=0;E[1:]=np.cumsum(dm*(1-f))
 return E+np.interp(edge/s,edge,C)

@njit
def fields(dm,r,M,Ms,a,travel,f,s,t):
 cm=dm*f;em=dm*(1-f);dest=s*r
 B=G*Ms/(r+a)+G*M*(potential_at(r,r,em)+potential_at(r,dest,cm)+t*potential_at(r,r,travel))
 Bs=G*Ms/(dest+a)+G*M*(potential_at(dest,r,em)+potential_at(dest,dest,cm)+t*potential_at(dest,r,travel))
 return B,Bs

@njit
def solve(dm,r,M,Ms,a,starcell,travel,s,v,model,t=0.,seed=0.,feedback=True):
 f=np.full(len(dm),seed);res=1.;converged=False
 for step in range(400):
  if feedback:B,Bs=fields(dm,r,M,Ms,a,travel,f,s,t)
  else:B=G*Ms/(r+a);Bs=G*Ms/(s*r+a)
  strength=np.maximum(Bs-B,0.) if model==1 else B
  crowd=(np.diff(profile(dm,f,s))+t*travel)*M/starcell if feedback else dm*M/starcell
  resistance=v*v*(1+crowd) if model==2 else np.full(len(dm),v*v)
  desired=strength/(strength+resistance);res=np.max(np.abs(desired-f))
  if res<1e-8:converged=True;break
  f=.65*f+.35*desired
 return f,profile(dm,f,s),converged,res,step+1

@njit
def trial_table(dm,r,M,Ms,a,starcell,travel,parameters,model):
 curves=np.empty((len(parameters),len(dm)+1));valid=np.empty(len(parameters),dtype=np.bool_)
 for i in range(len(parameters)):
  f,y,ok,res,n=solve(dm,r,M,Ms,a,starcell,travel,parameters[i,0],parameters[i,1],model)
  curves[i]=y;valid[i]=ok
 return curves,valid

def specific_grav_energy(inp,f,s):
 dm,r,M,Ms,a,starcell,travel=inp;cm=dm*f;em=dm*(1-f);dest=s*r
 ustar=-np.sum(em*G*Ms/(r+a)+cm*G*Ms/(dest+a))
 selfr=G*M*(potential_at(r,r,em)+potential_at(r,dest,cm));selfs=G*M*(potential_at(dest,r,em)+potential_at(dest,dest,cm))
 return float(ustar-.5*np.sum(em*selfr+cm*selfs))

# Independent direct spherical-shell potential check.
rtest=np.array([.2,1.,3.]);mtest=np.array([.1,.6,.3]);qtest=np.array([.1,.7,2.,10.])
assert np.max(np.abs(potential_at(qtest,rtest,mtest)-np.array([np.sum(mtest/np.maximum(q,rtest)) for q in qtest])))<1e-12
frozen=[];refits=[];maxzero=0.;maxrefine=0.;maxseed=0.
for old in prior:
 rows=[v for v in streams if v['geometry']==old['geometry'] and v['population']==old['population']];branch=-1 if old['source']=='distant' else 0;model=names.index(old['model'])
 out=[]
 for row in rows:
  target=np.interp(x,row['r_over_R'],row['target_cdf']);oi=next(v for v in old['rows'] if v['Name']==row['Name']);inp=setup(row,branch,512)
  zf,zy,zok,zres,zn=solve(*inp,*old['shared_parameters'],model,0.,0.,False)
  zero=float(np.max(np.abs(np.interp(x,np.linspace(0,1,513),zy)-oi['shared_cdf'])));maxzero=max(maxzero,zero);assert zero<.002
  versions=[]
  for mode,p in [('shared',old['shared_parameters']),('loo',oi['loo_parameters'])]:
   for ratio in [0.,.01,.1]:
    f,y,ok,res,n=solve(*inp,*p,model,ratio)
    f2,y2,ok2,res2,n2=solve(*inp,*p,model,ratio,1.)
    gap=float(np.max(np.abs(y-y2)));maxseed=max(maxseed,gap)
    yy=np.interp(x,np.linspace(0,1,513),y);total=yy+ratio*np.interp(x,np.linspace(0,1,513),np.r_[0.,np.cumsum(inp[-1])])
    assert abs(y[-1]-1)<1e-10 and np.min(np.diff(y))>=-1e-12 and np.min(f)>=0 and np.max(f)<=1
    work=specific_grav_energy(inp,np.zeros(512),p[0])-specific_grav_energy(inp,f,p[0]) if ratio==0 else None
    versions.append(dict(mode=mode,traveling_to_stored_ratio=ratio,converged=bool(ok),residual=float(res),iterations=int(n),alternate_seed_converged=bool(ok2),seed_max_cdf_difference=gap,compact_fraction=float(np.sum(inp[0]*f)),deposited_rms=float(np.sqrt(np.mean((yy-target)**2))),including_prescribed_travel_rms=float(np.sqrt(np.mean((total-target)**2))),specific_binding_energy_released_kms2=work,deposited_cdf=yy.tolist()))
  out.append(dict(Name=row['Name'],zero_feedback_max_cdf_error=zero,versions=versions))
 frozen.append(dict(geometry=old['geometry'],population=old['population'],source=old['source'],model=old['model'],previous_loo_rms=old['loo_rms'],rows=out))
 print('frozen',old['geometry'],old['population'],old['source'],old['model'],flush=True)

for geometry in ['companion_regular','standard_flat_FLRW']:
 rows=[v for v in streams if v['geometry']==geometry and v['population']=='Chabrier']
 for model,name in enumerate(names):
  curves=[];errors=[];invalid=[]
  for row in rows:
   yy,ok=trial_table(*setup(row,-1,256),params,model)
   y=np.array([np.interp(x,np.linspace(0,1,257),v) for v in yy]);truth=np.interp(x,row['r_over_R'],row['target_cdf'])
   curves.append(y);errors.append(np.where(ok,np.mean((y-truth)**2,axis=1),np.inf));invalid.append(int(np.sum(~ok)))
  errors=np.array(errors);best=int(np.argmin(errors.mean(axis=0)));per=[];loo=[]
  assert np.isfinite(errors[:,best]).all()
  for i,row in enumerate(rows):
   ix=int(np.argmin(np.delete(errors,i,axis=0).mean(axis=0)));assert np.isfinite(errors[i,ix]),'Failed omitted transfer';loo.append(errors[i,ix])
   inp=setup(row,-1,512);checks={}
   for mode,j in [('shared',best),('loo',ix)]:
    f,y,ok,res,n=solve(*inp,*params[j],model);f2,y2,ok2,res2,n2=solve(*inp,*params[j],model,0.,1.)
    yy=np.interp(x,np.linspace(0,1,513),y);drift=float(np.max(np.abs(yy-curves[i][j])));gap=float(np.max(np.abs(y-y2)));maxrefine=max(maxrefine,drift);maxseed=max(maxseed,gap)
    checks[mode]=dict(converged=bool(ok),residual=float(res),iterations=int(n),compact_fraction=float(np.sum(inp[0]*f)),refinement_max_cdf_change=drift,alternate_seed_converged=bool(ok2),seed_max_cdf_difference=gap,refined_cdf=yy.tolist(),specific_binding_energy_released_kms2=specific_grav_energy(inp,np.zeros(512),params[j,0])-specific_grav_energy(inp,f,params[j,0]))
   per.append(dict(Name=row['Name'],loo_parameters=params[ix].tolist(),loo_boundary=bool(ix//25 in [0,24] or ix%25 in [0,24]),checks=checks,target_cdf=np.interp(x,row['r_over_R'],row['target_cdf']).tolist()))
  old=next(v for v in prior if v['geometry']==geometry and v['population']=='Chabrier' and v['source']=='distant' and v['model']==name)
  result=dict(geometry=geometry,population='Chabrier',source='distant',model=name,shared_parameters=params[best].tolist(),shared_boundary=bool(best//25 in [0,24] or best%25 in [0,24]),invalid_grid_counts=invalid,shared_rms=float(np.sqrt(errors[:,best].mean())),loo_rms=float(np.sqrt(np.mean(loo))),old_no_feedback_loo_rms=old['loo_rms'],rows=per)
  refits.append(result);print('refit',geometry,name,'old/shared/LOO',*[round(100*result[k],3) for k in ['old_no_feedback_loo_rms','shared_rms','loo_rms']],flush=True)

out=dict(scope='Coupled deposited Newtonian potential and local exchange fixed point, not full time-dependent relativistic evolution',input_sha256={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in files},zero_feedback_max_cdf_error=maxzero,refit_max_refinement_cdf_change=maxrefine,max_seed_cdf_difference=maxseed,grid=x.tolist(),frozen=frozen,refits=refits)
(P/'feedback-results.json').write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8',newline='\n')
lines=['# Feedback from companion gravity into reservoir exchange','',
'This test updates the exchange-driving Newtonian potential using ordinary stars AND the current compact and extended companion deposits. It also updates the occupancy penalty from the current density. The one-third capture law, original finite inventory and baryonic gravity prescription are unchanged. This is a coupled snapshot closure, not a simulated source history or complete relativistic calculation.','',
'## What is coupled','',
'Known spherical potential: Phi(r)=-G[M(<r)/r+integral_r^infinity dM(rprime)/rprime]. Stellar Phi uses the original Hernquist proxy. The prior proposed rate fractions are recomputed in the combined field; solve f=F[Phi(f),rho(f)] with damped iteration. Compact packets move to s times their ORIGINAL radius, not repeatedly inward on each numerical pass. The fixed-point iteration has no assigned physical duration.','',
'Both populations conserve their combined mass-equivalent inventory. The baryonic density is held fixed; its orbits do not respond here. The rate laws and contraction distance remain phenomenological, and final support has not been solved. Thus this is more complete feedback than the initial-field test but still not everything responding dynamically to everything else.','',
'## Fixed parameters first, then shared refits','',
'Frozen tests preserve previous shared and omitted-target parameters for every geometry/population/source branch. Primary refits use Chabrier/distant inputs, separately for both geometries, on the declared 25x25 grid. Nonconvergent trial closures are excluded and their counts recorded. The targets have been previously inspected; five-to-one transfer is not blind evidence.','',
'| Geometry | Rule | Previous omitted error | Frozen with feedback | Refit with feedback, omitted error |','|---|---|---:|---:|---:|']
for v in refits:
 fr=next(f for f in frozen if f['geometry']==v['geometry'] and f['population']=='Chabrier' and f['source']=='distant' and f['model']==v['model'])
 arr=[next(z for z in h['versions'] if z['mode']=='loo' and z['traveling_to_stored_ratio']==0) for h in fr['rows']]
 score=np.sqrt(np.mean([z['deposited_rms']**2 for z in arr]))
 lines.append(f"| {v['geometry']} | {v['model']} | {100*v['old_no_feedback_loo_rms']:.2f} | {100*score:.2f} | {100*v['loo_rms']:.2f} |")
lines+=['','Errors are RMS cumulative-profile differences in percentage points, not velocities, lens residuals or significance. Refit grid spacing differs from the older grid; small changes must not be overinterpreted. Standard geometry is a comparison only. All frozen variations, failed trial counts, parameters, per-target curves and convergence diagnostics are in the JSON.','',
'## Traveling companions','',
'The instantaneous traveling inventory is unknown. At frozen parameters, separate scenarios use traveling/stored ratios 0, 0.01 and 0.1, with a spherical profile proportional to capture density divided by opacity. These are sensitivity assumptions, not measured supply or inferred cosmic ages. Its prescribed energy density contributes to the Newtonian driving potential; total enclosed output adds the traveling component explicitly. It is never normalized away. Relativistic pressure, directional stresses, photon gravity and a covariant field equation are not computed, so this does not determine the true gravity of radiation-like companions. The baseline stellar proxy already includes stellar mass; no separate photon inventory is inferred from it.','',
'## Energy and verification','',
'Potential energy before and after settling is evaluated for stars plus deposited self-gravity in the zero-travel case. The reported released specific binding energy is a required work/radiation channel, not extra stored mass or a completed energy budget. Changes in baryonic motion, support and radiation would have to close that ledger.','',f"Direct shell potential checks pass. Turning companion feedback off reproduces prior curves within {maxzero:.4g}. Refined selected profiles change by at most {maxrefine:.4g}; the largest difference between all-extended and all-compact seeds is {maxseed:.4g}. A converged fixed point is numerical consistency, not a dynamical stability proof. Nonconverged and alternate-seed flags are retained. No new stellar-motion/lensing likelihood fit or completed research goal is claimed."]
finding=['## Outcome','',
'With prior parameters frozen, deposited self-gravity increases retained-geometry omitted errors: binding 20.34 to 21.69 points, released binding 20.15 to 22.66, and occupancy penalty 21.85 to 25.63. Ignoring the deposits in the driving field was therefore a material approximation.','',
'Refitting recovers approximately the previous best performance: binding gives 20.04, released binding 20.59 and occupancy penalty 20.87 points. The earlier partial-migration benchmark was 20.01. The small binding difference is below the profile-refinement scale and is not evidence of an improvement. Standard-geometry occupancy feedback gives 22.93 points versus 23.47 previously, but this comparison is not adoption of expansion or a new lensing success.','',
'The shared retained-geometry binding fit uses s=0.1732 and v_ex=568.37 km/s, with compact fractions about 33.3-42.5%. All frozen runs and refit grid trials converged; all selected alternate seeds converged to essentially the same solution. This supports uniqueness within the tested seeds and parameter families, not a proof of global uniqueness or dynamical stability.','',
'Prescribed traveling energy is not the missing cure in these scenarios. At the old retained-geometry omitted parameters, increasing its ratio from zero to 0.1 changes deposited RMS from 21.69 to 21.86 for binding, 22.66 to 22.67 for released binding and 25.63 to 24.71 for occupancy. These compare deposit response only; additive traveling gravity is separately recorded and remains a relativistic modeling limitation.','',
'Feedback is now represented for the deposited component, but the hypothesis still lacks evolving baryons, physical arrival histories, angular/radiation stresses, and a demonstrated support-and-energy-loss mechanism. It has not become a complete all-energy gravitational model.','']
pos=lines.index('## Energy and verification');lines[pos:pos]=finding
(P/'feedback-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
selected=[v for v in refits if v['geometry']=='companion_regular']
for i,ax in enumerate(axes.flat):
 h=selected[0]['rows'][i];ax.plot(x,h['target_cdf'],'k',label='Fitted halo target')
 for v in selected:ax.plot(x,v['rows'][i]['checks']['loo']['refined_cdf'],label=v['model'])
 ax.set(title=h['Name'],xlabel='Radius / diagnostic region',ylabel='Enclosed / supplied inventory',xlim=(0,1),ylim=(0,1));ax.grid(alpha=.2)
axes.flat[0].legend(fontsize=7);fig.suptitle('Companion gravity feedback: omitted-galaxy predictions, retained geometry')
fig.savefig(P/'feedback.png',dpi=150)
