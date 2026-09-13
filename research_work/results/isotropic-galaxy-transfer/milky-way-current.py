"""Frozen Milky Way circular-speed transfer of exact-third and current feedback laws."""
from pathlib import Path
import json,hashlib,ast
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import gamma,gammainc
from numpy.polynomial.legendre import leggauss
from numba import njit
P=Path(__file__).resolve().parent;G=4.30091727003628e-6
paths=[P/'third-radiation-retention-results.json',P/'redistribution-compact-results.json',P/'migration-variants-results.json',P/'feedback-results.json',P.parent/'joint-galaxy-audit/milky-way-predictions.json',P.parent/'milky-way-capture/inputs.json',P/'feedback.py',P.parent/'joint-galaxy-audit/run.py',P.parent/'joint-galaxy-audit/gas.py']
cp=json.loads(paths[0].read_text())['models']['attenuated'];assert cp['q']==1/3
rotation=json.loads(paths[1].read_text())['models']['compact_partial']['parameters']
halo=next(v for v in json.loads(paths[2].read_text())['results'] if v['geometry']=='companion_regular' and v['population']=='Chabrier' and v['source']=='distant' and v['model']=='partial')['shared_parameters']
fb={v['model']:v['shared_parameters'] for v in json.loads(paths[3].read_text())['refits'] if v['geometry']=='companion_regular'}
old=json.loads(paths[4].read_text());observations=json.loads(paths[5].read_text())['eilers']['rows']
# Reuse only the pure spherical-shell potential, never execute another driver.
tree=ast.parse(paths[6].read_text());ns=dict(np=np,njit=njit)
exec(compile(ast.Module(body=[v for v in tree.body if isinstance(v,ast.FunctionDef) and v.name=='potential_at'],type_ignores=[]),str(paths[6]),'exec'),ns);potential_at=ns['potential_at']
components={'I':[(460*2.32e7,0,.3),(1700*2.32e7,5.3,.25),(1700*2.32e7,2.6,.8)],'II':[(1600*2.32e7,4.8,.25),(1700*2.32e7,2.,.8)]}
gas=[(53.1,7.,4.),(2180.,1.5,12.)]

def capture(rd,L,n=4096,order=96):
 r=np.r_[0.,np.geomspace(1e-5,1e5,n)];a=rd*cp['scale_to_disk'];xx=r/a;mu,w=leggauss(order)
 t=xx[:,None]*mu;B2=1+xx[:,None]**2*(1-mu*mu);B=np.sqrt(B2)
 tau=cp['k0_per_kpc']*a*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3));J=np.exp(-np.maximum(tau,0))@w/2
 X=(L/1e9)/rd**2;eta=X**(1/3)/(1+X**(1/3));C=2*cp['C_Msun_kpc3']*eta
 rho=C*J/(1+xx*xx)**2;m=4*np.pi*cumulative_trapezoid(r*r*rho,r,initial=0)
 T=np.pi*cp['k0_per_kpc']*a/2;area=np.pi*a*a*(T**(2/3)*gamma(1/3)*gammainc(1/3,T)-1+np.exp(-T));total=C/cp['k0_per_kpc']*area
 err=abs(m[-1]/total-1);assert err<.002
 return r,m,eta,err

def baryonic_binding(edge,variant,order):
 r=(edge[:-1]+edge[1:])/2;mu,w=leggauss(order)
 cyl=r[:,None]*np.sqrt(1-mu*mu);z=r[:,None]*mu;B=np.zeros(len(r))
 for M,a,b in components[variant]:B+=G*M*(1/np.sqrt(cyl*cyl+(a+np.sqrt(z*z+b*b))**2))@w/2
 # Monopole of the thin gas disks; no dark halo or speed-fitting input.
 safe=np.maximum(edge,1e-30);sigma=sum(S*1e6*np.exp(-hole/safe-safe/h) for S,h,hole in gas)
 mgas=2*np.pi*cumulative_trapezoid(edge*sigma,edge,initial=0)
 B+=G*potential_at(r,r,np.diff(mgas))
 return B

@njit
def feedback(edge,mass,B,s,v,released,seed):
 r=(edge[:-1]+edge[1:])/2;dm=np.diff(mass);f=np.full(len(r),seed);dest=s*r
 Bs=np.interp(dest,r,B);ok=False;res=1.
 for it in range(500):
  cm=dm*f;em=dm*(1-f)
  totalB=B+G*(potential_at(r,r,em)+potential_at(r,dest,cm))
  totalBs=Bs+G*(potential_at(dest,r,em)+potential_at(dest,dest,cm))
  strength=np.maximum(totalBs-totalB,0) if released else totalB
  p=strength/(strength+v*v);res=np.max(np.abs(p-f))
  if res<1e-8:ok=True;break
  f=.65*f+.35*p
 C=np.empty(len(edge));C[0]=0;C[1:]=np.cumsum(dm*f)
 E=mass-C;out=E+np.interp(edge/s,edge,C)
 return out,float(np.sum(dm*f)/mass[-1]),ok,res,it+1

def run(variant,rd,lf,n=4096,order=96):
 Mstar=sum(t[0] for t in components[variant]);L=Mstar/.5*lf
 edge,mass,eta,err=capture(rd,L,n,order);B=baryonic_binding(edge,variant,order)
 base=[v for v in old if v['baryons']==variant and v['model']=='baryons' and v['observable']=='vc_kms'];assert len(base)==38
 R=np.array([v['R_kpc'] for v in base]);vb=np.array([v['predicted'] for v in base]);obs=np.array([v['observed'] for v in base]);assert np.allclose(obs,[o['vc_kms'] for o in observations])
 curves={'ordinary_matter':np.zeros(len(edge)),'exact_third':mass}
 f,b=rotation;curves['rotation_trained_partial']=(1-f)*mass+f*np.interp(edge/np.exp(b),edge,mass)
 b,f=halo;curves['halo_trained_partial']=(1-f)*mass+f*np.interp(edge/np.exp(-b),edge,mass)
 diagnostics={}
 for model in ['binding','released_binding']:
  s,v=fb[model];y,f,ok,res,it=feedback(edge,mass,B,s,v,model=='released_binding',0.)
  y2,f2,ok2,res2,it2=feedback(edge,mass,B,s,v,model=='released_binding',1.)
  gap=float(np.max(np.abs(y-y2))/mass[-1]);assert ok and ok2 and gap<1e-6
  assert abs(y[-1]/mass[-1]-1)<1e-10 and np.min(np.diff(y))>=-1e-3
  curves['feedback_'+model]=y;diagnostics[model]=dict(compact_fraction=f,residual=res,iterations=it,seed_max_inventory_fraction_difference=gap)
 rows=[];scores=[]
 for model,ma in curves.items():
  pred=np.sqrt(vb*vb+G*np.interp(R,edge,ma)/R)
  for i,v in enumerate(base):rows.append(dict(model=model,R_kpc=float(R[i]),observed_kms=float(obs[i]),predicted_kms=float(pred[i]),ordinary_kms=float(vb[i]),residual_kms=float(pred[i]-obs[i]),split=v['split'],err_minus_kms=observations[i]['err_minus_kms'],err_plus_kms=observations[i]['err_plus_kms']))
  for subset in ['all','inner','outer']:
   use=np.ones(len(R),dtype=bool) if subset=='all' else np.array([v['split']==subset for v in base]);e=pred[use]-obs[use]
   scores.append(dict(model=model,subset=subset,n=int(use.sum()),RMSE_kms=float(np.sqrt(np.mean(e*e))),bias_kms=float(np.mean(e))))
 return dict(baryons=variant,Rd_kpc=rd,luminosity_proxy_factor=lf,L36_proxy_Lsun=L,eta=eta,total_deposit_Msun=float(mass[-1]),analytic_total_relative_error=err,feedback=diagnostics,scores=scores,rows=rows)

runs=[]
for variant in ['I','II']:
 for rd in [2.08,2.6,3.12]:
  for lf in [.5,1.,2.]:
   v=run(variant,rd,lf);runs.append(v);print(variant,rd,lf,[(s['model'],round(s['RMSE_kms'],2)) for s in v['scores'] if s['subset']=='all'],flush=True)
refinement=[]
for variant in ['I','II']:
 coarse=next(v for v in runs if v['baryons']==variant and v['Rd_kpc']==2.6 and v['luminosity_proxy_factor']==1.)
 fine=run(variant,2.6,1.,8192,192)
 delta=max(abs(a['predicted_kms']-b['predicted_kms']) for a,b in zip(coarse['rows'],fine['rows']));assert delta<.5
 refinement.append(dict(baryons=variant,max_speed_change_kms=delta,fine_scores=fine['scores']))
out=dict(scope='Frozen parameters predict inferred Milky Way circular speeds; not individual star orbits',input_sha256={str(f.relative_to(P.parents[2])):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},frozen_rotation_partial=rotation,frozen_halo_partial=halo,frozen_feedback_parameters=fb,runs=runs,refinement=refinement)
(P/'milky-way-current-results.json').write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8',newline='\n')
lines=['# Current model predictions for Milky Way circular speeds','',
'No parameters were fitted to Milky Way speeds. All capture, migration and exchange parameters are transferred from previous external-galaxy fits. The exact one-third law remains fixed. Both archived ordinary-matter baselines and every predeclared input sensitivity are retained.','',
'## Data and scope','',
'The 38 cached [Eilers et al.](https://arxiv.org/abs/1810.09466) circular-speed summaries cover about 5-25 kpc and are inferred from stellar kinematics using an axisymmetric Jeans analysis. They are not individual star velocities or direct acceleration readings. The source reports 2-5% systematic uncertainty near the Sun; the tabulated small formal errors do not capture every modeling uncertainty. Its separate dark-halo inference is not used here. The inner barred Galaxy is outside this comparison.','',
'Rd=2.6 kpc follows the earlier [Juric et al. star-count proxy](https://arxiv.org/abs/astro-ph/0510520). It is not a measured total 3.6-micron light scale. Needed luminosity is approximated by archived stellar mass / 0.5 Msun per Lsun. This nominal mass-to-light assumption is a conditional proxy, not a new photometric measurement. Predeclared sensitivity uses Rd=2.08,2.6,3.12 kpc and luminosity factors 0.5,1,2, without selecting a winning case.','',
'## Formula provenance','',
'Known gravity: v_c^2(R)=v_b^2(R)+G M_dep(<R)/R. The spherical deposited density and exact-third retention are the existing proposed capture model. Original ordinary disk/gas forces remain as archived. No dark halo is added.','',
'The 0.37% partial redistribution is trained on SPARC rotation curves; the much larger halo-trained partial fraction comes from normalized halo-profile fits. They are different calibrations. Current binding and released-binding parameters come from the six-lens retained-geometry feedback fits. None was selected using this Milky Way result.','',
'Feedback now uses a spherical angular average of the actual archived stellar Miyamoto-Nagai/Plummer potential plus the spherical monopole of its thin gas disks. It updates the compact and extended companion field. This is an explicitly spherical approximation for migration in a flattened galaxy; it does not infer three-dimensional bulge or vertical orbits. Gas and stellar baseline uncertainty is retained as two alternatives, not marginalized into a confidence interval.','',
'## Fiducial predictions','',
'RMS discrepancy across all 38 circular-speed bins, km/s:','',
'| Model | Ordinary baseline I | Ordinary baseline II |','|---|---:|---:|']
fid=[next(v for v in runs if v['baryons']==b and v['Rd_kpc']==2.6 and v['luminosity_proxy_factor']==1.) for b in ['I','II']]
for model in ['ordinary_matter','exact_third','rotation_trained_partial','halo_trained_partial','feedback_binding','feedback_released_binding']:
 vals=[next(s['RMSE_kms'] for s in v['scores'] if s['model']==model and s['subset']=='all') for v in fid]
 lines.append(f"| {model} | {vals[0]:.2f} | {vals[1]:.2f} |")
lines+=['','These are descriptive speed residuals, not chi-square significance or acceptance thresholds. The data were examined in earlier work, so frozen transfer is not a blind observation test.','',
'## Selected actual-versus-predicted rows','',
'Baseline I, fiducial inputs. Each selected radius is a real catalog row, not an interpolated measurement.','',
'| R kpc | Inferred circular speed | Ordinary matter | Exact third | Rotation-trained partial | Feedback binding |','|---|---:|---:|---:|---:|---:|']
for target in [5.27,8.19,12.25,15.22,20.27,24.82]:
 row=min([r for r in fid[0]['rows'] if r['model']=='exact_third'],key=lambda r:abs(r['R_kpc']-target));R=row['R_kpc']
 vals=[next(r['predicted_kms'] for r in fid[0]['rows'] if r['R_kpc']==R and r['model']==m) for m in ['ordinary_matter','exact_third','rotation_trained_partial','feedback_binding']]
 lines.append(f"| {R:.2f} | {row['observed_kms']:.2f} | "+' | '.join(f'{v:.2f}' for v in vals)+' |')
lines+=['','## Sensitivity and checks','',
'Every 18 baseline/scale/luminosity scenario is in the JSON with all 38 predictions, residuals and inner/outer scores. No parameter was adjusted to minimize these errors. Predicting individual speeds additionally requires orbital phases and velocity distributions; the radial curve cannot establish that.','',f"Doubling radial and angular resolution changes a fiducial predicted speed by at most {max(v['max_speed_change_kms'] for v in refinement):.4g} km/s. Analytic capture inventory checks pass. Feedback conserves deposited inventory and converges from opposite initial fractions. Attenuation paths extend to infinity; the outer radial integration is a numerical approximation, not a universe-size limit.",'',
'Absolute photon supply, formation work, support, source histories and traveling radiation stresses remain unclosed. The previous provisional vertical summaries are not promoted to independent constraints here. This run tests actual circular-speed predictions, without claiming a complete Milky Way model or completion of the six broader goals.']
finding=['## Outcome','',
'The fixed exact-third model substantially improves the Milky Way circular-speed curve relative to both ordinary baselines: RMS 6.76/10.41 km/s versus 52.57/62.34. The SPARC rotation-trained small redistribution changes this to 6.70/10.02. These are conditional predictions using the declared disk-scale and luminosity proxies, not fits to these speeds.','',
'The recent halo-trained large partial fraction overconcentrates the Milky Way and gives 35.57/28.56 km/s. The newer local feedback adjusts its mobile fraction to this galaxy instead of imposing the same approximately 36% split: binding gives about 14.2%/13.5% mobile, and released binding 9.6%/8.1%, for the full modeled inventory. Released-binding feedback improves baseline II to 5.75 km/s but worsens baseline I to 12.03. Therefore no variant dominates both ordinary-matter choices.','',
'Near the Sun at the actual 8.19 kpc row, inferred speed is 228.86 km/s. Under baseline I, ordinary matter predicts 188.43, exact third 226.60, and the rotation-trained partial variant 227.24. This is a useful radial prediction; it neither establishes individual stellar orbits nor proves the photon-to-gravity mechanism.','',
'The scale/luminosity sensitivity is material. For baseline I, exact-third RMS ranges roughly 6.6-18.5 km/s across the declared grid, while some feedback variants worsen substantially at the larger disk scale. These ranges are exploratory input sensitivity, not confidence intervals; the fiducial choice is retained, and no best-case proxy is selected after looking at the data. Further progress should independently constrain those light/matter inputs and test another kinematic sample, rather than adjusting them to force the curve to match.','']
pos=lines.index('## Sensitivity and checks');lines[pos:pos]=finding
(P/'milky-way-current-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
for ax,v in zip(axes,fid):
 obs=[r for r in v['rows'] if r['model']=='ordinary_matter'];R=[r['R_kpc'] for r in obs]
 ax.errorbar(R,[r['observed_kms'] for r in obs],yerr=[[r['err_minus_kms'] for r in obs],[r['err_plus_kms'] for r in obs]],fmt='k.',label='Inferred circular speed (formal errors)')
 for model in ['ordinary_matter','exact_third','rotation_trained_partial','feedback_binding','feedback_released_binding']:
  ax.plot(R,[r['predicted_kms'] for r in v['rows'] if r['model']==model],label=model)
 ax.set(title='Ordinary baseline '+v['baryons'],xlabel='Galactocentric radius (kpc)',ylabel='Circular speed (km/s)');ax.grid(alpha=.2)
axes[0].legend(fontsize=7);fig.suptitle('Milky Way: frozen external-galaxy parameters, no speed refit')
fig.savefig(P/'milky-way-current.png',dpi=150)
