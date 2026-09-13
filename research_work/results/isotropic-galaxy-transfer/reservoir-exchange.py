"""Frozen-environment two-reservoir exchange with local binding probabilities."""
from pathlib import Path
import json,hashlib
import numpy as np
P=Path(__file__).resolve().parent;G=4.30091e-6
files=['stream-capture-results.json','halo-deposition-map-results.json','third-retention-optics-results.json','migration-variants-results.json']
streams=json.loads((P/files[0]).read_text())['rows'];halos=json.loads((P/files[1]).read_text())['rows'];stars=json.loads((P/files[2]).read_text())['rows'];prior=json.loads((P/files[3]).read_text())['results']
x=np.linspace(0,1,401);ss=np.geomspace(.03,1.,49);vs=np.geomspace(10,3000,41)
params=np.array([(s,v) for s in ss for v in vs])
models=['binding','released_binding','occupancy_penalty']

def setup(row,branch,n):
 h=next(v for v in halos if all(v[k]==row[k] for k in ['Name','geometry','population']))
 st=next(v for v in stars if v['Name']==row['Name'] and v['model']=='attenuated_'+row['population'])
 M=next(v['original_deposit_enclosed_Msun'] for v in h['samples'] if v['aperture']=='5Re')
 edge=np.linspace(0,1,n+1);mid=(edge[:-1]+edge[1:])/2;R=row['R_kpc'];r=mid*R
 base=np.interp(edge,row['r_over_R'],row['case_cdfs'][branch]);dm=np.diff(base)
 a=R/5/1.8153;Ms=st['mass_Msun'];B=G*Ms/(r+a)
 vol=4*np.pi/3*np.diff((edge*R)**3);rho_input=M*dm/vol
 rho_star=Ms*a/(2*np.pi*r*(r+a)**3)
 return edge,dm,base,B,rho_input/rho_star,Ms,a,R,M

def evaluate(row,branch,model,parameters,n=2048):
 edge,dm,base,B,crowd,Ms,a,R,M=setup(row,branch,n);r=(edge[:-1]+edge[1:])/2*R
 curves=[];fractions=[]
 for s,v in parameters:
  strength=B if model!='released_binding' else np.maximum(G*Ms/(s*r+a)-B,0)
  resistance=v*v*(1+crowd) if model=='occupancy_penalty' else v*v
  f=strength/(strength+resistance)
  compact=dm*f;extended=dm-compact
  assert np.min(compact)>=0 and np.min(extended)>=-1e-14
  C=np.r_[0.,np.cumsum(compact)];E=np.r_[0.,np.cumsum(extended)]
  y=np.interp(x,edge,E)+np.interp(x/s,edge,C)
  assert abs(C[-1]+E[-1]-1)<1e-10 and abs(y[-1]-1)<1e-10
  assert np.min(np.diff(y))>=-1e-12
  assert np.min(y-np.interp(x,edge,base))>=-1e-10
  curves.append(y);fractions.append(C[-1])
 return np.array(curves),np.array(fractions),M

results=[];max_drift=0
for geometry in ['companion_regular','standard_flat_FLRW']:
 for population in ['Chabrier','Salpeter']:
  rows=[v for v in streams if v['geometry']==geometry and v['population']==population]
  for branch,label in [(-1,'distant'),(0,'near_1.1R')]:
   for model in models:
    curves=[];fractions=[];errors=[];baseerrors=[]
    for row in rows:
     y,f,_=evaluate(row,branch,model,params);truth=np.interp(x,row['r_over_R'],row['target_cdf']);base=np.interp(x,row['r_over_R'],row['case_cdfs'][branch])
     identity,_,_=evaluate(row,branch,model,[[1.,100.]])
     assert np.max(np.abs(identity[0]-base))<1e-4
     curves.append(y);fractions.append(f);errors.append(np.mean((y-truth)**2,axis=1));baseerrors.append(np.mean((base-truth)**2))
    errors=np.array(errors);best=int(np.argmin(errors.mean(axis=0)));loo=[];per=[]
    for i,row in enumerate(rows):
     ix=int(np.argmin(np.delete(errors,i,axis=0).mean(axis=0)));loo.append(errors[i,ix])
     fine,ff,M=evaluate(row,branch,model,params[[best,ix]],4096)
     drift=float(np.max(np.abs(fine-np.array([curves[i][best],curves[i][ix]]))));max_drift=max(max_drift,drift)
     assert drift<.002
     truth=np.interp(x,row['r_over_R'],row['target_cdf'])
     per.append(dict(Name=row['Name'],supplied_inventory_Msun=M,shared_compact_fraction=float(fractions[i][best]),loo_compact_fraction=float(fractions[i][ix]),loo_parameters=params[ix].tolist(),loo_boundary=bool(ix//41 in [0,48] or ix%41 in [0,40]),refinement_max_cdf_change=drift,baseline_max_error=float(np.max(np.abs(np.interp(x,row['r_over_R'],row['case_cdfs'][branch])-truth))),shared_max_error=float(np.max(np.abs(curves[i][best]-truth))),loo_max_error=float(np.max(np.abs(curves[i][ix]-truth))),shared_cdf=curves[i][best].tolist(),loo_cdf=curves[i][ix].tolist(),target_cdf=truth.tolist()))
    old=next(v for v in prior if v['geometry']==geometry and v['population']==population and v['source']==label and v['model']=='partial')
    result=dict(geometry=geometry,population=population,source=label,model=model,shared_parameters=params[best].tolist(),shared_boundary=bool(best//41 in [0,48] or best%41 in [0,40]),baseline_rms=float(np.sqrt(np.mean(baseerrors))),shared_rms=float(np.sqrt(errors[:,best].mean())),loo_rms=float(np.sqrt(np.mean(loo))),previous_partial_loo_rms=old['loo_rms'],rows=per)
    results.append(result);print(geometry,population,label,model,'base/shared/LOO',*[round(100*result[k],3) for k in ['baseline_rms','shared_rms','loo_rms']],result['shared_parameters'],flush=True)
out=dict(scope='Local two-state probability plus phenomenological inward transport; not self-consistent dynamical equilibrium',input_sha256={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in files},parameter_order=['radius_multiplier_s','exchange_velocity_scale_kms'],max_refinement_cdf_change=max_drift,grid=x.tolist(),results=results)
(P/'reservoir-exchange-results.json').write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8',newline='\n')
lines=['# Local exchange between compact and extended companion reservoirs','',
'Three shared exchange rules were executed with the one-third capture law and supplied inventory unchanged. Each initial shell is split into compact and extended states using its baryonic binding environment. The target halo is used only in fitting and scoring, never as a local driving field.','',
'## Formula provenance and meaning','',
'Known two-state rate balance: df/dt=k_on(1-f)-k_off*f, with fixed-environment stationary fraction f=k_on/(k_on+k_off). The proposed companion-specific ratios are:','',
'| Rule | k_on/k_off | Interpretation |','|---|---|---|',
'| Binding | B(r)/v_ex^2 | Strong binding favors compact states |',
'| Released binding | [B(s*r)-B(r)]/v_ex^2 | Available inward binding-energy change favors the transition |',
'| Occupancy penalty | B(r)/{v_ex^2[1+rho_input/rho_star]} | Larger relative companion density suppresses compact entry |','',
'B(r)=G M_star/(r+a_star) is the positive escape-binding scale of a known Hernquist stellar proxy. Stellar mass is from the ORIGINAL companion fit; a_star=Re/1.8153. The proxy is an approximation, not an independent measured three-dimensional matter field. The crowding factor is a proposed local penalty, not an established capacity law. We do not claim these ratios are unique new mathematics.','',
'The compact shell fraction moves to s*r; the rest remains at r. Thus the cumulative profile is F_final(r)=F_initial(r)-F_selected(r)+F_selected(r/s). This is known conservative transport. s still specifies a phenomenological movement distance; these tests do not independently derive wave support or a stopping radius. The states use the initial environment: no updating of occupation, self-gravity or detailed balance at the new destination has been solved. The arbitrary overall rate cancels; no relaxation time or cosmic age is inferred.','',
'## Shared and omitted-galaxy comparison','',
'Two parameters per rule are shared across six targets; separate five-target fits predict the omitted sixth. All targets have been inspected previously, so this is conditional cross-validation, not blind observational confirmation. RMS values below are cumulative-profile percentage points, not lensing/rotation errors. Standard geometry remains a benchmark only.','',
'| Geometry | Input | Exchange | Shared fit | Omitted transfer | Earlier partial transfer |','|---|---|---|---:|---:|---:|']
for r in results:
 if r['population']=='Chabrier':lines.append(f"| {r['geometry']} | {r['source']} | {r['model']} | {100*r['shared_rms']:.2f} | {100*r['loo_rms']:.2f} | {100*r['previous_partial_loo_rms']:.2f} |")
lines+=['','All population cases, compact fractions, shared/omitted parameters, boundary flags and individual curves are archived. No target-specific compact fraction was fitted directly. Source fields and mass proxies remain conditional inputs; normalized placement is not evidence for an adequate absolute photon supply or a successful projected lensing signal.','',
'## Verification and physical limits','',f"The compact and extended inventories remain nonnegative and sum to the original input. CDFs are monotone, inward ordering holds, and s=1 reproduces the captured profile. Doubling source cells changes cumulative fractions by at most {max_drift:.4g}. Parameter bounds were declared in reservoir-exchange-protocol.md and not expanded after results.",'',
'Mass-equivalent inventory conservation is not a full formation-energy ledger. The released gravitational energy requires a receiving channel, and rate ratios alone do not supply the work, momentum transfer or support needed after migration. The next physical completion would evolve exchange, transport, occupations and the potential together; no first-principles conversion or self-consistent exchange equilibrium has been derived here. No new stellar-motion or lensing fit is claimed.']
finding=['## Outcome','',
'Under retained geometry and distant input, the released-binding rule nearly matches the earlier partial-migration prediction: 20.15 versus 20.01 cumulative RMS points. Binding-only gives 20.34 and the occupancy penalty gives 21.85. None improves the primary retained-geometry benchmark. These differences are not significance estimates, and grid spacing and prior target exposure preclude treating small ranking differences as discoveries.','',
'The shared released-binding fit has s=0.1391 and v_ex=720.84 km/s. Its compact fractions range from about 30.9% to 37.6%, calculated from the initial stellar binding environment rather than fitted separately per target. This is a candidate local rationale for the earlier roughly 36% mobile fraction, not a first-principles derivation. Both models still choose a movement scale phenomenologically.','',
'It retains the same tension: concentrated targets improve, but J1112, J1621 and J1630 worsen relative to no migration. For example, J1112 maximum cumulative error rises from 5.98 to 41.73 points in released-binding omitted prediction, while J0037 falls from 62.25 to 36.44. A lower aggregate error does not mean all galaxies are reproduced.','',
'Under standard comparison geometry, occupancy suppression gives 23.47 points for distant input versus 23.55 without migration and 24.32 for earlier partial migration; the near-input result is 23.56. This small aggregate change is not a robust joint solution or a new lensing result. Standard geometry is not adopted as a premise for the hypothetical model.','',
'The useful outcome is a local exchange candidate with performance close to the best phenomenological split, using the same two shared parameters. Its missing physics is the simultaneous evolution of exchange, inward transport, occupation and self-gravity. These frozen initial-field rates do not yet establish that final reservoirs are in equilibrium or stable.','']
pos=lines.index('## Verification and physical limits');lines[pos:pos]=finding
(P/'reservoir-exchange-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
selected=[r for r in results if r['geometry']=='companion_regular' and r['population']=='Chabrier' and r['source']=='distant']
for i,ax in enumerate(axes.flat):
 h=selected[0]['rows'][i];ax.plot(x,h['target_cdf'],'k',label='Fitted halo target')
 for result in selected:ax.plot(x,result['rows'][i]['loo_cdf'],label=result['model'])
 ax.set(title=h['Name'],xlabel='Radius / diagnostic region',ylabel='Enclosed / supplied inventory',xlim=(0,1),ylim=(0,1));ax.grid(alpha=.2)
axes.flat[0].legend(fontsize=7);fig.suptitle('Local reservoir exchange: omitted-galaxy predictions, retained geometry')
fig.savefig(P/'reservoir-exchange.png',dpi=150)
