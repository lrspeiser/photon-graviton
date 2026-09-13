"""Variational wave support and circular-orbit support of fixed companion inventory."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.special import roots_genlaguerre,gamma,erf
from scipy.integrate import quad
P=Path(__file__).resolve().parent;G=4.30091e-6
files=['stream-capture-results.json','halo-deposition-map-results.json','third-retention-optics-results.json','migration-variants-results.json']
streams=json.loads((P/files[0]).read_text())['rows'];halos=json.loads((P/files[1]).read_text())['rows'];stars=json.loads((P/files[2]).read_text())['rows'];prior=json.loads((P/files[3]).read_text())
x=np.linspace(0,1,401);hs=np.logspace(1,5,41);ks=np.r_[0.,np.logspace(-8,-1,22)];lams=np.linspace(.05,1.5,117)

def inputs(row):
 h=next(v for v in halos if all(v[k]==row[k] for k in ['Name','geometry','population']))
 st=next(v for v in stars if v['Name']==row['Name'] and v['model']=='attenuated_'+row['population'])
 M=next(v['original_deposit_enclosed_Msun'] for v in h['samples'] if v['aperture']=='5Re')
 return M,st['mass_Msun'],row['R_kpc']/5/1.8153

def wave_radii(M,Ms,ast,params,n=64):
 z,w=roots_genlaguerre(n,.5);t=np.sqrt(z);w=w/gamma(1.5)
 h=params[:,0];K=params[:,1]
 def derivative(a):
  # a^4 d(e)/da is better scaled over the full radius range.
  return -1.5*h*h*a+G*M*a*a/np.sqrt(2*np.pi)+G*Ms*a**4*((t[None,:]/(a[:,None]*t+ast)**2)@w)-3*K*M/(2*np.pi)**1.5
 lo=np.full(len(params),np.log(1e-4));hi=np.full(len(params),np.log(1e5))
 dlo=derivative(np.exp(lo));dhi=derivative(np.exp(hi))
 for _ in range(65):
  mid=(lo+hi)/2;d=derivative(np.exp(mid));low=d<0;lo=np.where(low,mid,lo);hi=np.where(low,hi,mid)
 a=np.exp((lo+hi)/2);boundary=(dlo>=0)|(dhi<=0)
 return a,boundary

def wave_energy(a,M,Ms,ast,h,K,n=128):
 z,w=roots_genlaguerre(n,.5);w/=gamma(1.5)
 return 3*h*h/(4*a*a)-G*M/(np.sqrt(2*np.pi)*a)-G*Ms*np.sum(w/(a*np.sqrt(z)+ast))+K*M/((2*np.pi)**1.5*a**3)

def wave_cdf(r,a):
 z=np.asarray(r)/a
 return erf(z)-2/np.sqrt(np.pi)*z*np.exp(-z*z)

def particle(row,branch,lam,n=2048):
 M,Ms,ast=inputs(row);q=(np.arange(n)+.5)/n
 r0=np.interp(q,row['case_cdfs'][branch],np.array(row['r_over_R'])*row['R_kpc'])
 mb=lambda r:Ms*r*r/(r+ast)**2
 ell=lam*lam*r0*(mb(r0)+q*M)
 lo=np.zeros(n);hi=r0*max(1.,lam*lam)
 for _ in range(55):
  r=(lo+hi)/2;low=r*(mb(r)+q*M)<ell;lo=np.where(low,r,lo);hi=np.where(low,hi,r)
 r=(lo+hi)/2
 if not np.all(np.diff(r)>0):return np.full(len(x),np.nan),None,None
 balance=float(np.max(np.abs(r*(mb(r)+q*M)/ell-1)))
 density_derivative=M*np.gradient(q,r);star_derivative=2*Ms*r*ast/(r+ast)**3
 kappa2=G*((mb(r)+q*M)/r**3+(star_derivative+density_derivative)/r**2)
 assert np.min(kappa2)>0 and balance<1e-9
 return np.interp(x*row['R_kpc'],r,q,left=0,right=1),balance,float(kappa2.min())

results=[];max_wave_radius_drift=0;max_wave_energy_quad=0;max_particle_drift=0
for geometry in ['companion_regular','standard_flat_FLRW']:
 for population in ['Chabrier','Salpeter']:
  rows=[r for r in streams if r['geometry']==geometry and r['population']==population]
  for model in ['wave','repulsive_wave','circular_particles']:
   branches=[(-1,'distant'),(0,'near_1.1R')] if model=='circular_particles' else [(-1,'relaxed_no_source_shape')]
   for branch,label in branches:
    params=np.array([(h,0.) for h in hs]) if model=='wave' else np.array([(h,k) for h in hs for k in ks]) if model=='repulsive_wave' else lams[:,None]
    predictions=[];scores=[];extras=[];bases=[]
    for row in rows:
     M,Ms,ast=inputs(row);truth=np.interp(x,row['r_over_R'],row['target_cdf']);base=np.interp(x,row['r_over_R'],row['case_cdfs'][branch])
     if model!='circular_particles':
      a,bounds=wave_radii(M,Ms,ast,params)
      curves=np.array([wave_cdf(x*row['R_kpc'],v) for v in a]);extras.append((a,bounds))
     else:
      curves=np.array([particle(row,branch,float(p[0]),1024)[0] for p in params]);extras.append(None)
      unchanged=particle(row,branch,1.,4096)[0];assert np.max(np.abs(unchanged-base))<.001
     predictions.append(curves);scores.append(np.where(np.all(np.isfinite(curves),axis=1),np.mean((curves-truth)**2,axis=1),np.inf));bases.append(np.mean((base-truth)**2))
    scores=np.array(scores);best=int(np.argmin(scores.mean(axis=0)));per=[];loo=[]
    for i,row in enumerate(rows):
     ix=int(np.argmin(np.delete(scores,i,axis=0).mean(axis=0)));assert np.isfinite(scores[i,ix]), 'Omitted galaxy has no valid circular equilibrium for training-selected parameter';loo.append(scores[i,ix]);M,Ms,ast=inputs(row);diagnostics={}
     for mode,j in [('shared',best),('loo',ix)]:
      p=params[j]
      if model!='circular_particles':
       a=float(extras[i][0][j]);fine,bound=wave_radii(M,Ms,ast,p[None,:],128);af=float(fine[0]);drift=abs(a/af-1);max_wave_radius_drift=max(max_wave_radius_drift,drift)
       ee=wave_energy(af,M,Ms,ast,*p);delta=.002
       curvature=(wave_energy(af*np.exp(delta),M,Ms,ast,*p)-2*ee+wave_energy(af*np.exp(-delta),M,Ms,ast,*p))/delta**2
       starquad=quad(lambda u:4/np.sqrt(np.pi)*u*u*np.exp(-u*u)*(-G*Ms/(af*u+ast)),0,np.inf,epsabs=1e-6)[0]
       z,w=roots_genlaguerre(128,.5);starapprox=-G*Ms*np.sum(w/gamma(1.5)/(af*np.sqrt(z)+ast));err=abs(starapprox/starquad-1);max_wave_energy_quad=max(max_wave_energy_quad,err)
       diagnostics[mode]=dict(radius_kpc=a,refined_radius_kpc=af,radius_boundary=bool(bound[0]),radius_relative_refinement=drift,specific_energy=ee,gaussian_log_radius_curvature=curvature,gaussian_radial_minimum=bool(curvature>0 and not bound[0]),stellar_potential_quadrature_relative_error=err,wave_characteristic_speed_kms=float(np.sqrt(1.5)*p[0]/af),fraction_outside_R=float(1-wave_cdf(row['R_kpc'],af)))
      else:
       fine,balance,kmin=particle(row,branch,float(p[0]),4096);drift=float(np.max(np.abs(fine-predictions[i][j])));max_particle_drift=max(max_particle_drift,drift)
       diagnostics[mode]=dict(refinement_max_cdf_difference=drift,force_balance_relative_error=balance,min_epicyclic_frequency_squared=kmin,individual_circular_orbits_stable=bool(kmin>0),fraction_outside_R=float(1-fine[-1]))
     per.append(dict(Name=row['Name'],supplied_inventory_Msun=M,loo_parameters=params[ix].tolist(),loo_grid_boundary=bool(any(params[ix,k] in [params[:,k].min(),params[:,k].max()] for k in range(params.shape[1]) if np.ptp(params[:,k])>0)),baseline_max_error=float(np.max(np.abs(np.interp(x,row['r_over_R'],row['case_cdfs'][branch])-np.interp(x,row['r_over_R'],row['target_cdf'])))),shared_max_error=float(np.max(np.abs(predictions[i][best]-np.interp(x,row['r_over_R'],row['target_cdf'])))),loo_max_error=float(np.max(np.abs(predictions[i][ix]-np.interp(x,row['r_over_R'],row['target_cdf'])))),diagnostics=diagnostics,shared_cdf=predictions[i][best].tolist(),loo_cdf=predictions[i][ix].tolist(),target_cdf=np.interp(x,row['r_over_R'],row['target_cdf']).tolist()))
    old=next(v for v in prior['results'] if v['geometry']==geometry and v['population']==population and v['source']==('near_1.1R' if branch==0 else 'distant') and v['model']=='partial')
    result=dict(invalid_parameter_counts=[int(np.sum(~np.isfinite(v))) for v in scores],geometry=geometry,population=population,model=model,source=label,shared_parameters=params[best].tolist(),shared_grid_boundary=bool(any(params[best,k] in [params[:,k].min(),params[:,k].max()] for k in range(params.shape[1]) if np.ptp(params[:,k])>0)),baseline_rms=float(np.sqrt(np.mean(bases))),previous_partial_loo_rms=old['loo_rms'],shared_rms=float(np.sqrt(scores[:,best].mean())),loo_rms=float(np.sqrt(np.mean(loo))),rows=per)
    results.append(result);print(geometry,population,model,label,'base/shared/LOO',*[round(result[k]*100,3) for k in ['baseline_rms','shared_rms','loo_rms']],'params',result['shared_parameters'],flush=True)
out=dict(scope='Gaussian variational wave support and circular-orbit equilibria; not complete formation or stability proofs',input_sha256={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in files},wave_max_relative_radius_refinement=max_wave_radius_drift,wave_stellar_potential_max_quadrature_error=max_wave_energy_quad,particle_max_cdf_refinement=max_particle_drift,grid=x.tolist(),results=results)
(P/'wave-particle-results.json').write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8',newline='\n')
lines=['# Wave support versus orbiting particles','',
'The capture law and original deposited inventory remain fixed. These are equilibrium candidates for the captured component, using known physics with a proposed photon-companion origin. They do not establish that ordinary massless gravitons behave as massive scalar waves or as slow orbiting particles.','',
'## Equations and provenance','',
'Wave trial density: rho=M exp(-r^2/a^2)/(pi^1.5 a^3). Minimize E/M=3 h_eff^2/(4a^2)-GM/(sqrt(2pi)a)+<Phi_star>+K M/[(2pi)^1.5 a^3]. Here h_eff=hbar/m, and K>=0 corresponds to repulsive P=K rho^2 support. The Gaussian integrals, gradient-energy support and repulsive pressure are known structures. Our choice to put the entire captured inventory into this relaxed scalar trial state is a proposed model assumption. K=0 and K>=0 are tested separately.','',
'The [Hui et al. review](https://arxiv.org/abs/1610.08297) discusses massive scalar wave support and solitonic cores. We use local nonrelativistic support concepts, not its cosmological formation assumptions. Our Gaussian restriction is not a computed Schrodinger-Poisson ground state, and a whole halo need not occupy one ground state. No claim of a novel wave-support equation is made.','',
'Particle rule: assign j(q)=lambda sqrt[G r0(q)(M_star(<r0)+qM)] and solve j^2/G=r[M_star(<r)+qM]. The captured-profile quantile r0(q) sets the initial angular-momentum scale. Randomly oriented circular orbits give a spherical supported distribution. Circular force balance and angular-momentum conservation are known mechanics; the shared lambda and the photon-to-bound-particle interpretation are proposed inputs, not derived from incident waves. lambda=1 recovers the captured profile. This is an orbit construction, not an Eddington inversion.','',
'The [phase-space consistency analysis of Lacroix et al.](https://arxiv.org/abs/1805.02403) emphasizes that density and orbital distribution must be mutually consistent. Our positive circular-orbit construction addresses a restricted distribution; it does not derive a realistic orbital population or prove collective stability.','',
'## Comparison','',
'Supplied inventory comes from the ORIGINAL deposit inside diagnostic R=5 Re, not target NFW mass. A Hernquist stellar proxy uses the original companion-fit stellar mass and a_star=Re/1.8153. Material placed beyond R is retained in the model and not renormalized into the inner halo. Wave equilibrium forgets the initial stream shape, an assumed complete relaxation; particles retain the distant/near input through r0.','',
'Fit shared parameters on six targets; separately fit five and freeze for the sixth. Existing halo targets were already inspected. Errors are RMS cumulative-profile differences in percentage points, not observed velocity or lens residuals. Standard geometry remains a comparison only.','',
'| Geometry | Support | Source | Shared fit | Omitted-galaxy transfer | Previous partial transfer |','|---|---|---|---:|---:|---:|']
for r in results:
 if r['population']=='Chabrier':lines.append(f"| {r['geometry']} | {r['model']} | {r['source']} | {100*r['shared_rms']:.2f} | {100*r['loo_rms']:.2f} | {100*r['previous_partial_loo_rms']:.2f} |")
lines+=['','Both population proxies, all shared/omitted parameters, boundaries, individual errors, curves and support diagnostics are in the JSON. Changes of shape are not absolute-mass success. Our inferred NFW targets are restricted fits, some with poor stellar fits or parameter boundaries, not unique measured halos.','',
'## Support and numerical checks','',
'Wave energy curvature is checked against Gaussian expansion/contraction at the selected radius, including self-gravity and the stellar potential. Positive curvature is stability only in this one-radius trial family. It does not prove stability to arbitrary field perturbations, excited-state mixing, or nonlinear evolution. The scalar branch introduces an effective massive bound state and does not explain its production by photons.','',
'Particle force balance and positive epicyclic frequency are checked. These demonstrate stability of an individual circular orbit in the fixed positive-density potential, not collective halo stability. Angular momentum acquisition and disposal of energy during settling remain unmodeled.','',f"Refinement: wave radius changes by up to {max_wave_radius_drift:.4g} on doubling stellar-potential quadrature; direct stellar-potential integration differs by up to {max_wave_energy_quad:.4g}; quadrupling particle shells changes cumulative fractions by at most {max_particle_drift:.4g}. All particle force-balance and input-profile recovery checks pass. Boundaries and radius stability flags are preserved.",'',
'Full source history, relaxation times, photon conversion into bound states, support stresses for lensing, and a total energy ledger through formation remain open. No new rotation or lensing likelihood fit is claimed by an equilibrium shape comparison.']
findings=['## Outcome','',
'None of the three support candidates outperforms the earlier partial-migration prescription in retained-geometry omitted-target testing. Distant-input no-migration RMS was 26.74 points; the one-parameter wave model gives 39.08, the repulsive wave model 33.25, and circular particles 25.78. Earlier partial migration gave 20.01. Near-input particles improve their own baseline to 26.85 but also fail to beat the earlier 20.51 partial result. Standard-geometry transfer also fails to improve the distant no-migration benchmark.','',
'The wave trial states have positive energy curvature under Gaussian radial expansion/contraction: support can produce a preferred size, but that does not make it the correct halo size or shape. The repulsive shared fit selects the minimum h_eff allowed, so its extra pressure dominates rather than establishing a detected wave scale. Its fit does not transfer well. A Gaussian containing the entire inventory cannot represent arbitrary compact cores plus extended envelopes; these results do not exclude mixed/excited wave states or all wave-supported halos.','',
'The selected particle solutions satisfy circular force balance and positive individual-orbit epicyclic frequencies. Candidate parameters that reverse the ordering of mass shells were rejected as inconsistent with this ordered circular construction. Counts are recorded; rejection is not a claim that no more general crossing-orbit distribution exists. All shared and omitted selections have valid ordered solutions. Collective stability has not been established.','',
'What we gained is a restricted support calculation instead of a chosen stopping radius. What remains missing is a common formation and population rule that yields both a compact component and an extended one, while accounting for where binding energy and angular momentum go. No candidate earned a new motion/lensing fit by beating our existing profile benchmark.','']
pos=lines.index('## Support and numerical checks');lines[pos:pos]=findings
(P/'wave-particle-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
selected=[r for r in results if r['geometry']=='companion_regular' and r['population']=='Chabrier' and r['source']!='near_1.1R']
for i,ax in enumerate(axes.flat):
 h=selected[0]['rows'][i];ax.plot(x,h['target_cdf'],'k',label='Fitted halo target')
 for result in selected:ax.plot(x,result['rows'][i]['loo_cdf'],label=result['model'])
 ax.set(title=h['Name'],xlabel='Radius / diagnostic region',ylabel='Enclosed / supplied inventory',xlim=(0,1),ylim=(0,1));ax.grid(alpha=.2)
axes.flat[0].legend(fontsize=7)
fig.suptitle('Support candidates: omitted-galaxy predictions, retained geometry')
fig.savefig(P/'wave-particle.png',dpi=150)
