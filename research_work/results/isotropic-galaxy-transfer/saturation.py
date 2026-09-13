"""Saturable bound-state allocation with explicit unbound energy accounting."""
from pathlib import Path
import hashlib,json
import numpy as np
from numba import njit
P=Path(__file__).resolve().parent
files=['stream-capture-results.json','halo-deposition-map-results.json','third-retention-optics-results.json','migration-variants-results.json']
streams=json.loads((P/files[0]).read_text())['rows']
halos=json.loads((P/files[1]).read_text())['rows']
stars=json.loads((P/files[2]).read_text())['rows']
previous=json.loads((P/files[3]).read_text())
bs=np.linspace(0,np.log(100),33)
models=['constant','stellar_density','baryonic_gravity']

@njit
def allocate(m,cap,s,reverse=False):
    n=len(m);bound=np.zeros(n);blocked=0.
    for step in range(n):
        i=n-1-step if reverse else step
        amount=m[i];j0=min(i,int(s*(i+.5)))
        for j in range(j0,i+1):
            room=max(cap[j]-bound[j],0.)
            take=min(room,amount);bound[j]+=take;amount-=take
            if amount<1e-15:break
        blocked+=amount
    return bound,blocked

@njit
def sweep(m,basecap,params):
    out=np.zeros((len(params),len(m)+1));blocked=np.zeros(len(params))
    for k in range(len(params)):
        b=params[k,0];c=10**params[k,1]
        bound,left=allocate(m,basecap*c,np.exp(-b));blocked[k]=left
        out[k,1:]=np.cumsum(bound)
    return out,blocked

def setup(row,model,branch,n):
    h=next(h for h in halos if all(h[k]==row[k] for k in ['Name','geometry','population']))
    st=next(v for v in stars if v['Name']==row['Name'] and v['model']=='attenuated_'+row['population'])
    R=row['R_kpc'];edge=np.linspace(0,R,n+1);x=edge/R;mid=(edge[1:]+edge[:-1])/2
    M0=next(v['original_deposit_enclosed_Msun'] for v in h['samples'] if v['aperture']=='5Re')
    cdf=np.interp(x,row['r_over_R'],row['case_cdfs'][branch]);m=np.maximum(np.diff(cdf),0)
    vol=4*np.pi/3*np.diff(edge**3);a=(R/5)/1.8153;Ms=st['mass_Msun']
    if model=='constant':cap=vol/M0
    elif model=='stellar_density':cap=np.diff(Ms*edge**2/(edge+a)**2)/M0
    else:
        # G in kpc (km/s)^2/Msun; acceleration converted to SI.
        g=4.30091e-6*Ms/(mid+a)**2*1e6/3.085677581491367e19
        cap=vol*g/1e-10/M0
    truth=np.interp(x,row['r_over_R'],row['target_cdf'])
    return m,cap,truth,cdf,M0

results=[];largest_res=0;largest_order=0;worst_conservation=0
for geometry in ['companion_regular','standard_flat_FLRW']:
 for population in ['Chabrier','Salpeter']:
  rows=[r for r in streams if r['geometry']==geometry and r['population']==population]
  for branch,label in [(-1,'distant'),(0,'near_1.1R')]:
   for model in models:
    logs=np.linspace(-3,3,37) if model=='stellar_density' else np.linspace(3,12,37)
    params=np.array([(b,c) for b in bs for c in logs]);pred=[];errors=[];blocked=[];baseline=[]
    for row in rows:
        m,cap,truth,base,M0=setup(row,model,branch,160)
        y,z=sweep(m,cap,params)
        assert np.max(np.abs(y[:,-1]+z-1))<1e-10
        pred.append(y);blocked.append(z);errors.append(np.mean((y-truth)**2,axis=1));baseline.append(np.mean((base-truth)**2))
    errors=np.array(errors);best=int(np.argmin(errors.mean(axis=0)));per=[];loos=[]
    for i,row in enumerate(rows):
        ix=int(np.argmin(np.delete(errors,i,axis=0).mean(axis=0)));loos.append(errors[i,ix])
        m,cap,truth,base,M0=setup(row,model,branch,160)
        checks={};finecurves={}
        for kind,which in [('shared',best),('loo',ix)]:
            pp=params[which];fine_m,fine_cap,fine_truth,_,_=setup(row,model,branch,320)
            fine,free=allocate(fine_m,fine_cap*10**pp[1],np.exp(-pp[0]));rev,rf=allocate(fine_m,fine_cap*10**pp[1],np.exp(-pp[0]),True)
            yc=np.r_[0.,np.cumsum(fine)];yr=np.r_[0.,np.cumsum(rev)]
            conservation=abs(fine.sum()+free-1);worst_conservation=max(worst_conservation,conservation)
            assert conservation<1e-10 and np.all(fine<=fine_cap*10**pp[1]+1e-10) and np.min(fine)>=0
            drift=float(np.max(np.abs(pred[i][which]-yc[::2])));order=float(np.max(np.abs(yc-yr)))
            largest_res=max(largest_res,drift);largest_order=max(largest_order,order)
            unlimited,_=allocate(fine_m,np.full(len(fine_m),1e99),np.exp(-pp[0]))
            cap_effect=float(np.max(np.abs(yc-np.r_[0.,np.cumsum(unlimited)])))
            checks[kind]=dict(capacity_max_cdf_effect=cap_effect,resolution_max_cdf_change=drift,ordering_max_cdf_change=order,refined_blocked_fraction=float(free),reverse_order_blocked_fraction=float(rf),fine_rms=float(np.sqrt(np.mean((yc-fine_truth)**2))))
            finecurves[kind]=yc.tolist()
        per.append(dict(Name=row['Name'],supplied_original_inventory_Msun=M0,baseline_max_error=float(np.max(np.abs(base-truth))),shared_max_error=float(np.max(np.abs(pred[i][best]-truth))),loo_max_error=float(np.max(np.abs(pred[i][ix]-truth))),loo_parameters=params[ix].tolist(),loo_grid_boundary=bool(ix//37 in [0,32] or ix%37 in [0,36]),shared_blocked_fraction=float(blocked[i][best]),loo_blocked_fraction=float(blocked[i][ix]),checks=checks,refined_cdfs=finecurves))
    result=dict(geometry=geometry,population=population,source=label,model=model,shared_parameters=params[best].tolist(),shared_grid_boundary=bool(best//37 in [0,32] or best%37 in [0,36]),baseline_rms=float(np.sqrt(np.mean(baseline))),shared_rms=float(np.sqrt(errors[:,best].mean())),loo_rms=float(np.sqrt(np.mean(loos))),rows=per)
    old=next(v for v in previous['results'] if v['geometry']==geometry and v['population']==population and v['source']==label and v['model']=='partial')
    x=np.linspace(0,1,161)
    result['previous_partial_loo_rms_same_grid']=float(np.sqrt(np.mean([np.mean((np.interp(x,previous['grid'],v['loo_cdf'])-np.interp(x,previous['grid'],v['target_cdf']))**2) for v in old['rows']])))
    results.append(result)
    print(geometry,population,label,model,'baseline/shared/LOO pp',*[round(100*result[k],3) for k in ['baseline_rms','shared_rms','loo_rms']],flush=True)
out=dict(scope='Capacity-limited bound reservoir with traveling overflow; not a full gravity or lens fit',input_sha256={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in files},capacity_units='constant/gravity: log10 Msun/kpc^3; stellar: log10 dimensionless multiplier',parameter_order=['b','log10_capacity_amplitude'],largest_resolution_max_cdf_change=largest_res,largest_ordering_max_cdf_change=largest_order,max_inventory_conservation_error=worst_conservation,results=results)
(P/'saturation-results.json').write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8',newline='\n')
lines=['# Capacity-limited companion settling','',
'Three receiving-capacity rules were tested while keeping the one-third capture law fixed. Incoming packets fill available bound states inward of their original location; blocked packets try farther-out cells within their permitted migration interval. Any remainder stays in an explicitly counted traveling reservoir. No missing packet energy is renormalized into the bound profile.','',
'## Proposed rules and formula provenance','',
'Use mass-equivalent capacity rho_cap=u_cap/c^2. The E=mc^2 conversion and shell bookkeeping are known mathematics. These companion capacity rules are proposed, not derived physics:','',
'| Rule | Bound-state capacity density |','|---|---|',
'| Constant | rho_cap=C |','| Stellar density | rho_cap=zeta*rho_stars |','| Baryonic gravity | rho_cap=C*g_b/(1e-10 m/s^2) |','',
'Stellar density and acceleration are estimated with a known spherical Hernquist form, M_star(<r)=M_star*r^2/(r+a_star)^2 and a_star=Re/1.8153. This is a declared approximate proxy with the original companion-fit stellar mass, not a measured three-dimensional stellar distribution or an NFW-derived driving field. No halo target mass was used to set capacity or supplied inventory. The mass inputs themselves were fitted to prior stellar data, so they are not wholly independent measurements.','',
'For each source packet, the allowed interval is exp(-b)*r to r. The algorithm fills the innermost available cell, continues outward if full, and records unbound overflow. b is a phenomenological migration-rate-times-duration parameter; no universe age is imposed. This is an allocation rule, not a solved time-dependent drift equation or a gravitational force law. Source packets are processed inside-out; reversing order is reported as a sensitivity check.','',
'The supplied inventory is the original deposit mass inside diagnostic 5 Re, redistributed according to each conditional stream profile. R is not a measured halo boundary. Cumulative bound mass is divided by SUPPLIED mass, so its endpoint can be below one. The target is still a normalized halo profile; this does not verify the absolute halo mass, original photon supply or projected lensing. Traveling overflow also gravitates, but its field and stresses are not solved here; these comparisons concern the bound component only.','',
'## Shared fits and transfer','',
'Two parameters per rule were fitted on a declared grid. Each omitted-galaxy result fits the other five and freezes the rule for the omitted target. All targets were previously inspected. Numbers below are RMS cumulative-profile errors in percentage points, not velocity/lens errors or significance.','',
'| Geometry | Source | Capacity | No migration | Shared fit | Omitted-galaxy transfer |','|---|---|---|---:|---:|---:|']
for r in results:
 if r['population']=='Chabrier':lines.append(f"| {r['geometry']} | {r['source']} | {r['model']} | {100*r['baseline_rms']:.2f} | {100*r['shared_rms']:.2f} | {100*r['loo_rms']:.2f} |")
lines+=['','All population cases, parameter boundaries, per-target errors, original supplied inventories, blocked fractions, refined curves and sensitivity checks are in saturation-results.json. Grid bounds were not expanded after seeing results. A lower bound-component error achieved by leaving energy unbound is not by itself a successful total-gravity prediction.','',
'## Numerical and physical limits','',f"Packet accounting closes to {worst_conservation:.3g}; no bound cell exceeds its capacity. Doubling radial cells changes a selected cumulative curve by at most {largest_res:.4g}, and reversing packet order changes it by at most {largest_order:.4g}. These differences quantify unresolved numerical or assembly-order dependence rather than disappearing into a fit. The shared/omitted scores in the table use the declared 160-cell grid; refined scores are archived.",'',
'The bound-plus-traveling packet energy is conserved. Gravitational work released during settling, pressure/support, momentum exchange, self-gravity and time-dependent supply still require a physical completion. Capacity does not switch gravity off: filled bound cells and traveling energy remain gravitational sources. No research goal, new observational lens fit, or derivation of saturation is claimed complete.']
findings=['## Outcome','',
'None of these capacity prescriptions beats the previous partial-migration model on the retained-geometry omitted-galaxy comparison. Distant-input errors are about 24.64, 28.08 and 25.46 points for constant, stellar and gravitational capacities. The prior partial-migration value is about 20 points. Grid-matched comparisons to those archived predictions are included in the JSON; small changes from earlier numbers reflect radial sampling.','',
'The best distant-input shared fits choose high enough capacities to behave essentially like unsaturated contraction. Thus their improvement over no migration is not evidence for saturation. The near-source gravitational-capacity branch does actively redistribute packets and improves its own baseline to 23.62 points, but still does not outperform the previous partial model. These are already-inspected data and the capacity amplitudes are phenomenological.','',
'Some stellar-capacity omitted fits fail to bind about 28% of one target inventory. That energy remains in the traveling account; its gravitational field cannot be ignored in a physical model. Source ordering can change a selected bound cumulative profile by as much as 43.8 percentage points, much larger than the radial-grid effect (under 0.8 points). This is a serious assembly-history dependence of the allocation prescription, not numerical proof of a unique stable halo.','',
'The next physical requirement is a time-dependent transport and source-history rule allowing occupied states, incoming packets and traveling overflow to interact consistently. The present test gives no basis to claim that a universal saturation threshold has solved the halos, or to modify the one-third capture exponent.','']
pos=lines.index('## Numerical and physical limits');lines[pos:pos]=findings
(P/'saturation-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
selected=[r for r in results if r['geometry']=='companion_regular' and r['population']=='Chabrier' and r['source']=='distant']
for i,ax in enumerate(axes.flat):
 name=selected[0]['rows'][i]['Name'];row=next(r for r in streams if r['Name']==name and r['population']=='Chabrier' and r['geometry']=='companion_regular')
 x=np.linspace(0,1,321)
 ax.plot(x,np.interp(x,row['r_over_R'],row['target_cdf']),'k',label='Fitted halo target')
 ax.plot(x,np.interp(x,row['r_over_R'],row['case_cdfs'][-1]),color='gray',label='No migration')
 for result in selected:ax.plot(x,result['rows'][i]['refined_cdfs']['loo'],label=result['model'])
 ax.set(title=name,xlabel='Radius / diagnostic region',ylabel='Bound mass / supplied mass',xlim=(0,1),ylim=(0,1));ax.grid(alpha=.2)
axes.flat[0].legend(fontsize=7)
fig.suptitle('Saturable settling: omitted-galaxy predictions, retained geometry, distant input')
fig.savefig(P/'saturation.png',dpi=150)
