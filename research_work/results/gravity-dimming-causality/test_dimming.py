#!/usr/bin/env python3
"""JR-3: hidden-light diagnostics and an energy-conserving reversible exchange.

Run against the unchanged JR-2 complete package:
  python test_dimming.py --jr2 /path/to/Photon-Graviton-JR2 --output /fresh/path
No input files are changed. Fitted per-object transmissions are counterfactual,
not independent measurements, a universal causal law, or proof of conversion.
"""
from __future__ import annotations
import argparse, hashlib, io, json, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, solve_triangular
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import spearmanr

LOG10=np.log(10.)
MSUN=1.98847e30
LSUN=3.828e26
C=299792458.
YEAR=365.25*86400

def dump(p,x): p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def min_grid(fun,lo,hi,n=81):
    grid=np.linspace(lo,hi,n); val=np.array([fun(x) for x in grid]); choices=[(float(val[0]),lo),(float(val[-1]),hi)]
    for i in range(1,n-1):
        if val[i]<=val[i-1] and val[i]<=val[i+1]:
            res=minimize_scalar(fun,bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-10})
            choices.append((float(res.fun),float(res.x)))
    v,x=min(choices)
    return float(x),float(v)

def brief(rows):
    return dict(galaxies=len(rows),mean_galaxy_RMSE_kms=float(np.mean([r['RMSE_kms'] for r in rows])),
        median_fractional_RMS=float(np.median([r['fractional_RMS'] for r in rows])),
        mean_fractional_RMS=float(np.mean([r['fractional_RMS'] for r in rows])),
        raw_chi2=float(sum(r['chi2'] for r in rows)),
        within10=sum(r['fractional_RMS']<.1 for r in rows),within20=sum(r['fractional_RMS']<.2 for r in rows),
        positive_dimming=sum(r.get('net_dimming_fraction',0)>1e-5 for r in rows),
        median_net_dimming_fraction=float(np.median([r.get('net_dimming_fraction',0) for r in rows])))

def reversible():
    """Analytic, independent linear solve, exact exponential and ODE ledger."""
    fixtures=[]; max_stat=max_eq=max_balance=0.
    for kp in (.1,1.,10.):
      for km in (.1,1.,10.):
       for ke in (.1,1.,10.):
        for kc in (0.,.1,1.):
            M=np.array([[-ke-kp,km],[kp,-km-kc]])
            steady=np.linalg.solve(M,[-1.,0.])
            eg=1/(ke+kp*kc/(km+kc)); ec=kp*eg/(km+kc)
            err=float(np.max(abs(steady-[eg,ec]))/max(eg,ec)); max_stat=max(max_stat,err)
            balance=abs(ke*eg+kc*ec-1);max_balance=max(max_balance,balance)
            if kc==0:max_eq=max(max_eq,abs(ke*eg-1))
            fixtures.append(dict(k_plus=kp,k_minus=km,k_photon_escape=ke,k_companion_escape=kc,
                E_gamma=eg,E_chi=ec,L_out_over_L_in=ke*eg,companion_escape_power=kc*ec,
                maximum_linear_eigenvalue=float(np.max(np.linalg.eigvalsh((M+M.T)/2))) if False else float(np.max(np.real(np.linalg.eigvals(M)))),
                analytic_vs_linear_relative_error=err))
    evolved=[]; worst=0.
    for kp,km,ke,kc in [(1.,.2,1.,0.),(5.,.1,1.,0.),(1.,.2,1.,.2),(5.,.1,1.,.2)]:
        M=np.array([[-ke-kp,km],[kp,-km-kc]])
        timescale=1/abs(max(np.real(np.linalg.eigvals(M))));end=12*timescale
        # Eg, Ec, total emitted photons, total escaped companions, constant 1
        A=np.zeros((5,5));A[:2,:2]=M;A[0,4]=1;A[2,0]=ke;A[3,1]=kc
        y0=np.array([0.,0.,0.,0.,1.]);tt=np.linspace(0,end,81)
        exact=np.array([expm(A*t)@y0 for t in tt])
        sol=solve_ivp(lambda t,y:A@y,(0,end),y0,t_eval=tt,method='DOP853',rtol=2e-11,atol=2e-12)
        err=float(np.max(abs(sol.y.T-exact))/max(1,float(abs(exact).max())))
        closure=float(np.max(abs(sol.y[:4].sum(0)-tt))/max(1,end));worst=max(worst,err,closure)
        evolved.append(dict(k_plus=kp,k_minus=km,k_photon_escape=ke,k_companion_escape=kc,
            times=tt.tolist(),states=sol.y.T[:,:4].tolist(),photon_luminosity=(ke*sol.y[0]).tolist(),
            relative_exponential_ODE_error=err,relative_energy_ledger_error=closure,
            dimensionless_time=True,ode_success=bool(sol.success)))
    return dict(fixtures=fixtures,evolutions=evolved,
        checks=dict(stationary_linear_error=max_stat,stationary_energy_error=max_balance,
        zero_companion_escape_stationary_dimming=max_eq,evolution_max_error=worst),
        scope='Rate-limit finite energy system. L_out is all-direction bolometric escape, not one band/aperture.',
        identities=['At steady state L_in-L_out=k_chi_escape*E_chi.',
        'If k_chi_escape=0 and k_minus>0, steady L_out=L_in even with nonzero stored E_chi.',
        'Transient storage can dim. Spectral/directional redistribution can dim one measured channel.',
        'Energy cycling does not create additional net stored energy on each cycle.'])

def main():
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--jr2',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();out=args.output
    if out.exists():raise FileExistsError('Use a fresh output directory')
    out.mkdir(parents=True);start=time.monotonic()
    base=args.jr2;bundle=base/'original-JR1';repo=bundle/'repository'
    sys.path.insert(0,str(base/'code'));sys.path.insert(0,str(repo/'research_work/results/joint-response-iteration'))
    import inverse_properties as I
    import population_followup as P
    J=P.J
    frozen=json.loads((bundle/'evidence/R10_spheroid_stellar_population.json').read_text());p=frozen['parameters'];spec=frozen['specification']
    files=[q for q in repo.rglob('*') if q.is_file() and q.suffix!='.pyc']+[bundle/'evidence/R10_spheroid_stellar_population.json',base/'code/inverse_properties.py']
    before={str(q.relative_to(base)):hashlib.sha256(q.read_bytes()).hexdigest() for q in files}
    E=P.Experiment();cat=I.catalog(repo/J.INPUT_NAMES[1]);disks=[];replay=[]
    with zipfile.ZipFile(repo/J.INPUT_NAMES[0]) as z:
      for g in E.galaxies:
        raw=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(g['name']+'_rotmod.dat'))));raw=raw[raw[:,0]>0]
        disks.append(I.Disk(g,cat[g['name']],raw,J,P,p))
    def pred(d,lu,both=False):
        x=d.base.copy();x[2]=lu/LOG10
        if both:x[3]=lu/LOG10
        return d.forward(x)[0]
    def met(d,lu,both=False):
        y=pred(d,lu,both);r=y-d.g['observed'];f=r/d.g['observed']
        return dict(log_intrinsic_light_correction=float(lu),intrinsic_light_factor=float(np.exp(lu)),
            net_dimming_fraction=float(-np.expm1(-lu)),fractional_RMS=float(np.sqrt(np.mean(f*f))),
            RMSE_kms=float(np.sqrt(np.mean(r*r))),chi2=float(np.sum((r/d.g['error'])**2)),
            signed_fractional_mean=float(np.mean(f)),predicted_kms=y.tolist())
    rows=[]
    for i,d in enumerate(disks):
        base_m=met(d,0);old=next(r for r in frozen['metrics']['sparc_rows'] if r['name']==d.g['name'])
        replay.append(float(np.max(abs(np.array(base_m['predicted_kms'])-old['predicted_kms']))))
        cases=dict(baseline=base_m)
        for both in (False,True):
            tag='stars_and_gas' if both else 'stars'
            for f in (.25,.5):cases[tag+f'_dim_{int(f*100)}']=met(d,-np.log1p(-f),both)
            for obj in ('fractional','chi2'):
                def loss(lu):
                    r=(pred(d,lu,both)-d.g['observed'])/(d.g['observed'] if obj=='fractional' else d.g['error'])
                    return float(np.mean(r*r))
                ld,vd=min_grid(loss,0,LOG10)
                ls,vs=min_grid(loss,-LOG10,LOG10,n=161)
                for name,lu in [('dim_only',ld),('signed_exchange',ls)]:
                    m=met(d,lu,both);m['at_bound']=bool(min(abs(lu),abs(lu-LOG10))<1e-5) if name=='dim_only' else bool(abs(abs(lu)-LOG10)<1e-5)
                    cases[tag+'_'+name+'_'+obj]=m
        # Column of independent predicted velocities under stellar corrections.
        lower=pred(d,0);negative=0;max_drop=0.
        for lu in np.linspace(0,LOG10,101)[1:]:
            yy=pred(d,lu);diff=yy-lower;negative+=int((diff<-1e-9).sum());max_drop=min(max_drop,float(diff.min()));lower=yy
        g=d.g;A,rc,rt=P.source(g['Mstar'],g['Mgas'],g['Re'],p,False,g['Mbulge']/g['Mstar'])
        fs=10**p['logusph']*g['Mbulge']/(g['Mdisk']+10**p['logusph']*g['Mbulge'])
        q=p['q']+(p['q_sph']-p['q'])*fs;last=float(g['r'][-1]);mt=float(A*rt/J.G)
        ml=float(last**2*J.companion_force(last,A,rc,rt,q)/J.G)
        lc=float(d.c['L']*1e9)
        # B is an unknown conversion between catalog band solar units and bolometric solar units.
        tbase=mt*MSUN*C*C/(lc*LSUN)/YEAR; tl=ml*MSUN*C*C/(lc*LSUN)/YEAR
        ul=np.exp(cases['stars_dim_only_fractional']['log_intrinsic_light_correction'])
        budget=dict(Mchi_total_effective_Msun=mt,Mchi_within_last_radius_Msun=ml,
            last_radius_kpc=last,L36_catalog_solar_units=lc,
            normalized_energy_over_observed_luminosity_years=tbase,normalized_inner_energy_over_luminosity_years=tl,
            required_B_eta_T_years_at_best_dimming=None if ul-1<1e-7 else float(tbase/(ul-1)),
            required_B_eta_T_years_at_25pct_dimming=float(3*tbase),
            required_B_eta_T_years_at_50pct_dimming=float(tbase),
            scope='Required pre-existing R10 reservoir only; not extra gravity added on top. B, capture, history not measured.')
        rows.append(dict(name=g['name'],split=g['split'],catalog=d.c,original_outlier=base_m['fractional_RMS']>.2,
            cases=cases,stellar_correction_velocity_scan_negative_steps=negative,
            worst_velocity_scan_increment_kms=max_drop,energy=budget))
        if (i+1)%25==0:print('DISKS',i+1,flush=True)
    assert max(replay)<1e-8
    dump(out/'disk-results.json',rows)
    print('DISK SUMMARY',json.dumps({k:brief([r['cases'][k] for r in rows]) for k in rows[0]['cases']}),flush=True)
    # Self-consistent shared gravity-dependent band attenuation.
    ms=E.Ms;mg=E.Mg;re=E.Re;fs=E.Msph/E.Ms;u_s=10**p['logusph'];new_fs=u_s*fs/(1+(u_s-1)*fs)
    qs=p['q']+(p['q_sph']-p['q'])*new_fs
    vstar=np.array([np.interp(g['Re'],g['r'],10**p['logu']*(d.disk+u_s*d.bul)) for g,d in zip(E.galaxies,disks)])
    vgas=np.array([np.interp(g['Re'],g['r'],g['gas_v2']) for g in E.galaxies])
    def h(lu,index=None):
        # scalar or broadcasted arrays accepted for grid self-consistency search
        ids=slice(None) if index is None else index
        aa,cc,tt=P.source(ms[ids]*np.exp(lu),mg[ids],re[ids],p,False,fs[ids])
        chi=re[ids]*J.companion_force(re[ids],aa,cc,tt,qs[ids])
        bary=np.maximum(np.exp(lu)*vstar[ids]+vgas[ids],0)
        return chi/(chi+bary)
    h0=h(np.zeros(len(ms)));maxroots=0
    def solve_lu(tau,kind):
        nonlocal maxroots
        if tau==0:return np.zeros(len(ms))
        if kind=='constant':return np.full(len(ms),tau)
        grid=np.linspace(0,tau,66);f=grid[:,None]-tau*h(grid[:,None]);lus=[]
        for k in range(len(ms)):
            ii=np.flatnonzero((f[:-1,k]<=0)&(f[1:,k]>=0));maxroots=max(maxroots,len(ii))
            if not len(ii):raise ValueError('Missing opacity root')
            j=ii[0];lus.append(brentq(lambda t:t-tau*float(h(t,k)),grid[j],grid[j+1],xtol=1e-12))
        return np.array(lus)
    shared=[]
    for kind in ('constant','gravity_linked'):
        def loss(tau):
            lu=solve_lu(tau,kind)
            return float(np.mean([np.mean(((pred(d,lu[i])-d.g['observed'])/d.g['observed'])**2) for i,d in enumerate(disks) if d.g['split']=='train']))
        tau,cost=min_grid(loss,0,LOG10,n=41);lu=solve_lu(tau,kind);rr=[met(d,lu[i]) for i,d in enumerate(disks)]
        shared.append(dict(kind=kind,tau0=tau,training_mean_fractional_square=cost,
            maximum_self_consistency_error=float(np.max(abs(lu-tau*h(lu)))) if kind=='gravity_linked' else 0,
            splits={s:brief([r for d,r in zip(disks,rr) if d.g['split']==s]) for s in ('train','validation','test')},
            all=brief(rr),per_galaxy=[dict(name=d.g['name'],**r) for d,r in zip(disks,rr)]))
        print('SHARED',kind,tau,shared[-1]['all'],flush=True)
    dump(out/'shared-laws.json',shared)
    # Lensing: infer missing band light at fixed population shape and mass-to-light.
    lenses=[]
    for l in E.lenses:
        local=l.local_parameters(p,spec);aspec=dict(spec,stellar_nuisance=False)
        p0=p.copy();p0['logu']=local['logu'];p0['logusph']=0.;p0['beta']=local['beta']
        def lm(lu,beta=None,angle=True):
            pp=p0.copy();pp['logu']+=lu/LOG10
            if beta is not None:pp['beta']=beta
            v,f,_=l.prediction(pp,aspec);w=solve_triangular(l.chol,v-l.y,lower=True)
            return dict(intrinsic_light_factor=float(np.exp(lu)),net_dimming_fraction=float(-np.expm1(-lu)),beta=float(pp['beta']),
                stellar_fractional_RMS=float(np.sqrt(np.mean(((v-l.y)/l.y)**2))),stellar_chi2=float(w@w),
                deflection_error_at_observed=f,theta_observed=l.theta,
                theta_predicted=float(l.angle(pp,aspec)) if angle else None,predicted_vrms=v.tolist())
        baseline=lm(0);root=brentq(lambda lu:lm(lu,angle=False)['deflection_error_at_observed'],-LOG10,LOG10)
        beta,val=min_grid(lambda b:lm(root,b,False)['stellar_chi2'],-1,.35)
        star_lu,_=min_grid(lambda lu:lm(lu,angle=False)['stellar_chi2'],0,LOG10)
        joint_lu,_=min_grid(lambda lu:lm(lu,angle=False)['stellar_chi2']/len(l.y)+(lm(lu,angle=False)['deflection_error_at_observed']/.05)**2,0,LOG10)
        lenses.append(dict(name=l.name,cases=dict(baseline=baseline,ring_required=lm(root),ring_and_orbit=lm(root,beta),
            dim_for_stars=lm(star_lu),dim_joint=lm(joint_lu)),ring_requires_brightening=bool(root<0)))
    dump(out/'lens-results.json',lenses);print('LENSES DONE',flush=True)
    box=reversible();dump(out/'reversible-energy.json',box)
    cases=rows[0]['cases'].keys();outliers=[r for r in rows if r['original_outlier']]
    subgroups=dict(all=rows,original_outliers=outliers,
        overpredicted_outliers=[r for r in outliers if r['cases']['baseline']['signed_fractional_mean']>0],
        underpredicted_outliers=[r for r in outliers if r['cases']['baseline']['signed_fractional_mean']<=0])
    summaries={g:{k:brief([r['cases'][k] for r in rr]) for k in cases} for g,rr in subgroups.items()}
    counts={}
    for k in ('stars_dim_only_fractional','stars_signed_exchange_fractional','stars_and_gas_dim_only_fractional'):
        vals=[r['cases'][k] for r in rows]
        counts[k]=dict(no_exchange=sum(abs(v['net_dimming_fraction'])<1e-5 for v in vals),
            net_dimming=sum(v['net_dimming_fraction']>1e-5 for v in vals),
            net_brightening=sum(v['net_dimming_fraction']< -1e-5 for v in vals),
            upper_90pct_dimming=sum(v['net_dimming_fraction']>.9-1e-5 for v in vals),
            original_outliers_rescued=sum(r['original_outlier'] and r['cases'][k]['fractional_RMS']<.2 for r in rows))
    weights=[r['cases']['stars_signed_exchange_fractional']['net_dimming_fraction'] for r in rows]
    corr=spearmanr(weights,h0)
    energy={}
    for k in ('normalized_energy_over_observed_luminosity_years','normalized_inner_energy_over_luminosity_years'):
        v=[r['energy'][k] for r in rows];energy[k]=dict(min=float(np.min(v)),p10=float(np.percentile(v,10)),median=float(np.median(v)),p90=float(np.percentile(v,90)),max=float(np.max(v)))
    after={str(q.relative_to(base)):hashlib.sha256(q.read_bytes()).hexdigest() for q in files}
    assert before==after
    summary=dict(experiment='JR-3',protocol_commit='2841a36a9d443cb4ca6d513056e3def5b9d1ca08',baseline_commit='41ede5342a0ddabdff7dbea2a4f1951eceec65c6',
        scope='Frozen-R10 band-dimming diagnostic and separate reversible energy model, not a derived full conversion theory.',
        original_galaxies=149,original_rotation_rows=sum(len(d.g['r']) for d in disks),original_lenses=6,
        R10_velocity_replay_max_kms=max(replay),original_inputs_unchanged=True,
        disk_summaries=summaries,exchange_direction_counts=counts,
        shared_laws=[{k:v for k,v in row.items() if k!='per_galaxy'} for row in shared],
        shared_selfconsistent_max_roots=maxroots,stellar_scan_galaxies_with_decreases=sum(r['stellar_correction_velocity_scan_negative_steps']>0 for r in rows),
        signed_exchange_vs_predicted_companion_fraction=dict(spearman=float(corr.statistic),p_descriptive_only=float(corr.pvalue),warning='both inferred from same exposed data/model; not independent causal evidence'),
        energy_summary=energy,reversible_checks=box['checks'],lens_rows=lenses,
        elapsed_seconds=time.monotonic()-start,versions=dict(python=sys.version,numpy=np.__version__),
        limits=['Uniform photometric transmission at unchanged intrinsic population, not spectra or stellar-formation suppression.',
        'Independent opacity fits are not a causal universal law.',
        'Net gain requires a supplied companion population or intrinsic population change; signed fits do not derive it.',
        'Gas is inferred from line radiation; dimming a line does not destroy gas atoms.',
        'Stationary photon-escape result concerns bolometric all-direction energy, not a narrow band or aperture.',
        'Energy budget assumes effective mass gravitates as E/c^2; a different response needs an explicit new relation.',
        'No age, independent bolometric correction, source history, actual transmission or cluster model measured.'])
    dump(out/'summary.json',summary)
    dump(out/'input-manifest.json',dict(sha256=before,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print('FINAL',json.dumps({k:v for k,v in summary.items() if k not in ('disk_summaries','lens_rows')}),flush=True)

if __name__=='__main__':main()
