#!/usr/bin/env python3
"""JR-2 fixed R10 galaxy-property inversion. Preserves all input files.
Run: python inverse_properties.py --bundle /path/to/Photon-Graviton-JR1 --output /fresh/output
The output stores counterfactual properties, never corrected observations.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time, zipfile, io
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize_scalar, brentq
from scipy.linalg import solve_triangular
from scipy.stats import spearmanr


def dump(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def catalog(path):
    result={}
    for line in path.read_text().splitlines():
        f=line.split()
        if len(f)!=19: continue
        try:
            result[f[0]]=dict(type=int(f[1]),D=float(f[2]),eD=float(f[3]),method=int(f[4]),
                inc=float(f[5]),einc=float(f[6]),L=float(f[7]),eL=float(f[8]),Re=float(f[9]),
                SBe=float(f[10]),rd=float(f[11]),SBdisk=float(f[12]),MHI=float(f[13]),Q=int(f[17]))
        except ValueError: continue
    return result

class Disk:
    def __init__(self,g,c,raw,J,P,p):
        self.g,self.c,self.J,self.P,self.p=g,c,J,P,p
        self.bul=.7*raw[:,5]**2
        self.disk=.5*raw[:,4]**2
        self.proj=g['observed']*np.sin(np.deg2rad(c['inc']))
        self.errproj=g['error']*np.sin(np.deg2rad(c['inc']))
        self.base=np.array([0.,c['inc'],0.,0.])
        if min(c['eD'],c['einc'])<=0:raise ValueError('Missing property uncertainty')
    def forward(self,x):
        d,i,us,ug=10**x[0],x[1],10**x[2],10**x[3]
        g,p,J=self.g,self.p,self.J
        f=g['Mbulge']/g['Mstar']; u=10**p['logusph']
        A,rc,rt=self.P.source(g['Mstar']*d*d*us,g['Mgas']*d*d*ug,g['Re']*d,p,False,f)
        q=p['q']+(p['q_sph']-p['q'])*u*f/(1+(u-1)*f)
        rad=g['r']*d
        vb2=d*(ug*g['gas_v2']+10**p['logu']*us*(self.disk+u*self.bul))
        v=np.sqrt(np.maximum(vb2,0.)+rad*J.companion_force(rad,A,rc,rt,q))
        return v,v*np.sin(np.deg2rad(i))
    def residual(self,x,prior=False):
        v,proj=self.forward(x)
        r=(proj-self.proj)/self.errproj
        if prior:
            c=self.c
            r=np.r_[r,(10**x[0]-1)*c['D']/c['eD'],(x[1]-c['inc'])/c['einc'],x[2]/.15,x[3]/.10]
        return r
    def metrics(self,x):
        v,proj=self.forward(x);e=proj-self.proj;frac=e/self.proj
        c=self.c;g=self.g
        return dict(distance_factor=float(10**x[0]),distance_Mpc=float(c['D']*10**x[0]),
            distance_sigma=float((10**x[0]-1)*c['D']/c['eD']),inclination_deg=float(x[1]),
            inclination_shift_deg=float(x[1]-c['inc']),inclination_sigma=float((x[1]-c['inc'])/c['einc']),
            stellar_ML_factor=float(10**x[2]),gas_factor=float(10**x[3]),
            stellar_total_mass_factor=float(10**(2*x[0]+x[2])),gas_total_mass_factor=float(10**(2*x[0]+x[3])),
            fractional_RMS=float(np.sqrt(np.mean(frac*frac))),chi2=float(np.sum((e/self.errproj)**2)),
            reduced_chi2_N=float(np.mean((e/self.errproj)**2)),
            RMSE_catalog_equivalent_kms=float(np.sqrt(np.mean(e*e))/np.sin(np.deg2rad(c['inc']))),
            mean_signed_fractional_error=float(np.mean(frac)),predicted_projected_kms=proj.tolist(),
            predicted_intrinsic_kms=v.tolist())
    def fit(self,free,lo,hi,prior=False,starts=3):
        rng=np.random.default_rng(218+sum(map(ord,self.g['name']))+sum(free))
        lo=np.array(lo);hi=np.array(hi);x0=np.clip(self.base[free],lo+1e-8,hi-1e-8)
        def res(z):
            x=self.base.copy();x[free]=z;return self.residual(x,prior)
        zz=[x0]
        if len(free)>1:
            zz += [np.clip(x0+rng.normal(0,.12,len(free))*(hi-lo),lo+1e-7,hi-1e-7) for _ in range(starts-1)]
        fits=[least_squares(res,x,bounds=(lo,hi),max_nfev=350,ftol=1e-10,xtol=1e-10,gtol=1e-8,x_scale='jac') for x in zz]
        best=min(fits,key=lambda f:float(f.fun@f.fun));x=self.base.copy();x[free]=best.x
        out=self.metrics(x);out.update(objective=float(best.fun@best.fun),success=bool(best.success),
            starts=len(zz),bound_indices=[int(free[k]) for k in range(len(free)) if min(best.x[k]-lo[k],hi[k]-best.x[k])<1e-5],
            x=x.tolist(),nfev=int(best.nfev))
        return out


def disk_run(E,J,P,p,saved,out):
    cat=catalog(J.ROOT/J.INPUT_NAMES[1]);rows=[];ds=[];replay=[]
    with zipfile.ZipFile(J.ROOT/J.INPUT_NAMES[0]) as archive:
        for k,g in enumerate(E.galaxies):
            raw=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(g['name']+'_rotmod.dat'))));raw=raw[raw[:,0]>0]
            d=Disk(g,cat[g['name']],raw,J,P,p);ds.append(d);baseline=d.metrics(d.base)
            expected=next(t for t in saved['metrics']['sparc_rows'] if t['name']==g['name'])
            replay.append(float(np.max(np.abs(d.forward(d.base)[0]-expected['predicted_kms']))))
            c=d.c; flo=max(.1,1-2*c['eD']/c['D']);fhi=1+2*c['eD']/c['D']
            cases=dict(baseline=baseline,
                distance_only=d.fit([0],[np.log10(.2)],[np.log10(5)],starts=1),
                inclination_only=d.fit([1],[5.],[90.],starts=1),
                stellar_only=d.fit([2],[-1.],[1.],starts=1),
                gas_only=d.fit([3],[np.log10(.25)],[np.log10(4.)],starts=1),
                restricted_joint=d.fit([0,1,2,3],[np.log10(flo),max(5.,c['inc']-2*c['einc']),-.30,-.20],
                    [np.log10(fhi),min(90.,c['inc']+2*c['einc']),.30,.20],prior=True),
                wide_joint=d.fit([0,1,2,3],[np.log10(.2),5.,-1.,np.log10(.25)],
                    [np.log10(5.),90.,1.,np.log10(4.)],starts=4))
            row=dict(name=g['name'],split=g['split'],catalog=c,points=len(g['r']),
                input_Mstar=g['Mstar'],input_Mgas=g['Mgas'],bulge_fraction=g['Mbulge']/g['Mstar'],
                original_outlier=baseline['fractional_RMS']>.20,projected_velocity_proxy_kms=d.proj.tolist(),
                projected_velocity_error_proxy_kms=d.errproj.tolist(),catalog_radius_kpc=g['r'].tolist(),cases=cases)
            rows.append(row)
            if (k+1)%20==0:print('DISKS',k+1,flush=True)
    assert max(replay)<1e-8, max(replay)
    dump(out/'disk-results.json',rows)
    return ds,rows,dict(maximum_R10_velocity_replay_difference_kms=max(replay))


def lens_run(E,J,P,p,spec,saved,out):
    rows=[]
    for l in E.lenses:
        # Freeze all force coefficients. Offsets below describe galaxy properties.
        pp=l.local_parameters(p,spec);base_logu=pp['logu'];beta0=pp['beta']
        aspec=dict(spec,stellar_nuisance=False);p0=p.copy();p0['logu']=base_logu;p0['logusph']=0.;p0['beta']=beta0
        # The Lens subclass adds logusph=0; no previous object offset is applied twice.
        def predict(z):
            pd=p0.copy();pd['logu']=base_logu+z[0];pd['beta']=z[1]
            v,f,pars=l.prediction(pd,aspec)
            return v,f,pd
        v0,f0,_=predict([0,beta0]);old=next(a for a in saved['metrics']['lens_rows'] if a['name']==l.name)
        assert np.max(np.abs(v0-old['predicted_vrms_kms']))<1e-7
        def met(z):
            v,f,pd=predict(z);w=solve_triangular(l.chol,v-l.y,lower=True)
            kappa=-f;geom=1/(1+f);sigma_crit=J.CLIGHT**2/(4*np.pi*J.G*l.Dl*l.ratio)
            return dict(log_stellar_mass_shift_dex=float(z[0]),stellar_mass_factor=float(10**z[0]),
                mass_sigma_vs_published=float(z[0]/l.mass_log_error),stellar_mass_Msun=float(l.Mstar*10**(base_logu+z[0])),
                beta=float(z[1]),beta_shift=float(z[1]-beta0),
                stellar_fractional_RMS=float(np.sqrt(np.mean(((v-l.y)/l.y)**2))),
                stellar_chi2=float(w@w),stellar_RMSE_kms=float(np.sqrt(np.mean((v-l.y)**2))),
                theta_predicted_arcsec=float(l.angle(pd,aspec)),theta_observed_arcsec=l.theta,
                required_external_convergence=float(kappa),required_geometry_ratio_multiplier=float(geom),
                required_Dls_over_Ds=float(l.ratio*geom),required_geometry_ratio_physical=bool(l.ratio*geom<1),
                sheet_equivalent_projected_mass_inside_ring_Msun=float(kappa*sigma_crit*np.pi*l.b**2),
                predicted_vrms_kms=v.tolist())
        def starmass(z):return solve_triangular(l.chol,predict([z[0],beta0])[0]-l.y,lower=True)
        fstar=least_squares(starmass,[0],bounds=(-.30,.30),max_nfev=100)
        flo,fhi=predict([-.30,beta0])[1],predict([.30,beta0])[1]
        if flo*fhi<0:
            zm=brentq(lambda dm:predict([dm,beta0])[1],-.30,.30);ring=met([zm,beta0]);ring['ring_root_found']=True
        else:
            fm=minimize_scalar(lambda dm:predict([dm,beta0])[1]**2,bounds=(-.30,.30),method='bounded');ring=met([fm.x,beta0]);ring['ring_root_found']=False
        def fit(mode):
            def res(z):
                v,f,pd=predict(z);r=solve_triangular(l.chol,v-l.y,lower=True)
                if mode=='joint':r=np.r_[r,f/.05*np.sqrt(len(l.y))]
                return r
            starts=[[0,beta0],[.1,-.2],[-.1,.2]];best=min((least_squares(res,x,bounds=([-.30,-1.],[.30,.35]),max_nfev=250,
                ftol=1e-11,xtol=1e-11,gtol=1e-9) for x in starts),key=lambda f:f.fun@f.fun)
            m=met(best.x);m.update(success=bool(best.success),bound=bool(abs(best.x[0])>.29999 or best.x[1]<-.99999 or best.x[1]>.34999))
            return m
        rows.append(dict(name=l.name,previous_role='fit' if l.training else 'excluded_optimizer',
            observed_vrms_kms=l.y.tolist(),covariance=l.cov.tolist(),original_beta=beta0,
            original_mass_log_error_dex=l.mass_log_error,cases=dict(baseline=met([0,beta0]),
            mass_for_stars=met([fstar.x[0],beta0]),mass_for_ring=ring,
            mass_and_orbits_joint=fit('joint'),mass_and_orbits_stars_only=fit('stars'))))
        print('LENS',l.name,json.dumps({k:{x:v[x] for x in ['stellar_mass_factor','beta','stellar_fractional_RMS','theta_predicted_arcsec','required_external_convergence']} for k,v in rows[-1]['cases'].items()}),flush=True)
    dump(out/'lens-results.json',rows)
    return rows


def common_fit(ds,rows,out):
    # Common fractional distance, common inclination offset, common M/L and gas scale.
    # Equal-object standardized-square objective, no formula coefficient changes.
    lo=[np.log10(.5),-10.,-.3,-.2];hi=[np.log10(2.),10.,.3,.2]
    def residual(z,objects):
        values=[]
        for d in objects:
            x=np.array([z[0],np.clip(d.c['inc']+z[1],5.,90.),z[2],z[3]])
            values.append(d.residual(x)/np.sqrt(len(d.proj)))
        return np.concatenate(values)
    starts=[np.zeros(4),[.05,2,-.1,0],[-.1,-2,.1,0]]
    best=min((least_squares(residual,x,args=(ds,),bounds=(lo,hi),max_nfev=160) for x in starts),key=lambda f:f.fun@f.fun)
    metrics=[]
    for d in ds:
        x=np.array([best.x[0],np.clip(d.c['inc']+best.x[1],5.,90.),best.x[2],best.x[3]])
        metrics.append(dict(name=d.g['name'],**d.metrics(x)))
    result=dict(common_parameters=best.x.tolist(),distance_factor=float(10**best.x[0]),inclination_offset_deg=float(best.x[1]),
        stellar_ML_factor=float(10**best.x[2]),gas_factor=float(10**best.x[3]),objective=float(best.fun@best.fun),
        scope='In-sample common property-calibration diagnostic, not a shared correction validated on independent objects',rows=metrics)
    dump(out/'shared-calibration.json',result)
    return result


def summarize(disks,lenses,shared,replay):
    def agg(group,case):
        vals=[r['cases'][case] for r in group];return dict(galaxies=len(vals),
            mean_RMSE_catalog_equivalent_kms=float(np.mean([v['RMSE_catalog_equivalent_kms'] for v in vals])),
            median_fractional_RMS=float(np.median([v['fractional_RMS'] for v in vals])),
            mean_fractional_RMS=float(np.mean([v['fractional_RMS'] for v in vals])),
            raw_chi2=float(sum(v['chi2'] for v in vals)),within10pct=sum(v['fractional_RMS']<.10 for v in vals),
            within20pct=sum(v['fractional_RMS']<.20 for v in vals),mean_chi2_N_le2=sum(v['reduced_chi2_N']<=2 for v in vals),
            bound_hits=sum(bool(v.get('bound_indices')) for v in vals),
            median_distance_factor=float(np.median([v['distance_factor'] for v in vals])),
            median_inclination_shift_deg=float(np.median([v['inclination_shift_deg'] for v in vals])),
            median_stellar_ML_factor=float(np.median([v['stellar_ML_factor'] for v in vals])),
            median_gas_factor=float(np.median([v['gas_factor'] for v in vals])))
    cases=list(disks[0]['cases']);subgroups={'all':disks,'original_outliers':[r for r in disks if r['original_outlier']],
        'original_not_outliers':[r for r in disks if not r['original_outlier']]}
    summary={g:{c:agg(r,c) for c in cases} for g,r in subgroups.items()}
    # Descriptive direction statistics; condition on residual sign, not a causal test.
    directions={}
    for label,rows in [('initial_overprediction',[r for r in disks if r['cases']['baseline']['mean_signed_fractional_error']>0]),
                       ('initial_underprediction',[r for r in disks if r['cases']['baseline']['mean_signed_fractional_error']<=0])]:
        vals=[r['cases']['restricted_joint'] for r in rows]
        directions[label]=dict(galaxies=len(vals),distance_decreased=sum(v['distance_factor']<1 for v in vals),
            inclination_decreased=sum(v['inclination_shift_deg']<0 for v in vals),stellar_ML_decreased=sum(v['stellar_ML_factor']<1 for v in vals),
            gas_decreased=sum(v['gas_factor']<1 for v in vals),medians={k:float(np.median([v[k] for v in vals])) for k in
                ['distance_factor','inclination_shift_deg','stellar_ML_factor','gas_factor']})
    unresolved=[r for r in disks if r['cases']['restricted_joint']['fractional_RMS']>=.20]
    broad=[r for r in disks if r['cases']['wide_joint']['fractional_RMS']>=.20]
    correlations={}
    for attr in ['inc','SBe','D','Q','type']:
        x=[r['catalog'][attr] for r in disks];y=[r['cases']['restricted_joint']['fractional_RMS'] for r in disks]
        rr=spearmanr(x,y);correlations[attr]=dict(spearman=float(rr.statistic),p_descriptive=float(rr.pvalue))
    return dict(experiment='JR-2',frozen_model='JR-1 R10',replay=replay,disk_summary=summary,conditional_directions=directions,
        unresolved_restricted=[r['name'] for r in unresolved],unresolved_wide=[r['name'] for r in broad],
        descriptive_correlations=correlations,
        shared_calibration={k:v for k,v in shared.items() if k!='rows'},
        shared_performance=dict(mean_RMSE_catalog_equivalent_kms=float(np.mean([r['RMSE_catalog_equivalent_kms'] for r in shared['rows']])),
            median_fractional_RMS=float(np.median([r['fractional_RMS'] for r in shared['rows']])),within10pct=sum(r['fractional_RMS']<.1 for r in shared['rows']),within20pct=sum(r['fractional_RMS']<.2 for r in shared['rows'])),
        lens_results=lenses,limitations=['Counterfactual galaxy-property inference, not altered measurements or independent validation',
            'Restricted stellar/gas uncertainty scales are declared working assumptions',
            'Disk projected velocities reconstructed from catalog, not full 2D reanalysis',
            'Lens environment estimates fit one Einstein-radius condition algebraically, not a physical mass map',
            'No galaxy cluster forward model run; no timing, solar system, source-energy or microscopic result'])


def main():
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    out=args.output
    if (out/'summary.json').exists():raise FileExistsError('Use a fresh output')
    out.mkdir(parents=True,exist_ok=True)
    repo=args.bundle/'repository';code=repo/'research_work/results/joint-response-iteration';sys.path.insert(0,str(code))
    import population_followup as P
    J=P.J;saved=json.loads((args.bundle/'evidence/R10_spheroid_stellar_population.json').read_text());p=saved['parameters'];spec=saved['specification']
    inputs=[repo/k for k in J.INPUT_NAMES]+[code/'population_followup.py',code/'run.py',args.bundle/'evidence/R10_spheroid_stellar_population.json']
    hashes={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs};dump(out/'frozen-inputs.json',dict(sha256=hashes,parameters=p,specification=spec))
    t=time.monotonic();E=P.Experiment()
    ds,disks,replay=disk_run(E,J,P,p,saved,out)
    lenses=lens_run(E,J,P,p,spec,saved,out)
    shared=common_fit(ds,disks,out)
    result=summarize(disks,lenses,shared,replay);result['elapsed_seconds']=time.monotonic()-t
    result['unchanged_input_files_verified']=all(hashlib.sha256(f.read_bytes()).hexdigest()==hashes[str(f)] for f in inputs)
    assert result['unchanged_input_files_verified']
    dump(out/'summary.json',result)
    print('SUMMARY',json.dumps({k:v for k,v in result.items() if k!='lens_results'}),flush=True)
if __name__=='__main__':main()
