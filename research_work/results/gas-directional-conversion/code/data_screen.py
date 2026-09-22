"""JR-5: transfer the toy ray-conversion law to frozen SPARC/R10 inputs.

Observational numbers remain fixed. Gas profiles are source-only exponential
reconstructions, NOT measured three-dimensional HI maps. A conditional screen,
not a microscopic derivation, independent confirmation or lensing validation.
"""
from __future__ import annotations
import argparse,hashlib,io,json,time,zipfile
from pathlib import Path
import numpy as np
from scipy.special import i0e,i1e,k0e,k1e
from scipy.optimize import minimize_scalar,least_squares
from core import AxisModel,transfer,baseline,G


def gas_speed2(r,M,h):
    y=np.maximum(r/(2*h),1e-8)
    return 2*G*M/h*y*y*(i0e(y)*k0e(y)-i1e(y)*k1e(y))

def prepare(jr1,nr=256,nmu=32,lmax=12,height=.1,profiles=None):
    rec=json.loads((jr1/'evidence/R10_spheroid_stellar_population.json').read_text())
    p=rec['parameters'];repo=jr1/'repository';cat={}
    catalog=repo/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt'
    raw=repo/'temporal_candidate_audit/data/Rotmod_LTG.zip'
    for line in catalog.read_text().splitlines():
        f=line.split()
        if len(f)==19:
            try:cat[f[0]]=dict(rd=float(f[11]),Mgas=float(f[13])*1.33e9,RHI=float(f[14]),inc=float(f[5]),quality=int(f[17]))
            except ValueError:pass
    profiles={} if profiles is None else {d['name']:d for d in profiles}
    galaxies=[];source_records=[]
    with zipfile.ZipFile(raw) as z:
        for row in rec['metrics']['sparc_rows']:
            name=row['name'];c=cat[name]
            data=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(name+'_rotmod.dat'))));data=data[data[:,0]>0]
            r=np.array(row['r_kpc']);assert np.array_equal(r,data[:,0])
            M=c['Mgas'];target=data[:,3]*abs(data[:,3])
            if name in profiles:profile=profiles[name];h=profile['gas_scale_kpc']
            else:
                lo=np.log(.1*c['rd']);hi=np.log(30*c['rd'])
                denom=max(np.mean(target**2),1.)
                if M>0:
                    res=minimize_scalar(lambda lh:np.mean((gas_speed2(r,M,np.exp(lh))-target)**2)/denom,bounds=(lo,hi),method='bounded',options={'xatol':1e-10})
                    h=float(np.exp(res.x));relative=float(np.sqrt(res.fun));boundary=min(res.x-lo,hi-res.x)<1e-3
                else:h=c['rd'];relative=0.;boundary=False
                sigma0=M/(1.33*2*np.pi*h*h)/1e6
                rhimodel=h*np.log(sigma0) if sigma0>1 else 0.
                profile=dict(name=name,Mgas_Msun=M,gas_scale_kpc=h,gas_height_kpc=height*h,
                    gas_force_relative_RMS=relative,scale_bound=bool(boundary),catalog_RHI_kpc=c['RHI'],model_RHI_kpc=rhimodel,
                    negative_gas_force_rows=int((target<0).sum()),catalog_quality=c['quality'])
            source_records.append(profile)
            u=10**p['logusph'];mdisk=row['Mstar_input']-.0
            # q follows the original observed bulge contribution. No total Vobs.
            Mbul=float(np.max(.7*r*data[:,5]**2/G));fs=Mbul/row['Mstar_input']
            fs=u*fs/(1+(u-1)*fs);q=p['q']+(p['q_sph']-p['q'])*fs
            A,rc,rt=[row[k] for k in ('A_km2_s2','rc_kpc','rt_kpc')]
            m=AxisModel(A,rc,rt,q,row['Re_kpc'],nr,nmu,lmax)
            from core import gas_density
            m.gas=gas_density(m.mid,m.mu,M,h,height)
            g0=baseline(r,A,rc,rt,q)[0]
            baryon=np.maximum(data[:,3]*abs(data[:,3])+10**p['logu']*(.5*data[:,4]**2+u*.7*data[:,5]**2),0.)
            rebuilt=np.sqrt(baryon+r*g0)
            assert np.max(abs(rebuilt-np.array(row['predicted_kms'])))<1e-10
            galaxies.append(dict(name=name,split=row['split'],r=r,y=np.array(row['observed_kms']),err=np.array(row['error_kms']),
                baseline=rebuilt,gb_v2=baryon,model=m,profile=profile))
    return galaxies,source_records


def predict_one(g,alpha,beta,quadratic=True,monopole=False):
    m=g['model'];m.f=transfer(m.edges,m.gas,float(alpha),float(beta),quadratic)
    m.delta=m.rho0[:,None]*(2*m.f-1)
    m.coeff=(m.delta*m.w)@m.P.T*(2*m.ls+1)[None,:]
    m._potential()
    v2=g['gb_v2']+g['r']*m.radial(g['r'],0.,monopole)
    if np.any(v2<0):
        # Encode a negative inward circular force as a signed diagnostic speed;
        # do not hide it by clipping it to zero or deleting the row.
        return np.sign(v2)*np.sqrt(abs(v2))
    return np.sqrt(v2)


def score(galaxies,alpha,beta,quadratic=True,splits=None,monopole=False):
    rows=[]
    for g in galaxies:
        if splits is not None and g['split'] not in splits:continue
        v=predict_one(g,alpha,beta,quadratic,monopole)
        e=v-g['y'];fr=e/g['y'];old=(g['baseline']-g['y'])/g['y']
        rows.append(dict(name=g['name'],split=g['split'],points=len(v),r_kpc=g['r'].tolist(),observed_kms=g['y'].tolist(),
            error_kms=g['err'].tolist(),baseline_kms=g['baseline'].tolist(),predicted_kms=v.tolist(),
            RMS_kms=float(np.sqrt(np.mean(e*e))),fractional_RMS=float(np.sqrt(np.mean(fr*fr))),
            original_fractional_RMS=float(np.sqrt(np.mean(old*old))),original_mean_fractional_bias=float(np.mean(old)),
            chi2=float(np.sum((e/g['err'])**2)),nonpositive_circular_speed_count=int((v<=0).sum())))
    groups={}
    for split in ['train','validation','test','all']:
        rr=[r for r in rows if split=='all' or r['split']==split]
        if not rr:continue
        err=np.array([r['fractional_RMS'] for r in rr]);old=np.array([r['original_fractional_RMS'] for r in rr])
        groups[split]=dict(galaxies=len(rr),points=sum(r['points'] for r in rr),mean_RMS_kms=float(np.mean([r['RMS_kms'] for r in rr])),
            mean_fractional_square=float(np.mean(err*err)),median_fractional_RMS=float(np.median(err)),
            below_10pct=int((err<.1).sum()),below_20pct=int((err<.2).sum()),
            rescued_20pct=int(((old>=.2)&(err<.2)).sum()),spoiled_20pct=int(((old<.2)&(err>=.2)).sum()),
            raw_chi2=sum(r['chi2'] for r in rr),negative_rows=sum(r['nonpositive_circular_speed_count'] for r in rr))
    return dict(groups=groups,rows=rows)


def main():
    p=argparse.ArgumentParser(__doc__);p.add_argument('--jr1',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--max-evals',type=int,default=70);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    if (a.output/'data-summary.json').exists():raise FileExistsError('Preserve prior evidence')
    t0=time.monotonic();gs,profiles=prepare(a.jr1)
    (a.output/'gas-profiles.json').write_text(json.dumps(profiles,indent=2)+'\n')
    print('PREPARED',len(gs),'gas relative RMS median',np.median([p['gas_force_relative_RMS'] for p in profiles]),flush=True)
    train=[g for g in gs if g['split']=='train']
    n=len(train)
    specs=[dict(name='null',kind='null',quadratic=True),dict(name='forward_only',kind='forward',quadratic=True),
      dict(name='reverse_only',kind='reverse',quadratic=True),dict(name='both_linear',kind='both',quadratic=False),
      dict(name='forward_linear_reverse_quadratic',kind='both',quadratic=True)]
    records=[]
    for spec in specs:
        kind=spec['kind'];quadratic=spec['quadratic'];history=[]
        def unpack(x):
            z=10**np.asarray(x)
            return (z[0],0.) if kind=='forward' else (0.,z[0]) if kind=='reverse' else (z[0],z[1])
        def residual(x):
            alpha,beta=unpack(x)
            return np.concatenate([(predict_one(g,alpha,beta,quadratic)-g['y'])/g['y']/np.sqrt(len(g['y'])*n) for g in train])
        if kind=='null':alpha=beta=0.
        else:
            starts=[[-2.],[-.5],[.8]] if kind!='both' else [[-2.,-2.],[-.5,-.5],[.8,.8],[-1.,.8],[.8,-1.]]
            best=None
            for x in starts:
                tick=time.monotonic();res=least_squares(residual,x,bounds=(-6.,np.log10(30.)),max_nfev=a.max_evals,
                    ftol=2e-6,xtol=2e-6,gtol=2e-6,diff_step=1e-4)
                alpha,beta=unpack(res.x);objective=float(res.fun@res.fun)
                history.append(dict(start_log10=x,alpha=float(alpha),beta=float(beta),objective=objective,nfev=res.nfev,
                    success=bool(res.success),seconds=time.monotonic()-tick,
                    bound_hit=bool(np.min(np.minimum(res.x+6,np.log10(30.)-res.x))<1e-3)))
                if best is None or objective<best[0]:best=(objective,float(alpha),float(beta))
                print('FIT',spec['name'],alpha,beta,objective,res.nfev,flush=True)
            _,alpha,beta=best
        result=score(gs,alpha,beta,quadratic,splits=['train','validation'])
        record=dict(spec=spec,alpha=float(alpha),beta=float(beta),optimizer_starts=history,result=result)
        records.append(record)
        (a.output/(spec['name']+'.json')).write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        print('SCORES',spec['name'],json.dumps(result['groups']),flush=True)
    selected=min(records,key=lambda r:r['result']['groups']['validation']['mean_fractional_square'])
    decision=dict(selected=selected['spec']['name'],alpha=selected['alpha'],beta=selected['beta'],quadratic=selected['spec']['quadratic'],
        rule='29-role equal-galaxy mean fractional squared error, includes null; selection fixed before 31-role scoring',
        candidates=[dict(name=r['spec']['name'],alpha=r['alpha'],beta=r['beta'],groups=r['result']['groups']) for r in records])
    (a.output/'selection.json').write_text(json.dumps(decision,indent=2)+'\n')
    full=score(gs,selected['alpha'],selected['beta'],selected['spec']['quadratic'])
    (a.output/'selected-full.json').write_text(json.dumps(full,indent=2,allow_nan=False)+'\n')
    null=score(gs,0.,0.)
    (a.output/'null-full.json').write_text(json.dumps(null,indent=2)+'\n')
    # Also preserve the best genuinely gas-active model, even if null is chosen.
    active=min(records[1:],key=lambda r:r['result']['groups']['validation']['mean_fractional_square'])
    activefull=score(gs,active['alpha'],active['beta'],active['spec']['quadratic'])
    (a.output/'active-full.json').write_text(json.dumps(activefull,indent=2)+'\n')
    monopole=score(gs,active['alpha'],active['beta'],active['spec']['quadratic'],monopole=True)
    (a.output/'monopole-control.json').write_text(json.dumps(monopole,indent=2)+'\n')
    summary=dict(experiment='JR-5 real-data screen',selection=decision,selected_groups=full['groups'],baseline_groups=null['groups'],
        active_model=dict(name=active['spec']['name'],alpha=active['alpha'],beta=active['beta'],quadratic=active['spec']['quadratic'],groups=activefull['groups']),
        monopole_control_groups=monopole['groups'],seconds=time.monotonic()-t0,
        gas_reconstruction=dict(median_gas_force_relative_RMS=float(np.median([p['gas_force_relative_RMS'] for p in profiles])),
          above_20pct=int(sum(p['gas_force_relative_RMS']>.2 for p in profiles)),scale_bound_count=sum(p['scale_bound'] for p in profiles)),
        scope='Axisymmetric inferred gas profiles, fixed R10 source/capture scaffold; not resolved maps, fresh blind data or a microscopic energy-budget derivation')
    (a.output/'data-summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('FINAL',json.dumps(summary),flush=True)
if __name__=='__main__':main()
