"""Exact stationary endpoint for same-channel isotropic conservative re-emission."""
from pathlib import Path
import hashlib, io, json, zipfile
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

P=Path(__file__).resolve().parent;OLD=P.parent/'isotropic-galaxy-transfer'
BASE=P.parents[2]/'temporal_candidate_audit/data'
files=[OLD/'third-radiation-retention-results.json',OLD/'third-radiation-retention-predictions.json',P/'occupancy-transport-results.json']
reference=json.loads(files[0].read_text());saved=json.loads(files[1].read_text());previous=json.loads(files[2].read_text())
model=reference['models']['attenuated'];prior={r['galaxy']:r for r in previous['rows']}
source_C=model['source_C_before_retention_Msun_kpc3']
assert abs(source_C/(2*model['C_Msun_kpc3'])-1)<1e-12
for name,digest in reference['input_sha256'].items():
    assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
def eta(x):
    a=np.cbrt(x);return a/(1+a)
def mass_ratio(radii,a,n):
    nodes,weights=leggauss(n)
    x=radii[:,None]/a*(nodes+1)/2;xx=x[:,:,None];mu=nodes
    t=xx*mu;b2=1+xx*xx*(1-mu*mu);b=np.sqrt(b2)
    tau=model['k0_per_kpc']*a*(t/(2*b2*(b2+t*t))+(np.arctan(t/b)+np.pi/2)/(2*b**3))
    J=np.exp(-tau)@weights/2
    radial=x*x/(1+x*x)**2
    old=np.sum(radial*J*weights,axis=1)*radii/a/2
    numeric=np.sum(radial*weights,axis=1)*radii/a/2
    R=radii/a
    analytic=(np.arctan(R)-R/(1+R*R))/2
    assert max(abs(numeric/analytic-1))<1e-7
    return analytic/old

out=dict(scope='Conservative isotropic same-channel release in a stationary isotropic boundary bath: exact i=J=1 solution, no time-dependent retention or microscopic channel derived',
    input_sha256=reference['input_sha256'],model_input_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files},rows=[],scores=[])
cache=[]
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as archive:
    for row in saved:
        name=row['galaxy'];p=prior[name];a=p['a_kpc'];X=p['X']
        arr=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name+'_rotmod.dat'))))
        vb=arr[:,3]*abs(arr[:,3])+.5*arr[:,4]*abs(arr[:,4])+.7*arr[:,5]*abs(arr[:,5])
        good=np.isfinite(arr).all(axis=1)&(arr[:,0]>0)&(arr[:,1]>0)&(arr[:,2]>0)&(vb>0)
        arr,vb=arr[good],vb[good]
        assert arr[:,0].tolist()==row['R_kpc'] and arr[:,1].tolist()==row['observed_kms']
        pred=np.array(row['predicted_kms']);extra=pred**2-vb
        ratio64=mass_ratio(arr[:,0],a,64);ratio=mass_ratio(arr[:,0],a,128)
        v=np.sqrt(vb+extra*ratio);v64=np.sqrt(vb+extra*ratio64)
        error=float(max(abs(v-v64)));assert error<.01
        assert min(v-np.array(p['predicted_kms']))>-.1
        f=float(eta(X))
        out['rows'].append(dict(galaxy=name,split=row['split'],X=X,a_kpc=a,occupancy=f,
            total_deposit_mass_Msun=float(np.pi**2*source_C*a**3*f),
            R_kpc=arr[:,0].tolist(),predicted_kms=v.tolist(),extra_mass_ratio=ratio.tolist(),
            maximum_speed_change_from_reference_kms=float(max(abs(v-pred))),quadrature_speed_difference_kms=error))
        cache.append(dict(split=row['split'],vb=vb,extra=extra,ratio=ratio,obs=arr[:,1]))
for branch in ['reference','conservative_recycling']:
    def loss(amplitude,split):
        return float(np.mean([np.mean((np.sqrt(r['vb']+amplitude*r['extra']*(1 if branch=='reference' else r['ratio']))-r['obs'])**2) for r in cache if r['split']==split]))
    fit=minimize_scalar(lambda a:loss(a,'train'),bounds=(0,4),method='bounded',options={'xatol':1e-10})
    assert fit.success and 1e-7<fit.x<4-1e-7
    for mode,amp in [('frozen',1),('training_C_refit',float(fit.x))]:
        scores={s:float(np.sqrt(loss(amp,s))) for s in ['train','validation','test']}
        if branch=='reference' and mode=='frozen':
            for s,val in scores.items():
                assert min(abs(val-model[k][s]['RMSE_kms']) for k in ['scores','finer_scores'])<1e-8
        out['scores'].append(dict(branch=branch,mode=mode,C_multiplier=amp,RMSE_kms=scores))
assert len(saved)==149 and sum(len(r['R_kpc']) for r in saved)==3150
out['maximum_quadrature_speed_difference_kms']=max(r['quadrature_speed_difference_kms'] for r in out['rows'])
(P/'conservative-recycling-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['scores'],indent=2));print(out['maximum_quadrature_speed_difference_kms'])
