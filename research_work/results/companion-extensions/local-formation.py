"""Finite local filling in the prescribed full-opacity field: a conditional lower-field diagnostic."""
from pathlib import Path
import hashlib, io, json, zipfile
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

P=Path(__file__).resolve().parent; OLD=P.parent/'isotropic-galaxy-transfer'
BASE=P.parents[2]/'temporal_candidate_audit/data'
files=[OLD/'third-radiation-retention-results.json',OLD/'third-radiation-retention-predictions.json',P/'local-capacity-results.json']
ref=json.loads(files[0].read_text());saved=json.loads(files[1].read_text()); prior={r['galaxy']:r for r in json.loads(files[2].read_text())['rows']}
model=ref['models']['attenuated'];durations=[.01,.1,1.,10.,100.,1000.,None]
for name,digest in ref['input_sha256'].items():
    assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
def eta(x):
    q=np.cbrt(x);return q/(1+q)
def ratios(radii,a,X,n):
    nodes,weights=leggauss(n)
    r=radii[:,None]/a*(nodes+1)/2;x=r[:,:,None];mu=nodes
    t=x*mu;b2=1+x*x*(1-mu*mu);b=np.sqrt(b2)
    tau=model['k0_per_kpc']*a*(t/(2*b2*(b2+t*t))+(np.arctan(t/b)+np.pi/2)/(2*b**3))
    J=np.exp(-tau)@weights/2;Y=X*J
    radial=r*r/(1+r*r)**2
    old=np.sum(radial*J*eta(X)*weights,axis=1)
    tn,tw=leggauss(2*n);logs=6*np.log(10)*tn;theta=np.exp(logs)
    w=tw/(np.cosh(logs/3)+.5);w/=sum(w)
    rate=Y[:,:,None]+theta;feq=Y[:,:,None]/rate
    result={}
    last=np.zeros_like(Y)
    for u in durations:
        f=(feq*(1 if u is None else -np.expm1(-rate*u)))@w
        assert np.min(f-last)>-1e-12 and np.max(f)<1
        last=f
        result[str(u)]=np.sum(radial*f*weights,axis=1)/old
    return result

out=dict(scope='Empty initial local threshold storage, constant full-opacity incident field, shared dimensionless durations; no self-consistent time-dependent transfer or physical age',
    input_sha256=ref['input_sha256'],model_input_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files},scores=[],rows=[])
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
        low=ratios(arr[:,0],a,X,64);high=ratios(arr[:,0],a,X,128)
        error=max(float(max(abs(np.sqrt(vb+extra*high[str(u)])-np.sqrt(vb+extra*low[str(u)])))) for u in durations)
        assert error<.01,(name,error)
        cache.append(dict(galaxy=name,split=row['split'],vb=vb,extra=extra,ratios=high,obs=arr[:,1],quadrature_speed_difference_kms=error))
def loss(amp,split,key):
    return float(np.mean([np.mean((np.sqrt(r['vb']+amp*r['extra']*(1 if key=='reference' else r['ratios'][key]))-r['obs'])**2) for r in cache if r['split']==split]))
for key in ['reference']+[str(u) for u in durations]:
    hi=4.
    while loss(hi,'train',key)<loss(hi/2,'train',key):
        hi*=2
        assert hi<1e7
    fit=minimize_scalar(lambda a:loss(a,'train',key),bounds=(0,hi),method='bounded',options={'xatol':1e-10})
    assert fit.success and 1e-7<fit.x<hi-1e-7
    for mode,amp in [('frozen',1.),('training_C_refit',float(fit.x))]:
        score=dict(duration=key,mode=mode,C_multiplier=amp,RMSE_kms={s:float(np.sqrt(loss(amp,s,key))) for s in ['train','validation','test']})
        out['scores'].append(score)
        if mode=='training_C_refit':
            for r in cache:
                v=np.sqrt(r['vb']+amp*r['extra']*(1 if key=='reference' else r['ratios'][key]))
                out['rows'].append(dict(galaxy=r['galaxy'],split=r['split'],duration=key,RMSE_kms=float(np.sqrt(np.mean((v-r['obs'])**2)))))
selected=min([r for r in out['scores'] if r['mode']=='training_C_refit' and r['duration']!='reference'],key=lambda r:r['RMSE_kms']['train'])
out['training_selected_grid_duration']=selected
out['maximum_quadrature_speed_difference_kms']=max(r['quadrature_speed_difference_kms'] for r in cache)
assert len(saved)==149 and sum(len(r['R_kpc']) for r in saved)==3150
control=next(r for r in out['scores'] if r['duration']=='reference' and r['mode']=='frozen')
for split,v in control['RMSE_kms'].items():
    assert min(abs(v-model[k][split]['RMSE_kms']) for k in ['scores','finer_scores'])<1e-8
(P/'local-formation-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([r for r in out['scores'] if r['mode']=='training_C_refit'],indent=2));print(out['maximum_quadrature_speed_difference_kms'])
