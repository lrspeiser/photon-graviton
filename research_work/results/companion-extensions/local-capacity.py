"""Local equilibrium occupancy with fixed incident transport; exposed SPARC diagnostic."""
from pathlib import Path
import hashlib, io, json, zipfile
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

P=Path(__file__).resolve().parent
OLD=P.parent/'isotropic-galaxy-transfer'
BASE=P.parents[2]/'temporal_candidate_audit/data'
fitpath=OLD/'third-radiation-retention-results.json'
predpath=OLD/'third-radiation-retention-predictions.json'
xpath=P/'threshold-galaxies-results.json'
reference=json.loads(fitpath.read_text()); saved=json.loads(predpath.read_text())
model=reference['models']['attenuated']
xs={r['galaxy']:r['X'] for r in json.loads(xpath.read_text())['rows'] if r['half_width_decades']==0}
for name,digest in reference['input_sha256'].items():
    assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
rd={}
for line in (BASE/'SPARC_Lelli2016c.mrt').read_text().splitlines():
    f=line.split()
    if len(f)==19:
        try:rd[f[0]]=float(f[11])
        except ValueError:pass
def eta(x):
    a=np.cbrt(x)
    return a/(1+a)
def integrals(radii,a,X,n):
    nodes,weights=leggauss(n)
    r=radii[:,None]*(nodes[None,:]+1)/2
    x=r[:,:,None]/a
    mu=nodes[None,None,:]
    b2=1+x*x*(1-mu*mu); b=np.sqrt(b2); t=x*mu
    tau=model['k0_per_kpc']*a*(t/(2*b2*(b2+t*t))+(np.arctan(t/b)+np.pi/2)/(2*b**3))
    J=np.sum(np.exp(-np.maximum(tau,0))*weights[None,None,:],axis=2)/2
    assert np.all((J>0)&(J<=1+1e-12))
    radial=r*r/(1+(r/a)**2)**2
    old=np.sum(radial*J*eta(X)*weights[None,:],axis=1)*radii/2
    new=np.sum(radial*eta(X*J)*weights[None,:],axis=1)*radii/2
    assert np.all(new>=old*(1-1e-12))
    return new/old

out=dict(scope='Unattenuated capacity Cg with local eta(XJ); fixed reference J, no self-consistent occupancy feedback on transport; exposed inputs',
         input_sha256=reference['input_sha256'],model_input_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [fitpath,predpath,xpath]},rows=[],scores=[])
cache=[]
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as archive:
    for row in saved:
        name=row['galaxy'];arr=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name+'_rotmod.dat'))))
        vb=arr[:,3]*abs(arr[:,3])+.5*arr[:,4]*abs(arr[:,4])+.7*arr[:,5]*abs(arr[:,5])
        good=np.isfinite(arr).all(axis=1)&(arr[:,0]>0)&(arr[:,1]>0)&(arr[:,2]>0)&(vb>0)
        arr,vb=arr[good],vb[good]
        assert arr[:,0].tolist()==row['R_kpc'] and arr[:,1].tolist()==row['observed_kms']
        pred=np.array(row['predicted_kms']);extra=pred**2-vb
        assert min(extra)>-1e-8
        a=model['scale_to_disk']*rd[name];X=xs[name]
        ratios={n:integrals(arr[:,0],a,X,n) for n in [64,128]}
        coarse=np.sqrt(vb+extra*ratios[64]);v=np.sqrt(vb+extra*ratios[128])
        error=float(max(abs(v-coarse)));order=128
        if error>.001:
            ratio=integrals(arr[:,0],a,X,256)
            refined=np.sqrt(vb+extra*ratio);error=float(max(abs(refined-v)));v=refined;order=256
        else:ratio=ratios[128]
        assert error<.01,(name,error)
        cache.append(dict(split=row['split'],vb=vb,extra=extra,ratio=ratio,obs=arr[:,1]))
        out['rows'].append(dict(galaxy=name,split=row['split'],X=X,a_kpc=a,
            R_kpc=arr[:,0].tolist(),observed_kms=arr[:,1].tolist(),predicted_kms=v.tolist(),
            extra_gravity_ratio=ratio.tolist(),RMSE_kms=float(np.sqrt(np.mean((v-arr[:,1])**2))),
            max_speed_change_kms=float(max(abs(v-pred))),quadrature_order=order,quadrature_speed_difference_kms=error))
for branch in ['reference','local']:
    def loss(amplitude,split):
        values=[]
        for r in cache:
            if r['split']!=split:continue
            ratio=1 if branch=='reference' else r['ratio']
            v=np.sqrt(r['vb']+amplitude*r['extra']*ratio)
            values.append(np.mean((v-r['obs'])**2))
        return float(np.mean(values))
    fit=minimize_scalar(lambda a:loss(a,'train'),bounds=(0,4),method='bounded',options={'xatol':1e-10})
    assert fit.success and 1e-7<fit.x<4-1e-7
    for mode,amplitude in [('frozen',1),('training_C_refit',float(fit.x))]:
        scores={split:float(np.sqrt(loss(amplitude,split))) for split in ['train','validation','test']}
        if branch=='reference' and mode=='frozen':
            for split,value in scores.items():
                assert min(abs(value-model[k][split]['RMSE_kms']) for k in ['scores','finer_scores'])<1e-8
        out['scores'].append(dict(branch=branch,mode=mode,C_multiplier=amplitude,RMSE_kms=scores))
assert len(saved)==149 and sum(len(r['R_kpc']) for r in saved)==3150
out['maximum_quadrature_speed_difference_kms']=max(r['quadrature_speed_difference_kms'] for r in out['rows'])
out['maximum_speed_change_kms']=max(r['max_speed_change_kms'] for r in out['rows'])
(P/'local-capacity-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['scores'],indent=2));print(out['maximum_quadrature_speed_difference_kms'],out['maximum_speed_change_kms'])
