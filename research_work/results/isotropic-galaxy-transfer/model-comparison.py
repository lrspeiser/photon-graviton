"""Matched observed-rotation comparison with explicit fitting budgets."""
from pathlib import Path
import io,json,zipfile,hashlib
import numpy as np
from scipy.optimize import minimize,minimize_scalar
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'temporal_candidate_audit/data'
G=4.30091727003628e-6;kpc=3.085677581491367e19
reference=json.loads((HERE/'third-radiation-retention-results.json').read_text())
for name,digest in reference['input_sha256'].items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
splits=json.loads((BASE/'sparc_frozen.json').read_text())['split'];catalog={}
for line in (BASE/'SPARC_Lelli2016c.mrt').read_text().splitlines():
    f=line.split()
    if len(f)==19:
        try:catalog[f[0]]=(float(f[11]),float(f[7]))
        except ValueError:pass
saved={r['galaxy']:r for r in json.loads((HERE/'third-radiation-retention-predictions.json').read_text())}
data=[]
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as archive:
    for split,names in splits.items():
        for name in names:
            arr=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name+'_rotmod.dat'))))
            vb=arr[:,3]*abs(arr[:,3])+.5*arr[:,4]*abs(arr[:,4])+.7*arr[:,5]*abs(arr[:,5])
            good=np.isfinite(arr).all(axis=1)&(arr[:,0]>0)&(arr[:,1]>0)&(arr[:,2]>0)&(vb>0)
            arr,vb=arr[good],vb[good];rd,L=catalog[name]
            assert arr[:,0].tolist()==saved[name]['R_kpc'] and arr[:,1].tolist()==saved[name]['observed_kms']
            data.append(dict(name=name,split=split,r=arr[:,0],y=arr[:,1],err=arr[:,2],vb=vb,rd=rd,X=L/rd**2,comp=np.array(saved[name]['predicted_kms'])))
assert len(data)==149 and sum(len(d['y']) for d in data)==3150
train=data[:89]
def mond(d,a0):
    gb=d['vb']/d['r'];a=a0*kpc/1e6
    return np.sqrt(d['r']*(.5*gb+np.sqrt(.25*gb**2+a*gb)))
def nfw(d,rho,rs):
    x=d['r']/rs;term=np.log1p(x)-x/(1+x)
    tiny=x<1e-4;term[tiny]=x[tiny]**2/2-2*x[tiny]**3/3+3*x[tiny]**4/4
    return np.sqrt(d['vb']+G*4*np.pi*rho*rs**3*term/d['r'])
def global_nfw(d,p):return nfw(d,10**p[0]*d['X']**p[2],10**p[1]*d['rd'])
def error(d,v,ix=slice(None)):return float(np.mean(np.log10(v[ix]/d['y'][ix])**2))
fitmond=minimize_scalar(lambda p:np.mean([error(d,mond(d,10**p)) for d in train]),bounds=(-12,-8),method='bounded',options={'xatol':1e-10})
assert fitmond.success
bounds=[(-2,12),(-1,2),(-2,2)]
opts=[minimize(lambda p:np.mean([error(d,global_nfw(d,p)) for d in train]),start,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-12,'maxiter':500}) for start in [[7,0,0],[6,1,.5],[8,-.5,-.5]]]
assert all(o.success for o in opts)
best=min(opts,key=lambda o:o.fun)
def limit_nfw(d,p):return np.sqrt(d['vb']+10**p[0]*d['rd']*d['X']**p[1]*d['r'])
limitfits=[minimize(lambda p:np.mean([error(d,limit_nfw(d,p)) for d in train]),start,method='L-BFGS-B',bounds=[(-8,8),(-2,2)],options={'ftol':1e-12,'maxiter':500}) for start in [[2,0],[3,.5],[1,-.5]]]
assert all(o.success for o in limitfits)
limitbest=min(limitfits,key=lambda o:o.fun)
def boundary(p,b):return any(min(abs(v-lo),abs(v-hi))<1e-5 for v,(lo,hi) in zip(p,b))
models={'baryons':lambda d:np.sqrt(d['vb']),'companion_third':lambda d:d['comp'],'MOND_simple_fitted':lambda d:mond(d,10**fitmond.x),'MOND_simple_fixed':lambda d:mond(d,1.2e-10),'NFW_shared_scaling':lambda d:global_nfw(d,best.x)}
models['NFW_shared_large_scale_limit']=lambda d:limit_nfw(d,limitbest.x)
rows=[];summary=[]
for name,fn in models.items():
    for d in data:
        v=fn(d);rows.append(dict(model=name,galaxy=d['name'],split=d['split'],R_kpc=d['r'].tolist(),observed_kms=d['y'].tolist(),measurement_error_kms=d['err'].tolist(),predicted_kms=v.tolist()))
    for split in splits:
        dd=[d for d in data if d['split']==split]
        summary.append(dict(model=name,split=split,n=len(dd),RMSE_kms=float(np.sqrt(np.mean([np.mean((fn(d)-d['y'])**2) for d in dd]))),log_RMS=float(np.sqrt(np.mean([error(d,fn(d)) for d in dd]))),mean_diagonal_chi2_per_point_per_galaxy=float(np.mean([np.mean(((fn(d)-d['y'])/d['err'])**2) for d in dd]))))
adaptive=[];outer=[]
for d in data[89:]:
    n=int(np.ceil(.6*len(d['y'])));idx=slice(0,n);b=[(-2,12),(-2,3)]
    fits=[minimize(lambda p:error(d,nfw(d,10**p[0],10**p[1]),idx),s,method='L-BFGS-B',bounds=b,options={'ftol':1e-12,'maxiter':500}) for s in [[7,0],[6,1],[8,-1]]]
    good=[f for f in fits if f.success];assert good
    fit=min(good,key=lambda o:o.fun);v=nfw(d,10**fit.x[0],10**fit.x[1])
    adaptive.append(dict(galaxy=d['name'],split=d['split'],n_inner=n,rho_s_Msun_kpc3=float(10**fit.x[0]),r_s_kpc=float(10**fit.x[1]),boundary=boundary(fit.x,b),outer_observed_kms=d['y'][n:].tolist(),outer_predicted_kms=v[n:].tolist()))
    for name,fn in dict(models,NFW_target_inner_fit=lambda _:v).items():
        pred=fn(d)[n:];obs=d['y'][n:]
        outer.append(dict(model=name,split=d['split'],galaxy=d['name'],mse=float(np.mean((pred-obs)**2)),log_mse=float(np.mean(np.log10(pred/obs)**2))))
outer_summary=[]
for name in list(models)+['NFW_target_inner_fit']:
    for split in ['validation','test']:
        rr=[r for r in outer if r['model']==name and r['split']==split]
        outer_summary.append(dict(model=name,split=split,RMSE_kms=float(np.sqrt(np.mean([r['mse'] for r in rr]))),log_RMS=float(np.sqrt(np.mean([r['log_mse'] for r in rr])))))
out=dict(status='Selected rotation prescriptions; exposed partitions; no full-theory ranking',input_sha256=reference['input_sha256'],MOND_a0_SI=float(10**fitmond.x),MOND_boundary=bool(abs(fitmond.x+12)<1e-5 or abs(fitmond.x+8)<1e-5),NFW_shared_parameters=best.x.tolist(),NFW_boundary=boundary(best.x,bounds),NFW_starts=[dict(parameters=o.x.tolist(),loss=float(o.fun),success=bool(o.success)) for o in opts],global_parameter_counts=dict(baryons=0,companion_third=3,MOND_simple_fitted=1,MOND_simple_fixed=0,NFW_shared_scaling=3),summary=summary,outer_summary=outer_summary,adaptive_NFW_boundary_count=sum(r['boundary'] for r in adaptive),adaptive_NFW=adaptive)
out['NFW_large_scale_limit_parameters']=limitbest.x.tolist()
out['global_parameter_counts']['NFW_shared_large_scale_limit']=2
# The simple MOND result must satisfy mu(g/a0)*g=g_b, mu(x)=x/(1+x).
for d in data:
    gb=d['vb']/d['r'];g=mond(d,10**fitmond.x)**2/d['r'];a0=10**fitmond.x*kpc/1e6
    assert np.max(abs((g*g/(g+a0))/gb-1))<1e-12
out['MOND_implicit_equation_check']=True
(HERE/'model-comparison-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
(HERE/'model-comparison-predictions.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print('a0',out['MOND_a0_SI'],'NFW',out['NFW_shared_parameters'],'boundaries',out['NFW_boundary'],out['adaptive_NFW_boundary_count'])
for r in summary:print(r['model'],r['split'],r['RMSE_kms'],r['log_RMS'])
print('OUTER');print(json.dumps(outer_summary,indent=2))

