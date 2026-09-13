"""One postulated companion-feedback law tested on exposed brightness and redshift."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];data=ROOT/'shared_interaction_test/data'
source=data/'Pantheon+SH0ES.dat';covsource=data/'Pantheon+SH0ES_STAT+SYS.cov'
df=pd.read_csv(source,sep=r'\s+');raw=np.loadtxt(covsource);N=int(raw[0]);C=raw[1:].reshape(N,N);C=(C+C.T)/2
cal=np.flatnonzero(df.IS_CALIBRATOR.to_numpy()==1)
ev=np.flatnonzero((df.IS_CALIBRATOR.to_numpy()==0)&(df.zHD.to_numpy()>=.1))
Cc=C[np.ix_(cal,cal)];Ce=C[np.ix_(ev,ev)];Cec=C[np.ix_(ev,cal)]
w=cho_solve(cho_factor(Cc),np.ones(len(cal)));w/=w.sum()
mu=df.CEPH_DIST.to_numpy()[cal];dc=10**((mu-25)/5)
Mbase=float(w@(df.m_b_corr.to_numpy()[cal]-mu))
var=float(w@Cc@w);cross=Cec@w;V=Ce+var-cross[:,None]-cross[None,:]
z=df.zHD.to_numpy()[ev];zh=df.zHEL.to_numpy()[ev];obs=df.m_b_corr.to_numpy()[ev];A=np.log1p(z)
alpha=.0002488993286382367;c=299792.458
def integrated(a,chi):return a+chi*(np.expm1(a)-a)
def shift_at_distance(distance,chi):
    target=alpha*np.asarray(distance);lo=np.zeros_like(target);hi=target.copy()
    for _ in range(52):
        mid=(lo+hi)/2;below=integrated(mid,chi)<target
        lo=np.where(below,mid,lo);hi=np.where(below,hi,mid)
    ans=(lo+hi)/2
    assert np.max(abs(integrated(ans,chi)-target))<1e-12
    return ans
def predict(chi):
    ac=shift_at_distance(dc,chi)
    M=Mbase-5*float(w@ac)/np.log(10)
    D=integrated(A,chi)/alpha
    return M+25+5*np.log10(D*(1+zh)),M
tr=np.flatnonzero(z<.3);te=np.flatnonzero(z>=.3)
ftr=cho_factor(V[np.ix_(tr,tr)]);fte=cho_factor(V[np.ix_(te,te)])
def objective(chi):
    pred,_=predict(chi);r=obs[tr]-pred[tr]
    return float(r@cho_solve(ftr,r))
fit=minimize_scalar(objective,bounds=(0,10),method='bounded',options={'xatol':1e-10})
assert fit.success
chi=min([0.,float(fit.x),10.],key=objective)
base,_=predict(0);pred,M=predict(chi);r0=obs-base;r=obs-pred
bins=[]
for lo,hi in [(.1,.3),(.3,.6),(.6,1),(1,3)]:
    ix=np.flatnonzero((z>=lo)&(z<hi));sub=V[np.ix_(ix,ix)]
    bw=cho_solve(cho_factor(sub),np.ones(len(ix)));bw/=bw.sum()
    bins.append(dict(z_min=lo,z_max=hi,count=len(ix),baseline_mean=float(bw@r0[ix]),feedback_mean=float(bw@r[ix])))
gs=ROOT/'redshift_paper/all_164_groups.csv';gal=pd.read_csv(gs)
d=gal.catalog_distance_mpc.to_numpy();gz=gal.observed_cmb_z.to_numpy()
gp=np.expm1(shift_at_distance(d,chi));g0=np.expm1(alpha*d)
grows=[]
for split in ['train','validation','test']:
    ix=gal.split.to_numpy()==split
    grows.append(dict(historical_exposed_split=split,count=int(ix.sum()),
        baseline_rmse_kms=float(np.sqrt(np.mean((c*(g0[ix]-gz[ix]))**2))),
        feedback_rmse_kms=float(np.sqrt(np.mean((c*(gp[ix]-gz[ix]))**2)))))
out=dict(status='Phenomenological feedback transfer on reused standardized observations',alpha0_per_mpc=alpha,
    chi=chi,boundary_hit=chi in [0.,10.],M_calibrator=M,train_count=len(tr),transfer_count=len(te),
    train_baseline_chi2=objective(0),train_fitted_chi2=objective(chi),
    transfer_baseline_chi2=float(r0[te]@cho_solve(fte,r0[te])),
    transfer_feedback_chi2=float(r[te]@cho_solve(fte,r[te])),
    scores_use_same_covariance_without_fitted_parameter_uncertainty=True,bins=bins,galaxy_results=grows,
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,covsource,gs,HERE/'feedback-protocol.md']})
(HERE/'feedback-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
pd.DataFrame(dict(CID=df.iloc[ev].CID.to_numpy(),zHD=z,observed_magnitude=obs,predicted_magnitude=pred,
    residual=r,role=np.where(z<.3,'feedback_fit','frozen_transfer_reused'))).to_csv(HERE/'feedback-brightness-predictions.csv',index=False,lineterminator='\n')
gal['feedback_predicted_z']=gp;gal['feedback_residual_kms']=c*(gp-gz)
gal.to_csv(HERE/'feedback-galaxy-predictions.csv',index=False,lineterminator='\n')
print(json.dumps(out,indent=2))
