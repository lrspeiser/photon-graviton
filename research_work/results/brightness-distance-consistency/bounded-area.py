"""Fit one bounded area response and freeze for farther supernovae."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize_scalar
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
REGULAR='--regular' in sys.argv
basepath=ROOT/'shared_interaction_test/data'
source=basepath/'Pantheon+SH0ES.dat';covsource=basepath/'Pantheon+SH0ES_STAT+SYS.cov'
old=json.loads((HERE/'opacity-results.json').read_text())
for p in [source,covsource]:assert hashlib.sha256(p.read_bytes()).hexdigest()==old['hashes'][p.name]
df=pd.read_csv(source,sep=r'\s+');raw=np.loadtxt(covsource);N=int(raw[0]);cov=raw[1:].reshape(N,N)
assert len(df)==N and np.max(abs(cov-cov.T))<1e-7
cov=(cov+cov.T)/2
cal=np.flatnonzero(df.IS_CALIBRATOR.to_numpy()==1)
ev=np.flatnonzero((df.IS_CALIBRATOR.to_numpy()==0)&(df.zHD.to_numpy()>=.1))
Cc=cov[np.ix_(cal,cal)];Ce=cov[np.ix_(ev,ev)];Cec=cov[np.ix_(ev,cal)]
w=cho_solve(cho_factor(Cc),np.ones(len(cal)));w/=w.sum()
alpha=old['alpha_per_mpc'];mu=df.CEPH_DIST.to_numpy()[cal];Dcal=10**((mu-25)/5)
Kcal=2.5*alpha*Dcal/np.log(10);fcal=-np.expm1(-alpha*Dcal)
M0=float(w@(df.m_b_corr.to_numpy()[cal]-mu-2*Kcal))
Mvar=float(w@Cc@w);cross=Cec@w;V=Ce+Mvar-cross[:,None]-cross[None,:]
z=df.zHD.to_numpy()[ev];zh=df.zHEL.to_numpy()[ev];distance=np.log1p(z)/alpha
baseline=M0+25+5*np.log10(distance*(1+zh));f=1-1/(1+zh);obs=df.m_b_corr.to_numpy()[ev]
train=np.flatnonzero(z<.3);test=np.flatnonzero(z>=.3)
assert len(train)==466 and len(test)==494 and len(cal)==77
fac=cho_factor(V[np.ix_(train,train)]);testfac=cho_factor(V[np.ix_(test,test)])
def area(value,fraction):return 1+fraction/(1+value*fraction) if REGULAR else 1+value*fraction
def pred(eta):return baseline+2.5*np.log10(area(eta,f))-float(w@(2.5*np.log10(area(eta,fcal))))
def score(res,factor):return float(res@cho_solve(factor,res))
def loss(eta):return score((obs-pred(eta))[train],fac)
upper=100 if REGULAR else 10
opt=minimize_scalar(loss,bounds=(0,upper),method='bounded',options={'xatol':1e-10})
eta=min([0.,float(upper),float(opt.x)],key=loss);prediction=pred(eta)
baseline_score=score((obs-baseline)[train],fac)
assert abs(baseline_score-old['train_baseline_chi2'])<1e-7
assert abs(score((obs-baseline)[test],testfac)-old['transfer_baseline_chi2'])<1e-7
bins=[]
for lo,hi in [(.1,.3),(.3,.6),(.6,1),(1,3)]:
    ix=np.flatnonzero((z>=lo)&(z<hi));vc=V[np.ix_(ix,ix)]
    bw=cho_solve(cho_factor(vc),np.ones(len(ix)));bw/=bw.sum()
    bins.append(dict(z_min=lo,z_max=hi,n=len(ix),baseline_mean_residual=float(bw@(obs-baseline)[ix]),bounded_area_mean_residual=float(bw@(obs-prediction)[ix])))
if REGULAR:
    for row in bins:row['regular_area_mean_residual']=row.pop('bounded_area_mean_residual')
out=dict(status='Exposed standardized-data fit; postulated optical response, not a derived metric',boundary=eta in [0,upper],optimizer_success=bool(opt.success),alpha_per_mpc=alpha,b=1,train_baseline_score=baseline_score,train_revised_score=loss(eta),farther_baseline_score=score((obs-baseline)[test],testfac),farther_revised_score=score((obs-prediction)[test],testfac),M_calibrated=M0-float(w@(2.5*np.log10(area(eta,fcal)))),area_ratio_at_z1=area(eta,.5),maximum_area_ratio=area(eta,1.),bins=bins)
out['q' if REGULAR else 'eta']=eta
if REGULAR:out['observer_optical_focusing_per_Mpc2']=alpha**2*(3*eta+13/4)
rows=[dict(CID=str(df.iloc[i].CID),zHD=float(z[j]),zHEL=float(zh[j]),role='training' if z[j]<.3 else 'frozen_transfer_exposed',observed_magnitude=float(obs[j]),baseline_prediction=float(baseline[j]),revised_prediction=float(prediction[j])) for j,i in enumerate(ev)]
for fn,obj in [('bounded-area-results.json',out),('bounded-area-predictions.json',rows)]:
    (HERE/(fn.replace('bounded-area','regular-area') if REGULAR else fn)).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
