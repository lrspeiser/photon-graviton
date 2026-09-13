"""Fit one optional photon-removal coefficient, then transfer across redshift groups."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
data=ROOT/'shared_interaction_test/data'
source=data/'Pantheon+SH0ES.dat'
covsource=data/'Pantheon+SH0ES_STAT+SYS.cov'
df=pd.read_csv(source,sep=r'\s+')
raw=np.loadtxt(covsource);N=int(raw[0]);C=raw[1:].reshape(N,N)
assert len(df)==N and np.max(abs(C-C.T))<1e-7
C=(C+C.T)/2
cal=np.flatnonzero(df.IS_CALIBRATOR.to_numpy()==1)
ev=np.flatnonzero((df.IS_CALIBRATOR.to_numpy()==0)&(df.zHD.to_numpy()>=.1))
Cc=C[np.ix_(cal,cal)];Ce=C[np.ix_(ev,ev)];Cec=C[np.ix_(ev,cal)]
w=cho_solve(cho_factor(Cc),np.ones(len(cal)));w/=w.sum()
alpha=.0002488993286382367
mu=df.CEPH_DIST.to_numpy()[cal]
Dcal=10**((mu-25)/5)
Kcal=2.5*alpha*Dcal/np.log(10)
M0=float(w@(df.m_b_corr.to_numpy()[cal]-mu-2*Kcal))
Kbar=float(w@Kcal)
Mvar=float(w@Cc@w);cross=Cec@w
V=Ce+Mvar-cross[:,None]-cross[None,:]
z=df.zHD.to_numpy()[ev];zh=df.zHEL.to_numpy()[ev]
distance=np.log1p(z)/alpha
base=M0+25+5*np.log10(distance*(1+zh))
X=2.5*np.log10(1+zh)-Kbar
obs=df.m_b_corr.to_numpy()[ev];r0=obs-base
train=np.flatnonzero(z<.3);test=np.flatnonzero(z>=.3)
Vtt=V[np.ix_(train,train)];fac=cho_factor(Vtt)
wx=cho_solve(fac,X[train]);den=float(X[train]@wx);fitw=wx/den
epsilon_raw=float(fitw@r0[train]);epsilon=max(0.,epsilon_raw)
assert epsilon_raw>0, 'Boundary fit requires a different uncertainty treatment'
variance=float(fitw@Vtt@fitw)
pred=base+epsilon*X;res=obs-pred
Vhh=V[np.ix_(test,test)];q=V[np.ix_(test,train)]@fitw
Vpred=Vhh+variance*np.outer(X[test],X[test])-np.outer(q,X[test])-np.outer(X[test],q)
def score(r,cov):return float(r@cho_solve(cho_factor(cov),r))
bins=[]
for lo,hi in [(.1,.3),(.3,.6),(.6,1),(1,3)]:
    ix=np.flatnonzero((z>=lo)&(z<hi));vc=V[np.ix_(ix,ix)]
    bw=cho_solve(cho_factor(vc),np.ones(len(ix)));bw/=bw.sum()
    bins.append(dict(z_min=lo,z_max=hi,count=len(ix),
        baseline_gls_mean_residual=float(bw@r0[ix]),
        revised_gls_mean_residual=float(bw@res[ix]),
        note='Correlated descriptive bin means; no independent-bin likelihood'))
out=dict(status='Retrospective standardized-magnitude revision, not a blind or raw-data test',
    alpha_per_mpc=alpha,b_fixed=1,calibrator_rows=len(cal),calibrator_unique_SN=int(df.iloc[cal].CID.nunique()),
    training_rows=len(train),transfer_test_rows=len(test),
    epsilon=epsilon,epsilon_formal_sigma=float(np.sqrt(variance)),
    opacity_per_mpc=epsilon*alpha,M_at_zero_opacity=M0,
    M_at_fitted_opacity=M0-epsilon*Kbar,
    train_baseline_chi2=score(r0[train],Vtt),train_revised_chi2=score(res[train],Vtt),
    transfer_baseline_chi2=score(r0[test],Vhh),
    transfer_revised_chi2_with_fit_cross_covariance=score(res[test],Vpred),
    transfer_revised_chi2_on_baseline_covariance=score(res[test],Vhh),
    bins=bins,survival_fraction_at_z1=2**(-epsilon),
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,covsource,HERE/'opacity-protocol.md']})
(HERE/'opacity-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
pd.DataFrame(dict(CID=df.iloc[ev].CID.to_numpy(),zHD=z,zHEL=zh,
    role=np.where(z<.3,'opacity_fit','frozen_transfer_reused'),observed_magnitude=obs,
    baseline_prediction=base,revised_prediction=pred,baseline_residual=r0,revised_residual=res)).to_csv(HERE/'opacity-predictions.csv',index=False,lineterminator='\n')
print(json.dumps(out,indent=2))
