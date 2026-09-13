"""Fit the constant rate using low-z brightness and transfer to other observables."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
data=ROOT/'shared_interaction_test/data'
source=data/'Pantheon+SH0ES.dat';covsource=data/'Pantheon+SH0ES_STAT+SYS.cov'
df=pd.read_csv(source,sep=r'\s+')
raw=np.loadtxt(covsource);N=int(raw[0]);C=raw[1:].reshape(N,N);C=(C+C.T)/2
cal=np.flatnonzero(df.IS_CALIBRATOR.to_numpy()==1)
ev=np.flatnonzero((df.IS_CALIBRATOR.to_numpy()==0)&(df.zHD.to_numpy()>=.1))
Cc=C[np.ix_(cal,cal)];Ce=C[np.ix_(ev,ev)];Cec=C[np.ix_(ev,cal)]
w=cho_solve(cho_factor(Cc),np.ones(len(cal)));w/=w.sum()
mu=df.CEPH_DIST.to_numpy()[cal];dc=10**((mu-25)/5);dbar=float(w@dc)
Mbase=float(w@(df.m_b_corr.to_numpy()[cal]-mu))
var=float(w@Cc@w);cross=Cec@w;V=Ce+var-cross[:,None]-cross[None,:]
z=df.zHD.to_numpy()[ev];zh=df.zHEL.to_numpy()[ev];obs=df.m_b_corr.to_numpy()[ev]
alpha0=.0002488993286382367;c=299792.458
def prediction(alpha):
    M=Mbase-5*alpha*dbar/np.log(10)
    return M+25+5*np.log10((1+zh)*np.log1p(z)/alpha)
base=prediction(alpha0);r0=obs-base
tr=np.flatnonzero(z<.3);te=np.flatnonzero(z>=.3)
vf=cho_factor(V[np.ix_(tr,tr)]);ones=np.ones(len(tr))
weights=cho_solve(vf,ones);weights/=weights.sum()
offset=float(weights@r0[tr]);sigma=float(np.sqrt(weights@V[np.ix_(tr,tr)]@weights))
def shift(alpha):return -5*(alpha-alpha0)*dbar/np.log(10)-5*np.log10(alpha/alpha0)
def invert(target):return brentq(lambda alpha:shift(alpha)-target,20/c,120/c,xtol=1e-16)
alpha=invert(offset);interval=sorted([c*invert(offset-sigma),c*invert(offset+sigma)])
pred=prediction(alpha);r=obs-pred
assert np.max(abs((pred-base)-offset))<1e-10
def score(res,ix):return float(res[ix]@cho_solve(cho_factor(V[np.ix_(ix,ix)]),res[ix]))
bins=[];B=[]
for lo,hi in [(.1,.3),(.3,.6),(.6,1),(1,3)]:
    ix=np.flatnonzero((z>=lo)&(z<hi));sub=V[np.ix_(ix,ix)]
    bw=cho_solve(cho_factor(sub),np.ones(len(ix)));bw/=bw.sum()
    full=np.zeros(len(ev));full[ix]=bw;B.append(full)
    bins.append(dict(z_min=lo,z_max=hi,count=len(ix),baseline_mean=float(bw@r0[ix]),revised_mean=float(bw@r[ix])))
B=np.array(B);bc=B@V@B.T;means=B@r
contrasts=[dict(bin_index=i,minus_first=float(means[i]-means[0]),
                formal_sigma=float(np.sqrt(bc[i,i]+bc[0,0]-2*bc[i,0]))) for i in range(1,4)]
gal_source=ROOT/'redshift_paper/all_164_groups.csv'
gal=pd.read_csv(gal_source);distance=gal.catalog_distance_mpc.to_numpy();gz=gal.observed_cmb_z.to_numpy()
gres0=c*(np.expm1(alpha0*distance)-gz);gres=c*(np.expm1(alpha*distance)-gz)
gs=[]
for split in ['train','validation','test']:
    ix=gal.split.to_numpy()==split
    gs.append(dict(historical_exposed_split=split,count=int(ix.sum()),
        original_rmse_kms=float(np.sqrt(np.mean(gres0[ix]**2))),
        brightness_rate_rmse_kms=float(np.sqrt(np.mean(gres[ix]**2))),
        brightness_rate_bias_kms=float(np.mean(gres[ix]))))
out=dict(status='Retrospective cross-observable comparison, no new blind outcomes',
    fitted_alpha_per_mpc=alpha,c_alpha_kms_mpc=c*alpha,formal_delta_chi2_one_c_alpha=interval,
    fitted_magnitude_offset=offset,formal_offset_sigma=sigma,
    original_c_alpha=c*alpha0,nearby_training_bootstrap95_c_alpha=[72.41722973559975,76.88738437794862],
    brightness_train_count=len(tr),brightness_transfer_count=len(te),
    train_baseline_chi2=score(r0,tr),train_fitted_chi2=score(r,tr),
    transfer_baseline_chi2=score(r0,te),transfer_frozen_point_chi2=score(r,te),
    scores_use_same_covariance_not_fitted_parameter_prediction_bands=True,
    bins=bins,invariant_bin_contrasts=contrasts,galaxy_results=gs,
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,covsource,gal_source,HERE/'rate-protocol.md']})
(HERE/'rate-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
pd.DataFrame(dict(CID=df.iloc[ev].CID.to_numpy(),zHD=z,observed_magnitude=obs,predicted_magnitude=pred,
    role=np.where(z<.3,'rate_fit','frozen_transfer_reused'),residual=r)).to_csv(HERE/'rate-brightness-predictions.csv',index=False,lineterminator='\n')
gal['brightness_rate_predicted_z']=np.expm1(alpha*distance);gal['brightness_rate_residual_kms']=gres
gal.to_csv(HERE/'rate-galaxy-predictions.csv',index=False,lineterminator='\n')
print(json.dumps(out,indent=2))
