"""Fixed-form, fixed-kappa Pantheon+ brightness test with Cepheid-only calibration.
Full released covariance including repeated-SN and calibrator correlations.
No high-redshift data used to fit kappa or primary absolute magnitude.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.integrate import quad

P=Path(__file__).resolve().parent; D=P/'data'
cfg=json.loads((P/'protocol.json').read_text())
df=pd.read_csv(D/'Pantheon+SH0ES.dat',sep=r'\s+')
raw=np.loadtxt(D/'Pantheon+SH0ES_STAT+SYS.cov')
N=int(raw[0]);C=raw[1:].reshape(N,N); assert len(df)==N
asym=float(np.max(abs(C-C.T))); assert asym<1e-7 # released decimal-rounding asymmetry
C=(C+C.T)/2
cal=np.flatnonzero(df.IS_CALIBRATOR.to_numpy()==1)
hi=np.flatnonzero((df.zHD.to_numpy()>=.1)&(df.IS_CALIBRATOR.to_numpy()==0))
Cc=C[np.ix_(cal,cal)];Ch=C[np.ix_(hi,hi)];Hc=C[np.ix_(hi,cal)]
fac=cho_factor(Cc);one=np.ones(len(cal));w=cho_solve(fac,one);w/=sum(w)
dcal=df.m_b_corr.to_numpy()[cal]-df.CEPH_DIST.to_numpy()[cal]
M=float(w@dcal);var=float(w@Cc@w)
cross=Hc@w
V=Ch+var-cross[:,None]-cross[None,:]
vf=cho_factor(V);ones=np.ones(len(hi));vi1=cho_solve(vf,ones)
# 1 parsec defined by 648000/pi AU; Julian year and c exact.
pc_m=149597870700.*648000/np.pi
mpc_in_mly=pc_m/(299792458.*31557600.)
kap=cfg['kappa_per_Mly']*mpc_in_mly
H0=299792.458*kap
z=df.zHD.to_numpy()[hi]; zh=df.zHEL.to_numpy()[hi]
DL=(1+zh)*np.log1p(z)/kap
DLhd=(1+z)*np.log1p(z)/kap
DLref=np.array([(1+h)*quad(lambda zz:1/np.sqrt(.3*(1+zz)**3+.7),0,zz,epsabs=1e-11)[0]/kap for zz,h in zip(z,zh)])
obs=df.m_b_corr.to_numpy()[hi]
def stats(dl):
    pred=M+25+5*np.log10(dl);r=obs-pred
    off=float(vi1@r/(ones@vi1));sig=float((ones@vi1)**-.5)
    chi=float(r@cho_solve(vf,r));sh=r-off
    return dict(chi_squared=chi,degrees_of_freedom=len(hi),
        rmse_mag=float(np.sqrt(np.mean(r*r))),mean_residual_mag=float(np.mean(r)),
        fitted_shape_only_offset_mag=off,offset_sigma_mag=sig,
        shape_only_chi_squared=float(sh@cho_solve(vf,sh)),shape_only_dof=len(hi)-1),r,pred
s,r,pred=stats(DL);sr,rr,pr=stats(DLref);sh,_,_=stats(DLhd)
edges=[.1,.3,.6,1.,3.];b=[];binweights=[]
for l,u in zip(edges[:-1],edges[1:]):
    j=np.flatnonzero((z>=l)&(z<u));sub=V[np.ix_(j,j)];f=cho_factor(sub);o=np.ones(len(j));ww=cho_solve(f,o);ww/=sum(ww)
    mean=float(ww@r[j]);sigma=float(np.sqrt(ww@sub@ww))
    bw=np.zeros(len(hi));bw[j]=ww;binweights.append(bw)
    b.append(dict(z_min=l,z_max=u,light_curves=len(j),unique_SN=int(df.iloc[hi[j]].CID.nunique()),
        mean_residual_mag=mean,mean_sigma_mag=sigma,
        implied_distance_ratio=10**(mean/5),
        mean_reference_residual_mag=float(ww@rr[j]),
        caveat='Bin means are correlated; GLS weights may be negative'))
W=np.array(binweights);bcov=W@V@W.T;bm=W@r
bf=cho_factor(bcov);bo=np.ones(len(b));bw=cho_solve(bf,bo);bw/=sum(bw)
bmean=float(bw@bm);bchi=float((bm-bmean)@cho_solve(bf,bm-bmean))
contrasts=[]
for j in range(1,len(b)):
    se=float(np.sqrt(bcov[j,j]+bcov[0,0]-2*bcov[j,0]))
    contrasts.append(dict(bin_index=j,minus_first_bin_mag=float(bm[j]-bm[0]),sigma_mag=se))

# Independent conditional-Gaussian prediction using calibrator residuals.
# This is a labeled covariance sensitivity check, not a change to primary protocol.
B=cho_solve(fac,Hc.T).T
q=np.ones(len(hi))-B@one
Vcond=Ch-B@Hc.T+var*np.outer(q,q)
rcond=r-B@(dcal-M)
condchi=float(rcond@cho_solve(cho_factor(Vcond),rcond))

out=dict(protocol=cfg,rows_total=N,calibrator_light_curves=len(cal),
    calibrator_unique_SN=int(df.iloc[cal].CID.nunique()),
    higher_redshift_light_curves=len(hi),higher_redshift_unique_SN=int(df.iloc[hi].CID.nunique()),
    z_range=[float(min(z)),float(max(z))],kappa_per_Mpc=kap,H0_equivalent_km_s_Mpc=H0,
    M_calibrator_only=M,M_calibrator_sigma=var**.5,
    calibrator_chi_squared=float((dcal-M)@cho_solve(fac,dcal-M)),calibrator_dof=len(cal)-1,
    fixed_temporal=s,fixed_reference=sr,redshift_bins=b,
    bin_trend=dict(covariance_mag_squared=bcov.tolist(),constant_bin_mean_chi_squared=bchi,
        degrees_of_freedom=3,contrasts_to_first_bin=contrasts),
    frame_sensitivity=dict(max_prediction_difference_mag=float(max(abs(5*np.log10(DLhd/DL)))),
        fixed_temporal_using_zHD_only=sh),
    conditional_calibration_sensitivity_chi_squared=condchi,
    covariance_max_asymmetry=asym,
    caveats=['Standardized published magnitudes and covariance inherit SALT2, selection, host, dust and fiducial-cosmology assumptions.',
      'Primary absolute magnitude calibrated only on Cepheid hosts; these distances have their own ladder assumptions.',
      'Some SN objects may overlap earlier duration studies; different observable, not a wholly independent sky sample.',
      'No full light-curve refit or new joint cosmological inference.',
      'Fixed reference is not best-fit LCDM and its absolute residual includes the deliberately fixed kappa scale.'])
(P/'brightness_results.json').write_text(json.dumps(out,indent=2)+'\n')
pd.DataFrame(dict(row_index=hi,CID=df.iloc[hi].CID.to_numpy(),survey=df.iloc[hi].IDSURVEY.to_numpy(),zHD=z,zHEL=zh,
    m_b_corr=obs,predicted_fixed_temporal=pred,residual_fixed_temporal=r,
    predicted_fixed_reference=pr,residual_fixed_reference=rr)).to_csv(P/'brightness_predictions.csv',index=False)
print(json.dumps(out,indent=2))
