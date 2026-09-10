"""Post-result stellar normalization diagnostic; see repair-amendment.md."""
import json
from pathlib import Path
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares

H=Path(__file__).resolve().parent; G=4.30091727003628e-6
data=json.loads((H.parent/'milky-way-capture/inputs.json').read_text())
initial=json.loads((H/'results.json').read_text()); gc=json.loads((H/'gas-cache.json').read_text())
er=data['eilers']['rows']; br=data['bovy']['rows']
R=np.array([x['R_kpc'] for x in er]); v=np.array([x['vc_kms'] for x in er]); rz=np.array([x['R_kpc'] for x in br]); z=1.1; radialz=np.sqrt(rz*rz+z*z)
k=np.array([x['Kz_over_2piG_Msun_pc2'] for x in br]); ke=np.array([x['Kz_error'] for x in br])
gasr=PchipInterpolator(gc['R_grid'],gc['radial_force'])(R); gasz=np.array(gc['vertical_force'])
sr=np.sqrt(np.array([.5*(x['err_minus_kms']+x['err_plus_kms']) for x in er])**2+(.03*v)**2)
cal=R<=15; zcal=np.arange(len(br))%2==0

@lru_cache(None)
def grid(high):
    rr=np.r_[0,np.geomspace(.0001,100,1799 if high else 899)]
    mu,w=leggauss(192 if high else 96); return rr,mu,w

def mass(theta,high=False):
    rr,mu,w=grid(high); r=rr[:,None]; beta=10**theta[1]; rs=10**(-.3)
    b2=r*r*(1-mu*mu); B=np.sqrt(rs*rs+b2)
    tau=beta*rs*rs/B*(np.arctan(r*mu/B)+np.arctan(np.sqrt(np.maximum(0,10000-b2))/B))
    rho=beta/(1+(rr/rs)**2)*.5*np.sum(w*np.exp(-tau),axis=1)
    m=cumulative_trapezoid(4*np.pi*rr*rr*rho,rr,initial=0)
    return PchipInterpolator(rr,m/m[-1]*10**theta[0])

out={}; rows=[]
for variant in ['I','II']:
    # Recover frozen stellar forces by subtracting gas from the original baryons-only predictions.
    allrows=json.loads((H/'milky-way-predictions.json').read_text())
    vr=np.array([x['predicted'] for x in allrows if x['baryons']==variant and x['model']=='baryons' and x['observable']=='vc_kms'])
    kz=np.array([x['predicted'] for x in allrows if x['baryons']==variant and x['model']=='baryons' and x['observable']!='vc_kms'])
    stellarR=vr*vr/R-gasr; stellarZ=kz*2*np.pi*G*1e6-gasz
    def prediction(p,high=False):
        M=mass(p,high)
        vp=np.sqrt(R*(p[2]*stellarR+gasr)+G*M(R)/R)
        kp=(p[2]*stellarZ+gasz+G*M(radialz)*z/radialz**3)/(2*np.pi*G*1e6)
        return vp,kp
    baseline=initial['milky_way'][variant]['capture_fits']['well']['log_parameters']
    check=prediction([baseline[0],baseline[1],1],True)
    prev=np.array([x['predicted'] for x in allrows if x['baryons']==variant and x['model']=='capture_well' and x['observable']=='vc_kms'])
    assert np.max(abs(check[0]-prev))<1e-9
    def residual(p):
        vp,kp=prediction(p)
        return np.r_[(vp[cal]-v[cal])/sr[cal],(kp[zcal]-k[zcal])/ke[zcal],(p[2]-1)/.2]
    fits=[least_squares(residual,[baseline[0],baseline[1],s],bounds=([9,-5,.5],[13,1,1.5]),xtol=1e-10,ftol=1e-10,gtol=1e-9,max_nfev=200) for s in [.7,1,1.3]]
    f=min(fits,key=lambda x:np.sum(x.fun*x.fun)); vp,kp=prediction(f.x,True); vl,kl=prediction(f.x)
    rms=lambda a:float(np.sqrt(np.mean(a*a)))
    out[variant]=dict(log_M100=float(f.x[0]),log_beta0=float(f.x[1]),stellar_mass_scale=float(f.x[2]),rs_fixed_kpc=float(10**(-.3)),optimizer_success=bool(f.success),
        radial_inner_RMSE_kms=rms(vp[cal]-v[cal]),radial_outer_RMSE_kms=rms(vp[~cal]-v[~cal]),
        vertical_calibration_diagonal_score=float(np.sum(((kp[zcal]-k[zcal])/ke[zcal])**2)),vertical_reserved_diagonal_score=float(np.sum(((kp[~zcal]-k[~zcal])/ke[~zcal])**2)),
        vertical_all_diagonal_score=float(np.sum(((kp-k)/ke)**2)),vertical_calibration_n=int(zcal.sum()),vertical_reserved_n=int((~zcal).sum()),
        max_resolution_change=float(max(np.max(abs(vp-vl)),np.max(abs(kp-kl)))))
    for i in range(len(R)): rows.append(dict(baryons=variant,observable='vc_kms',split='inner' if cal[i] else 'outer',R_kpc=float(R[i]),observed=float(v[i]),predicted=float(vp[i])))
    for i in range(len(rz)): rows.append(dict(baryons=variant,observable='Kz_over_2piG',split='calibration' if zcal[i] else 'reserved',R_kpc=float(rz[i]),observed=float(k[i]),predicted=float(kp[i])))
for name,value in [('repair-results.json',out),('repair-predictions.json',rows)]:
    (H/name).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))
