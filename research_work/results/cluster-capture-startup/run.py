"""Conditional finite-flight absorption; no astronomical inputs or fitted parameters."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
OUT=Path(__file__).resolve().parent

def volume(k,t,n):
    x,w=roots_legendre(n); r=(x+1)/2; wr=w/2
    mu,wm=roots_legendre(2*n)
    ell=r[:,None]*mu+np.sqrt(1-r[:,None]**2+r[:,None]**2*mu**2)
    q=np.exp(-k*ell)*np.maximum(t-ell,0)
    return float(2*k/t*np.sum(wr*r*r*(q@wm)))

def chord(k,t):
    def terms(y,which):
        length=2*y; m=min(length,t)
        if which==0:
            value=quad(lambda s:k*np.exp(-k*s)*(t-s),0,m,epsabs=1e-12)[0]
        elif which==1:
            value=-np.expm1(-k*m)/k
        else:
            value=max(t-length,0)*np.exp(-k*length)
        return 2*y*value/t
    return [quad(lambda y:terms(y,j),0,1,points=[t/2] if t<2 else None,epsabs=1e-11)[0] for j in range(3)]

rows=[]
for k in [.1,1,10]:
    steady=quad(lambda y:2*y*(-np.expm1(-2*k*y)),0,1)[0]
    for t in [.1,.5,1,2,10]:
        deposit,travel,escape=chord(k,t)
        coarse=volume(k,t,256); fine=volume(k,t,512)
        assert abs(deposit+travel+escape-1)<1e-10
        assert abs(fine-deposit)<2e-5
        assert abs(fine-coarse)<2e-5
        assert deposit<=steady+1e-10
        rows.append(dict(optical_depth=k,duration_in_R_over_c=t,deposited_fraction=deposit,
            traveling_fraction=travel,escaped_fraction=escape,steady_capture_fraction=steady,
            deposit_over_steady_estimate=deposit/steady,volume_chord_difference=abs(fine-deposit),
            refinement_difference=abs(fine-coarse),ledger_error=abs(deposit+travel+escape-1)))
result=dict(scope='Dimensionless sphere R=c=1, fixed opacity, boundary bath switched on at t=0, initially empty interior. Permanent stationary deposits. No observed cluster fit, support dynamics, focusing or photon production.',rows=rows)
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(rows,indent=2))
