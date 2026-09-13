"""Conditional Coma shear-shape fit on inner bins and frozen outer prediction."""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar
HERE=Path(__file__).resolve().parent
raw=json.loads((HERE.parent/'cluster-observation-readiness/kubo-figure-data.json').read_text())
rr=raw['rows'];R=np.array([d['published_radius_h_inverse_Mpc'] for d in rr]);R=R/R[0]
y=np.array([d['shear_t'] for d in rr]);err=np.array([d['plotted_sigma_t'] for d in rr])
mu,w=leggauss(96);t,wt=leggauss(160);theta=(t+1)*np.pi/4;wt=wt*np.pi/4
results=[]
for T in [0.,100.]:
    x=np.geomspace(1e-7,1e8,16001)
    if T:
        u=x[:,None]*mu;B2=1+x[:,None]**2*(1-mu**2);B=np.sqrt(B2)
        tau=(2*T/np.pi)*(u/(2*B2*(B2+u*u))+(np.arctan(u/B)+np.pi/2)/(2*B**3))
        J=.5*np.sum(np.exp(-np.maximum(tau,0))*w,axis=1)
    else:J=np.ones_like(x)
    rho=J/(1+x*x)**2
    logrho=PchipInterpolator(np.log(x),np.log(np.maximum(rho,1e-300)))
    # Project along z=sqrt(1+b^2) tan(theta), resolving the central core too.
    b=np.geomspace(1e-5,1e5,6001)
    scale=np.sqrt(1+b*b)[:,None]
    rad=np.sqrt(b[:,None]**2+(scale*np.tan(theta))**2)
    dens=np.exp(logrho(np.log(rad)))
    Sigma=2*np.sum(dens*scale/np.cos(theta)**2*wt,axis=1)
    integ=Sigma[0]*b[0]**2/2+cumulative_trapezoid(Sigma*b,b,initial=0)
    delta=2*integ/b**2-Sigma
    F=PchipInterpolator(np.log(b),delta)
    def evaluate(loga):
        f=F(np.log(R/10**loga));amp=max(0.,float(np.sum(f[:3]*y[:3]/err[:3]**2)/np.sum((f[:3]/err[:3])**2)))
        pred=amp*f;loss=float(np.sum(((pred[:3]-y[:3])/err[:3])**2))
        return loss,amp,pred
    grid=np.linspace(-2,2,161);losses=[evaluate(v)[0] for v in grid];i=int(np.argmin(losses))
    opt=minimize_scalar(lambda a:evaluate(a)[0],bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),method='bounded',options={'xatol':1e-8})
    candidates=[(evaluate(v)[0],v) for v in [float(opt.x),-2.,2.]]
    _,loga=min(candidates);loss,amp,pred=evaluate(loga)
    results.append(dict(model='transparent' if T==0 else 'strong_interception',T=T,scale_in_first_radius_units=10**loga,positive_amplitude=amp,scale_boundary=abs(abs(loga)-2)<1e-5,training_diagonal_sum=loss,outer_diagonal_sum=float(np.sum(((pred[3:]-y[3:])/err[3:])**2)),predicted_shear=pred.tolist()))
    if T==0:
        # Exact projection for rho=(1+r^2)^-2.
        assert np.max(abs(Sigma/(np.pi/2*(1+b*b)**-1.5)-1))<1e-3
f=R**-2;amp=max(0.,float(np.sum(f[:3]*y[:3]/err[:3]**2)/np.sum((f[:3]/err[:3])**2)));pred=amp*f
results.append(dict(model='point_mass_shape',positive_amplitude=amp,training_diagonal_sum=float(np.sum(((pred[:3]-y[:3])/err[:3])**2)),outer_diagonal_sum=float(np.sum(((pred[3:]-y[3:])/err[3:])**2)),predicted_shear=pred.tolist()))
out=dict(scope='Shape-only exploratory transfer; not full cluster mass or geometry',observed=rr,radius_ratios=R.tolist(),models=results)
(HERE/'coma-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(results,indent=2))
