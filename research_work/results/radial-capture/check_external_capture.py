from pathlib import Path
import os
import json
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import cumulative_trapezoid,quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',Path(__file__).resolve().parents[3]/'research_work/generated'))/'radial-capture'
mu,w=roots_legendre(512)
x=np.linspace(0,1,4001)
records=[]
fig,ax=plt.subplots(figsize=(6.5,4.5))
for tau in [.001,.1,1,10]:
    # Uniform isotropic incident specific intensity, pure absorption, unit R.
    length=x[:,None]*mu+np.sqrt(1-x[:,None]**2+x[:,None]**2*mu**2)
    angular=np.exp(-tau*length)@w
    shell=2*tau*x*x*angular
    deposited=cumulative_trapezoid(shell,x,initial=0)
    # Independent incoming chord integral: impact parameter area measure 2b db.
    expected=quad(lambda b:2*b*(-np.expm1(-2*tau*np.sqrt(1-b*b))),0,1,epsabs=1e-12)[0]
    relative_error=abs(deposited[-1]/expected-1)
    assert relative_error<1e-4
    v=np.sqrt(deposited[1:]/x[1:]/deposited[-1])
    ax.plot(x[1:],v,label=f'capture depth = {tau:g}')
    records.append(dict(capture_depth=tau,absorbed_fraction=float(expected),radial_integral_relative_error=float(relative_error),half_radius_velocity_ratio=float(np.sqrt(deposited[2000]/.5/deposited[-1])),edge_to_center_deposition_density_ratio=float(angular[-1]/angular[0])))
ax.plot(x,x,'k--',alpha=.5,label='uniform deposit density')
ax.set(xlabel='Radius / outer radius R',ylabel='Deposited speed contribution / value at R',title='Uniform external illumination: an outer rise, not a flat profile')
ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.text(.5,.01,'Conditional spherical absorption model; no scattering or gravitational focusing.',ha='center',fontsize=8)
fig.tight_layout(rect=(0,.035,1,1));fig.savefig(OUT/'external-capture.png',dpi=170)
(OUT/'external-capture-results.json').write_text(json.dumps(dict(scope='Uniform external isotropic companion intensity; constant absorption in sphere; permanent stationary deposits; Newtonian deposited contribution only',checks=records),indent=2),encoding='utf-8')
print(json.dumps(records,indent=2))
