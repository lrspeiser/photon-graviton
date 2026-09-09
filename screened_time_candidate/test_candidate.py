"""Synthetic homogeneous-patch tests; not a galaxy or CMB data fit."""
import numpy as np
from scipy.linalg import eigh
import json
from pathlib import Path
p=Path(__file__).resolve().parent
k=7.731496595524618e-5 # per million light years
u=np.linspace(0,1,20001) # gamma*t
A=10.; eps=1e-4
results={}
for label,M in [('void',0.),('dense',1000.)]:
    B=A*M
    mat=np.array([[A*A+B*B,-A*A*np.sqrt(eps)],[-A*A*np.sqrt(eps),eps*A*A]])
    w2,Q=eigh(mat); w2=np.maximum(w2,0); w=np.sqrt(w2)
    s=1/(1+M*M)
    v=Q.T@np.array([s,1/np.sqrt(eps)])
    # sin(wu)/w handles the zero eigenfrequency continuously.
    pos=Q@(v[:,None]*u[None,:]*np.sinc(w[:,None]*u[None,:]/np.pi))
    vel=Q@(v[:,None]*np.cos(w[:,None]*u[None,:]))
    chi,psi=pos[0],pos[1]*np.sqrt(eps)
    dc,dp=vel[0],vel[1]*np.sqrt(eps)
    energy=.5*dc**2+.5*dp**2/eps+.5*A*A*(chi-psi)**2+.5*B*B*chi**2
    results[label]={'chi_at_u1':float(chi[-1]),'psi_at_u1':float(psi[-1]),'chi_over_psi_at_u1':float(chi[-1]/psi[-1]),'psi_dot_at_u1':float(dp[-1]),'maximum_relative_energy_error':float(np.max(np.abs(energy/energy[0]-1))),'minimum_mode_frequency_squared':float(w2[0])}
results['redshift_comparison_same_local_slope']=[{'R_Mly':R,'z_original_exponential':float(np.expm1(k*R)),'z_homogeneous_rolling_field':k*R} for R in [1,100,1000,10000]]
results['gamma_per_year']=k/1e6
results['one_over_gamma_Gyr']=1/(k*1000)
results['caveats']=['Dimensionless mass and inertia parameters are illustrative, not estimated from data.','Fixed-density homogeneous patches omit boundaries, matter dynamics, photon backreaction and gravity.','The two-field subsystem conserves energy; the full light-plus-field model is not yet closed.','Homogeneous rolling-field redshift is not a derived global prediction of the screened model.']
(p/'results.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
