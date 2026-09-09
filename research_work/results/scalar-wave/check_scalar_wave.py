"""Closed scalar plus electromagnetic-mode truncation of Z(chi)F^2 theory.

Natural units; no charges/absorbers; flat prescribed metric. Two quadratures
represent one traveling Fourier mode. Not a halo or cosmological solution.
"""
from pathlib import Path
import os
import json,shutil
import numpy as np
from scipy.integrate import solve_ivp
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',Path(__file__).resolve().parents[3]/'research_work/generated'))/'scalar-wave'
OUT.mkdir(parents=True,exist_ok=True)
def run(g,m,k,A0,rtol):
    t=np.linspace(0,400,8001)
    def rhs(_,y):
        chi,v,a,b,ad,bd=y;Z=1+g*chi
        if Z<=0:raise ValueError('Photon kinetic coefficient crossed zero')
        return [v,.5*g*(ad*ad+bd*bd-k*k*(a*a+b*b))-m*m*chi,ad,bd,-g*v/Z*ad-k*k*a,-g*v/Z*bd-k*k*b]
    Z0=1+g
    # Leading adiabatic amplitude derivative initially zero because chi_dot=0.
    sol=solve_ivp(rhs,(0,400),[1,0,A0,0,0,-k*A0],t_eval=t,method='DOP853',rtol=rtol,atol=rtol*1e-3)
    assert sol.success
    chi,v,a,b,ad,bd=sol.y;Z=1+g*chi
    em=.5*Z*(ad**2+bd**2+k*k*(a*a+b*b))
    scalar=.5*(v*v+m*m*chi*chi)
    total=em+scalar
    # Exact phase derivative of complex A=a+ib, no finite difference required.
    freq=-(a*bd-b*ad)/(a*a+b*b)
    # Wronskian conserved by the mode equation, distinct from instantaneous EM energy.
    wronskian=Z*(a*bd-b*ad)
    return dict(time=t,chi=chi,Z=Z,em=em,scalar=scalar,total=total,frequency=freq,wronskian=wronskian)
cases=[]
for g,m in [(0,.03),(.2,.03),(.2,.3)]:
    coarse=run(g,m,1,1e-3,1e-9);fine=run(g,m,1,1e-3,1e-11)
    energy_error=float(np.max(abs(fine['total']/fine['total'][0]-1)))
    wronskian_error=float(np.max(abs(fine['wronskian']/fine['wronskian'][0]-1)))
    frequency_convergence=float(np.max(abs(coarse['frequency']-fine['frequency'])))
    assert energy_error<1e-8 and wronskian_error<1e-8 and frequency_convergence<1e-6
    if g==0:assert np.max(abs(fine['frequency']-1))<1e-10
    cases.append(dict(g=g,m_over_k=m,minimum_Z=float(fine['Z'].min()),maximum_Z=float(fine['Z'].max()),max_total_energy_relative_drift=energy_error,max_wronskian_relative_drift=wronskian_error,max_frequency_refinement_difference=frequency_convergence,frequency_over_k_range=[float(fine['frequency'].min()),float(fine['frequency'].max())],photon_mode_energy_ratio_range=[float((fine['em']/fine['em'][0]).min()),float((fine['em']/fine['em'][0]).max())],final_photon_mode_energy_ratio=float(fine['em'][-1]/fine['em'][0]),final_scalar_energy_change=float(fine['scalar'][-1]-fine['scalar'][0])))
    if m==.03 and g==.2:
        sample={key:value[::80].tolist() for key,value in fine.items()}
result=dict(scope=__doc__,cases=cases,sampled_adiabatic_case=sample,conclusion='Action-specific mode behavior, not a derived companion-particle source or irreversible redshift law; field energy and Wronskian checked separately')
(OUT/'scalar-wave-checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
shutil.copy2(__file__,OUT/Path(__file__).name)
print(json.dumps(cases,indent=2))
