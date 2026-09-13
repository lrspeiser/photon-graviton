"""Homogeneous shell supply with competing receivers and exact energy ledger."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm

OUT=Path(__file__).resolve().parent

def ray_companion(x,b):
    if b==1: return x*np.exp(-x)
    return (np.exp(-x)-np.exp(-b*x))/(b-1)

def bath(x,b):
    photon=-np.expm1(-x)
    if b==0: companion=x-photon
    elif b==1: companion=-np.expm1(-x)-x*np.exp(-x)
    else: companion=(photon-(-np.expm1(-b*x))/b)/(b-1)
    deposited=x-photon-companion
    return np.array([photon,companion,deposited])

rows=[]
for b in [0.,.01,.1,1.,10.]:
    for x in [.1,1.,10.,100.,1000.]:
        u=bath(x,b)
        generator=np.array([[-1,0,0,1],[1,-b,0,0],[0,b,0,0],[0,0,0,0]],float)
        independent=(expm(generator*x)@np.array([0.,0.,0.,1.]))[:3]
        shell=quad(lambda s:ray_companion(s,b),0,x,epsabs=1e-10,epsrel=1e-10)[0]
        err=max(abs(u-independent))/max(x,1)
        assert err<1e-10 and abs(shell-u[1])/max(x,1)<1e-10
        assert min(u)>-1e-10 and abs(sum(u)-x)<1e-10*max(1,x)
        rows.append(dict(beta_over_alpha=b,alpha_c_T=x,photon_energy=float(u[0]),companion_energy=float(u[1]),deposited_energy=float(u[2]),
            energy_unit='luminosity_density/(alpha*c)',instantaneous_capture_power_over_emitted_power=float(b*u[1]),
            matrix_relative_error=float(err),shell_integral_absolute_error=float(abs(shell-u[1]))))

capture=[]
for tau in [.1,1.,10.]:
    analytic=1-(1-(1+2*tau)*np.exp(-2*tau))/(2*tau*tau)
    numeric=quad(lambda y:2*y*(-np.expm1(-2*tau*y)),0,1,epsabs=1e-13)[0]
    assert abs(analytic-numeric)<1e-12
    capture.append(dict(kappa_R=tau,absorbed_fraction=analytic,cross_section_over_R_squared=float(np.pi*analytic)))

# Illustrative populations: area competition fixes allocation, independently of total emissivity.
number_density=np.array([1e-5,2e-6]);radii=np.array([1.,2.]);tau=np.array([1.,10.])
fractions=1-(1-(1+2*tau)*np.exp(-2*tau))/(2*tau*tau)
sigma=np.pi*radii**2*fractions;beta=float(number_density@sigma)
per_receiver=sigma/beta;allocation=number_density*per_receiver
assert abs(sum(allocation)-1)<1e-14
result=dict(scope='Analytic homogeneous comparison, not measured cosmic luminosity density or actual capture population. No fixed age or size. No spectral/time-stretch solution or lensing fit.',
    rows=rows,capture_cross_sections=capture,
    illustrative_populations=dict(number_density_per_Mpc3=number_density.tolist(),radius_Mpc=radii.tolist(),effective_cross_section_Mpc2=sigma.tolist(),
        beta_per_Mpc=beta,mean_free_path_Mpc=1/beta,steady_power_per_receiver_over_luminosity_density_Mpc3=per_receiver.tolist(),fraction_of_all_emission_allocated=allocation.tolist()))
(OUT/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(dict(cases=len(rows),max_matrix_error=max(r['matrix_relative_error'] for r in rows),capture=capture,populations=result['illustrative_populations']),indent=2))
