"""Restricted unscreened model: derivations, published constraint, FIRAS diagnostic.
Run with numpy/scipy. No new blind data test or spacecraft residual fit.
"""
from pathlib import Path
import json, math, shutil
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import least_squares

P=Path(__file__).resolve().parent
kappa=0.000077315
gamma=kappa/1e6 # per Julian year, c = one light-year/year at reference epoch
T0=2.72548
branches=[]
for name,a,d in [('fixed atoms, changing magnetic response',0,0),
                 ('equal electric and magnetic response, fixed mass',1,0),
                 ('equal responses, compensating mass',1,2),
                 ('equal responses, fixed atomic length',1,1)]:
    p=2*a-d
    branches.append(dict(name=name,a=a,b=2-a,mass_exponent=d,clock_p=p,
        measured_redshift_exponent=1-p,alpha_exponent=1-a,
        bohr_length_exponent=a-d,alpha_drift_per_year=(1-a)*gamma))

# Published Lange et al. 2021 summary; conditional standard clock sensitivities.
mu,sigma=1e-18,1.1e-18
lo,hi=(mu-1.96*sigma)/gamma,(mu+1.96*sigma)/gamma

# Homogeneous radiation + field Hamiltonian, dimensionless K=1, V=0.
# H = v^2/2 + U0 exp(-chi). Radiation lost energy feeds field motion.
energy=[]
for U0 in [0.001,0.1]:
    sol=solve_ivp(lambda t,y:[y[1],U0*np.exp(-y[0])],[0,5],[0,0.2],
                  method='DOP853',rtol=1e-12,atol=1e-14,dense_output=True)
    chi,v=sol.sol(np.linspace(0,5,1001))
    rad=U0*np.exp(-chi); field=0.5*v*v; H=rad+field
    energy.append(dict(initial_radiation=U0,final_stretch=float(np.exp(chi[-1])),
        radiation_energy_change=float(rad[-1]-rad[0]),
        field_energy_change=float(field[-1]-field[0]),
        max_relative_total_energy_error=float(np.max(np.abs(H/H[0]-1))),
        final_field_velocity=float(v[-1])))

# Planck occupation preservation under omega -> omega/S.
x=np.geomspace(0.01,100,1000)
occupation=[]
for S in [2,10,3000/T0]:
    before=1/np.expm1(x)
    after=1/np.expm1((x/S)/(1/S))
    occupation.append(dict(stretch=S,max_relative_occupation_error=float(np.max(abs(after/before-1))),
        temperature_ratio=1/S,energy_density_ratio_fixed_coordinate_volume=1/S,
        photon_density_ratio_fixed_coordinate_volume=1))

# Real FIRAS residuals: fit one temperature and supplied Galaxy template.
if not (P/'firas.txt').exists():
    shutil.copyfile(P.parent/'temporal_candidate_audit/data/firas.txt',P/'firas.txt')
f=np.loadtxt(P/'firas.txt'); nu=f[:,0]*100*299792458.
res=f[:,2]/1000; sig=f[:,3]/1000; gal=f[:,4]/1000
def B(T):
    return 2*6.62607015e-34*nu**3/299792458.**2/np.expm1(6.62607015e-34*nu/(1.380649e-23*T))*1e20
def residual(v):
    return (B(2.725+v[0])-B(2.725)+v[1]*gal-res)/sig
fit=least_squares(residual,[0,0],xtol=1e-12,ftol=1e-12,gtol=1e-10)
Tout=2.725+fit.x[0]
degeneracy=[]
for source in [10,300,3000,6000]:
    R=math.log(source/Tout)/kappa
    predicted=source*math.exp(-kappa*R)
    degeneracy.append(dict(source_temperature_K=source,distance_Mly=R,
        predicted_temperature_K=predicted))
out=dict(status='Restricted model calculations and published-data diagnostics; not validation of time theory',
    kappa_per_Mly=kappa,gamma_per_year=gamma,branches=branches,
    clock_constraint=dict(source='https://arxiv.org/abs/2010.06620',
        alpha_drift_mean=mu,alpha_drift_sigma=sigma,
        conditional_95_interval_on_1_minus_a=[lo,hi],
        caveat='Assumes field rate gamma, this electromagnetic model, and applicable published sensitivity coefficients'),
    conservative_synthetic_energy=energy,planck_occupation=occupation,
    firas=dict(source='https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html',
        channels=len(f),fitted_reference_temperature_K=float(Tout),
        fitted_galaxy_amplitude=float(fit.x[1]),diagonal_chi_squared=float(sum(residual(fit.x)**2)),
        illustrative_source_distance_degeneracy=degeneracy,
        caveat='Reused public residual table; diagonal errors only; not an absolute recalibration or a new time-field detection'),
    example_3000K=dict(stretch=3000/T0,z=3000/T0-1,
        exponential_distance_Mly=math.log(3000/T0)/kappa,
        rolling_distance_Mly=(3000/T0-1)/kappa,
        caveat='3000 K is an arbitrary example, not an assumed Big Bang or inferred source temperature'))
(P/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
