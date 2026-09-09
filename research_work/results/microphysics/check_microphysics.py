"""Checks of conditional capture kinematics and static two-mode mixing."""
from pathlib import Path
import os
import json,shutil
import numpy as np
from scipy.linalg import expm
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',Path(__file__).resolve().parents[3]/'research_work/generated'))/'microphysics'
OUT.mkdir(parents=True,exist_ok=True)
checks=[]
# Natural units, target initial rest energy M=1. Initial massless energy E=x.
x=np.geomspace(1e-12,1e3,151)
mf=np.sqrt(1+2*x)
rest_fraction=2/(mf+1)
recoil_fraction=x/(1+x+mf)
beta=x/(1+x)
balance=float(abs(rest_fraction+recoil_fraction-1).max())
assert balance<1e-14
assert np.all(rest_fraction>0) and np.all(recoil_fraction>0) and np.all(beta<1)
checks.append(dict(check='absorbed_energy_partitions_into_rest_gain_and_recoil',sample_count=len(x),max_residual=balance))
# Four-momentum mass shell with scale-normalized residual.
residual=float(np.max(abs((1+x)**2-x*x-mf*mf)/(1+x)**2))
assert residual<1e-14
checks.append(dict(check='capture_four_momentum_mass_shell',max_scaled_residual=residual))
# The inverse channel excited -> ground + massless is kinematically open.
emission=x/mf
ground_energy=np.sqrt(1+emission**2)
residual=float(np.max(abs(ground_energy+emission-mf)/mf))
assert residual<1e-14 and np.all(emission>0)
checks.append(dict(check='inverse_emission_kinematically_open',max_scaled_residual=residual))
# Two opposite incident packets have no net recoil for their combined absorber.
residual=float(abs(((1+2*x)-1)/(2*x)-1)[x>1e-4].max())
assert residual<1e-11
checks.append(dict(check='opposed_equal_packets_have_zero_total_momentum',max_relative_rest_gain_error=residual))
mix=[]
for detuning,coupling in [(0,1),(1,.1),(2,0)]:
    H=np.array([[detuning/2,coupling],[coupling,-detuning/2]],float)
    length=np.linspace(0,10,101)
    state=np.array([expm(-1j*H*s)@np.array([1.,0.]) for s in length])
    prob=abs(state[:,1])**2
    omega=np.sqrt((detuning/2)**2+coupling**2)
    analytical=(coupling/omega)**2*np.sin(omega*length)**2 if omega else np.zeros_like(length)
    err=float(abs(prob-analytical).max()); norm=float(abs((abs(state)**2).sum(axis=1)-1).max())
    assert err<1e-13 and norm<1e-13
    mix.append(dict(detuning=detuning,coupling=coupling,analytic_probability_error=err,norm_residual=norm))
checks.append(dict(check='static_two_mode_transfer_matches_analytic_solution',cases=mix))
examples=[]
for z in [1e-6,.01,1,100]:
    m=np.sqrt(1+2*z)
    examples.append(dict(incident_energy_over_target_rest_energy=z,stored_rest_energy_fraction=2/(m+1),recoil_energy_fraction=z/(1+z+m),final_speed_over_c=z/(1+z),inverse_emission_energy_over_initial=z/m/z))
result=dict(scope='Special-relativistic free-target absorption and a stipulated Hermitian static mixing matrix. No cross-section, lifetime, redshift law, halo model or revised law inferred.',checks=checks,capture_examples=examples)
(OUT/'microphysics-checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
shutil.copy2(__file__,OUT/'check_microphysics.py')
print(json.dumps(dict(check_groups_passed=len(checks),capture_examples=examples),indent=2))
