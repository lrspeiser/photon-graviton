"""Finite-site protected-state kinetics with thermal backflow and energy counters."""
from pathlib import Path
import json,math
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
P=Path(__file__).resolve().parent

def matrix(A,B,kd,ku,eps,g):
 Q=np.array([[-A,B+g,0.],[A,-B-g-kd,ku],[0.,kd,-ku]])
 M=np.zeros((6,6));M[:3,:3]=Q
 M[3,:3]=[A,-B,0] # Net incoming photon energy, Ea=1.
 M[4,:3]=[0,eps*kd,-eps*ku] # Net energy to thermal bath.
 M[5,:3]=[0,g,0] # Radiative loss from bright state.
 return M

if __name__=='__main__':
 out=dict(scope='Alternative single-excitation site, not the prior unlimited bosonic pair mode; rates illustrative except thermal detailed balance',energy_unit='bright excitation energy Ea=1',time_unit='inverse illumination rate A=1',rows=[])
 A=1.;B=1.01**3;g=.001
 for eps in [.1,.5,.9]:
  for T in [.01,.1,1.]:
   for k in [.1,1.,10.]:
    n=1/math.expm1(eps/T);kd=k*(n+1);ku=k*n
    M=matrix(A,B,kd,ku,eps,g);initial=np.array([1.,0,0,0,0,0]);after=expm(M*20)@initial
    sol=solve_ivp(lambda t,y:M@y,[0,20],initial,rtol=1e-10,atol=1e-12)
    err=float(np.max(abs(sol.y[:,-1]-after)));assert err<1e-8
    dark=expm(matrix(0.,0.,kd,ku,eps,g)*1000)@after
    for state in [after,dark]:
     stored=state[1]+(1-eps)*state[2]
     assert abs(state[:3].sum()-1)<1e-10 and min(state[:3])>=-1e-12
     assert abs(state[3]-stored-state[4]-state[5])<1e-9
    f=math.exp(-eps/T);pb=1/(1+f*(1+(B+g)/A));pa=f*pb;p0=(B+g)/A*pa
    steady=np.array([p0,pa,pb]);assert np.max(abs(M[:3,:3]@steady))<1e-12
    total=g+kd+ku;slow=2*g*ku/(total+math.sqrt(total*total-4*g*ku))
    out['rows'].append(dict(relaxation_energy_fraction=eps,bath_temperature_over_Ea=T,relaxation_scale=k,up_over_down=ku/kd,population_after_pump=after[:3].tolist(),stored_energy_after_pump=float(after[1]+(1-eps)*after[2]),net_input_energy=float(after[3]),net_bath_energy_after_pump=float(after[4]),radiated_energy_after_pump=float(after[5]),population_after_dark=dark[:3].tolist(),stored_energy_after_dark=float(dark[1]+(1-eps)*dark[2]),steady_protected_probability=pb,steady_stored_energy=pa+(1-eps)*pb,dark_slow_decay_rate=slow,dark_slow_mode_half_time=math.log(2)/slow,independent_ODE_max_difference=err))
 (P/'protected-state-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps([r for r in out['rows'] if r['relaxation_scale']==1 and r['bath_temperature_over_Ea']==.1],indent=2))
