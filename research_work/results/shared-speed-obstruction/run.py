"""Conditional Hamiltonian obstruction and explicit proposed loophole checks."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
HERE=Path(__file__).resolve().parent

def main():
    momenta=[1.,2.,4.];rows=[]
    for nf in [1.01,1.3,2.]:
        a=(nf-1)/6
        for p in momenta:
            for target in [None,2.]:
                b=0 if target is None else target*(1-1/nf)
                ef=p/nf+b
                work=quad(lambda t:(-p+(0 if target is None else target))*a/(1+a*t)**2,0,6,epsabs=1e-12)[0]
                assert abs(ef-p-work)<1e-12
                rows.append({'final_n':nf,'momentum':p,'offset_tuned_momentum':target,'initial_energy':p,'final_energy':ef,'energy_fraction_retained':ef/p,'integrated_explicit_time_work':work,'final_speed':1/nf})
    # A momentum-dependent correction cancels the energy drift, but also the changed speed.
    n=1.3;p=2.;h=1e-5
    compensated=lambda p:p/n+p*(1-1/n)
    speed=(compensated(p+h)-compensated(p-h))/(2*h)
    assert abs(speed-1)<1e-10 and abs(speed-1/n)>.2
    # Common clock conversion removes the homogeneous frequency shift for both sectors.
    clock=[]
    for nf in [1.01,1.3,2.]:
        q=1/nf;reference=1/nf;local=reference/q
        assert local==1
        clock.append({'n':nf,'q':q,'photon_reference_energy_ratio':reference,'companion_reference_energy_ratio':reference,'photon_local_energy_ratio':local,'companion_local_energy_ratio':local})
    result={'conditional_assumptions':['Differentiable canonical ray Hamiltonian H_c(t,x,p).','One canonical momentum, with no extra evolving internal state.','Energy being retained is this Hamiltonian in a common reference standard.','Equal nondispersive speed dH_c/dp=c0/n(t,x) for a continuum of positive momenta.','Retention means dH_c/dt=0 along every free ray, not just one selected energy or endpoint.'],
      'derived_equations':['H_c=c0*p/n+B(t,x)','dH_c/dt=partial_t H_c=-c0*p*n_t/n^2+B_t','Subtracting two distinct p values forces n_t=0 if both energies are retained everywhere.'],
      'offset_cases':rows,'momentum_dependent_compensation':{'n':n,'restored_energy_speed':speed,'required_shared_speed':1/n},'common_clock_controls':clock,
      'scope':'Conditional mathematical restriction, not a general no-go theorem for companion theories, quantum interactions, internal energy storage or modified dynamics.'}
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('18 offset cases, 3 common-clock controls and momentum-dependent compensation checked.')

if __name__=='__main__':main()
