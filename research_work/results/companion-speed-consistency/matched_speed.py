"""Small-wave check of a proposed positive Hamiltonian, not a full optical fit."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent

def potential(n,dx):return np.sum(np.diff(np.r_[np.log(n),np.log(n[0])])**2)/(2*dx)
def force(n,dx):
    log=np.log(n)
    return (np.roll(log,-1)-2*log+np.roll(log,1))/(dx*dx*n)

def case(N,cells):
    L=8.;dx=L/cells;x=np.arange(cells)*dx;k=2*np.pi/L;mode=np.cos(k*x)
    initial=np.r_[N+1e-5*mode,np.zeros(cells)]
    discrete_omega=2*np.sin(k*dx/2)/(dx*N)
    def rhs(t,y):
        n,p=y[:cells],y[cells:]
        if min(n)<=0:raise ValueError('n must remain positive')
        return np.r_[p,force(n,dx)]
    end=.35*2*np.pi/discrete_omega
    sol=solve_ivp(rhs,[0,end],initial,dense_output=True,method='DOP853',rtol=1e-11,atol=1e-14,max_step=.5*dx*N)
    assert sol.success
    def amplitude(t):return (sol.sol(t)[:cells]-N)@mode*2/cells
    zero=brentq(amplitude,.2*2*np.pi/discrete_omega,.3*2*np.pi/discrete_omega,xtol=1e-10)
    measured=np.pi/(2*zero);energies=[]
    for t in np.linspace(0,end,101):
        y=sol.sol(t);energies.append(dx*np.sum(y[cells:]**2)/2+potential(y[:cells],dx))
    err=float(np.max(abs(np.array(energies)-energies[0]))/energies[0])
    assert err<1e-5 and abs(measured/discrete_omega-1)<1e-5
    return {'background_n':N,'cells':cells,'measured_phase_speed':float(measured/k),'continuum_photon_speed':1/N,'discrete_predicted_phase_speed':float(discrete_omega/k),'frequency_relative_error':float(measured/discrete_omega-1),'relative_energy_error':err}

# Independent numerical gradient of the discrete energy checks the nonlinear force.
rng=np.random.default_rng(329);n=1+.1*rng.random(32);dx=.25;epsilon=1e-6
numeric=[]
for i in range(len(n)):
    e=np.zeros(len(n));e[i]=epsilon
    numeric.append(-(potential(n+e,dx)-potential(n-e,dx))/(2*epsilon*dx))
gradient_error=float(max(abs(np.array(numeric)-force(n,dx))))
assert gradient_error<1e-7
rows=[case(N,cells) for N in [1.,1.2,2.] for cells in [128,256]]
result={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'proposed_gradient_energy':'K*c0^2/2 * integral (partial_x ln n)^2 dx','canonical_wave_equation':'n_tt = c0^2/n * partial_xx ln n, without photons; add the Hamiltonian photon source when coupled.',
 'linearized_speed':'c0/N about constant n=N>0','finite_difference_energy_gradient_max_error':gradient_error,'runs':rows,
 'scope':'Tests positive-energy variational consistency and small-wave speed only. No finite photon conversion, detector clocks, lossless transport in evolving backgrounds, capture or gravity prediction.'}
(HERE/'matched-speed-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
