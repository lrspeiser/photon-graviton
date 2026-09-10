"""Characteristic solution of the bounded chiral EM/reservoir Hamiltonian."""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def profile(x):
    u=np.clip(x/4,0,1);inside=(x>0)&(x<4)
    W=np.where(inside,np.sin(np.pi*u)**4,0.)
    Wx=np.where(inside,np.pi*np.sin(np.pi*u)**3*np.cos(np.pi*u),0.)
    return W,Wx

def simulate(a,b,N,tight=False):
    labels=np.linspace(-2,-1,N)
    F0=2*np.sin(np.pi*(labels+2))**2*np.cos(16*np.pi*(labels+2))
    e0=F0**2/2
    initial=np.array([labels,np.ones(N),F0,np.zeros(N)])
    def rhs(t,z):
        x,J,F,P=z.reshape(4,N);W,Wx=profile(x)
        n=1+(b+a*t)*W;v=1/n;vx=-(b+a*t)*Wx/n**2;vs=-a*W/n**2
        return np.array([v,vx*J,-vx*F,-vx*P-vs*F**2/2]).ravel()
    times=np.linspace(0,24,481)
    sol=solve_ivp(rhs,[0,24],initial.ravel(),method='DOP853',rtol=2e-12 if tight else 2e-10,
                  atol=2e-13 if tight else 2e-11,dense_output=True,t_eval=times,max_step=.04 if tight else .08)
    assert sol.success
    states=sol.y.T.reshape(len(times),4,N)
    x,J,F,P=states.transpose(1,0,2)
    W,_=profile(x);v=1/(1+(b+a*times[:,None])*W)
    ledger=J*(v*F**2/2+P)
    error=float(np.max(abs(ledger-e0)))
    amplitude_error=float(np.max(abs(J*F-F0)))
    assert error<1e-7 and amplitude_error<1e-7 and J.min()>0 and x[-1].min()>8
    photon=float(simpson(J[-1]*F[-1]**2/2,x=labels))
    reservoir=float(simpson(J[-1]*P[-1],x=labels))
    start=float(simpson(e0,x=labels))
    probes=[]
    for j in [0,N//2,N-1]:
        arrive=brentq(lambda t:sol.sol(t).reshape(4,N)[0,j]-8,0,24,xtol=1e-12)
        yy=sol.sol(arrive).reshape(4,N)
        probes.append(dict(initial_x=float(labels[j]),arrival_time=arrive,
            arrival_J=float(yy[1,j]),observed_frequency_ratio=float(1/yy[1,j])))
    # Boundary arrival derivative versus the independent spatial Jacobian.
    mid=N//2
    def arrival(j):return brentq(lambda t:sol.sol(t).reshape(4,N)[0,j]-8,0,24,xtol=1e-12)
    derivative=-(arrival(mid+1)-arrival(mid-1))/(labels[mid+1]-labels[mid-1])
    result=dict(a=a,static_amplitude=b,rays=N,tight=tight,
        initial_energy=start,photon_energy_fraction=photon/start,reservoir_gain_fraction=reservoir/start,
        maximum_characteristic_energy_error=error,maximum_wave_amplitude_identity_error=amplitude_error,
        minimum_reservoir_density=float(P.min()),minimum_J=float(J.min()),
        final_packet_position_range=[float(x[-1].min()),float(x[-1].max())],
        final_characteristics_in_conversion_region=bool(np.any((x[-1]>0)&(x[-1]<4))),
        arrival_probes=probes,central_arrival_interval_stretch=derivative,
        central_timing_vs_frequency_error=abs(derivative-probes[1]['arrival_J']))
    return result

def main():
    rows=[simulate(0,.4,129),simulate(.1,0,129),simulate(.1,0,257),simulate(.1,0,257,True)]
    assert abs(rows[0]['photon_energy_fraction']-1)<1e-8 and abs(rows[0]['reservoir_gain_fraction'])<1e-12
    assert abs(rows[-1]['photon_energy_fraction']-rows[-2]['photon_energy_fraction'])<1e-8
    assert rows[-1]['central_timing_vs_frequency_error']<1e-6
    exact_stretch=float(np.exp(.1*1.5))
    assert abs(rows[-1]['photon_energy_fraction']-1/exact_stretch)<1e-8
    assert abs(rows[-1]['central_arrival_interval_stretch']-exact_stretch)<1e-8
    files=[Path(__file__),HERE.parent/'directional-wave-reservoir/run.py',HERE.parent/'clock-ruler-completion/report.md']
    result=dict(scope='Bounded prescribed-profile EM/reservoir characteristic test; not a derived gravitational environment.',
        profile='W(x)=sin^4(pi*x/4) for 0<x<4, zero outside; n=1+(b+a*s)W.',
        clock_branch='s=t; P initially zero; positive/right-moving sector.',rows=rows,
        exact_integrated_profile=1.5,exact_dynamic_stretch=exact_stretch,
        endpoint_energy_flux_over_canonical_momentum_density=exact_stretch,
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        endpoint_atomic_coupling_derived=False,environment_momentum_accounted=False,
        total_matter_and_gravity_action_completed=False,holdouts_opened=False)
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
