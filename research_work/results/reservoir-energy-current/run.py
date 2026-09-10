"""Local transport test of the synchronized ideal packet-reservoir branch."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=HERE.parent/'internal-companion-reservoir/run.py'
A=.1;E0=.1;T=8.;L=12.

def density_average(edges,shift):
    # Unit-integral compact initial number density 2*cos^2(pi*(x-1.5)).
    y=np.clip(edges-shift,1,2)
    primitive=y+np.sin(2*np.pi*(y-1.5))/(2*np.pi)
    return np.diff(primitive)/np.diff(edges)

def simulate(nx):
    edges=np.linspace(0,L,nx+1);dx=L/nx
    number=density_average(edges,0)
    photon=number.copy();seed=E0*number;moving=np.zeros(nx);stationary=np.zeros(nx)
    initial=float(dx*np.sum(photon+seed));t=0.;outgoing=0.;crossing_moving=0.;crossing_static=0.
    face=int(round(5/dx));assert abs(edges[face]-5)<1e-12
    maximum_global_error=0.;maximum_control_error=0.;steps=0
    static_outgoing=0.
    while t<T:
        dt=min(.4*dx,T-t);mid=t+dt/2
        for start,stop in [(t,mid)]:
            updated=photon*(1+A*start)/(1+A*stop)
            transfer=photon-updated;photon=updated
            moving+=transfer;stationary+=transfer
        v=1/(1+A*mid)
        def advect(u):
            flux=np.r_[0.,v*u]
            return u-dt/dx*np.diff(flux),flux
        photon,fg=advect(photon);seed,fs=advect(seed);moving,fr=advect(moving)
        number,_=advect(number)
        outgoing+=dt*(fg[-1]+fs[-1]+fr[-1])
        static_outgoing+=dt*(fg[-1]+fs[-1])
        crossing_moving+=dt*(fg[face]+fs[face]+fr[face])
        crossing_static+=dt*(fg[face]+fs[face])
        updated=photon*(1+A*mid)/(1+A*(t+dt))
        transfer=photon-updated;photon=updated
        moving+=transfer;stationary+=transfer
        t+=dt;steps+=1
        energies=[dx*np.sum(photon+seed+moving)+outgoing,
                  dx*np.sum(photon+seed+stationary)+static_outgoing]
        maximum_global_error=max(maximum_global_error,max(abs(e-initial) for e in energies))
        control=[dx*np.sum((photon+seed+moving)[:face])+crossing_moving,
                 dx*np.sum((photon+seed+stationary)[:face])+crossing_static]
        maximum_control_error=max(maximum_control_error,max(abs(e-initial) for e in control))
    displacement=np.log1p(A*T)/A
    exact=density_average(edges,displacement)
    assert min(photon.min(),seed.min(),moving.min(),stationary.min())>=-1e-14
    assert maximum_global_error<1e-11 and maximum_control_error<1e-11
    # Uniform-time exact transfer fractions must hold in the moving branch.
    np.testing.assert_allclose(photon,number/(1+A*T),atol=1e-11)
    np.testing.assert_allclose(moving,number*(1-1/(1+A*T)),atol=1e-11)
    result=dict(cells=nx,steps=steps,dx=dx,initial_total_energy=initial,
        final_photon_energy=float(dx*photon.sum()),final_seed_energy=float(dx*seed.sum()),
        final_moving_reservoir_gain=float(dx*moving.sum()),
        final_stationary_reservoir_gain=float(dx*stationary.sum()),
        moving_gain_left_of_x5=float(dx*moving[:face].sum()),
        stationary_gain_left_of_x5=float(dx*stationary[:face].sum()),
        moving_total_energy_crossed_x5=float(crossing_moving),
        stationary_total_energy_crossed_x5=float(crossing_static),
        maximum_global_energy_error=maximum_global_error,
        maximum_control_volume_error=maximum_control_error,
        number_density_L1_error=float(dx*np.sum(abs(number-exact))),
        exact_packet_support=[1+displacement,2+displacement])
    exact_static=quad(lambda x:2*np.cos(np.pi*(x-1.5))**2*(1-np.exp(-A*(5-x))),1,2,
                      epsabs=1e-13,epsrel=1e-13)[0]
    result['exact_stationary_gain_left_of_x5']=exact_static
    result['stationary_control_error']=abs(result['stationary_gain_left_of_x5']-exact_static)
    return result

def main():
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),SOURCE]}
    runs=[simulate(n) for n in [480,960,1920]]
    assert all(runs[i+1]['number_density_L1_error']<runs[i]['number_density_L1_error'] for i in range(2))
    assert all(runs[i+1]['stationary_control_error']<runs[i]['stationary_control_error'] for i in range(2))
    result=dict(scope='Synthetic local continuity test for the stipulated ideal packet branch; no electromagnetic field theory or observational fit.',
        parameters=dict(a=A,seed_energy_per_packet=E0,initial_propagation_energy_per_packet=1,time=T,domain_length=L),
        units='Dimensionless diagnostic; no astronomical calibration.',runs=runs,source_hashes=hashes,
        seed_energy_included=True,creation_mechanism_derived=False,
        relativistic_stress_tensor_derived=False,holdouts_opened=False)
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(runs[-1],indent=2))

if __name__=='__main__':main()
