"""One-dimensional EM chiral sectors with postulated advected clock reservoirs."""
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
L=20.;T=4.

def derivative(y,k):return np.fft.ifft(1j*k*np.fft.fft(y,axis=-1),axis=-1).real

def simulate(N,a):
    x=np.arange(N)*L/N;k=2*np.pi*np.fft.fftfreq(N,d=L/N)
    initial=np.array([np.exp(-((x-6)/.65)**2)*np.cos(2*np.pi*40*x/L),
                      .7*np.exp(-((x-14)/.8)**2)*np.cos(2*np.pi*60*x/L)])
    y=np.vstack([initial,np.zeros((2,N))]);dx=L/N
    energy0=float(dx*np.sum(initial**2/2));dt0=.2*dx;t=0.;max_error=0.;phases=[];times=[]
    def rhs(time,state):
        v=1/(1+a*time);dy=derivative(state,k)
        return np.array([-v*dy[0],v*dy[1],-v*dy[2]+a*v*v*state[0]**2/2,
                         v*dy[3]+a*v*v*state[1]**2/2])
    while t<T:
        dt=min(dt0,T-t)
        b1=rhs(t,y);b2=rhs(t+dt/2,y+dt*b1/2);b3=rhs(t+dt/2,y+dt*b2/2);b4=rhs(t+dt,y+dt*b3)
        y+=dt*(b1+2*b2+2*b3+b4)/6;t+=dt
        v=1/(1+a*t);energy=dx*np.sum(v*(y[0]**2+y[1]**2)/2+y[2]+y[3])
        max_error=max(max_error,abs(energy-energy0)/energy0)
        phases.append([np.angle(np.fft.fft(y[0])[40]),np.angle(np.fft.fft(y[1])[60])]);times.append(t)
    shift=np.log1p(a*T)/a if a else T;v=1/(1+a*T)
    exact=np.array([np.fft.ifft(np.fft.fft(initial[j])*np.exp((-1 if j==0 else 1)*1j*k*shift)).real for j in range(2)])
    reservoir_exact=(1-v)*exact**2/2
    # Each acquired-energy profile follows its own EM direction.
    current_error=float(np.max(abs(y[2:]-reservoir_exact)))
    phase=np.unwrap(np.array(phases),axis=0)
    omega=np.abs((phase[-1]-phase[-3])/(times[-1]-times[-3]))
    reference_time=(times[-1]+times[-3])/2
    expected=np.array([2*np.pi*40/L,2*np.pi*60/L])/(1+a*reference_time)
    result=dict(cells=N,a=a,time=T,initial_energy=energy0,
        propagation_energy_fraction=float(dx*np.sum(v*(y[0]**2+y[1]**2)/2)/energy0),
        acquired_reservoir_fraction=float(dx*np.sum(y[2]+y[3])/energy0),
        maximum_fractional_energy_error=float(max_error),
        amplitude_max_error=float(np.max(abs(y[:2]-exact))),
        directional_reservoir_max_error=current_error,
        reservoir_minimum=float(y[2:].min()),
        measured_frequency_over_initial=(omega/np.array([2*np.pi*40/L,2*np.pi*60/L])).tolist(),
        frequency_fractional_error=float(np.max(abs(omega/expected-1))),
        exact_displacement=shift)
    return result

def local_identity(N=1024):
    x=np.arange(N)*L/N;k=2*np.pi*np.fft.fftfreq(N,d=L/N);a=.1
    D=.4+.1*np.sin(2*np.pi*x/L);B=.2*np.cos(4*np.pi*x/L)
    s=np.array([.2+.03*np.sin(2*np.pi*x/L),.25+.02*np.cos(4*np.pi*x/L)])
    P=np.array([.3+.04*np.cos(2*np.pi*x/L),.2+.03*np.sin(4*np.pi*x/L)])
    sign=np.array([1.,-1.])[:,None];v=1/(1+a*s);vp=-a*v*v;sx=derivative(s,k)
    e=np.array([(D+B)**2/4,(D-B)**2/4])
    st=1-sign*v*sx;Pt=-sign*v*derivative(P,k)-vp*e
    Hd=.5*(v[0]*(D+B)+v[1]*(D-B));Hb=.5*(v[0]*(D+B)-v[1]*(D-B))
    Dt=-derivative(Hb,k);Bt=-derivative(Hd,k)
    et=np.array([.5*(D+B)*(Dt+Bt),.5*(D-B)*(Dt-Bt)])
    ht=np.sum(vp*st*e+v*et+Pt*st-P*sign*(vp*st*sx+v*derivative(st,k)),axis=0)
    J=np.sum(sign*v*(v*e+P*st),axis=0)
    residual=float(np.max(abs(ht+derivative(J,k))))
    assert residual<1e-10 and np.min(st)>0
    return residual

def main():
    runs=[simulate(N,a) for a in [0.,.1] for N in [512,1024,2048]]
    for a in [0.,.1]:
        r=[s for s in runs if s['a']==a]
        assert r[-1]['maximum_fractional_energy_error']<1e-6
        assert r[-1]['amplitude_max_error']<r[0]['amplitude_max_error']
        assert r[-1]['frequency_fractional_error']<1e-5
        assert r[-1]['directional_reservoir_max_error']<1e-6
    timing=[]
    for distance in [1.,3.,5.]:
        a=.1;factor=np.exp(a*distance)
        for birth in [0.,.5,1.]:
            def arrive(tb):
                return brentq(lambda ta:quad(lambda t:1/(1+a*t),tb,ta,epsabs=1e-12)[0]-distance,tb,tb+20,xtol=1e-13)
            first,second=arrive(birth),arrive(birth+.2)
            error=max(abs((second-first)/.2-factor),abs((1+a*first)/(1+a*birth)-factor))
            assert error<1e-10
            timing.append(dict(distance=distance,birth_time=birth,arrival_time=first,
                interval_ratio=(second-first)/.2,frequency_ratio=(1+a*birth)/(1+a*first),
                exact_stretch=factor,maximum_error=error))
    result=dict(scope='Hypothetical one-dimensional local Hamiltonian wave/reservoir candidate; no astronomical fit.',
        runs=runs,inhomogeneous_local_energy_identity_residual=local_identity(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        local_ordinary_light_speed_and_clocks_proved=False,three_dimensional_covariant_theory=False,
        quantum_stability_proved=False,holdouts_opened=False,boundary_ray_timing_controls=timing)
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
