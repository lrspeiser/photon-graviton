"""1D periodic Maxwell-scalar Hamiltonian dynamics; dimensionless mechanism test."""
from pathlib import Path
import json
import numpy as np
H=Path(__file__).resolve().parent
L=100.;g=.3;mass=.2;T=45.

def run(n,dt,two):
    dx=L/n;x=(np.arange(n)-n/2)*dx
    k=2*np.pi*np.fft.fftfreq(n,d=dx);k[n//2]=0
    def derivative(a):return np.fft.ifft(1j*k*np.fft.fft(a)).real
    first=np.exp(-.5*((x+15)/2)**2)
    second=np.exp(-.5*((x-15)/2)**2) if two else np.zeros(n)
    A=first+second;P=-derivative(first)+derivative(second)
    y=np.array([A,P,np.zeros(n),np.zeros(n)])
    def rhs(y):
        A,P,phi,pi=y;Z=np.exp(g*phi);Ax=derivative(A)
        phix=derivative(phi)
        return np.array([P/Z,derivative(Z*Ax),pi,
             derivative(phix)-mass*mass*phi+.5*g*(P*P/Z-Z*Ax*Ax)])
    def energy(y):
        A,P,phi,pi=y;Z=np.exp(g*phi)
        em=dx*np.sum(.5*(P*P/Z+Z*derivative(A)**2))
        scalar=dx*np.sum(.5*(pi*pi+derivative(phi)**2+mass*mass*phi*phi))
        return float(em),float(scalar)
    em0,sc0=energy(y);e0=em0+sc0
    def momentum(y):return float(-dx*np.sum(y[1]*derivative(y[0])+y[3]*derivative(y[2])))
    p0=momentum(y)
    max_drift=0.;max_momentum_drift=0.;records=[]
    steps=round(T/dt)
    for step in range(steps):
        k1=rhs(y);k2=rhs(y+dt*k1/2);k3=rhs(y+dt*k2/2);k4=rhs(y+dt*k3)
        y+=dt/6*(k1+2*k2+2*k3+k4)
        if step%100==99 or step==steps-1:
            em,sc=energy(y);max_drift=max(max_drift,abs((em+sc)/e0-1))
            max_momentum_drift=max(max_momentum_drift,abs(momentum(y)-p0)/e0)
            records.append({'time':(step+1)*dt,'electromagnetic_energy':em,'scalar_energy':sc,
                            'scalar_fraction_of_initial_energy':sc/e0})
    em,sc=energy(y)
    result={'n':n,'dt':dt,'counterpropagating_pulses':two,'initial_energy':e0,
        'final_electromagnetic_energy':em,'final_scalar_energy':sc,
        'final_scalar_fraction':sc/e0,'max_total_energy_relative_drift':max_drift,
        'initial_momentum':p0,'max_momentum_drift_over_initial_energy':max_momentum_drift,
        'max_abs_final_scalar_amplitude':float(np.max(np.abs(y[2]))),
        'maximum_scalar_fraction_during_run':max(r['scalar_fraction_of_initial_energy'] for r in records)}
    assert max_drift<1e-6 and max_momentum_drift<1e-6
    if not two:assert sc/e0<1e-20
    return result,records

cases=[];histories=[]
for two in [False,True]:
    for n,dt in [(512,.02),(1024,.01)]:
        result,history=run(n,dt,two);cases.append(result)
        histories.append({'n':n,'dt':dt,'counterpropagating_pulses':two,'history':history})
        print(json.dumps(result),flush=True)
collision=[r for r in cases if r['counterpropagating_pulses']]
refinement=abs(collision[0]['final_scalar_fraction']/collision[1]['final_scalar_fraction']-1)
assert refinement<1e-4 and collision[1]['final_scalar_fraction']>1e-6
summary={'classification':'Dimensionless field interaction demonstration; not astrophysical conversion, deposition or redshift validation',
    'parameters':{'box_length':L,'coupling_g':g,'scalar_mass':mass,'final_time':T,'pulse_width':2.,'initial_centers':[-15,15]},
    'Z':'exp(g phi), positive completion with same linear phi F^2 vertex',
    'cases':cases,'collision_fraction_relative_refinement_change':refinement,
    'identity_of_phi_with_lens_massive_source_established':False}
for file,obj in [('results.json',summary),('energy-histories.json',histories)]:
    (H/file).write_text(json.dumps(obj,indent=2)+'\n',newline='\n')
