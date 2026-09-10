"""Carrier-resolved extension of the closed gauge-kinetic pulse calculation.

Known Hamiltonian dynamics; chosen dimensionless initial conditions, no data fit.
"""
from pathlib import Path
import json
import numpy as np

HERE = Path(__file__).resolve().parent
L, T, WIDTH, G, MASS = 240., 80., 3., .3, .2

def run(n, dt, carrier, scalar_mass=MASS):
    dx = L/n
    x = (np.arange(n)-n/2)*dx
    k = 2*np.pi*np.fft.fftfreq(n, d=dx)
    kd = k.copy(); kd[n//2] = 0
    def deriv(a):
        return np.fft.ifft(1j*kd*np.fft.fft(a)).real
    def pulse(center):
        a = np.exp(-.5*((x-center)/WIDTH)**2)*np.cos(carrier*(x-center))
        return a / np.sqrt(2*dx*np.sum(deriv(a)**2))
    left, right = pulse(-15), pulse(15)
    y = np.array([left+right, -deriv(left)+deriv(right), np.zeros(n), np.zeros(n)])
    def rhs(y):
        a, p, phi, pi = y
        z = np.exp(G*phi); ax = deriv(a)
        return np.array([p/z, deriv(z*ax), pi,
            deriv(deriv(phi))-scalar_mass**2*phi+.5*G*(p*p/z-z*ax*ax)])
    def energies(y):
        a,p,phi,pi = y; z=np.exp(G*phi)
        return np.array([dx*np.sum(.5*(p*p/z+z*deriv(a)**2)),
            dx*np.sum(.5*(pi*pi+deriv(phi)**2+scalar_mass**2*phi**2))])
    initial = energies(y).sum()
    max_drift = 0.
    for step in range(round(T/dt)):
        a=rhs(y); b=rhs(y+dt*a/2); c=rhs(y+dt*b/2); d=rhs(y+dt*c)
        y += dt*(a+2*b+2*c+d)/6
        if step % 200 == 199:
            max_drift=max(max_drift, abs(energies(y).sum()/initial-1))
    final = energies(y)
    # Right outgoing wave, at x=65. Window excludes the opposite packet.
    # Compare to the exact freely translated initial packet with the same window.
    center = -15+T
    window = .5*(1-np.tanh((np.abs(x-center)-20)/2))
    free = np.fft.ifft(np.fft.fft(left)*np.exp(-1j*k*T)).real
    reference = -deriv(free)*window
    outgoing = .5*(y[1]/np.exp(G*y[2])-deriv(y[0]))*window
    positive = k>0
    freq=k[positive]
    def spectrum(field):
        return np.abs(np.fft.fft(field)[positive])**2
    p0,p1=spectrum(reference),spectrum(outgoing)
    norm0,norm1=p0.sum(),p1.sum()
    mean0=float(np.sum(freq*p0)/norm0); mean1=float(np.sum(freq*p1)/norm1)
    pdf0,pdf1=p0/norm0,p1/norm1
    # Wave-action proxy integral P(k)/k, meaningful for asymptotically free EM.
    action_ratio=float(np.sum(p1/freq)/np.sum(p0/freq))
    z_deviation=float(np.sqrt(np.sum(outgoing**2*(np.exp(G*y[2])-1)**2)/np.sum(outgoing**2)))
    result=dict(n=n,dt=dt,carrier=carrier,scalar_mass=scalar_mass,initial_energy=float(initial),
        final_scalar_fraction=float(final[1]/initial),
        max_total_energy_relative_drift=max_drift,
        outgoing_field_power_ratio=float(norm1/norm0),
        outgoing_wave_action_proxy_ratio=action_ratio,
        initial_mean_wavenumber=mean0,final_mean_wavenumber=mean1,
        mean_wavenumber_ratio=mean1/mean0,
        normalized_spectrum_total_variation=float(.5*np.abs(pdf1-pdf0).sum()),
        field_weighted_kinetic_coefficient_departure=z_deviation,
        spectrum_scope='Windowed right-moving classical field diagnostic, not a photon count or detector model')
    assert max_drift<2e-5
    return result,dict(carrier=carrier,scalar_mass=scalar_mass,k=freq.tolist(),
        initial_normalized_power=pdf0.tolist(),final_normalized_power=pdf1.tolist())

if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--resonant', action='store_true')
    args=parser.parse_args()
    if args.resonant:
        runs=[run(n,dt,1.,2.) for n,dt in [(2048,.01),(4096,.005)]]
        (HERE/'resonant-results.json').write_text(json.dumps(
            dict(cases=[row[0] for row in runs],spectrum=runs[-1][1]),indent=2)+'\n',newline='\n')
        print(json.dumps([row[0] for row in runs],indent=2),flush=True)
        raise SystemExit(0)
    cases=[]; spectra=[]
    for carrier in [1.,2.,4.]:
        for n,dt in [(2048,.01),(4096,.005)]:
            result,spectrum=run(n,dt,carrier)
            cases.append(result)
            if n==4096: spectra.append(spectrum)
            print(json.dumps(result),flush=True)
    result=dict(classification='Dimensionless carrier-resolved mechanism diagnostic; no observational fitting',
        parameters=dict(L=L,T=T,width=WIDTH,g=G,scalar_mass=MASS,
                        packet_energy=.5,initial_centers=[-15,15]),cases=cases)
    for name,obj in [('results.json',result),('spectra.json',spectra)]:
        (HERE/name).write_text(json.dumps(obj,indent=2)+'\n',newline='\n')
