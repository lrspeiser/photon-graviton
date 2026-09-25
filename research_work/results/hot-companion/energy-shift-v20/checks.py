#!/usr/bin/env python3
"""Bounded follow-up to the structured-medium energy-shift model.

Exact diagonalization of an effective fixed-excitation exchange Hamiltonian.
NOT a moving many-body field simulation, not an astronomical fit, and not a
proof of physical state preparation or a completed gravitational law.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

ROOT=Path(__file__).resolve().parent
C=0.04/(4*np.pi)
KAPPA=1.3537825795626388
OMEGA=1.2

def exchange(positions: np.ndarray, occupancy: int):
    n=len(positions)
    states=[s for s in range(1<<n) if s.bit_count()==occupancy]
    idx={s:i for i,s in enumerate(states)}
    H=np.zeros((len(states),len(states)))
    pairs=[]
    for i in range(n):
        for j in range(i+1,n):
            delta=positions[i]-positions[j]
            r=float(np.linalg.norm(delta))
            if r<=0: raise ValueError('Coincident positions')
            coupling=-C*np.exp(-KAPPA*r)/r
            rows=[]; cols=[]
            for a,s in enumerate(states):
                if ((s>>i)&1) != ((s>>j)&1):
                    b=idx[s^(1<<i)^(1<<j)]
                    H[a,b]+=coupling
                    rows.append(a); cols.append(b)
            pairs.append((i,j,delta,r,np.array(rows,int),np.array(cols,int)))
    vals,vecs=eigh(H)
    gs=vecs[:,0]
    F=np.zeros_like(positions)
    corr=[]
    for i,j,delta,r,rows,cols in pairs:
        expectation=float(np.dot(gs[rows],gs[cols])) if len(rows) else 0.
        derivative=C*np.exp(-KAPPA*r)*(1/r**2+KAPPA/r)
        fi=-derivative*expectation*delta/r
        F[i]+=fi; F[j]-=fi
        corr.append(expectation)
    # In a fixed-excitation sector sigma_minus changes N, so its expectation
    # vanishes. Also verify using an explicitly embedded full-Hilbert vector.
    full=np.zeros(1<<n,dtype=complex)
    full[np.array(states)]=gs
    coherent=[]
    for i in range(n):
        val=0j
        for s in range(1<<n):
            if (s>>i)&1: val+=np.conj(full[s^(1<<i)])*full[s]
        coherent.append(float(abs(val)))
    return dict(exchange_energy=float(vals[0]),total_energy=float(OMEGA*occupancy+vals[0]),
                corr_min=float(min(corr)),corr_max=float(max(corr)),
                rms_force=float(np.sqrt(np.mean(np.sum(F**2,axis=1)))),
                net_force=float(np.linalg.norm(F.sum(axis=0))),
                max_one_point_coherence=float(max(coherent))), F

def energy_only(M,r,kind):
    intensity=M/r**2
    if kind=='linear': return intensity
    if kind=='root': return np.sqrt(intensity)
    if kind=='log': return np.log(intensity)
    raise ValueError(kind)

def local_force(M,r,kind):
    # Magnitudes only; sign depends on actual microscopic coupling.
    if kind=='linear': return 2*M/r**3
    if kind=='root': return np.sqrt(M)/r**2
    if kind=='log': return 2/r
    raise ValueError(kind)

def main():
    rng=np.random.default_rng(250925)
    positions=np.array([[x,y,z] for x in (-.5,.5) for y in (-.5,.5) for z in (-.5,.5)],float)
    positions+=rng.normal(0,.035,positions.shape)
    rows=[]
    max_fd=0.
    for n in range(9):
        row,F=exchange(positions,n)
        row['excitations']=n
        row['fraction_excited']=n/8
        # independent positional energy derivative for every body's x coordinate
        for i in range(8):
            h=1e-5
            pp=positions.copy(); pm=positions.copy()
            pp[i,0]+=h; pm[i,0]-=h
            ep=exchange(pp,n)[0]['exchange_energy']
            em=exchange(pm,n)[0]['exchange_energy']
            fd=-(ep-em)/(2*h)
            max_fd=max(max_fd,abs(fd-F[i,0]))
        rows.append(row)
    reference=rows[1]['rms_force']
    for row in rows:
        row['rms_force_vs_one_excitation']=row['rms_force']/reference
    scalings=[]
    max_energy_fd=0.
    for kind in ('linear','root','log'):
        pM=math.log(local_force(4,2,kind)/local_force(1,2,kind))/math.log(4)
        pr=-math.log(local_force(1,4,kind)/local_force(1,2,kind))/math.log(2)
        h=1e-5; M=2.; r=3.
        numeric=abs((energy_only(M,r+h,kind)-energy_only(M,r-h,kind))/(2*h))
        max_energy_fd=max(max_energy_fd,abs(numeric-local_force(M,r,kind)))
        scalings.append(dict(local_energy=kind,mass_exponent=pM,distance_exponent=pr))
    # Positive fixed-weight Yukawa sum: r^2 |F| is nonincreasing.
    kvals=np.geomspace(1e-5,1e3,501)
    weights=rng.lognormal(0,2,len(kvals)); weights/=sum(weights)
    radii=np.geomspace(1e-3,1e3,201)
    x=kvals[:,None]*radii[None,:]
    forces=(weights[:,None]*np.exp(-x)*(1+x)/radii[None,:]**2).sum(axis=0)
    exponents=-np.diff(np.log(forces))/np.diff(np.log(radii))
    out=dict(scope='Exact fixed-sector equilibrium exchange check; no motion or bath dynamics integrated.',
             repository_reviewed='67075c2c59dbe3c7cf55af4eb1d9747bcf000ffb',
             parameters=dict(C=C,kappa=KAPPA,omega=OMEGA,positions=positions.tolist()),
             finite_population=rows,
             validation=dict(max_absolute_force_energy_derivative_error=max_fd,
                             max_net_force=max(r['net_force'] for r in rows),
                             max_one_point_coherence=max(r['max_one_point_coherence'] for r in rows),
                             energy_function_derivative_error=max_energy_fd,
                             min_positive_mixture_exponent=float(min(exponents))),
             local_energy_scalings=scalings,
             target=dict(mass_exponent=.5,distance_exponent=1.))
    dest=ROOT/'results.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['validation'],indent=2))
    for r in rows: print(r['excitations'],r['exchange_energy'],r['rms_force_vs_one_excitation'],r['corr_min'],r['corr_max'])
    print(json.dumps(scalings,indent=2))

if __name__=='__main__': main()
