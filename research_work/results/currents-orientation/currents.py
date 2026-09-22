#!/usr/bin/env python3
"""JR-6: 3D, source-driven companion-current toy. Not galaxy validation.

Usage: python code/currents.py --output results-new
The angular Fourier reduction solves a three-dimensional PDE, not a 2D source
sheet. Momentum, entrainment and a physical conversion interaction remain open.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass, asdict, replace
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.sparse import bmat, diags, eye, kron
from scipy.sparse.linalg import spsolve, expm_multiply

@dataclass(frozen=True)
class Setup:
    nr: int = 24
    nz: int = 32
    rmax: float = 6.
    zmax: float = 3.
    dp: float = .5
    dc: float = .05
    loss_p: float = .5
    loss_c: float = .1
    omega_s: float = .6
    omega_c: float = .6
    differential: bool = False
    epsilon: float = .6
    softening: float = .2

class Grid:
    def __init__(self, s: Setup):
        self.s = s
        self.rf = np.linspace(0, s.rmax, s.nr+1)
        self.r = (self.rf[:-1]+self.rf[1:])/2
        self.z = (np.arange(s.nz)+.5)*2*s.zmax/s.nz-s.zmax
        self.dr, self.dz = s.rmax/s.nr, 2*s.zmax/s.nz
        self.wr = np.diff(self.rf**2)/2
        self.weight = np.repeat(self.wr*self.dz, s.nz)
        self.R, self.Z = np.meshgrid(self.r, self.z, indexing='ij')
        self.size = s.nr*s.nz
        h = np.exp(-self.R/2)/np.cosh(self.Z/.35)**2
        h = h/(h+.15)
        self.kp = (.5*np.exp(2*h)).ravel()
        self.km = (.5*np.exp(4*h*h)).ravel()
        self.J0 = np.exp(-self.R**2/(2*1.5**2)-self.Z**2/(2*.35**2)).ravel()
        self.J0 /= 2*np.pi*(self.J0@self.weight)
        self.J2 = s.epsilon*(self.R**2/(self.R**2+.8**2)).ravel()*self.J0
        omega = s.omega_c/(1+self.R/2) if s.differential else np.full_like(self.R,s.omega_c)
        self.omega = omega.ravel()
    def lap(self, m: int):
        n = self.s.nr
        upper = self.rf[1:-1]/self.dr/self.wr[:-1]
        lower = self.rf[1:-1]/self.dr/self.wr[1:]
        center = -np.r_[upper,0.] - np.r_[0.,lower] - m*m/self.r**2
        lr = diags([lower,center,upper],[-1,0,1],shape=(n,n),format='csr')
        nz = self.s.nz
        center_z = -2*np.ones(nz);center_z[[0,-1]]=-1
        lz = diags([np.ones(nz-1),center_z,np.ones(nz-1)],[-1,0,1],format='csr')/self.dz**2
        return (kron(lr,eye(nz))+kron(eye(n),lz)).tocsr()
    def matrix(self,m:int):
        s = self.s;L = self.lap(m)
        a = -s.dp*L + diags(self.kp+s.loss_p-1j*m*s.omega_s)
        d = -s.dc*L + diags(self.km+s.loss_c+1j*m*(self.omega-s.omega_s))
        return bmat([[a,diags(-self.km)],[diags(-self.kp),d]],format='csc')

def solve(s:Setup):
    g=Grid(s);answers=[];res=[]
    for m,j in [(0,g.J0),(2,g.J2)]:
        A=g.matrix(m);rhs=np.r_[j,np.zeros(g.size)];y=spsolve(A,rhs)
        answers.append(y);res.append(float(np.linalg.norm(A@y-rhs)/max(np.linalg.norm(rhs),1e-300)))
    p0,c0=np.split(answers[0].real,2);p2,c2=np.split(answers[1],2)
    ep=2*np.pi*(p0@g.weight);ec=2*np.pi*(c0@g.weight)
    total_r2=2*np.pi*((c0*g.R.ravel()**2)@g.weight)
    quadrupole=np.pi*((np.conj(c2)*g.R.ravel()**2)@g.weight)
    normalized_q=abs(quadrupole)/total_r2
    conv=2*np.pi*((g.kp*p0-g.km*c0)@g.weight)
    diagnostics=dict(energy_P=float(ep),energy_C=float(ec),outgoing_power=float(s.loss_p*ep+s.loss_c*ec),
        balance_error=float(abs(s.loss_p*ep+s.loss_c*ec-1)),
        reaction_C_balance_error=float(abs(conv-s.loss_c*ec)),
        linear_residual_max=max(res),min_P_over_all_phi=float(np.min(p0-np.abs(p2))),
        min_C_over_all_phi=float(np.min(c0-np.abs(c2))),
        normalized_C_quadrupole=float(normalized_q),C_major_axis_degrees=float(np.angle(quadrupole)*90/np.pi),
        channel_azimuth_mean_independent_of_flow=True)
    return g,(p0,c0,p2,c2),diagnostics

def readout(g:Grid,fields,nphi:int=64):
    """Quadrature of one fixed Plummer potential; no lensing gain."""
    p0,c0,p2,c2=fields;s=g.s
    phi=(np.arange(nphi)+.5)*2*np.pi/nphi
    shape=(s.nr,s.nz,nphi)
    def volume(u0,u2):
        return u0.reshape(s.nr,s.nz,1)+np.real(u2.reshape(s.nr,s.nz,1)*np.exp(2j*phi))
    C=volume(c0,c2);P=volume(p0,p2)
    x=np.broadcast_to(g.R[...,None]*np.cos(phi),shape).ravel()
    y=np.broadcast_to(g.R[...,None]*np.sin(phi),shape).ravel()
    z=np.broadcast_to(g.Z[...,None],shape).ravel()
    w=np.broadcast_to(g.weight.reshape(s.nr,s.nz,1)*2*np.pi/nphi,shape).ravel()
    masses=np.vstack([C.ravel()*w,(P+C).ravel()*w])
    angles=np.linspace(0,2*np.pi,32,endpoint=False);radius=2.5
    rows=[]
    for a in angles:
        ca,sa=np.cos(a),np.sin(a);b=np.array([radius*ca,radius*sa,0.])
        dx,dy,dz=x-b[0],y-b[1],z
        inv=(dx*dx+dy*dy+dz*dz+s.softening**2)**-1.5
        acc=np.column_stack([masses@(dx*inv),masses@(dy*inv),masses@(dz*inv)])
        inv2=4/(dx*dx+dy*dy+s.softening**2)
        bend=np.column_stack([masses@(dx*inv2),masses@(dy*inv2)])
        row=dict(azimuth_degrees=float(a*180/np.pi))
        for i,key in enumerate(('C','both')):
            row[key]=dict(inward_force=float(-acc[i,0]*ca-acc[i,1]*sa),
                tangential_force=float(-acc[i,0]*sa+acc[i,1]*ca),
                inward_deflection=float(-bend[i,0]*ca-bend[i,1]*sa),
                tangential_deflection=float(-bend[i,0]*sa+bend[i,1]*ca))
        rows.append(row)
    stats={}
    for key in ('C','both'):
        fr=np.array([r[key]['inward_force'] for r in rows]);al=np.array([r[key]['inward_deflection'] for r in rows]);ft=np.array([r[key]['tangential_force'] for r in rows])
        stats[key]=dict(mean_inward_force=float(fr.mean()),force_semirange_over_mean=float(np.ptp(fr)/(2*fr.mean())),
            lens_semirange_over_mean=float(np.ptp(al)/(2*al.mean())),max_tangential_over_mean_radial=float(np.max(abs(ft))/fr.mean()))
    return dict(radius=radius,azimuths=rows,statistics=stats,nphi=nphi)

def scalar_controls():
    # Exact linear convective-memory response of a spatial harmonic.
    # Not a prediction of actual stellar or planet spin coupling.
    ratios=[0,.1,1,3,10,100]
    memory=[dict(m_times_relative_frequency_tau=x,amplitude=1/np.sqrt(1+x*x),phase_lag_degrees=np.degrees(np.arctan(x))) for x in ratios]
    # Distribution over independent emission-axis directions; angular rms from
    # sum of STF spin-axis tensors. No relation to coherent wave phase claimed.
    rng=np.random.default_rng(20260921);alignment=[]
    for N in (10,100,1000,10000):
        samples=[]
        for _ in range(100):
            n=rng.normal(size=(N,3));n/=np.linalg.norm(n,axis=1)[:,None]
            Q=n.T@n/N-np.eye(3)/3
            samples.append(float(np.linalg.norm(Q)/np.sqrt(2/3)))
        alignment.append(dict(independent_axes=N,mean_relative_anisotropy=float(np.mean(samples)),RMS_relative_anisotropy=float(np.sqrt(np.mean(np.array(samples)**2))),analytic_RMS=1/np.sqrt(N)))
    return dict(memory_transfer=memory,axis_averaging=alignment,full_alignment_relative_anisotropy=1.,
       spin_sign_invariance='n*n^T is exactly unchanged by n -> -n; spin sense needs a current or other odd-in-spin coupling')

def transient_control():
    s=Setup(nr=8,nz=10);g=Grid(s);A=g.matrix(2);b=np.r_[g.J2,np.zeros(g.size)];periodic=spsolve(A,b)
    # In the source-rotating frame dz/dt=-A z+b. The transient from zero is
    # z(T)=periodic-exp(-A*T)periodic. Propagate an augmented affine system
    # independently rather than using that solution expression.
    B=bmat([[-A,b[:,None]],[None,None]],format='csc') if False else None
    from scipy.sparse import csr_matrix
    M=bmat([[-A,csr_matrix(b[:,None])],[csr_matrix((1,len(b))),csr_matrix((1,1))]],format='csc')
    T=80.;initial=np.r_[np.zeros(len(b),complex),1.+0j]
    y=expm_multiply(M*T,initial,traceA=(-A.diagonal().sum())*T)
    return dict(nr=8,nz=10,time=T,relative_distance_to_periodic=float(np.linalg.norm(y[:-1]-periodic)/np.linalg.norm(periodic)),
                affine_constant_error=float(abs(y[-1]-1)),description='independent augmented matrix-exponential evolution from empty m=2 fields')

def main():
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    out=args.output
    if out.exists():raise FileExistsError('Use a new output directory')
    out.mkdir(parents=True);start=time.monotonic();base=Setup()
    cases=dict(aligned=base,stationary_companion=replace(base,omega_c=0),opposed=replace(base,omega_c=-.6),
        shear=replace(base,differential=True),both_reversed=replace(base,omega_s=-.6,omega_c=-.6),
        axisymmetric_aligned=replace(base,epsilon=0),axisymmetric_opposed=replace(base,epsilon=0,omega_c=-.6),
        shorter_residence=replace(base,loss_c=1),longer_residence=replace(base,loss_c=.02))
    output={}; stored={}
    for label,s in cases.items():
        g,f,d=solve(s);r=readout(g,f);stored[label]=(g,f)
        np.savez_compressed(out/(label+'.npz'),R=g.R,Z=g.Z,P0=f[0],C0=f[1],P2=f[2],C2=f[3],weight=g.weight)
        output[label]=dict(setup=asdict(s),diagnostics=d,readout=r)
        print(label,json.dumps(d),flush=True)
    refinements={}
    for nr,nz,np_ in ((48,64,128),(72,96,192)):
        for label in ('aligned','opposed'):
            s=replace(cases[label],nr=nr,nz=nz);g,f,d=solve(s);r=readout(g,f,np_)
            key=f'{label}_{nr}';refinements[key]=dict(setup=asdict(s),diagnostics=d,readout=r)
            np.savez_compressed(out/(key+'.npz'),R=g.R,Z=g.Z,P0=f[0],C0=f[1],P2=f[2],C2=f[3],weight=g.weight)
            print(key,json.dumps(d),flush=True)
    # Expanded boundary at the fine radial and vertical step sizes.
    s=replace(base,nr=64,nz=96,rmax=8,zmax=4.5);g,f,d=solve(s)
    refinements['aligned_boundary']=dict(setup=asdict(s),diagnostics=d,readout=readout(g,f,128))
    comparisons={}
    for label in ('aligned','opposed'):
        a=refinements[label+'_48'];b=refinements[label+'_72'];vals={}
        for channel in ('C','both'):
            for obs in ('inward_force','inward_deflection'):
                x=np.array([v[channel][obs] for v in a['readout']['azimuths']]);y=np.array([v[channel][obs] for v in b['readout']['azimuths']])
                vals[channel+'_'+obs+'_max_relative_change']=float(np.max(abs(x/y-1)))
        vals['quadrupole_relative_change']=float(a['diagnostics']['normalized_C_quadrupole']/b['diagnostics']['normalized_C_quadrupole']-1)
        comparisons[label]=vals
    fields=stored['aligned'][1];reverse=stored['both_reversed'][1]
    controls=dict(axisymmetric_flow_flip_max_difference=float(max(np.max(abs(a-b)) for a,b in zip(stored['axisymmetric_aligned'][1],stored['axisymmetric_opposed'][1]))),
        both_flow_reversal_conjugation_max_difference=float(max(np.max(abs(a-np.conj(b))) for a,b in zip(fields,reverse))),
        primary_total_C_spread=float(np.ptp([output[k]['diagnostics']['energy_C'] for k in ('aligned','opposed','stationary_companion','shear')])),
        transient=transient_control())
    summary=dict(experiment='JR-6',scope='3D periodic source-current solution; no new galaxy fitting or full gravitational backreaction',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=output,refinements=refinements,
        comparisons=comparisons,controls=controls,analytic_and_ensemble_controls=scalar_controls(),seconds=time.monotonic()-start)
    (out/'results.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('CONTROLS',json.dumps(controls),flush=True);print('REFINEMENT',json.dumps(comparisons),flush=True)
if __name__=='__main__':main()
