#!/usr/bin/env python3
"""Shared Gaussian one-loop response of the completed area/stress candidate.

Quartic partners are retained in the exact Fock definition and small Fock test;
they are NOT included in this one-loop Gaussian response. All species, mirrors,
slab bands and inherited coefficients enter both source calculations identically.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import numpy as np
import area_stress_model as m
sys.path.insert(0,str(m.ROOT/'validation'))
from check_shared_matter_ir import bubble, contact


def lift_spectrum(L,width,charge):
    momenta=2*np.pi*(np.arange(L)+.5)/L-np.pi
    p=np.stack(np.meshgrid(momenta,momenta,momenta,indexing='ij'),axis=-1).reshape(-1,3)
    gamma=[np.kron(np.eye(width),g) for g in m.GAMMA]
    beta=np.kron(np.eye(width),m.BETA)
    H=np.broadcast_to(m.slab.slab_hamiltonian((0,0,0),m.slab.wilson_mass(charge),width),(len(p),4*width,4*width)).copy()
    for i in range(3):
        H+=np.sin(p[:,i])[:,None,None]*gamma[i]+(1-np.cos(p[:,i]))[:,None,None]*beta
    residual=float(np.max(abs(H[0]-m.slab.slab_hamiltonian(p[0],m.slab.wilson_mass(charge),width))))
    if residual>1.e-12:
        raise ArithmeticError('Bloch source reconstruction mismatch')
    wall=m.slab.slab_hamiltonian((0,0,0),m.slab.wilson_mass(charge),width)
    O=wall+3*beta
    C=.5*O@O+3*np.eye(len(O))
    H=m.X_ROOT*H+m.X_ROOT**2*C
    E,V=np.linalg.eigh(H)
    return p,H,E,V,gamma,beta,O,C,residual


def shear_bloch(p,width,charge,s):
    """Independent nonlinear definition used to verify source derivatives."""
    wall=m.slab.slab_hamiltonian((0,0,0),m.slab.wilson_mass(charge),width)
    gamma=[np.kron(np.eye(width),g) for g in m.GAMMA];beta=np.kron(np.eye(width),m.BETA)
    a=np.exp(s*np.array([0.,1.,-1.])/np.sqrt(2.))
    H=np.broadcast_to(wall,(len(p),*wall.shape)).copy()
    for i in range(3):
        H+=a[i]*np.sin(p[:,i])[:,None,None]*gamma[i]+a[i]**2*(1-np.cos(p[:,i]))[:,None,None]*beta
    O=wall+np.sum(a*a)*beta
    C=.5*O@O+.5*np.sum(a*a+a**4)*np.eye(len(O))
    return m.X_ROOT*H+m.X_ROOT**2*C


def one_volume(L,width,max_mode,charge_groups=None):
    names=('photon','frame')
    sums={n:{'photon_static':0.,'photon_temporal_Euclidean':0.,'frame_static':0.,'frame_temporal_Euclidean':0.,
             'frame_bubble':0.,'frame_contact_kinetic':0.,'frame_contact_regulator':0.,'frame_contact_completion':0.,
             'longitudinal_photon_static':0.} for n in range(max_mode+1)}
    source_errors=[];gaps=[];Ward=0.
    for charge,multiplicity in (Counter(abs(q) for q in m.slab.CHARGES) if charge_groups is None else charge_groups).items():
        p,H,E,V,gamma,beta,O,C,err=lift_spectrum(L,width,charge)
        half=2*width
        gaps.append({'abs_charge':charge,'indirect_filled_band_gap':float(np.min(E[:,half])-np.max(E[:,half-1])),
                     'maximum_occupied_energy':float(np.max(E[:,half-1])),
                     'minimum_empty_energy':float(np.min(E[:,half]))})
        assert gaps[-1]['indirect_filled_band_gap']>0
        sy=np.sin(p[:,1])[:,None,None]*gamma[1];sz=np.sin(p[:,2])[:,None,None]*gamma[2]
        wy=(1-np.cos(p[:,1]))[:,None,None]*beta;wz=(1-np.cos(p[:,2]))[:,None,None]*beta
        Jg=m.X_ROOT*((sy-sz)/np.sqrt(2.)+np.sqrt(2.)*(wy-wz))
        Ckin=m.X_ROOT*.5*(sy+sz)
        Creg=m.X_ROOT*2*(wy+wz)
        Ccomp=np.broadcast_to(m.X_ROOT**2*(2*(O@beta+beta@O)+10*np.eye(len(O))),H.shape)
        Cg=Ckin+Creg+Ccomp
        JA=m.X_ROOT*charge*(np.cos(p[:,1])[:,None,None]*gamma[1]+np.sin(p[:,1])[:,None,None]*beta)
        CA=m.X_ROOT*charge**2*(-sy+np.cos(p[:,1])[:,None,None]*beta)
        CL=m.X_ROOT*charge**2*(-np.sin(p[:,0])[:,None,None]*gamma[0]+np.cos(p[:,0])[:,None,None]*beta)
        cp=contact(V,CA,half);cl=contact(V,CL,half)
        contacts=[contact(V,X,half) for X in (Ckin,Creg,Ccomp)]
        step=2.e-4;smallp=p[:2];H0=shear_bloch(smallp,width,charge,0)
        Hp=shear_bloch(smallp,width,charge,step);Hm=shear_bloch(smallp,width,charge,-step)
        source_errors.append({'abs_charge':charge,'first':float(np.max(abs((Hp-Hm)/(2*step)-Jg[:2]))),
                              'second':float(np.max(abs((Hp+Hm-2*H0)/step**2-Cg[:2])))})
        assert source_errors[-1]['first']<1.e-6 and source_errors[-1]['second']<2.e-6
        for n in range(max_mode+1):
            shift=np.roll(np.arange(L**3).reshape(L,L,L),-n,axis=0).ravel();q=2*np.pi*n/L
            ps,pt=bubble(E,V,E[shift],V[shift],JA,half)
            gs,gt=bubble(E,V,E[shift],V[shift],Jg,half)
            row=sums[n]
            row['photon_static']+=multiplicity*(ps+cp)
            row['photon_temporal_Euclidean']+=multiplicity*pt
            row['frame_static']+=multiplicity*(gs+sum(contacts))
            row['frame_temporal_Euclidean']+=multiplicity*gt
            row['frame_bubble']+=multiplicity*gs
            for key,val in zip(('frame_contact_kinetic','frame_contact_regulator','frame_contact_completion'),contacts):
                row[key]+=multiplicity*val
            if n:
                JL=m.X_ROOT*charge*(np.cos(p[:,0]+q/2)[:,None,None]*gamma[0]+np.sin(p[:,0]+q/2)[:,None,None]*beta)
                Ward=max(Ward,float(np.max(abs(2*np.sin(q/2)*JL-charge*(H[shift]-H)))))
                ls,_=bubble(E,V,E[shift],V[shift],JL,half)
                row['longitudinal_photon_static']+=multiplicity*(ls+cl)
        del H,E,V,sy,sz,wy,wz,Jg,Cg,JA,CA,CL,Ckin,Creg,Ccomp
    rows=[{'n':n,'q':2*np.pi*n/L,'qhat':2*np.sin(np.pi*n/L),**sums[n]} for n in range(max_mode+1)]
    fits={}
    for name in names:
        qh=np.asarray([r['qhat'] for r in rows[1:]])
        raw=np.asarray([r[name+'_static'] for r in rows[1:]])
        intercept=rows[0][name+'_static'];Y=raw-intercept
        design=np.column_stack((qh**2,qh**4));coef,*_=np.linalg.lstsq(design,Y,rcond=None)
        narrow,*_=np.linalg.lstsq(design[:-1],Y[:-1],rcond=None)
        fits[name]={'raw_uniform_response':intercept,'temporal_Euclidean_coefficient':rows[0][name+'_temporal_Euclidean'],
                    'spatial_q_squared':float(coef[0]),'spatial_q_fourth':float(coef[1]),
                    'window_shift_in_q_squared':float(narrow[0]-coef[0]),
                    'heldout_last_mode_residual':float(abs(design[-1]@narrow-Y[-1])),
                    'constant_subtracted_from_physical_H':False}
    long=max(abs(r['longitudinal_photon_static']) for r in rows[1:])
    assert Ward<1.e-10 and long<1.e-8
    common_gap=min(g['minimum_empty_energy'] for g in gaps)-max(g['maximum_occupied_energy'] for g in gaps)
    assert common_gap>0
    return {'L':L,'width':width,'spatial_dimension':3,'source_derivative_errors':source_errors,
            'band_filling_controls':gaps,'common_filling_gap':common_gap,
            'operator_photon_Ward_residual':Ward,'longitudinal_response_residual':long,
            'raw_response_rows':rows,'source_derivative_fits':fits,
            'area_energy_on_uniform_flat_connection':0.,
            'physical_pole_fit_permitted':False}


def coordinate_stress_symbol(L,width=2):
    """Geometric coordinate-stress source vs half-density generator.

    X=x+epsilon xi, xi=exp(iqx); nu=1+epsilon delta, e_xx=1/(1+epsilon delta),
    delta=i sin(q) xi. Independent Omega=I is the flat pulled-back spin frame.
    Keep the canonical Wilson and completion sources. Failure of this Gaussian
    candidate is not a no-go for every possible interacting completion.
    """
    q=2*np.pi/L;p=np.full(3,np.pi/L);pp=p+np.array([q,0,0])
    ga=[np.kron(np.eye(width),g) for g in m.GAMMA];beta=np.kron(np.eye(width),m.BETA)
    wall=m.slab.slab_hamiltonian((0,0,0),m.slab.wilson_mass(1),width)
    O=wall+3*beta;C=.5*O@O+3*np.eye(len(O));K=(-1j*ga[0]-beta)/2
    delta=1j*np.sin(q)
    d=np.exp(1j*p[0])-1;dp=np.exp(1j*pp[0])-1
    w=1-np.cos(p[0]);wp=1-np.cos(pp[0])
    Jkin=-.5*(np.sin(p[0])+np.sin(pp[0]))*ga[0]
    Jwil=-.5*(wp+w+np.conj(dp)*d)*beta
    def edge_derivative(q):
        return .25j*(1+np.exp(1j*q))*ga[0]+(.75+.25*np.exp(1j*q))*beta
    A=edge_derivative(q);Ad=edge_derivative(-q).conj().T
    z=-(1.5+.5*np.exp(-1j*q))
    Jcomp=.5*z*(O@beta+beta@O)+A@K.conj().T+K@Ad+np.exp(-1j*q)*(Ad@K+K.conj().T@A)
    H0=m.slab.slab_hamiltonian(p,m.slab.wilson_mass(1),width)
    H1=m.slab.slab_hamiltonian(pp,m.slab.wilson_mass(1),width)
    comm=1j*.5*(np.sin(p[0])+np.sin(pp[0]))*(H0-H1)
    source=delta*(Jkin+Jwil)
    lifted_source=m.X_ROOT*source+m.X_ROOT**2*delta*Jcomp
    discrepancy=lifted_source-m.X_ROOT*comm
    return {'L':L,'q':q,'external_p':p.tolist(),
            'effective_geometry_Ward_error':float(np.linalg.norm(source-comm,2)),
            'completed_Gaussian_Ward_error':float(np.linalg.norm(discrepancy,2)),
            'completed_Ward_error_over_q':float(np.linalg.norm(discrepancy,2)/q),
            'completion_source_over_q':float(np.linalg.norm(m.X_ROOT**2*delta*Jcomp,2)/q),
            'derived_small_q_limit':float(m.X_ROOT**2*np.linalg.norm(O@beta+beta@O+3*np.eye(len(O)),2)),
            'internal_spin_or_U1_Ward':False}


def finite_graph_source_control():
    """Validate the independent coordinate-source symbol on actual graph blocks."""
    L,width=3,1;sites,nb=m.lattice(L);N=L**3
    q=2*np.pi/L;phase=q*np.asarray(sites)[:,0];delta=np.sin(q)*np.cos(phase)
    links=np.broadcast_to(np.eye(4),(N,3,4,4)).copy();step=1.e-5
    Hs=[]
    for s in (-step,step):
        e=np.repeat(np.eye(3)[None],N,axis=0);e[:,0,0]=1/(1+s*delta)
        nu,gamma,metric,_=m.geometry(e)
        Hs.append(m.matter_matrix(L,width,1,nu,gamma,metric,links)['completed'])
    derivative=(Hs[1]-Hs[0])/(2*step)
    p=np.full(3,np.pi/L);pp=p+np.array([q,0,0])
    v=np.exp(1j*np.asarray(sites)@p)/np.sqrt(N);vp=np.exp(1j*np.asarray(sites)@pp)/np.sqrt(N)
    measured=np.kron(vp.conj()[None,:],np.eye(4))@derivative@np.kron(v[:,None],np.eye(4))
    ga=m.GAMMA;beta=m.BETA;wall=m.slab.slab_hamiltonian((0,0,0),m.slab.wilson_mass(1),1);O=wall+3*beta;K=(-1j*ga[0]-beta)/2
    Jkin=-.5*(np.sin(p[0])+np.sin(pp[0]))*ga[0]
    Jwil=-.5*((1-np.cos(p[0]))+(1-np.cos(pp[0]))+np.conj(np.exp(1j*pp[0])-1)*(np.exp(1j*p[0])-1))*beta
    A=.25j*(1+np.exp(1j*q))*ga[0]+(.75+.25*np.exp(1j*q))*beta
    Ad=(.25j*(1+np.exp(-1j*q))*ga[0]+(.75+.25*np.exp(-1j*q))*beta).conj().T
    Jcomp=.5*(-1.5-.5*np.exp(-1j*q))*(O@beta+beta@O)+A@K.conj().T+K@Ad+np.exp(-1j*q)*(Ad@K+K.conj().T@A)
    predicted=.5*np.sin(q)*(m.X_ROOT*(Jkin+Jwil)+m.X_ROOT**2*Jcomp)
    residual=float(np.max(abs(predicted-measured)))
    assert residual<1.e-8
    return {'L':L,'width':width,'coordinate_source_graph_vs_Bloch_residual':residual}


def run(quick=False,volumes=None):
    sizes=(8,12) if quick else (12,16,24)
    width=2 if quick else 6
    volumes=[one_volume(L,width,3 if quick else 4) for L in sizes] if volumes is None else volumes
    if [v['L'] for v in volumes]!=list(sizes) or any(v['width']!=width for v in volumes):
        raise ValueError('Volume checkpoints do not match the declared run')
    ward=[coordinate_stress_symbol(L) for L in (16,32,64,128,256)]
    assert ward[-1]['completed_Ward_error_over_q']>1.e-3
    assert abs(ward[-1]['completed_Ward_error_over_q']/ward[-1]['derived_small_q_limit']-1)<5.e-4
    assert ward[-1]['effective_geometry_Ward_error']<ward[0]['effective_geometry_Ward_error']
    return {'execution_pass':True,'mode':'quick' if quick else 'full',
            'charges':list(m.slab.CHARGES),'x':m.X_ROOT,'separate_sector_adjustments':0,
            'volumes':volumes,'coordinate_stress_Ward':ward,
            'independent_graph_source_control':finite_graph_source_control(),
            'full_Fock_completion_contains_quartic_terms':True,
            'quartic_interaction_corrections_in_large_volume_loop':False,
            'common_kinetic_photon_frame_matching_established':False,
            'physical_common_cone_test_pass':False,
            'decision':'The specified Gaussian completed matter sector has a nonvanishing leading coordinate-stress Ward defect. Exact internal Spin(3)/U(1) covariance does not repair spatial diffeomorphism symmetry. Do not publish a physical photon/gravity speed fit or subtract the response by hand.',
            'claim_boundary':'This failure concerns the declared candidate and Gaussian source. Full interacting Fock, frame, connection and continuum mechanisms remain untested.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--quick',action='store_true');p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run(a.quick)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(json.dumps(r,indent=2))
