#!/usr/bin/env python3
"""Classical quadratic 3D area-reference kernel, without a TT insertion.

Both Fourier quadratures of six symmetric-frame components and nine connection
components are retained. A quadratic Fourier block is not an exact many-body
sector of the finite-spin quantum model.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import null_space, expm
import area_stress_model as m
sys.path.insert(0,str(m.ROOT/'microscopic'))
from check_constrained_local_reduction import H_BASIS, gauge_map, scalar_constraint_vector, KINETIC_METRIC
U_INHERITED=0.006789913753821769


def quadratic_data(L, mode=(1,0,0)):
    k=2*np.pi*np.asarray(mode,dtype=float)/L
    # Exact quadrature of the degree-two plane-wave polynomial, not a chain.
    phase=2*np.pi*np.arange(8)/8
    h_basis=np.asarray(H_BASIS)
    B0=m.geometry(np.eye(3)[None])[3][0]
    deltaB=[]
    for h in h_basis:
        dg=-.5*np.einsum('ai,abc->ibc',h,m.GAMMA)
        db=np.empty((3,3,4,4),dtype=complex)
        for i in range(3):
            for j in range(3):
                db[i,j]=.5*np.trace(h)*B0[i,j]+(dg[i]@m.GAMMA[j]+m.GAMMA[i]@dg[j]-dg[j]@m.GAMMA[i]-m.GAMMA[j]@dg[i])/(2j)
        deltaB.append(db)
    deltaB=np.asarray(deltaB)
    def energy(v):
        h=v[:12].reshape(2,6);c=v[12:].reshape(2,3,3)
        db=(np.einsum('a,aijbc->ijbc',h[0],deltaB)[None]*np.cos(phase)[:,None,None,None,None]
            +np.einsum('a,aijbc->ijbc',h[1],deltaB)[None]*np.sin(phase)[:,None,None,None,None])
        def connection(i,offset=0.):
            f=phase+k[i]/2+offset
            return (np.sin(f)[:,None,None]*np.einsum('a,abc->bc',c[0,i],m.SPIN)
                    +np.cos(f)[:,None,None]*np.einsum('a,abc->bc',c[1,i],m.SPIN))
        E=0.
        for i in range(3):
            for j in range(i+1,3):
                parts=[connection(i),connection(j,k[i]),-connection(i,k[j]),-connection(j)]
                D=sum(parts);S2=np.zeros_like(D)
                for a in range(4):
                    for b in range(a+1,4):
                        S2+=.5j*(parts[a]@parts[b]-parts[b]@parts[a])
                e=-m.X_ROOT*(db[:,i,j]@D+B0[i,j]@S2)+.5*m.X_ROOT**2*(B0[i,j]@B0[i,j])@(D@D)
                E+=np.trace(e,axis1=-2,axis2=-1).real.mean()/4
        return float(E)
    Hess=np.zeros((30,30));eye=np.eye(30);diagonal=np.asarray([energy(v) for v in eye])
    for a in range(30):
        Hess[a,a]=2*diagonal[a]
        for b in range(a):
            Hess[a,b]=Hess[b,a]=energy(eye[a]+eye[b])-diagonal[a]-diagonal[b]
    A=Hess[12:,12:];B=Hess[:12,12:]
    condition=float(np.linalg.cond(A))
    rank=int(np.linalg.matrix_rank(A,tol=1.e-10))
    if rank!=18:
        raise ArithmeticError(f'connection Hessian rank {rank}; do not hide its null space')
    reduced=Hess[:12,:12]-B@np.linalg.solve(A,B.T)
    K=2*reduced  # volume average cos^2=1/2
    kg=2*np.sin(k/2)
    G1=gauge_map(kg);S1=scalar_constraint_vector(kg)
    G=np.zeros((12,6));G[:6,:3]=G1;G[6:,3:]=G1
    S=np.zeros((2,12));S[0,:6]=S1;S[1,6:]=S1
    gauge_error=float(np.linalg.norm(K@G)/max(np.linalg.norm(K)*np.linalg.norm(G),1.e-30))
    constraints=np.zeros((8,24));constraints[:6,12:]=G.T;constraints[6:,:12]=S
    physical=null_space(constraints)
    H=np.zeros((24,24));H[:12,:12]=K
    H[12:18,12:18]=U_INHERITED*KINETIC_METRIC
    H[18:,18:]=U_INHERITED*KINETIC_METRIC
    symplectic=np.block([[np.zeros((12,12)),np.eye(12)],[-np.eye(12),np.zeros((12,12))]])
    leakage=float(np.linalg.norm(constraints@symplectic@H@physical)/max(np.linalg.norm(H),1.e-30))
    # Fresh generator from the coframe/Bianchi relation, including placement.
    d=kg*np.exp(-.5j*k)
    hb=np.asarray(H_BASIS)
    Gc=np.column_stack([np.einsum('aij,ij->a',hb,np.outer(d,e)+np.outer(e,d)) for e in np.eye(3)])
    def realify(A):
        return np.block([[A.real,A.imag],[-A.imag,A.real]])
    Gback=realify(Gc)
    fresh_gauge=float(np.linalg.norm(K@Gback)/max(np.linalg.norm(K)*np.linalg.norm(Gback),1.e-30))
    # Diagonal momenta shift one link to a common vertex; no fitted coefficient.
    phases=[]
    for tensor in hb:
        i,j=np.argwhere(abs(tensor)>1.e-8)[0]
        phases.append(np.exp(.5j*(k[i]+k[j])))
    D=realify(np.diag(phases))
    U0=np.kron(np.eye(2),KINETIC_METRIC)
    Umatch=D.T@U0@D
    Sreal=np.zeros((2,12));Sreal[0,:6]=S1;Sreal[1,6:]=S1
    Smatch=Sreal@D
    Cmatch=np.zeros((8,24));Cmatch[:6,12:]=Gback.T;Cmatch[6:,:12]=Smatch
    Hmatch=H.copy();Hmatch[12:,12:]=U_INHERITED*Umatch
    Nmatch=null_space(Cmatch)
    compatible_leakage=float(np.linalg.norm(Cmatch@symplectic@Hmatch@Nmatch)/max(np.linalg.norm(Hmatch),1.e-30))
    bracket_error=float(np.linalg.norm(Cmatch@symplectic@Cmatch.T))
    # Reduce the actual constraints; no target dispersion selects polarizations.
    reduced_basis=null_space(np.vstack((Gback.T,Smatch)))
    Kr=reduced_basis.T@K@reduced_basis
    Ur=reduced_basis.T@Umatch@reduced_basis
    if np.linalg.norm(Ur-np.eye(len(Ur)))>1.e-10:
        raise ArithmeticError('Reduced kinetic metric is not identity; solve a generalized eigenproblem')
    spring=np.linalg.eigvalsh(Kr)
    gaps=np.sqrt(np.maximum(np.linalg.eigvalsh(U_INHERITED*Ur@Kr),0))
    initial=Nmatch@np.arange(1.,Nmatch.shape[1]+1);initial/=np.linalg.norm(initial)
    ev=expm(.7*symplectic@Hmatch)@initial
    preservation=float(np.max(abs(Cmatch@ev)))
    assert fresh_gauge<1.e-11 and compatible_leakage<1.e-10 and bracket_error<1.e-10
    assert reduced_basis.shape[1]==4 and np.min(spring)>0 and preservation<1.e-10
    placement={
        'derived_backward_vector_gauge_residual':fresh_gauge,
        'local_trace_placement':'trace at common vertex = sum_i pi_ii(x+hat_i)',
        'new_adjustable_coefficients':0,
        'first_class_bracket_error':bracket_error,
        'constraint_evolution_relative_residual':compatible_leakage,
        'unprojected_evolution_constraint_error':preservation,
        'real_physical_configuration_modes':int(reduced_basis.shape[1]),
        'polarizations_per_signed_momentum':int(reduced_basis.shape[1]//2),
        'positive_spring_eigenvalues':spring.tolist(),
        'quadratic_gravity_gaps':gaps.tolist(),
        'effective_gravity_speed_at_this_momentum':(gaps/np.linalg.norm(kg)).tolist(),
        'quantum_gravity_pole_proven':False,
        'assumption':'Inherited DeWitt kinetic coefficients with local placement matched to the newly derived coframe gauge generator; not a full quantum derivation of the kinetic rule.'}
    kinetic_spectrum=np.linalg.eigvalsh(K)
    threshold=1.e-8*max(np.linalg.norm(K),1.e-15)
    inertia={'negative':int(np.sum(kinetic_spectrum < -threshold)),
             'zero':int(np.sum(abs(kinetic_spectrum)<=threshold)),
             'positive':int(np.sum(kinetic_spectrum>threshold))}
    return {'L':L,'mode_vector':list(mode),'k':k.tolist(),'khat':kg.tolist(),'connection_rank':rank,'connection_condition_number':condition,
            'connection_eigenvalues':np.linalg.eigvalsh(A).tolist(),
            'reduced_frame_Hessian':K.tolist(),'frame_Hessian_eigenvalues':kinetic_spectrum.tolist(),
            'frame_Hessian_inertia':inertia,
            'vector_gauge_relative_residual':gauge_error,
            'linear_constraint_evolution_relative_residual':leakage,
            'placement_matched_kinetic':placement,
            'TT_projector_inserted':False}, energy, Hess


def exact_energy(v,L,mode=(1,0,0)):
    phase=2*np.pi*np.arange(L)/L
    hc=np.einsum('a,aij->ij',v[:6],np.asarray(H_BASIS))
    hs=np.einsum('a,aij->ij',v[6:12],np.asarray(H_BASIS))
    frames=np.asarray([expm(-.5*(np.cos(p)*hc+np.sin(p)*hs)) for p in phase])
    _,_,_,area=m.geometry(frames)
    c=v[12:].reshape(2,3,3);links=np.empty((L,3,4,4),dtype=complex)
    for n in range(L):
        for i in range(3):
            p=phase[n]+np.pi*mode[i]/L
            vector=np.sin(p)*c[0,i]+np.cos(p)*c[1,i]
            links[n,i]=expm(1j*np.einsum('a,aij->ij',vector,m.SPIN))
    total=0.
    for n in range(L):
        for i in range(3):
            for j in range(i+1,3):
                ni=(n+mode[i])%L;nj=(n+mode[j])%L
                W=links[n,i]@links[ni,j]@links[nj,i].conj().T@links[n,j].conj().T
                total+=m.area_energy(area[n,i,j],W)/L
    return total


def run(quick=False):
    rows=[]
    for L in ((8,16,24) if quick else (8,12,16,24,32,48,64,96,128)):
        for mode in ((1,0,0),(1,1,0),(1,2,1)):
            row,_,_=quadratic_data(L,mode);rows.append(row)
    fits=[]
    for mode in ((1,0,0),(1,1,0),(1,2,1)):
        subset=[r for r in rows if r['mode_vector']==list(mode) and r['L']>= (8 if quick else 32)]
        q2=np.asarray([np.sum(np.asarray(r['khat'])**2) for r in subset])
        for index in (0,2):
            speeds=np.asarray([r['placement_matched_kinetic']['effective_gravity_speed_at_this_momentum'][index] for r in subset])
            design=np.column_stack((np.ones(len(q2)),q2,q2*q2)) if len(q2)>=4 else np.column_stack((np.ones(len(q2)),q2))
            coef,*_=np.linalg.lstsq(design,speeds**2,rcond=None)
            narrow,*_=np.linalg.lstsq(design[1:],speeds[1:]**2,rcond=None)
            fits.append({'mode_vector':list(mode),'polarization':index//2,
                         'extrapolated_classical_c_squared':float(coef[0]),
                         'c_squared_over_inherited_U_times_x':float(coef[0]/(U_INHERITED*m.X_ROOT)),
                         'fit_window_shift_in_c_squared':float(narrow[0]-coef[0]),
                         'physical_quantum_common_cone':False})
    row,E,Hess=quadratic_data(8)
    rng=np.random.default_rng(2309);v=rng.normal(size=30);v/=np.linalg.norm(v)
    eps=np.array([.01,.02,.04,.08]);errors=[]
    for s in eps:
        measured=(exact_energy(s*v,8)+exact_energy(-s*v,8))/2
        errors.append(abs(measured-s*s*E(v)))
    power=float(np.polyfit(np.log(eps),np.log(errors),1)[0])
    assert 3.8<power<4.2
    assert max(abs(E(v)-.5*v@Hess@v),abs(exact_energy(np.zeros(30),8)))<1.e-12
    for row in rows:
        if row['L']!=8:
            row.pop('reduced_frame_Hessian')
            row.pop('connection_eigenvalues')
    return {'execution_pass':True,'mode':'quick' if quick else 'full','area_kernel':rows,
            'exact_move_Taylor_error_power':power,
            'common_x':m.X_ROOT,'frame_kinetic_U_inherited':U_INHERITED,
            'classical_low_momentum_fits':fits,
            'old_unshifted_constraint_gate_pass':bool(all(r['linear_constraint_evolution_relative_residual']<1.e-10 for r in rows)),
            'placement_matched_quadratic_constraint_gate_pass':bool(all(r['placement_matched_kinetic']['constraint_evolution_relative_residual']<1.e-10 for r in rows)),
            'full_nonlinear_constraints_constructed':False,
            'photon_3D_quantum_pole_available':False,
            'physical_speed_matching_performed':False,
            'claim_boundary':'Classical quadratic kernel derived from the specified area/reference move. Background auxiliary elimination is not a full finite quantum-gravity reduction.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--quick',action='store_true');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=run(a.quick);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(json.dumps(r,indent=2))
