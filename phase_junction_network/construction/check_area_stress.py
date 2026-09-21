#!/usr/bin/env python3
"""Executable construction gates; internal covariance is NOT gravity closure."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import block_diag, expm
import area_stress_model as m


def maxabs(a):
    return float(np.max(abs(a)))


def filled_curvature(H,J,C):
    E,U=np.linalg.eigh(H);half=len(E)//2
    V=U.conj().T@J@U
    bubble=-2*np.sum(abs(V[half:,:half])**2/(E[half:,None]-E[None,:half]))
    contact=np.trace(U[:,:half].conj().T@C@U[:,:half]).real
    return {'bubble':float(bubble),'contact':float(contact),'total':float(bubble+contact)}


def finite_area(samples):
    rng=np.random.default_rng(6202026)
    rows=[]
    for _ in range(samples):
        A=rng.normal(size=(3,3)); e=expm((A+A.T)*.1)[None]
        _,_,_,B=m.geometry(e); B=B[0,1,2]
        v=rng.normal(size=3);T=np.einsum('a,aij->ij',v,m.SPIN)
        W=expm(1j*T*.1)
        move=m.area_move(B,W);H=m.completed(move)
        Q=np.eye(4)-2*m.X_ROOT*move
        factor=(Q.conj().T@Q+Q@Q.conj().T-2*np.eye(4))/4
        sine=(W-W.conj().T)/(2j); cosine=(W+W.conj().T)/2
        derived=(-m.X_ROOT*np.trace(B@sine)+m.X_ROOT**2*np.trace(B@B@(np.eye(4)-cosine))).real/4
        rotation=expm(1j*np.einsum('a,aij->ij',rng.normal(size=3),m.SPIN))
        slope=-m.X_ROOT*np.trace(B@T).real/4
        step=1.e-5
        measured=(m.area_energy(B,expm(1j*step*T))-m.area_energy(B,expm(-1j*step*T)))/(2*step)
        phase=expm(1j*T*.05)
        # An equal-weight two-path coherent block contains (W-I)/2 exactly.
        had=np.kron(np.asarray([[1,1],[1,-1]])/np.sqrt(2),np.eye(4))
        interferometer=had@block_diag(W,np.eye(4))@had
        path_error=maxabs(-1j*B@interferometer[:4,4:]-move)
        flat=m.area_energy(B,np.eye(4))
        no_reference=m.completed(B/(2j))
        row={'hermiticity':maxabs(H-H.conj().T),'factorization':maxabs(H-factor),
             'trace_inventory':abs(np.trace(H).real/4-derived),
             'covariance':maxabs(m.completed(m.area_move(rotation@B@rotation.conj().T,rotation@W@rotation.conj().T))-rotation@H@rotation.conj().T),
             'orientation_reversal':abs(m.area_energy(B,W)-m.area_energy(-B,W.conj().T)),
             'flat_energy':abs(flat), 'linear_curvature_error':abs(measured-slope),
             'two_path_block_error':path_error,
             'lower_eigenvalue':float(np.linalg.eigvalsh(H)[0]),
             'no_reference_flat_energy':float(np.trace(no_reference).real/4)}
        assert all(row[k]<2.e-9 for k in ('hermiticity','factorization','trace_inventory','covariance','orientation_reversal','flat_energy','linear_curvature_error','two_path_block_error'))
        assert row['lower_eigenvalue'] >= -.5-1.e-12
        assert row['no_reference_flat_energy']>0
        rows.append(row)
    maxima={k:max(r[k] for r in rows) for k in rows[0]}
    maxima['lower_eigenvalue']=min(r['lower_eigenvalue'] for r in rows)
    # Explicit finite locked coordinate sector, not a continuous-spin claim.
    p=3; shift=np.roll(np.eye(p),1,axis=0);Z=np.diag(np.exp(2j*np.pi*np.arange(p)/p))
    Zh=np.kron(Z,np.eye(p));Xc=np.kron(np.eye(p),shift);Zc=np.kron(np.eye(p),Z)
    relative=Zc@Zh.conj().T
    lock=p/(4*np.pi)*(4*np.eye(p*p)-Xc-Xc.T-relative-relative.conj().T)
    Bspin=m.geometry(np.eye(3)[None])[3][0,1,2]
    b=np.diag([.8,1.,1.2]);Bfull=np.kron(np.kron(b,np.eye(p*p)),Bspin)
    W=block_diag(*[expm(1j*theta*m.SPIN[0]) for theta in (-.2,0.,.2)])
    Wfull=np.kron(np.eye(p*p),W)
    HB=m.completed(m.area_move(Bfull,Wfull))
    Lfull=np.kron(lock,np.eye(p*4))
    lockcomm=maxabs(HB@Lfull-Lfull@HB)
    assert lockcomm<1.e-12
    E,V=np.linalg.eigh(Lfull)
    P=V[:,np.isclose(E,E[0],atol=1.e-10)]
    rng=np.random.default_rng(62);init=P@rng.normal(size=P.shape[1]);init/=np.linalg.norm(init)
    evolved=expm(-1j*7*(HB+Lfull))@init
    leakage=float(np.linalg.norm(evolved-P@(P.conj().T@evolved)))
    assert leakage<1.e-11
    return {'samples':samples,'maxima':maxima,'finite_coordinate_dimension':len(HB),
            'lock_commutator':lockcomm,'unprojected_lock_leakage':leakage,
            'new_architecture':'equal-weight coherent loop/reference channel with fixed quadrature phase',
            'linear_term':'-x normalized_trace(B F)',
            'retained_completion':'x^2 normalized_trace[B^2 (I-Re W)]',
            'microscopic_unique_or_symmetry_protected':False,
            'finite_graviton_phase_established':False}


def graph_test(quick):
    L=3;width=1 if quick else 2;N=L**3;block=4*width
    rng=np.random.default_rng(882)
    frames=np.asarray([expm(.07*(A+A.T)) for A in rng.normal(size=(N,3,3))])
    nu,gamma,metric,area=m.geometry(frames)
    links=np.asarray([[expm(1j*np.einsum('a,aij->ij',v,m.SPIN)) for v in line]
                       for line in rng.normal(scale=.05,size=(N,3,3))])
    phases=rng.normal(scale=.05,size=(N,3))
    sites,nb=m.lattice(L)
    data=m.matter_matrix(L,width,-1,nu,gamma,metric,links,phases)
    rotations=np.asarray([expm(1j*np.einsum('a,aij->ij',v,m.SPIN)) for v in rng.normal(scale=.3,size=(N,3))])
    gp,ap,lp=m.transform_sources(gamma,area,links,nb,rotations)
    transformed=m.matter_matrix(L,width,-1,nu,gp,metric,lp,phases)
    U=block_diag(*[np.kron(np.eye(width),S) for S in rotations])
    covariance={key:maxabs(transformed[key]-U@data[key]@U.conj().T) for key in ('target','completed','kinetic','wilson','completion')}
    assert max(covariance.values())<3.e-11
    area_cov=abs(m.plaquette_energy(L,area,links)-m.plaquette_energy(L,ap,lp))
    assert area_cov<1.e-11
    lam=rng.normal(scale=.3,size=N)
    gauge=phases+lam[:,None]-lam[nb]
    charge=-1
    G=np.diag(np.repeat(np.exp(1j*charge*lam),block))
    gauged=m.matter_matrix(L,width,charge,nu,gamma,metric,links,gauge)
    gauge_error=maxabs(gauged['completed']-G@data['completed']@G.conj().T)
    assert gauge_error<2.e-11
    # Local spin source variation tests include links and Clifford fields.
    gen=block_diag(*[np.kron(np.eye(width),sum(v[a]*m.SPIN[a] for a in range(3))) for v in rng.normal(size=(N,3))])
    H=data['completed'];J=1j*(gen@H-H@gen);C=-(gen@(gen@H-H@gen)-(gen@H-H@gen)@gen)
    ward=filled_curvature(H,J,C)
    assert abs(ward['total'])<1.e-9
    # Directly differentiate geometry + links, not a conjugated output matrix.
    gsmall=[gen[a*block:a*block+4,a*block:a*block+4] for a in range(N)]
    step=2.e-4;Hs=[]
    for s in (-step,step):
        Rs=np.asarray([expm(1j*s*j) for j in gsmall])
        gg,aa,ll=m.transform_sources(gamma,area,links,nb,Rs)
        Hs.append(m.matter_matrix(L,width,charge,nu,gg,metric,ll,phases)['completed'])
    d_error=maxabs((Hs[1]-Hs[0])/(2*step)-J)
    c_error=maxabs((Hs[1]+Hs[0]-2*H)/step**2-C)
    assert d_error<1.e-6 and c_error<1.e-5
    # Omit spin-link response while varying the frame: required failure.
    bad=[]
    for s in (-step,step):
        Rs=np.asarray([expm(1j*s*j) for j in gsmall])
        gg,_,_=m.transform_sources(gamma,area,links,nb,Rs)
        bad.append(m.matter_matrix(L,width,charge,nu,gg,metric,links,phases)['completed'])
    missing_link=maxabs((bad[1]-bad[0])/(2*step)-J)
    assert missing_link>1.e-3
    # Fresh internal-generator bracket, not diffeomorphism algebra.
    algebra=0.
    for a in range(3):
        b=(a+1)%3;c=(a+2)%3
        algebra=max(algebra,maxabs(m.SPIN[a]@m.SPIN[b]-m.SPIN[b]@m.SPIN[a]-1j*m.SPIN[c]))
    assert algebra<1.e-14
    lap=rng.uniform(.8,1.2,size=N);HN=m.lapse_source(H,lap,block)
    init=rng.normal(size=len(H))+1j*rng.normal(size=len(H));init/=np.linalg.norm(init)
    evolved=expm(-1j*.7*HN)@init
    conservation=max(abs(np.linalg.norm(evolved)-1),abs(np.vdot(evolved,HN@evolved)-np.vdot(init,HN@init)))
    assert conservation<1.e-11
    # Coordinate-unitary orbit is a control, not a geometric gravity claim.
    K=(data['shifts'][0]-data['shifts'][0].conj().T)/(2j)
    xi=np.diag(np.repeat(np.sin(2*np.pi*np.asarray(sites)[:,0]/L),block))
    D=(xi@K+K@xi)/2
    JD=1j*(D@H-H@D);CD=-(D@(D@H-H@D)-(D@H-H@D)@D)
    coordinate_orbit=filled_curvature(H,JD,CD)
    assert abs(coordinate_orbit['total'])<1.e-9
    hermiticity=maxabs(H-H.conj().T)
    return {'L':L,'width':width,'matrix_dimension':len(H),'spin_covariance':covariance,
            'area_scalar_covariance':area_cov,'charge_covariance':gauge_error,
            'spin_second_variation_Ward':ward,'spin_first_derivative_error':d_error,
            'spin_second_derivative_error':c_error,'missing_spin_link_negative_control':missing_link,
            'internal_generator_algebra':algebra,'unprojected_energy_and_norm_error':float(conservation),
            'coordinate_unitary_orbit_control':coordinate_orbit,
            'hermiticity':hermiticity,
            'gravitational_constraint_status':'NOT_ESTABLISHED: internal/background covariance and unitary orbit are not nonlinear scalar/vector constraint closure'}


def flat_recovery(quick):
    L=3;width=1 if quick else 2;N=L**3
    nu,gamma,metric,area=m.geometry(np.repeat(np.eye(3)[None],N,axis=0))
    links=np.broadcast_to(np.eye(4),(N,3,4,4)).copy()
    data=m.matter_matrix(L,width,1,nu,gamma,metric,links)
    O=m.slab.slab_hamiltonian((0,0,0),m.slab.wilson_mass(1),width)+3*np.kron(np.eye(width),m.BETA)
    C=.5*O@O+3*np.eye(4*width)
    target_bands=[];completed_bands=[]
    for n in np.ndindex(L,L,L):
        p=2*np.pi*(np.asarray(n)+.5)/L
        H=m.slab.slab_hamiltonian(p,m.slab.wilson_mass(1),width)
        target_bands.extend(np.linalg.eigvalsh(H))
        completed_bands.extend(np.linalg.eigvalsh(m.X_ROOT*H+m.X_ROOT**2*C))
    errors={'effective_flat_spectrum':maxabs(np.linalg.eigvalsh(data['target'])-np.sort(target_bands)),
            'completed_flat_spectrum':maxabs(np.linalg.eigvalsh(data['completed'])-np.sort(completed_bands))}
    assert max(errors.values())<1.e-10
    return errors


def fock_test():
    # Four-component spinors on two sites: 8 modes and 256 Fock states.
    def target(s):
        a,b=np.exp(s/np.sqrt(2.)),np.exp(-s/np.sqrt(2.))
        O=m.slab.slab_hamiltonian((0,0,0),m.slab.wilson_mass(1),1)+(1+a*a+b*b)*m.BETA
        K=(-1j*a*m.GAMMA[1]-a*a*m.BETA)/2
        return np.block([[O,K],[K.conj().T,O]])
    raw=target(0.);model=m.full_fock_completion(raw,4)
    H=model['full'];Q=model['quartic_fock']
    one=np.asarray([1<<i for i in range(8)])
    assert maxabs(H[one][:,one]-model['onebody'])<1.e-12
    assert maxabs(Q[one][:,one])<1.e-12
    assert np.linalg.norm(Q)>1.e-3
    small=m.full_fock_completion(raw,4,m.X_ROOT/2)
    assert maxabs(small['quartic_fock']-Q/4)<1.e-12
    Jlocal=block_diag(m.SPIN[0],-.7*m.SPIN[1])
    J=m.bilinear_fock(Jlocal)
    D=1j*(J@H-H@J);C=-(J@(J@H-H@J)-(J@H-H@J)@J)
    step=1.e-4
    matrices=[]
    for s in (-step,step):
        U=expm(1j*s*Jlocal)
        matrices.append(m.full_fock_completion(U@raw@U.conj().T,4)['full'])
    first=maxabs((matrices[1]-matrices[0])/(2*step)-D)
    second=maxabs((matrices[1]+matrices[0]-2*H)/step**2-C)
    assert first<1.e-7 and second<1.e-6
    # Geometry stress must also differentiate the quartic part.
    fp=m.full_fock_completion(target(step),4)
    fm=m.full_fock_completion(target(-step),4)
    stress_quartic=(fp['quartic_fock']-fm['quartic_fock'])/(2*step)
    assert np.linalg.norm(stress_quartic)>1.e-3
    init=np.arange(1,257,dtype=float);init/=np.linalg.norm(init)
    evolved=expm(-1j*.3*H)@init
    error=float(abs(np.vdot(evolved,H@evolved)-np.vdot(init,H@init)))
    assert error<1.e-11
    return {'fermion_modes':8,'Fock_dimension':256,'quartic_operator_norm':float(np.linalg.norm(Q,2)),
            'quartic_frame_stress_norm':float(np.linalg.norm(stress_quartic,2)),
            'one_particle_reduction_residual':maxabs(H[one][:,one]-model['onebody']),
            'full_Fock_spin_first_Ward_error':first,'full_Fock_spin_second_Ward_error':second,
            'unprojected_interacting_energy_error':error,
            'loop_boundary':'Large-volume response is the Gaussian one-loop contribution of this candidate; quartic corrections are not included in that response.'}


def run(quick=False):
    return {'execution_pass':True,'mode':'quick' if quick else 'full',
            'area':finite_area(8 if quick else 40),'graph':graph_test(quick),
            'flat_recovery':flat_recovery(quick),'full_Fock_completion':fock_test(),'x':m.X_ROOT,
            'separate_sector_adjustments':0,'historical_results_modified':False,
            'full_unified_theory_claimed':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--quick',action='store_true');p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=run(a.quick);a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(json.dumps(r,indent=2))
