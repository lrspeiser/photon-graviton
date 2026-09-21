#!/usr/bin/env python3
"""Two explicit coherent-mediator constructions, rather than fitted corrections.

1. Exactly eliminate a four-stage gapped ring to determine all oriented return
   amplitudes. Preserve its absolute transfer gain and frequency dependence.
2. Replace the site-pinched matter completion with a LOCAL doubled fermion
   Hamiltonian. This is a different microscopic ansatz, not an equivalence to
   the former interacting Fock completion. It removes the known constant
   simultaneous-low-momentum stress defect without adding a counterterm.

Neither construction is a derivation of Newton's force law or full nonlinear
spacetime constraints. The ring transfer is not silently completed a second time
and relabeled an exact low-energy gravitational energy.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import numpy as np
import scipy
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'construction'))
import area_stress_model as old
from check_geometry_response import coordinate_stress_symbol

X = old.X_ROOT
PARENT = '0e8fc4087f8a676ed872daa970acc38756c78af9'


def adj(a):
    return a.conj().T


def norm(a):
    return float(np.linalg.norm(a, 2))


def ring_parameters(x=X):
    if not 0 < x < .5:
        raise ValueError('Declared gapped small-coupling branch requires 0<x<1/2')
    d = 1 + 2*x*x
    alpha = (d + np.sqrt(d*d-4*x*x))/2
    r = x/alpha
    rho = r**4
    return d, alpha, r, rho


def ring_shift(W):
    """Four elementary spatial transports; S^4 has the plaquette holonomy.

    Internal gauge chosen so the three first link transports are identity.
    Arbitrary four-link unitary products are related by local frame rotations.
    """
    n = len(W)
    if W.shape != (n,n) or norm(adj(W)@W-np.eye(n)) > 1e-10:
        raise ValueError('Unitary holonomy required')
    S = np.zeros((4*n,4*n), complex)
    for j in range(3):
        S[(j+1)*n:(j+2)*n,j*n:(j+1)*n] = np.eye(n)
    S[:n,3*n:] = W
    return S


def ring_dirac(W):
    """LOCAL Hermitian dilation whose square is the existing completed ring.

    Q=(1+2x^2)I-x(S+S†)=alpha(I-rS)†(I-rS).
    K=[0 A; A† 0], A=sqrt(alpha)(I-rS).
    Introducing this chiral partner is an explicit architecture choice. Its
    couplings are algebraically fixed by Q, not fitted to a gravitational speed.
    """
    d, alpha, r, rho = ring_parameters()
    S = ring_shift(W)
    A = np.sqrt(alpha)*(np.eye(len(S))-r*S)
    Z = np.zeros_like(A)
    K = np.block([[Z,A],[adj(A),Z]])
    Q = d*np.eye(len(S))-X*(S+adj(S))
    return K,Q,A


def return_block(W):
    """Lower-left zero-energy resolvent of K, node 0 to node 0."""
    _,alpha,_,rho = ring_parameters()
    return np.linalg.solve(np.eye(len(W))-rho*W,np.eye(len(W)))/np.sqrt(alpha)


def normalized_shape(W):
    """Normalization ONLY for reporting shape; raw gain is not changed."""
    _,_,_,rho=ring_parameters()
    I=np.eye(len(W))
    return (1-rho)*(W-I)@np.linalg.inv(I-rho*W)


def ring_bridge(W):
    """Two low ports, physical/reference rings, equal splitter amplitudes.

    All external port couplings are x/sqrt(2). The opposite sign of the reference
    return implements destructive interference at W=I. No vacuum spectrum is
    subtracted after the fact. The entire microscopic Hermitian block is kept.
    """
    n=len(W); K,_,_=ring_dirac(W); Kref,_,_=ring_dirac(np.eye(n))
    N=len(K); Ktot=np.block([[K,np.zeros_like(K)],[np.zeros_like(K),Kref]])
    V=np.zeros((2*n,2*N),complex)
    # left low port sees upper chiral node 0; right sees lower chiral node 0
    V[:n,:n]=X/np.sqrt(2)*np.eye(n)
    V[:n,N:N+n]=X/np.sqrt(2)*np.eye(n)
    V[n:,4*n:5*n]=X/np.sqrt(2)*np.eye(n)
    V[n:,N+4*n:N+5*n]=-X/np.sqrt(2)*np.eye(n)
    H=np.block([[np.zeros((2*n,2*n)),V],[adj(V),Ktot]])
    return H,V,Ktot


def solve_ring(quick):
    d,alpha,r,rho=ring_parameters()
    rng=np.random.default_rng(261101)
    errs={name:0. for name in ('dilation_square','return_resolvent','reference_bridge',
                                'energy_dependent_pole','induced_residue','path_sum')}
    minimum_gap=10.; largest_residue=0.
    for _ in range(3 if quick else 10):
        raw=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); generator=(raw+adj(raw))/2
        W=expm(1j*.6*generator)
        K,Q,A=ring_dirac(W); n=len(W)
        target=np.block([[Q,np.zeros_like(Q)],[np.zeros_like(Q),Q]])
        errs['dilation_square']=max(errs['dilation_square'],norm(K@K-target))
        R=np.linalg.inv(K)
        errs['return_resolvent']=max(errs['return_resolvent'],norm(R[4*n:5*n,:n]-return_block(W)))
        minimum_gap=min(minimum_gap,float(np.min(abs(np.linalg.eigvalsh(K)))))
        H,V,heavy=ring_bridge(W)
        Heff=-V@np.linalg.solve(heavy,adj(V))
        difference=return_block(W)-return_block(np.eye(n))
        predicted=np.block([[np.zeros((n,n)),-.5*X*X*adj(difference)],
                             [-.5*X*X*difference,np.zeros((n,n))]])
        errs['reference_bridge']=max(errs['reference_bridge'],norm(Heff-predicted))
        # The original norm of a retained low component is not one. Keep Z!
        Z=np.eye(2*n)+V@np.linalg.solve(heavy@heavy,adj(V))
        eps=1e-5
        def effective_kernel(z):
            return z*np.eye(2*n)+V@np.linalg.solve(heavy-z*np.eye(len(heavy)),adj(V))
        dz=(effective_kernel(eps)-effective_kernel(-eps))/(2*eps)
        errs['induced_residue']=max(errs['induced_residue'],norm(dz-Z))
        largest_residue=max(largest_residue,norm(Z-np.eye(2*n)))
        eigenvalues,eigenvectors=np.linalg.eigh(H)
        for idx in np.argsort(abs(eigenvalues))[:2*n]:
            e=float(eigenvalues[idx]); v=eigenvectors[:2*n,idx]
            errs['energy_dependent_pole']=max(errs['energy_dependent_pole'],float(np.linalg.norm(effective_kernel(e)@v)))
        P=np.zeros_like(W); power=np.eye(n)
        for winding in range(1,9):
            power=power@W
            P+=(1-rho)**2*rho**(winding-1)*(power-np.eye(n))
        errs['path_sum']=max(errs['path_sum'],norm(P-normalized_shape(W)))
    H0,V0,K0=ring_bridge(np.eye(2))
    assert norm(-V0@np.linalg.solve(K0,adj(V0)))<1e-14
    assert max(errs.values()) < 2e-9
    # Analytic Taylor coefficients, independently fit the odd measured transfer.
    # Use extended precision series to avoid tiny coefficients hiding in a fit.
    m2=(1+rho)/(1-rho)
    m3=(1+4*rho+rho*rho)/(1-rho)**2
    R4=(1+10*rho+rho*rho)/(1-rho)**2
    theta=np.linspace(.01,.09,14)
    signal=(1-rho)**2*np.sin(theta)/(1+rho*rho-2*rho*np.cos(theta))
    coefficients=np.linalg.lstsq(np.column_stack([theta,theta**3,theta**5,theta**7]),signal,rcond=None)[0]
    measured_m3=-6*coefficients[1]/coefficients[0]
    assert abs(measured_m3-m3)<1e-8
    raw_slope=.5*X*X*rho/(np.sqrt(alpha)*(1-rho)**2)
    Z0=np.eye(4)+V0@np.linalg.solve(K0@K0,adj(V0))
    assert norm(Z0-np.trace(Z0).real/4*np.eye(4))<1e-12
    z0=float(np.trace(Z0).real/4)
    # Obtain the ACTUAL small low-band slope, not the static Schur slope alone.
    phase_rows=[]
    for phase in (.02,.05,.1,.2):
        W=np.diag([np.exp(1j*phase),np.exp(-1j*phase)])
        H,_,_=ring_bridge(W)
        low=sorted(np.linalg.eigvalsh(H),key=abs)[:4]
        phase_rows.append({'phase':phase,'low_band_positive_energy':float(max(low)),
                           'static_transfer_magnitude':float(.5*X*X*norm(return_block(W)-return_block(np.eye(2)))),
                           'linear_dynamical_estimate':raw_slope/z0*phase})
    return {'x':X,'primitive_stages':4,'completed_ring_onsite':d,
            'alpha':alpha,'one_step_decay_r':r,'one_winding_return_rho':rho,
            'a_n':'(1-rho)^2 rho^(n-1); derived oriented return shape, not fitted path amplitudes',
            'normalized_m2':m2,'normalized_m3':m3,'normalized_R4_transfer_power':R4,
            'R4_definition':'Minus the fourth over second theta derivative of |P|^2 at zero; not a gravitational field equation',
            'm3_from_independent_transfer_fit':float(measured_m3),
            'normalized_oriented_curvature_response':'(1-rho)^2 sin(theta)/(1+rho^2-2rho cos(theta))',
            'raw_static_bridge_slope':raw_slope,'flat_induced_kinetic_residue':z0,
            'small_phase_dynamical_low_band_slope':raw_slope/z0,
            'auxiliary_gap_lower_bound':float(np.sqrt(alpha)*(1-r)),
            'minimum_sampled_auxiliary_gap':minimum_gap,
            'largest_induced_residue_change':largest_residue,
            'identity_residuals':errs,'physical_low_band_examples':phase_rows,
            'new_free_continuous_coefficients':0,
            'explicit_architecture_choices':['chiral Hermitian dilation of the completed four-step ring',
              'two opposite-chirality low ports with equal fixed-amplitude physical/reference channels'],
            'not_claimed':'Transfer law is not yet the full area Hamiltonian or an acceleration-distance law. Completing this effective amplitude again would define another model; it is not done here.'}


def coherent_parent(h):
    A=X*h; I=np.eye(len(h))
    return np.block([[A,A],[A,I+A]])


def light_function(e):
    """Stable evaluation of a+(1-sqrt(1+4a²))/2."""
    a=X*np.asarray(e)
    return a-2*a*a/(1+np.sqrt(1+4*a*a))


def heavy_function(e):
    a=X*np.asarray(e)
    return a+(1+np.sqrt(1+4*a*a))/2


def function_matrix(h):
    e,u=np.linalg.eigh(h)
    return (u*light_function(e))@adj(u)


def light_frechet(h1,J,h0):
    """Off-diagonal source of f(h), with coincident-eigenvalue limits."""
    e1,u1=np.linalg.eigh(h1); e0,u0=np.linalg.eigh(h0)
    delta=e1[:,None]-e0[None,:]
    fdelta=light_function(e1)[:,None]-light_function(e0)[None,:]
    divided=np.zeros_like(delta)
    mask=abs(delta)>1e-10
    np.divide(fdelta,delta,out=divided,where=mask)
    mid=(e1[:,None]+e0[None,:])/2
    a=X*mid
    divided[~mask]=(X*(1-2*a/np.sqrt(1+4*a*a)))[~mask]
    return u1@(divided*(adj(u1)@J@u0))@adj(u0)


def coordinate_vertices(L,p=None,width=2,charge=1):
    q=2*np.pi/L
    p=np.full(3,np.pi/L) if p is None else np.asarray(p,float)
    pp=p+np.array([q,0,0])
    gamma=np.kron(np.eye(width),old.GAMMA[0]); beta=np.kron(np.eye(width),old.BETA)
    d=np.exp(1j*p[0])-1; dp=np.exp(1j*pp[0])-1
    Jk=-.5*(np.sin(p[0])+np.sin(pp[0]))*gamma
    Jw=-.5*((1-np.cos(p[0]))+(1-np.cos(pp[0]))+np.conj(dp)*d)*beta
    J=1j*np.sin(q)*(Jk+Jw)
    h0=old.slab.slab_hamiltonian(p,old.slab.wilson_mass(charge),width)
    h1=old.slab.slab_hamiltonian(pp,old.slab.wilson_mass(charge),width)
    D=.5*(np.sin(p[0])+np.sin(pp[0]))
    return q,h0,h1,J,D


def solve_matter(quick):
    rows=[]
    for L in ((16,64,256) if quick else (16,32,64,128,256,512,1024)):
        q,h0,h1,J,D=coordinate_vertices(L)
        old_error=coordinate_stress_symbol(L)['completed_Ward_error_over_q']
        R0=J-1j*D*(h0-h1)
        Jparent=X*np.block([[J,J],[J,J]])
        R=Jparent-1j*D*(coherent_parent(h0)-coherent_parent(h1))
        target_R=X*np.block([[R0,R0],[R0,R0]])
        assert norm(R-target_R)<1e-14
        light_J=light_frechet(h1,J,h0)
        light_R=light_J-1j*D*(function_matrix(h0)-function_matrix(h1))
        rows.append({'L':L,'q':q,'old_site_pinched_defect_over_q':old_error,
                     'coherent_parent_defect_over_q':norm(R)/q,
                     'exact_parent_error_bound_over_q':2*X*norm(R0)/q,
                     'light_band_static_operator_defect_over_q':norm(light_R)/q})
    power=float(np.polyfit(np.log([r['q'] for r in rows][-3:]),
                           np.log([r['coherent_parent_defect_over_q'] for r in rows][-3:]),1)[0])
    assert 2.95<power<3.05
    # Essential distinction: fixing the low-p artifact does not make the UV
    # geometric regulator diffeomorphism-covariant or complete loop integrals.
    q,h0,h1,J,D=coordinate_vertices(4096,p=[1.1,.7,.4])
    UV=X*np.block([[J,J],[J,J]])-1j*D*(coherent_parent(h0)-coherent_parent(h1))
    finite_uv=norm(UV)/q
    assert finite_uv>.001
    # Local real-space geometric reconstruction, not response-by-conjugation.
    L=3;width=1;sites,_=old.lattice(L);N=L**3
    qs=2*np.pi/L;phase=qs*np.asarray(sites)[:,0];delta=np.sin(qs)*np.cos(phase)
    links=np.broadcast_to(np.eye(4),(N,3,4,4)).copy()
    parents=[];targets=[];step=1e-5
    for s in (-step,0.,step):
        frame=np.broadcast_to(np.eye(3),(N,3,3)).copy();frame[:,0,0]=1/(1+s*delta)
        volume,gamma,metric,_=old.geometry(frame)
        h=old.matter_matrix(L,width,1,volume,gamma,metric,links)['target']
        targets.append(h); parents.append(coherent_parent(h))
    p=np.full(3,np.pi/L);pp=p+np.array([qs,0,0]);v=np.exp(1j*np.asarray(sites)@p)/np.sqrt(N);vp=np.exp(1j*np.asarray(sites)@pp)/np.sqrt(N)
    right=np.kron(np.eye(2),np.kron(v[:,None],np.eye(4)))
    left=np.kron(np.eye(2),np.kron(vp.conj()[None,:],np.eye(4)))
    measured=left@((parents[2]-parents[0])/(2*step))@right
    _,_,_,J,D=coordinate_vertices(L,width=1)
    predicted=X*np.block([[J,J],[J,J]])/(2j)
    graph_error=norm(measured-predicted)
    assert graph_error<1e-8
    e=np.linalg.eigvalsh(targets[1]); actual=np.linalg.eigvalsh(parents[1]);expected=np.sort(np.r_[light_function(e),heavy_function(e)])
    spectrum_error=float(np.max(abs(actual-expected)))
    assert spectrum_error<1e-11 and np.min(heavy_function(e))>.5
    # A local bilinear parent has an exact, local Fock lift. It is NOT the old
    # quartic interacting completed-Fock model, which remains unmodified.
    hsmall=np.array([[.13,.2+.05j],[.2-.05j,-.11]])
    hparent=coherent_parent(hsmall)
    Hf=old.bilinear_fock(hparent)
    single=np.linalg.eigvalsh(hparent)
    sums=sorted(sum(single[i] for i in range(4) if mask>>i&1) for mask in range(16))
    fock_error=float(np.max(abs(np.linalg.eigvalsh(Hf)-sums)))
    assert fock_error<1e-12
    return {'microscopic_parent':'H_joint=[[x h,x h],[x h,I+x h]]; h contains the complete existing geometric kinetic/Wilson/wall operator',
            'effective_light_band':'f(h)=x h+[I-sqrt(I+4 x^2 h^2)]/2',
            'exact_inverse_propagator':'K(z)=z-xh-x^2 h(z-I-xh)^(-1)h',
            'series':'f(h)=x h-x^2 h^2+x^4 h^4-2x^6 h^6+...',
            'source_rule':'delta f=f^[1](h_left,h_right) delta h; heavy-resolvent and residue terms must be retained for physical poles',
            'new_sector_specific_adjustments':0,
            'spatial_support':'The microscopic doubled parent has the same spatial hopping range as h. Eliminating the gapped partner generates longer-range effective paths; the uneliminated parent remains local.',
            'all_matter_species_receive_same_map':True,
            'new_architecture_not_old_model_equivalence':True,
            'low_p_stress_rows':rows,'defect_over_q_small_q_power':power,
            'fixed_high_internal_momentum_defect_over_q':finite_uv,
            'full_microscopic_coordinate_Ward_identity_solved':False,
            'complete_interacting_matter_loop_calculated':False,
            'new_low_momentum_massless_partner':False,
            'heavy_branch_gap_bound':.5,
            'graph_source_error':graph_error,'finite_graph_spectrum_error':spectrum_error,
            'full_small_fock_spectrum_error':fock_error,
            'physical_gravity_formula_derived':False}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--quick',action='store_true')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    original=ROOT/'microscopic/fermionic_multipair_vacuum_results.json'
    before=original.read_bytes()
    report={'source_parent':PARENT,'mode':'quick' if args.quick else 'full',
            'ring_solution':solve_ring(args.quick),'coherent_matter_solution':solve_matter(args.quick),
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'archived_multipair_sha256':hashlib.sha256(before).hexdigest(),
            'result':'Explicit path-weight solution and elimination of the known constant simultaneous-low-momentum source artifact, not a complete gravity theory',
            'source_sha256':{str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in (Path(__file__),ROOT/'construction/area_stress_model.py',ROOT/'construction/check_geometry_response.py',ROOT/'microscopic/chiral_matter_model.py')},
            'historical_results_modified':False}
    assert original.read_bytes()==before
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
