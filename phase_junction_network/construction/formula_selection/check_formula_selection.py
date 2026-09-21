#!/usr/bin/env python3
"""Bounded test of whether transported projectors select a unique gravity formula.

All alternatives are candidate local operators, not proved nonlinear gravity
solutions. A covariance identity is not a spacetime-constraint derivation.
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

CONSTRUCTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CONSTRUCTION))
import area_stress_model as model
from check_area_kernel import quadratic_data

PARENT = 'b6f55a632e707ede4b872cd0092c85f485b459c5'
PROTOCOL = 'e72c53296622ea25494f88d344c6b8a1a84cd87d'
X = model.X_ROOT


def norm(a):
    return float(np.linalg.norm(a, 2))


def projectors(n, block):
    if n % block:
        raise ValueError('Block size must divide dimension')
    return [np.diag(np.repeat(np.arange(n//block) == a, block).astype(float))
            for a in range(n//block)]


def dephase(h, ps):
    return sum(p @ h @ p for p in ps)


def completion(h, ps):
    d = dephase(h, ps)
    return dephase(h @ h, ps) - .5*d @ d


def fixed_derivative(h, dh, ps):
    d, dd = dephase(h, ps), dephase(dh, ps)
    return dephase(h @ dh + dh @ h, ps) - .5*(d @ dd + dd @ d)


def projector_derivative(h, ps, dps):
    def delta(a):
        return sum(dp @ a @ p + p @ a @ dp for p, dp in zip(ps, dps))
    d, dd = dephase(h, ps), delta(h)
    return delta(h @ h) - .5*(d @ dd + dd @ d)


def projector_tests(seed, samples):
    rng = np.random.default_rng(seed)
    maximum = dict(covariance=0., finite_difference=0., block_recovery=0.)
    controls = []
    for _ in range(samples):
        n, block = 8, 4
        r = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
        h = (r+r.conj().T)/8
        r = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
        d = (r+r.conj().T)/8
        ps = projectors(n,block)
        dh = 1j*(d @ h-h @ d)
        dps = [1j*(d @ p-p @ d) for p in ps]
        c = completion(h,ps)
        fixed = fixed_derivative(h,dh,ps)
        moving = projector_derivative(h,ps,dps)
        target = 1j*(d @ c-c @ d)
        maximum['covariance'] = max(maximum['covariance'],norm(fixed+moving-target))
        step = 1.e-4
        def fd(step):
            ss = [expm(1j*t*d) for t in (-step,step)]
            cs = [completion(s @ h @ s.conj().T, [s @ p @ s.conj().T for p in ps]) for s in ss]
            return (cs[1]-cs[0])/(2*step)
        measured = (4*fd(step/2)-fd(step))/3
        maximum['finite_difference'] = max(maximum['finite_difference'],norm(measured-fixed-moving))
        expected = model.full_fock_completion(h,block)['onebody']
        maximum['block_recovery'] = max(maximum['block_recovery'],norm(expected-X*h-X**2*c))
        controls.append(norm(fixed-target))
    assert max(maximum.values()) < 2.e-9
    assert min(controls) > 1.e-3
    return {'samples':samples, 'max_errors':maximum,
            'minimum_fixed_projector_negative_control':min(controls),
            'formula':'delta C = D_H C[delta H] + D_P C[delta P]',
            'coefficient_of_projector_term':1,
            'meaning':'One is fixed by differentiation, not by fitting G or a pole speed.',
            'physical_coordinate_symmetry_proven':False}


def fock_control():
    width = 1
    o = model.slab.slab_hamiltonian((0,0,0),model.slab.wilson_mass(1),width)+3*model.BETA
    k = (-1j*model.GAMMA[0]-model.BETA)/2
    h = np.block([[o,k],[k.conj().T,o]])
    d = np.kron(np.array([[0,-1j],[1j,0]]),np.eye(4))*.37
    ps = projectors(8,4)
    moves = [-ps[0]@h@ps[0]/2, -ps[1]@h@ps[1]/2, -ps[0]@h@ps[1]]
    dh = 1j*(d@h-h@d)
    dps = [1j*(d@p-p@d) for p in ps]
    dms = []
    for a in range(2):
        p,dp=ps[a],dps[a]
        dms.append(-(dp@h@p+p@dh@p+p@h@dp)/2)
    dms.append(-(dps[0]@h@ps[1]+ps[0]@dh@ps[1]+ps[0]@h@dps[1]))
    df = model.bilinear_fock(d)
    hs, ds = [], []
    for move, dm in zip(moves,dms):
        b,db = model.bilinear_fock(move),model.bilinear_fock(dm)
        hs.append(model.completed(b))
        ds.append(-X*(db+db.conj().T)+X**2*(db.conj().T@b+b.conj().T@db+db@b.conj().T+b@db.conj().T))
    full,derivative = sum(hs),sum(ds)
    exact = 1j*(df@full-full@df)
    generated = model.full_fock_completion(h,4)
    residual=norm(derivative-exact)
    assert residual<1.e-10 and norm(full-generated['full'])<1.e-12
    return {'dimension':256,'derivative_covariance_error':residual,
            'quartic_remainder_norm':norm(generated['quartic_fock']),
            'quartic_stress_norm':norm(1j*(df@generated['quartic_fock']-generated['quartic_fock']@df)),
            'interacting_spacetime_constraints_proven':False}


def stress_terms(L,width=2):
    q=2*np.pi/L; p=np.full(3,np.pi/L); pp=p+np.array([q,0,0])
    gamma=np.kron(np.eye(width),model.GAMMA[0]); beta=np.kron(np.eye(width),model.BETA)
    wall=model.slab.slab_hamiltonian((0,0,0),model.slab.wilson_mass(1),width)
    o=wall+3*beta; anti=o@beta+beta@o
    k=(-1j*gamma-beta)/2
    z=-(1.5+.5*np.exp(-1j*q))
    def a(q):return .25j*(1+np.exp(1j*q))*gamma+(.75+.25*np.exp(1j*q))*beta
    aa,ad=a(q),a(-q).conj().T
    jcomp=.5*z*anti+aa@k.conj().T+k@ad+np.exp(-1j*q)*(ad@k+k.conj().T@aa)
    geo_c=1j*np.sin(q)*jcomp
    # Independent finite Fourier integration of i[D,H]. The polynomial has
    # degree at most two, so sixteen quadrature points integrate exactly.
    t=2*np.pi*np.arange(16)/16
    def hh(t):return o+np.sin(t)[:,None,None]*gamma-np.cos(t)[:,None,None]*beta
    comm=1j*.5*(np.sin(t)+np.sin(t+q))[:,None,None]*(hh(t)-hh(t+q))
    jj={r:np.mean(np.exp(-1j*t*r)[:,None,None]*comm,axis=0) for r in (0,1,-1)}
    derivative_c=.5*(o@jj[0]+jj[0]@o)
    derivative_c+=jj[1]@k.conj().T+k@(np.exp(1j*q)*jj[-1])+jj[-1]@k+k.conj().T@(np.exp(-1j*q)*jj[1])
    analytic_c=-.25j*np.sin(q)*anti
    assert norm(derivative_c-analytic_c)<1.e-12
    # C is spatially constant, so [D,C]=0 and the P derivative is -D_H C[i[D,H]].
    projector=-derivative_c
    d0=np.exp(1j*p[0])-1; d1=np.exp(1j*pp[0])-1
    jkin=-.5*(np.sin(p[0])+np.sin(pp[0]))*gamma
    jwil=-.5*((1-np.cos(p[0]))+(1-np.cos(pp[0]))+np.conj(d1)*d0)*beta
    geo=1j*np.sin(q)*(jkin+jwil)
    hp=model.slab.slab_hamiltonian(p,model.slab.wilson_mass(1),width)
    hpp=model.slab.slab_hamiltonian(pp,model.slab.wilson_mass(1),width)
    target=1j*.5*(np.sin(p[0])+np.sin(pp[0]))*(hp-hpp)
    before=X*(geo-target)+X**2*geo_c
    after=X*(geo-target)+X**2*(geo_c+projector)
    return {'L':L,'q':q,'before_error_over_q':norm(before)/q,'after_error_over_q':norm(after)/q,
            'original_limit':X**2*norm(anti+3*np.eye(len(o))),
            'transported_projector_limit':X**2*norm(.75*anti+3*np.eye(len(o))),
            'fourier_derivative_error':norm(derivative_c-analytic_c)}


def actual_graph_control():
    L=3; sites,_=model.lattice(L); n=L**3; block=4
    frames=np.repeat(np.eye(3)[None],n,axis=0)
    nu,gamma,g,_=model.geometry(frames)
    links=np.broadcast_to(np.eye(4),(n,3,4,4)).copy()
    result=model.matter_matrix(L,1,1,nu,gamma,g,links)
    h=result['target']; ps=projectors(len(h),block)
    coords=np.asarray(sites); q=2*np.pi/L
    xi=np.diag(np.repeat(np.cos(q*coords[:,0]),block))
    shift=result['shifts'][0]; k=(shift-shift.conj().T)/(2j)
    d=(xi@k+k@xi)/2
    dh=1j*(d@h-h@d); dps=[1j*(d@p-p@d) for p in ps]
    c=completion(h,ps)
    cp=projector_derivative(h,ps,dps)
    ch=fixed_derivative(h,dh,ps)
    target=1j*(d@c-c@d)
    assert norm(c-result['completion'])<1.e-11
    assert norm(cp+ch-target)<1.e-10
    p0=np.full(3,np.pi/L); p1=p0+np.array([q,0,0])
    v0=np.kron((np.exp(1j*coords@p0)/np.sqrt(n))[:,None],np.eye(4))
    v1=np.kron((np.exp(1j*coords@p1)/np.sqrt(n))[:,None],np.eye(4))
    o=model.slab.slab_hamiltonian((0,0,0),model.slab.wilson_mass(1),1)+3*model.BETA
    analytic=.125j*np.sin(q)*(o@model.BETA+model.BETA@o) # cos source has amplitude 1/2
    error=norm(v1.conj().T@cp@v0-analytic)
    assert error<1.e-10
    return {'dimension':len(h),'projector_stress_graph_vs_symbol_error':error,
            'simultaneous_orbit_derivative_error':norm(cp+ch-target),
            'physical_geometric_source_not_replaced_by_conjugation':True}


def polynomial(w,coefficients):
    if not coefficients or any(n<1 or int(n)!=n for n in coefficients):
        raise ValueError('Positive integer windings required')
    return sum(a*(np.linalg.matrix_power(w,int(n))-np.eye(len(w))) for n,a in coefficients.items())


def poly_energy(b,w,coefficients):
    p=polynomial(w,coefficients)
    return float(np.trace(model.completed(b@p/(2j))).real/len(b))


def formula_tests(seed,samples):
    rules=[{'name':'single_loop','coefficients':{1:1.}},
           {'name':'double_loop_normalized','coefficients':{2:.5}},
           {'name':'equal_derivative_weight_mixture','coefficients':{1:.5,2:.25}}]
    rng=np.random.default_rng(seed)
    maximum=dict(trace_formula=0.,covariance=0.,orientation=0.,hermiticity=0.,flat_energy=0.)
    minimum=0.
    for _ in range(samples):
        frames=expm(rng.normal(size=(3,3))*.1)
        b=model.geometry(frames[None])[3][0,0,1]
        r=rng.normal(size=3); r/=np.linalg.norm(r)
        f=np.einsum('a,aij->ij',r,model.SPIN)*.37
        w=expm(1j*f)
        s=expm(1j*np.einsum('a,aij->ij',rng.normal(size=3),model.SPIN))
        for rule in rules:
            co=rule['coefficients']; p=polynomial(w,co)
            h=model.completed(b@p/(2j)); energy=poly_energy(b,w,co)
            formula=np.trace(-X*b@((p-p.conj().T)/(2j))+.5*X**2*(b@b)@(p@p.conj().T)).real/len(b)
            maximum['trace_formula']=max(maximum['trace_formula'],abs(energy-formula))
            maximum['covariance']=max(maximum['covariance'],abs(energy-poly_energy(s@b@s.conj().T,s@w@s.conj().T,co)))
            maximum['orientation']=max(maximum['orientation'],abs(energy-poly_energy(-b,w.conj().T,co)))
            maximum['hermiticity']=max(maximum['hermiticity'],norm(h-h.conj().T))
            maximum['flat_energy']=max(maximum['flat_energy'],abs(poly_energy(b,np.eye(len(b)),co)))
            minimum=min(minimum,float(np.min(np.linalg.eigvalsh(h))))
    assert max(maximum.values())<1.e-12 and minimum>=-.5-1.e-12
    b=.8; sz=np.diag([1.,-1.]); values=[]
    for rule in rules:
        co=rule['coefficients']; moments={j:sum(a*n**j for n,a in co.items()) for j in (1,2,3)}
        m1,m2,m3=[moments[j] for j in (1,2,3)]
        assert abs(m1-1)<1.e-15
        d1,d2,d3,d4=-X*b*m1,X**2*b*b*m1*m1,X*b*m3,X**2*b*b*(3*m2*m2-4*m1*m3)
        odd_estimates=[];even_estimates=[]
        for theta in (.04,.02,.01):
            plus=poly_energy(b*sz,expm(1j*theta*sz),co)
            minus=poly_energy(b*sz,expm(-1j*theta*sz),co)
            odd_estimates.append(((plus-minus)/2-d1*theta)/theta**3)
            even_estimates.append(((plus+minus)/2-d2*theta**2/2)/theta**4)
        odd=(4*odd_estimates[-1]-odd_estimates[-2])/3
        even=(4*even_estimates[-1]-even_estimates[-2])/3
        assert abs(odd-d3/6)<1.e-8 and abs(even-d4/24)<1.e-8
        r3,r4=-d3/d1,-d4/d2
        assert abs((r4-r3)-3*(m3/m1-(m2/m1)**2))<1.e-12
        values.append({'name':rule['name'],'path_coefficients':co,'moments':moments,
                       'curvature_derivatives_at_zero':[d1,d2,d3,d4],
                       'R3_minus_E3_over_E1':r3,'R4_minus_E4_over_E2':r4,
                       'nonlinear_weight_variance_identity':r4-r3,
                       'cubic_derivative_numerical_error':abs(6*odd-d3),
                       'quartic_derivative_numerical_error':abs(24*even-d4),
                       'energy_at_phase_0p4':poly_energy(b*sz,expm(.4j*sz),co)})
    assert values[0]['R3_minus_E3_over_E1']!=values[1]['R3_minus_E3_over_E1']
    return {'samples':samples,'operator_errors':maximum,'minimum_move_eigenvalue':minimum,'rules':values,
            'new_sector_adjustments':0,
            'rule_class':'Finite-support coherent repeated-plaquette polynomials with real path weights; not a classification of all theories.'}


def directional_energy(v,L,mode,co):
    from check_area_kernel import H_BASIS
    phase=2*np.pi*np.arange(L)/L
    hc=np.einsum('a,aij->ij',v[:6],np.asarray(H_BASIS)); hs=np.einsum('a,aij->ij',v[6:12],np.asarray(H_BASIS))
    frames=np.asarray([expm(-.5*(np.cos(p)*hc+np.sin(p)*hs)) for p in phase])
    area=model.geometry(frames)[3]
    c=v[12:].reshape(2,3,3);links=np.empty((L,3,4,4),dtype=complex)
    for n in range(L):
        for i in range(3):
            p=phase[n]+np.pi*mode[i]/L
            links[n,i]=expm(1j*np.einsum('a,aij->ij',np.sin(p)*c[0,i]+np.cos(p)*c[1,i],model.SPIN))
    total=0.
    for n in range(L):
        for i in range(3):
            for j in range(i+1,3):
                ni,nj=(n+mode[i])%L,(n+mode[j])%L
                w=links[n,i]@links[ni,j]@links[nj,i].conj().T@links[n,j].conj().T
                total+=poly_energy(area[n,i,j],w,co)/L
    return total


def same_quadratic_kernel(seed,samples):
    rng=np.random.default_rng(seed);L=12;mode=(1,2,1)
    result,_,hessian=quadratic_data(L,mode)
    errors=[]
    for _ in range(samples):
        v=rng.normal(size=30);v/=np.linalg.norm(v)
        target=float(v@hessian@v)
        for co in ({1:1.},{2:.5},{1:.5,2:.25}):
            fd=[]
            for epsilon in (.004,.002):
                fd.append((directional_energy(epsilon*v,L,mode,co)+directional_energy(-epsilon*v,L,mode,co))/epsilon**2)
            estimate=(4*fd[1]-fd[0])/3
            errors.append(abs(estimate-target))
    assert max(errors)<2.e-8
    return {'L':L,'mode':list(mode),'random_directions':samples,
            'maximum_shared_Hessian_directional_error':max(errors),
            'linear_polarizations':result['placement_matched_kinetic']['polarizations_per_signed_momentum'],
            'proof':'Im P(exp(iF))=F+O(F^3), P P_dagger=F^2+O(F^4) when sum n*a_n=1; the traced area energy is identical through total field order two.',
            'full_nonlinear_constraints_shared':False}


def run(quick=False):
    original=CONSTRUCTION.parent/'microscopic/fermionic_multipair_vacuum_results.json'
    original_hash=hashlib.sha256(original.read_bytes()).hexdigest()
    assert original_hash=='9fc66061e31fe53c4b8a016c5bda9eda8a5f00e0d7efa044314a443a1b5d9400'
    rows=[stress_terms(L) for L in ((32,128,512) if quick else (16,32,64,128,256,512,1024))]
    assert abs(rows[-1]['after_error_over_q']/rows[-1]['transported_projector_limit']-1)<1.e-4
    report={'execution_pass':True,'mode':'quick' if quick else 'full','source_parent':PARENT,'protocol_commit':PROTOCOL,
            'projector_derivative':projector_tests(260921,2 if quick else 6),
            'full_Fock_control':fock_control(),'actual_graph_control':actual_graph_control(),
            'coordinate_stress_refinement':rows,
            'area_formula':formula_tests(210926,4 if quick else 32),
            'common_quadratic_kernel':same_quadratic_kernel(92126,1 if quick else 4),
            'physical_status':{'unique_gravity_formula_established':False,
                               'conditional_single_loop_curvature_formula_derived':True,
                               'projector_correction_alone_removes_leading_physical_stress_defect':False,
                               'Newtonian_force_or_G_derived':False,
                               'nonlinear_constraint_closure_or_common_quantum_cone_proven':False},
            'decision':'STOP consistency-only expansion. Transporting the projectors alone retains a leading regulator-stress defect, and equal weak-field data admit different nonlinear path moments. A microscopic derivation of path weights and full regulator transport is required for formula selection.',
            'archived_multipair_sha256':original_hash,
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}}
    report['source_sha256']={str(p.relative_to(CONSTRUCTION.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                            [Path(__file__),CONSTRUCTION/'area_stress_model.py',CONSTRUCTION/'check_area_kernel.py']}
    assert hashlib.sha256(original.read_bytes()).hexdigest()==original_hash
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=run(args.quick)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'execution_pass':report['execution_pass'],'physical_status':report['physical_status'],
                      'last_stress_row':report['coordinate_stress_refinement'][-1],
                      'quadratic':report['common_quadratic_kernel']},indent=2))
