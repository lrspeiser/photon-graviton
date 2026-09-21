#!/usr/bin/env python3
"""Necessary local nonlinear gates; do not relabel these as full gravity closure.

Derive a longer-range energy current from actual completed moves, derive the
connection commutator from holonomy multiplication, audit local volume gauge
fixing, and test candidate nonlinear spatial/scalar brackets without projection.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import expm, null_space

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'microscopic'))
from all_sector_move_common import X_ROOT, completed_hamiltonian
from check_constrained_local_reduction import (
    H_BASIS, TRACE_VECTOR, gauge_map, scalar_constraint_vector, REDUCED_POISSON,
    KINETIC_METRIC, fierz_pauli_matrix)


def completed_move_energy_algebra():
    e = np.eye(3)
    m01, m12 = np.outer(e[0], e[1]), np.outer(e[1], e[2])
    h01 = completed_hamiltonian(np.zeros((3,3)), [m01], X_ROOT)
    h12 = completed_hamiltonian(np.zeros((3,3)), [m12], X_ROOT)
    j01, j12 = m01-m01.T, m12-m12.T
    j02 = np.outer(e[0], e[2])-np.outer(e[2], e[0])
    commutator = h01@h12-h12@h01
    derived = X_ROOT**2*j02-X_ROOT**3*(j01+j12)
    residual = float(np.max(abs(commutator-derived)))
    assert residual < 1.e-14
    # Nearest-neighbor currents alone cannot represent the induced 0->2 term.
    outside = commutator+X_ROOT**3*(j01+j12)
    H = h01+h12
    state = np.asarray([1., 1j, 0.5], dtype=complex)
    state /= np.linalg.norm(state)
    energy = float(np.vdot(state, H@state).real)
    conservation = []
    for t in (0., 1., 5., 20.):
        evolved = expm(-1j*t*H)@state  # no projection at any time
        conservation.append(max(abs(np.vdot(evolved, H@evolved).real-energy),
                                abs(np.vdot(evolved, evolved).real-1)))
    assert max(conservation) < 1.e-12
    return {'energy_density_definition': 'Each complete h_M, including its diagonal partner',
            'derived_commutator': '[h01,h12]=x^2 J02-x^3(J01+J12)',
            'coefficient_x': X_ROOT, 'identity_residual': residual,
            'nearest_neighbor_current_closure_residual': float(np.linalg.norm(outside)),
            'unprojected_total_energy_and_norm_error': float(max(conservation)),
            'meaning': 'A necessary new composite current is derived, not fitted; '
                       'these transport densities have not been identified with full gravitational constraints'}


def holonomy_audit():
    rng = np.random.default_rng(26092026)
    A, B = rng.normal(size=(2,3,3))
    A, B = 1j*(A+A.T)/2, 1j*(B+B.T)/2
    A /= np.linalg.norm(A)
    B /= np.linalg.norm(B)
    comm = A@B-B@A
    epsilons = np.asarray([.02, .03, .05, .08])
    full, missing, square = [], [], []
    for epsilon in epsilons:
        W = expm(epsilon*A)@expm(epsilon*B)@expm(-epsilon*A)@expm(-epsilon*B)
        full.append(float(np.linalg.norm(W-np.eye(3)-epsilon**2*comm)))
        missing.append(float(np.linalg.norm(W-np.eye(3))))
        H = completed_hamiltonian(np.zeros((3,3)), [W], X_ROOT)
        constant = (-2*X_ROOT+2*X_ROOT**2)*np.eye(3)
        square.append(float(np.linalg.norm(H-constant)))
    powers = {name: float(np.polyfit(np.log(epsilons), np.log(values), 1)[0])
              for name, values in [('BCH_remainder',full),('missing_commutator',missing),('completed_unitary_loop',square)]}
    assert 2.9 < powers['BCH_remainder'] < 3.1
    assert 1.9 < powers['missing_commutator'] < 2.1
    assert 3.8 < powers['completed_unitary_loop'] < 4.2
    return {'derived_curvature_term': 'W=I+epsilon^2[A,B]+O(epsilon^3) for a constant connection',
            'powers': powers,
            'decision': 'The Hermitian completed unitary loop starts at curvature squared. '
                        'A first-order frame/area insertion is still needed; this test does not derive it.'}


def cubic_kinetic_bootstrap():
    """Solve affine spatial-density covariance, not a preinserted GR action.

    The nonlinear tensor-density transformation is an explicit assumption to be
    matched to the finite junction algebra. Under that assumption the four
    coefficients are uniquely fixed by the already selected quadratic metric.
    """
    rng = np.random.default_rng(5742)
    design, targets = [], []
    for _ in range(160):
        pi, B = rng.normal(size=(2,3,3))
        pi = (pi+pi.T)/2
        dh = B+B.T
        dpi = -B.T@pi-pi@B+np.trace(B)*pi
        T2 = np.trace(pi@pi)-.5*np.trace(pi)**2
        dT2 = np.sum((2*pi-np.trace(pi)*np.eye(3))*dpi)
        design.append([np.trace(dh)*np.trace(pi@pi),
                       np.trace(dh)*np.trace(pi)**2,
                       np.trace(dh@pi@pi), np.trace(pi)*np.trace(dh@pi)])
        targets.append(np.trace(B)*T2-dT2)
    design, targets = np.asarray(design),np.asarray(targets)
    coefficients,_,rank,_ = np.linalg.lstsq(design[:80],targets[:80],rcond=None)
    heldout = float(np.max(abs(design[80:]@coefficients-targets[80:])))
    assert rank == 4 and heldout < 1.e-11
    assert np.max(abs(coefficients-[-.5,.25,2.,-1.])) < 1.e-12
    return {'basis':['tr(h) tr(pi^2)','tr(h) (tr pi)^2','tr(h pi^2)','tr(pi) tr(h pi)'],
            'coefficients':coefficients.tolist(), 'rank':int(rank),
            'heldout_covariance_residual':heldout,
            'derived_T3':'2 tr(h pi^2) - tr(pi) tr(h pi) - 1/2 tr(h) tr(pi^2) + 1/4 tr(h) (tr pi)^2',
            'assumption':'g transforms as a covariant spatial metric; pi as a weight-one contravariant density',
            'microscopic_junction_matching_established':False}


def local_volume_rank():
    rows = []
    for vector in ((1,0,0),(1,2,3),(2,-1,1)):
        k = np.asarray(vector, dtype=float)
        C = np.zeros((6,12))
        C[:3,6:] = gauge_map(k).T
        C[3,:6] = scalar_constraint_vector(k)
        C[4,:6], C[5,6:] = TRACE_VECTOR, TRACE_VECTOR
        rank = int(np.linalg.matrix_rank(C))
        bracket = C@REDUCED_POISSON@C.T
        second_class = int(np.linalg.matrix_rank(bracket))
        first_class = rank-second_class
        configurations = int((12-2*first_class-second_class)/2)
        assert (rank, second_class, first_class, configurations) == (6,4,2,2)
        # Dirac-bracket evolution incorporates the second-class conditions in
        # the generator itself. Initial data are constrained once, never reprojected.
        omega = REDUCED_POISSON
        omega_d = omega-omega@C.T@np.linalg.pinv(bracket)@C@omega
        Hessian = np.zeros((12,12))
        Hessian[:6,:6] = fierz_pauli_matrix(k)
        Hessian[6:,6:] = KINETIC_METRIC
        initial = null_space(C)@np.arange(1.,7.)[:null_space(C).shape[1]]
        initial /= np.linalg.norm(initial)
        energy0 = initial@Hessian@initial/2
        errors, energy_errors = [], []
        for time in (0.,.1,1.,3.):
            evolved = expm(time*omega_d@Hessian)@initial
            errors.append(float(np.max(abs(C@evolved))))
            energy_errors.append(float(abs(evolved@Hessian@evolved/2-energy0)))
        assert max(errors) < 1.e-10 and max(energy_errors) < 1.e-10
        rows.append({'k': list(vector), 'constraint_rank': rank,
                     'linear_unprojected_constraint_error':max(errors),
                     'linear_unprojected_energy_error':max(energy_errors),
                     'second_class_constraints': second_class,
                     'remaining_first_class_constraints': first_class,
                     'physical_configurations': configurations,
                     'bracket_matrix': bracket.tolist()})
    C0 = np.zeros((2,12))
    C0[0,:6], C0[1,6:] = TRACE_VECTOR, TRACE_VECTOR
    rank0 = int(np.linalg.matrix_rank(C0))
    r0 = int(np.linalg.matrix_rank(C0@REDUCED_POISSON@C0.T))
    homogeneous_count = int((12-2*(rank0-r0)-r0)/2)
    assert homogeneous_count == 5
    return {'nonzero_momentum': rows,
            'homogeneous_configurations_after_global_pair': homogeneous_count,
            'decision': 'Local determinant/trace conditions gauge-fix two old gauge directions. '
                        'They cannot be appended as an independent second-class pair while '
                        'claiming all four old constraints remain first class.'}


def derivative(field, axis, spectral=False):
    L = field.shape[axis]
    if not spectral:
        return (np.roll(field,-1,axis)-np.roll(field,1,axis))/(4*np.pi/L)
    k = np.fft.fftfreq(L, 1./L)
    shape = [1]*field.ndim
    shape[axis] = L
    return np.fft.ifft(1j*k.reshape(shape)*np.fft.fft(field,axis=axis),axis=axis).real


def advect(vector, field, spectral):
    return sum(vector[...,a].reshape(vector.shape[:3]+(1,)*(field.ndim-3))
               *derivative(field,a,spectral) for a in range(3))


def lie(vector, metric, spectral):
    D = np.stack([derivative(vector,a,spectral) for a in range(3)], axis=-2)
    return (advect(vector,metric,spectral)+np.einsum('...ik,...kj->...ij',D,metric)
            +np.einsum('...ik,...jk->...ij',metric,D))


def bracket_test(L, spectral=False):
    x,y,z = 2*np.pi*np.indices((L,L,L))/L
    xi = np.stack((np.sin(y)+.2*np.cos(z)+.2*np.sin(x),np.sin(z)+.3*np.cos(x)+.1*np.sin(y),np.sin(x)+.15*np.sin(z)),axis=-1)
    eta = np.stack((np.cos(z),np.cos(x)+.1*np.sin(y),np.cos(y)),axis=-1)
    h = np.zeros((L,L,L,3,3))
    h[...,0,0],h[...,1,1],h[...,2,2] = .07*np.sin(x+y),.05*np.cos(y+z),.03*np.sin(z+x)
    h[...,0,1]=h[...,1,0]=.02*np.cos(x-z)
    h[...,1,2]=h[...,2,1]=.01*np.sin(x+y)
    g = h+np.eye(3)
    comm = lie(xi,lie(eta,g,spectral),spectral)-lie(eta,lie(xi,g,spectral),spectral)
    zeta = advect(xi,eta,spectral)-advect(eta,xi,spectral)
    target = lie(zeta,g,spectral)
    relative = float(np.linalg.norm(comm-target)/np.linalg.norm(target))
    # Quadratic self-energy from the existing reduced first-order frame sector.
    pi = h*2
    pi[...,0,0] += .1*np.cos(x)
    pi[...,0,1] += .1*np.sin(x)*np.cos(y)
    pi[...,1,0] = pi[...,0,1]
    N, M = 1+.2*np.sin(x)+.3*np.cos(y), 1+.4*np.sin(2*x)+.2*np.cos(z)
    gradN = np.stack([derivative(N,a,spectral) for a in range(3)],axis=-1)
    gradM = np.stack([derivative(M,a,spectral) for a in range(3)],axis=-1)
    HessN = np.stack([derivative(gradN,a,spectral) for a in range(3)],axis=-2)
    HessM = np.stack([derivative(gradM,a,spectral) for a in range(3)],axis=-2)
    trpi = np.trace(pi,axis1=-2,axis2=-1)
    kinetic_gradient = pi-.5*trpi[...,None,None]*np.eye(3)
    AN = HessN-np.trace(HessN,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)
    AM = HessM-np.trace(HessM,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)
    direct = float(np.mean(np.sum((M[...,None,None]*AN-N[...,None,None]*AM)*kinetic_gradient,axis=(-2,-1))))
    divpi = sum(derivative(pi[...,a,:],a,spectral) for a in range(3))
    shift = N[...,None]*gradM-M[...,None]*gradN
    expected = float(np.mean(np.sum(divpi*shift,axis=-1)))
    # Include both kinetic and gradient energy; eliminated auxiliary connection
    # energy is already represented by this Fierz-Pauli form, not counted twice.
    dh = np.stack([derivative(h,a,spectral) for a in range(3)],axis=-3)
    divh = sum(derivative(h[...,a,:],a,spectral) for a in range(3))
    traceh = np.trace(h,axis1=-2,axis2=-1)
    gradtrace = np.stack([derivative(traceh,a,spectral) for a in range(3)],axis=-1)
    T = .5*(np.sum(pi*pi,axis=(-2,-1))-.5*trpi**2)
    V = .5*(np.sum(dh*dh,axis=(-3,-2,-1))-2*np.sum(divh*divh,axis=-1)
             +2*np.sum(divh*gradtrace,axis=-1)-np.sum(gradtrace*gradtrace,axis=-1))
    volume_variation = np.trace(lie(xi,np.broadcast_to(np.eye(3),g.shape),spectral),axis1=-2,axis2=-1)
    return {'L': L, 'derivative': 'spectral smooth-mode control' if spectral else 'local central difference',
            'vector_bracket_relative_residual': relative,
            'scalar_self_energy_bracket_direct': direct,
            'scalar_self_energy_bracket_target': expected,
            'scalar_self_energy_bracket_error': abs(direct-expected),
            'quadratic_frame_self_energy_mean': float(np.mean(T+V)),
            'local_volume_change_under_old_vector_generator_max': float(np.max(abs(volume_variation)))}


def run(quick=False):
    sizes = (12,16,24) if quick else (12,16,24,32,48)
    rows = [bracket_test(L) for L in sizes]
    control = bracket_test(16,True)
    assert control['vector_bracket_relative_residual'] < 1.e-12
    assert control['scalar_self_energy_bracket_error'] < 1.e-12
    assert rows[-1]['vector_bracket_relative_residual'] < rows[0]['vector_bracket_relative_residual']
    assert rows[0]['vector_bracket_relative_residual'] > .001
    assert rows[0]['local_volume_change_under_old_vector_generator_max'] > .1
    assert rows[0]['scalar_self_energy_bracket_error'] > 1.e-6
    power = float(np.polyfit(np.log(sizes), np.log([r['vector_bracket_relative_residual'] for r in rows]),1)[0])
    return {'audit_pass': True, 'full_nonlinear_gravity_closed': False,
            'mode': 'quick' if quick else 'full',
            'completed_move_energy_algebra': completed_move_energy_algebra(),
            'holonomy': holonomy_audit(), 'local_volume': local_volume_rank(),
            'cubic_kinetic_bootstrap': cubic_kinetic_bootstrap(),
            'candidate_bracket_refinement': rows, 'smooth_spectral_control': control,
            'local_bracket_error_power': power,
            'scalar_source_candidate': 'S(h)+epsilon*[T_frame(pi)+V_FP(h)]; '
                                       'at first order {C[N],C[M]} gives the momentum source in the smooth control',
            'not_established': ['junction derivation of the full nonlinear frame/area vertex',
                                'closed finite-lattice gravitational constraint algebra',
                                'unprojected preservation of all nonlinear gravitational constraints',
                                'nonlinear positive physical spectrum or vacuum cancellation'],
            'decision': 'Retain the derived commutator/current terms. The naive local central-difference '
                        'nonlinear extension fails exact closure; the volume pair changes constraint '
                        'classification. Neither defect is repaired by manually projecting evolved states.'}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=run(args.quick)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))
