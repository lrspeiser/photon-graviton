#!/usr/bin/env python3
"""One free-slab fermion loop for photon and frame sources on the same 3D grid.

This is a source-response calculation for the existing domain-wall Hamiltonian,
not the interacting finite-spin gauge determinant or a completed quantum-gravity
model. All wall, mirror and doubler states of the chosen slab are retained.
No constant response is removed from the physical kernel. Differencing the
kernel is used only to *measure* derivative coefficients.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'microscopic'))
import chiral_matter_model as matter


def build_spectrum(L, width, charge):
    momenta = 2*np.pi*(np.arange(L)+0.5)/L-np.pi  # antiperiodic matter; no zero-mode ambiguity
    p = np.stack(np.meshgrid(momenta, momenta, momenta, indexing='ij'), axis=-1).reshape(-1, 3)
    identity = np.eye(width)
    gammas = [np.kron(identity, g) for g in matter.SPATIAL_GAMMAS]
    gamma5 = np.kron(identity, matter.GAMMA5)
    base = matter.slab_hamiltonian((0., 0., 0.), matter.wilson_mass(charge), width)
    H = np.broadcast_to(base, (len(p), *base.shape)).copy()
    for axis in range(3):
        H += np.sin(p[:, axis])[:, None, None]*gammas[axis]
        H += (1-np.cos(p[:, axis]))[:, None, None]*gamma5
    source_residual = float(np.max(np.abs(H[0]-matter.slab_hamiltonian(
        p[0], matter.wilson_mass(charge), width))))
    assert source_residual < 1.e-12
    energies, vectors = np.linalg.eigh(H)
    half = 2*width
    assert np.all(energies[:, :half] < 0) and np.all(energies[:, half:] > 0)
    return p, H, energies, vectors, gammas, gamma5, source_residual


def bubble(energies, vectors, shifted_e, shifted_v, vertex, half, batch=256):
    static, temporal = 0., 0.
    for start in range(0, len(energies), batch):
        sl = slice(start, start+batch)
        V = vertex[sl] if vertex.ndim == 3 else vertex
        transformed = shifted_v[sl].conj().swapaxes(-2, -1) @ V @ vectors[sl]
        w1 = abs(transformed[:, half:, :half])**2
        w2 = abs(transformed[:, :half, half:])**2
        d1 = shifted_e[sl, half:, None]-energies[sl, None, :half]
        d2 = energies[sl, None, half:]-shifted_e[sl, :half, None]
        static -= float(np.sum(w1/d1)+np.sum(w2/d2))
        temporal += float(np.sum(w1/d1**3)+np.sum(w2/d2**3))
    return static/len(energies), temporal/len(energies)


def contact(vectors, second, half):
    total = 0.
    for start in range(0, len(vectors), 256):
        sl = slice(start, start+256)
        occ = vectors[sl, :, :half]
        total += float(np.einsum('pij,pij->', occ.conj(), second[sl] @ occ).real)
    return total/len(vectors)


def one_volume(L, width, max_mode):
    modes = range(max_mode+1)
    totals = {n: {'photon_static': 0., 'frame_static': 0.,
                  'photon_temporal': 0., 'frame_temporal': 0.,
                  'longitudinal_static': 0., 'frame_volume_preserving_static': 0.} for n in modes}
    ward_residual, source_residual = 0., 0.
    charge_groups = Counter(abs(q) for q in matter.CHARGES)
    for charge, multiplicity in charge_groups.items():
        p, H, E, vectors, gamma, g5, residual = build_spectrum(L, width, charge)
        source_residual = max(source_residual, residual)
        half = 2*width
        photon = charge*(np.cos(p[:, 1])[:, None, None]*gamma[1]
                         + np.sin(p[:, 1])[:, None, None]*g5)
        frame = (np.sin(p[:, 1])[:, None, None]*gamma[1]
                 - np.sin(p[:, 2])[:, None, None]*gamma[2])/np.sqrt(2.)
        photon_second = charge**2*(-np.sin(p[:, 1])[:, None, None]*gamma[1]
                                   +np.cos(p[:, 1])[:, None, None]*g5)
        longitudinal_second = charge**2*(-np.sin(p[:, 0])[:, None, None]*gamma[0]
                                         +np.cos(p[:, 0])[:, None, None]*g5)
        frame_second = .5*(np.sin(p[:, 1])[:, None, None]*gamma[1]
                             +np.sin(p[:, 2])[:, None, None]*gamma[2])
        cf = contact(vectors, frame_second, half)
        ct = contact(vectors, photon_second, half)
        cl = contact(vectors, longitudinal_second, half)
        for n in modes:
            shift = np.roll(np.arange(L**3).reshape(L,L,L), -n, axis=0).ravel()
            Eq, Vq = E[shift], vectors[shift]
            ps, pt = bubble(E, vectors, Eq, Vq, photon, half)
            fs, ft = bubble(E, vectors, Eq, Vq, frame, half)
            values = totals[n]
            values['photon_static'] += multiplicity*(ct+ps)
            values['frame_static'] += multiplicity*fs  # e enters the actual slab linearly: H_ee=0
            values['frame_volume_preserving_static'] += multiplicity*(cf+fs)
            values['photon_temporal'] += multiplicity*pt
            values['frame_temporal'] += multiplicity*ft
            q = 2*np.pi*n/L
            longitudinal = charge*(np.cos(p[:, 0]+q/2)[:, None, None]*gamma[0]
                                    +np.sin(p[:, 0]+q/2)[:, None, None]*g5)
            if n:
                identity_error = 2*np.sin(q/2)*longitudinal-charge*(H[shift]-H)
                ward_residual = max(ward_residual, float(np.max(abs(identity_error))))
                ls, _ = bubble(E, vectors, Eq, Vq, longitudinal, half)
                values['longitudinal_static'] += multiplicity*(cl+ls)
        del H, E, vectors, photon, frame, photon_second, longitudinal_second, frame_second
    rows = [{'momentum_index': n, 'momentum': 2*np.pi*n/L,
             'lattice_momentum': 2*np.sin(np.pi*n/L), **totals[n]} for n in modes]
    coefficients = {}
    for sector in ('photon', 'frame'):
        qhat = np.asarray([r['lattice_momentum'] for r in rows[1:]])
        differences = np.asarray([r[sector+'_static']-rows[0][sector+'_static'] for r in rows[1:]])
        design = np.column_stack((qhat**2, qhat**4))
        c, _, rank, _ = np.linalg.lstsq(design, differences, rcond=None)
        assert rank == 2
        narrow, *_ = np.linalg.lstsq(design[:3], differences[:3], rcond=None)
        coefficients[sector] = {
            'masslike_static_response_NOT_subtracted': rows[0][sector+'_static'],
            'd_Pi_d_omega_squared_at_zero': rows[0][sector+'_temporal'],
            'spatial_q_squared_coefficient': float(c[0]),
            'spatial_q_fourth_coefficient': float(c[1]),
            'three_mode_spatial_coefficient': float(narrow[0]),
            'window_change_in_spatial_coefficient': float(narrow[0]-c[0]),
            'maximum_fit_residual': float(np.max(abs(design@c-differences))),
            'fourth_mode_heldout_residual': float(abs(design[-1]@narrow-differences[-1])),
        }
    longitudinal = max(abs(r['longitudinal_static']) for r in rows[1:])
    assert ward_residual < 1.e-10
    assert longitudinal < 2.e-9
    # A nonzero masslike response is retained, not hidden behind an IR fit.
    return {'L': L, 'spatial_dimension': 3, 'slab_width': width,
            'source_H_reproduction_residual': source_residual,
            'operator_gauge_Ward_residual': ward_residual,
            'maximum_longitudinal_static_response': longitudinal,
            'momentum_rows': rows, 'infrared_response_coefficients': coefficients,
            'volume_preserving_frame_response_at_zero': rows[0]['frame_volume_preserving_static']}


def uniform_background_control():
    """Independent ground-energy finite differences validate contact + bubble."""
    totals = {name: [0., 0.] for name in ('photon', 'frame_linear', 'frame_volume_preserving')}
    for charge, multiplicity in Counter(abs(q) for q in matter.CHARGES).items():
        p,H,E,U,gamma,g5,_ = build_spectrum(4,4,charge)
        half = 8
        Sy = np.sin(p[:,1])[:,None,None]*gamma[1]
        Sz = np.sin(p[:,2])[:,None,None]*gamma[2]
        photon = charge*(np.cos(p[:,1])[:,None,None]*gamma[1]+np.sin(p[:,1])[:,None,None]*g5)
        photon_second = charge**2*(-Sy+np.cos(p[:,1])[:,None,None]*g5)
        frame = (Sy-Sz)/np.sqrt(2.)
        bs,_ = bubble(E,U,E,U,photon,half)
        fs,_ = bubble(E,U,E,U,frame,half)
        expected = {'photon':bs+contact(U,photon_second,half),
                    'frame_linear':fs,
                    'frame_volume_preserving':fs+contact(U,.5*(Sy+Sz),half)}
        E0 = float(np.mean(np.sum(E[:,:half],axis=1)))
        for name in totals:
            estimates = []
            for step in (.002,.001):
                perturbed = []
                for amplitude in (-step,step):
                    if name == 'photon':
                        Hy = H + (np.sin(p[:,1]+charge*amplitude)-np.sin(p[:,1]))[:,None,None]*gamma[1]
                        Hy += (np.cos(p[:,1])-np.cos(p[:,1]+charge*amplitude))[:,None,None]*g5
                    elif name == 'frame_linear':
                        Hy = H+amplitude*frame
                    else:
                        Hy = H+np.expm1(amplitude/np.sqrt(2.))*Sy+np.expm1(-amplitude/np.sqrt(2.))*Sz
                    perturbed.append(float(np.mean(np.sum(np.linalg.eigvalsh(Hy)[:,:half],axis=1))))
                estimates.append((sum(perturbed)-2*E0)/step**2)
            extrapolated = (4*estimates[1]-estimates[0])/3
            totals[name][0] += multiplicity*expected[name]
            totals[name][1] += multiplicity*extrapolated
    report={name:{'bubble_plus_contact':v[0], 'ground_energy_second_derivative':v[1],
                  'absolute_difference':abs(v[0]-v[1])} for name,v in totals.items()}
    assert max(v['absolute_difference'] for v in report.values()) < 2.e-4
    return report


def run(quick=False):
    sizes, width = ((8, 12), 4) if quick else ((12, 16, 24), 6)
    volumes = [one_volume(L, width, 4) for L in sizes]
    return {'audit_pass': True, 'mode': 'quick' if quick else 'full',
            'matter_content': list(matter.CHARGES),
            'normalization': 'Unrescaled existing slab_hamiltonian, used identically for photon and frame vertices',
            'photon_vertex': 'Peierls derivative of BOTH kinetic and Wilson hopping',
            'frame_vertex': 'Exact existing e_yy=-e_zz derivative, divided by sqrt(2)',
            'new_counterterms': 0, 'mirrors_and_doublers_retained': True,
            'volumes': volumes,
            'uniform_background_control': uniform_background_control(),
            'infrared_gate': {
                'measured_several_momenta_and_volumes': True,
                'matched_bare_kinetic_normalizations_available': False,
                'frame_stress_Ward_identity_established': False,
                'physical_speed_comparison_permitted': False,
                'finite_momentum_shift_relabelled_as_speed': False,
                'volume_preserving_frame_contact_included_as_control': True,
            },
            'common_cone_established': False,
            'minimal_frame_source_is_massless': bool(all(abs(v['infrared_response_coefficients']['frame']['masslike_static_response_NOT_subtracted']) < 1.e-8 for v in volumes)),
            'blocked_by': [
                'The existing slab frame derivative has an uncancelled zero-momentum shear response even on the determinant-one exponential path; its derived contact is retained and no compensating counterterm was added',
                'No microscopic map fixes this slab normalization relative to the interacting photon and dressed-frame kinetic residues',
                'No stress-energy Ward identity or nonlinear frame coupling is supplied by the existing slab',
                'The photon finite-spin 3D phase is not established by the reduced chains'],
            'claim_boundary': 'A complete free-slab fermion loop with identical matter in both responses, not a self-consistent interacting gauge/frame pole calculation or mirror completion'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run(args.quick)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps(report, indent=2, sort_keys=True))
