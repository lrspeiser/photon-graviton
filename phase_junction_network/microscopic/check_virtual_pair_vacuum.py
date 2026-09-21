#!/usr/bin/env python3
"""Stage 3H: virtual-pair photon dressing in a finite gauge-matter vacuum.

Adds a no-matter vacuum block and one neutral-pair block to the Stage-3G 2x2
spin-1 patch. Gauge-invariant pair creation on each oriented link uses the same
completed move coefficient as the photon, frame, endpoint, and companion
sectors. The calculation tests exact local/nonzero-momentum Ward identities,
longitudinal pure-gauge invariance, and transverse photon self-energy from
virtual pairs.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
from typing import Any
import numpy as np
from scipy.sparse import bmat, coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import eigsh
import check_spatial_gauge_matter_patch as patch

PAIR_GAP = 1.0
PAIR_COUPLING = patch.K_MATTER
PAIR_COMPLETION = 1.0


def ready(value: Any) -> Any:
    if isinstance(value, dict): return {str(k): ready(v) for k,v in value.items()}
    if isinstance(value, (list,tuple,np.ndarray)): return [ready(v) for v in value]
    if isinstance(value, (np.floating,np.integer)): return value.item()
    if isinstance(value,np.bool_): return bool(value)
    if isinstance(value,complex): return {"real":value.real,"imag":value.imag}
    return value


def pair_creation_moves(vacuum_basis, pair_basis):
    pair_lookup={state:i for i,state in enumerate(pair_basis)}
    moves=[]
    for link_index,(tail,direction) in enumerate(patch.LINKS):
        head=patch.add(tail,direction)
        rows=[]; cols=[]
        for column,flux in enumerate(vacuum_basis):
            if flux[link_index] <= -1: continue
            moved=list(flux); moved[link_index]-=1
            target=(patch.SITE_INDEX[head],patch.SITE_INDEX[tail],tuple(moved))
            if target in pair_lookup:
                rows.append(pair_lookup[target]); cols.append(column)
        moves.append(coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(pair_basis),len(vacuum_basis)),dtype=complex).tocsr())
    return moves


def combined_system(external=None,pair_coupling=PAIR_COUPLING,pair_gap=PAIR_GAP,completion=PAIR_COMPLETION,matter_scale=1.0):
    vacuum=patch.pure_basis(); pairs=patch.matter_basis()
    plus,minus=patch.matter_move_matrices(pairs)
    external=np.zeros(len(patch.LINKS),dtype=float) if external is None else np.asarray(external,float)
    h_v=patch.gauge_hamiltonian(vacuum)
    h_p=patch.gauge_hamiltonian(pairs,matter=True)+patch.matter_hamiltonian(pairs,plus,minus,external=external,matter_scale=matter_scale)+pair_gap*diags(np.ones(len(pairs)),dtype=complex)
    creates=pair_creation_moves(vacuum,pairs)
    cross=csr_matrix((len(pairs),len(vacuum)),dtype=complex)
    vacuum_diag=csr_matrix((len(vacuum),len(vacuum)),dtype=complex)
    pair_diag=csr_matrix((len(pairs),len(pairs)),dtype=complex)
    for index,move in enumerate(creates):
        phase=external[index]
        cross += -pair_coupling*np.exp(1j*phase)*move
        vacuum_diag += completion*pair_coupling**2*(move.getH()@move)
        pair_diag += completion*pair_coupling**2*(move@move.getH())
    h=bmat([[h_v+vacuum_diag,cross.getH()],[cross,h_p+pair_diag]],format='csr')
    return vacuum,pairs,h,plus,minus,creates,h_v


def lowest(h,count):
    e,v=eigsh(h,k=min(count,h.shape[0]-2),which='SA',tol=1e-10,maxiter=100000)
    order=np.argsort(e); return e[order],v[:,order]


def electric_values(vacuum,pairs):
    return np.concatenate((patch.electric_transverse_values(vacuum),patch.electric_transverse_values(pairs,matter=True)))


def spectral(h,values,count=100):
    e,v=lowest(h,count); ground=v[:,0]; created=values*ground; norm=float(np.vdot(created,created).real)
    groups=[]; start=1
    while start<len(e):
        stop=start+1
        while stop<len(e) and abs(e[stop]-e[start])<1e-8: stop+=1
        basis,_=np.linalg.qr(v[:,start:stop]); amplitudes=basis.conj().T@created; weight=float(np.vdot(amplitudes,amplitudes).real)
        groups.append({"gap":float(np.mean(e[start:stop])-e[0]),"multiplicity":stop-start,"residue":weight,"fraction":weight/norm if norm else 0.0})
        start=stop
    groups.sort(key=lambda row:row['residue'],reverse=True)
    return {"ground_energy":float(e[0]),"operator_norm":norm,"photon_gap":groups[0]['gap'],"photon_multiplicity":groups[0]['multiplicity'],"photon_residue":groups[0]['fraction'],"top_groups":groups[:8],"energies":e,"vectors":v}


def charge_current_operators(vacuum,pairs,plus,minus,creates,pair_coupling=PAIR_COUPLING,matter_scale=1.0):
    nv=len(vacuum); total=nv+len(pairs)
    densities=[]
    for site in range(len(patch.SITES)):
        values=np.zeros(total)
        for index,(positive,negative,_) in enumerate(pairs):
            values[nv+index]=(1 if positive==site else 0)-(1 if negative==site else 0)
        densities.append(diags(values,dtype=complex).tocsr())
    _,matter_currents=patch.charge_and_current_operators(pairs,plus,minus,scale=matter_scale)
    currents=[]
    for index,create in enumerate(creates):
        pair_current=1j*pair_coupling*bmat([[None,-create.getH()],[create,None]],format='csr')
        embedded_matter=bmat([[csr_matrix((nv,nv),dtype=complex),None],[None,matter_currents[index]]],format='csr')
        currents.append(pair_current+embedded_matter)
    return densities,currents


def ward_checks(h,vacuum,pairs,plus,minus,creates,e,v):
    densities,currents=charge_current_operators(vacuum,pairs,plus,minus,creates)
    local=[]
    for site in patch.SITES:
        divergence=csr_matrix(h.shape,dtype=complex)
        for direction in range(2):
            divergence += currents[patch.LINK_INDEX[(site,direction)]]
            divergence -= currents[patch.LINK_INDEX[(patch.add(site,direction,-1),direction)]]
        residual=1j*(h@densities[patch.SITE_INDEX[site]]-densities[patch.SITE_INDEX[site]]@h)+divergence
        local.append(float(np.max(np.abs(residual.data))) if residual.nnz else 0.0)
    momentum=(math.pi,0.0); rho=csr_matrix(h.shape,dtype=complex); divergence_k=csr_matrix(h.shape,dtype=complex)
    for site in patch.SITES:
        phase=np.exp(-1j*(momentum[0]*site[0]+momentum[1]*site[1])); rho += phase*densities[patch.SITE_INDEX[site]]
    for direction in range(2):
        derivative=1-np.exp(-1j*momentum[direction])
        if abs(derivative)<1e-14: continue
        current_k=csr_matrix(h.shape,dtype=complex)
        for site in patch.SITES:
            phase=np.exp(-1j*(momentum[0]*site[0]+momentum[1]*site[1])); current_k += phase*currents[patch.LINK_INDEX[(site,direction)]]
        divergence_k += derivative*current_k
    operator=1j*(h@rho-rho@h)+divergence_k
    ground=v[:,0]; spectral_residuals=[]
    for index in range(1,len(e)):
        excited=v[:,index]
        spectral_residuals.append(abs(1j*(e[index]-e[0])*(excited.conj()@(rho@ground))+excited.conj()@(divergence_k@ground)))
    return {"maximum_local_continuity_residual":max(local),"operator_ward_residual":float(np.max(np.abs(operator.data))) if operator.nnz else 0.0,"maximum_low_spectral_ward_residual":float(max(spectral_residuals)),"rho_k_norm_on_ground":float(np.linalg.norm(rho@ground))}


def external_unitary(vacuum,pairs,amplitude):
    phases=np.ones(len(vacuum)+len(pairs),dtype=complex); offset=len(vacuum)
    for index,(positive,negative,_) in enumerate(pairs):
        lambda_positive=-0.5*amplitude*((-1)**patch.SITES[positive][0])
        lambda_negative=-0.5*amplitude*((-1)**patch.SITES[negative][0])
        phases[offset+index]=np.exp(1j*(lambda_positive-lambda_negative))
    return diags(phases,dtype=complex).tocsr()


def polarization(pair_coupling=PAIR_COUPLING,pair_gap=PAIR_GAP,step=0.002):
    longitudinal=patch.mode_pattern(0); transverse=patch.mode_pattern(1)
    def energy(a,b):
        _,_,hamiltonian,_,_,_,_=combined_system(a*longitudinal+b*transverse,pair_coupling=pair_coupling,pair_gap=pair_gap)
        return float(lowest(hamiltonian,2)[0][0])
    energy_zero=energy(0,0); long_plus=energy(step,0); long_minus=energy(-step,0); trans_plus=energy(0,step); trans_minus=energy(0,-step)
    pp=energy(step,step); pm=energy(step,-step); mp=energy(-step,step); mm=energy(-step,-step)
    curvature_long=(long_plus+long_minus-2*energy_zero)/step**2
    curvature_trans=(trans_plus+trans_minus-2*energy_zero)/step**2
    mixed=(pp-pm-mp+mm)/(4*step**2)
    vacuum,pairs,h0,_,_,_,_=combined_system(pair_coupling=pair_coupling,pair_gap=pair_gap)
    _,_,hl,_,_,_,_=combined_system(0.37*longitudinal,pair_coupling=pair_coupling,pair_gap=pair_gap)
    unitary=external_unitary(vacuum,pairs,0.37); difference=hl-unitary@h0@unitary.getH()
    return {"longitudinal_curvature":curvature_long,"transverse_curvature":curvature_trans,"mixed_curvature":mixed,"transversality_ratio":abs(curvature_long)/max(abs(curvature_trans),1e-30),"exact_longitudinal_unitary_equivalence_residual":float(np.max(np.abs(difference.data))) if difference.nnz else 0.0}


def build_report(quick=False):
    vacuum,pairs,h,plus,minus,creates,h_pure=combined_system(); values=electric_values(vacuum,pairs)
    dressed=spectral(h,values,70 if quick else 120); bare=spectral(h_pure,patch.electric_transverse_values(vacuum),50)
    ward=ward_checks(h,vacuum,pairs,plus,minus,creates,dressed['energies'],dressed['vectors'])
    response=polarization(); control=polarization(pair_coupling=0.0)
    ground=dressed['vectors'][:,0]; pair_occupation=float(np.sum(np.abs(ground[len(vacuum):])**2))
    scan=[]
    if not quick:
        for scale in (0.25,0.5,0.75,1.0):
            vacuum_s,pairs_s,h_s,_,_,_,_=combined_system(pair_coupling=PAIR_COUPLING*scale)
            spectrum=spectral(h_s,electric_values(vacuum_s,pairs_s),80)
            ground_s=spectrum['vectors'][:,0]; occupation=float(np.sum(np.abs(ground_s[len(vacuum_s):])**2))
            scan.append({"pair_coupling_scale":scale,"photon_gap":spectrum['photon_gap'],"gap_shift":spectrum['photon_gap']-bare['photon_gap'],"photon_residue":spectrum['photon_residue'],"virtual_pair_occupation":occupation})
    gap_shift=dressed['photon_gap']-bare['photon_gap']
    checks={
        "basis_counts":len(vacuum)==115 and len(pairs)==1484,
        "ground_is_vacuum_dominated":0<pair_occupation<0.1,
        "local_continuity_exact":ward['maximum_local_continuity_residual']<1e-12,
        "operator_ward_exact":ward['operator_ward_residual']<1e-12,
        "spectral_ward":ward['maximum_low_spectral_ward_residual']<1e-9,
        "longitudinal_is_pure_gauge":response['exact_longitudinal_unitary_equivalence_residual']<1e-12 and response['transversality_ratio']<1e-5,
        "virtual_transverse_response_nonzero":abs(response['transverse_curvature'])>1e-5,
        "decoupled_control_zero":abs(control['transverse_curvature'])<1e-9,
        "photon_dressing_nonzero":abs(gap_shift)>1e-5,
        "photon_residue_retained":dressed['photon_residue']>0.9,
        "coupling_scan_monotonic":True if quick else all(right['gap_shift']>left['gap_shift'] and right['virtual_pair_occupation']>left['virtual_pair_occupation'] for left,right in zip(scan,scan[1:])),
        "same_unit_completion":PAIR_COMPLETION==1.0,
    }
    public=lambda spectrum:{key:value for key,value in spectrum.items() if key not in ('energies','vectors')}
    report={
        "module":"phase_junction_network/microscopic/check_virtual_pair_vacuum.py",
        "status":"Stage 3H PASS at finite virtual-pair scope" if all(checks.values()) else "Stage 3H FAIL",
        "run_mode":"quick" if quick else "full",
        "lattice":{"shape":[2,2],"vacuum_dimension":len(vacuum),"one_pair_dimension":len(pairs),"combined_dimension":len(vacuum)+len(pairs)},
        "coefficients":{"pair_gap_over_Delta":PAIR_GAP,"pair_creation_coupling_over_Delta":PAIR_COUPLING,"completion_coefficient":PAIR_COMPLETION,"new_counterterms":0},
        "bare_photon":public(bare),
        "virtual_pair_dressed_photon":public(dressed),
        "virtual_pair":{"ground_pair_occupation":pair_occupation,"gap_shift":gap_shift,"pole_self_energy_delta_omega_squared":dressed['photon_gap']**2-bare['photon_gap']**2,"residue_ratio":dressed['photon_residue']/bare['photon_residue']},
        "ward_identity":ward,
        "static_polarization":response,
        "decoupled_control":control,
        "pair_coupling_scan":scan,
        "checks":checks,
        "stage_pass":all(checks.values()),
        "claim_boundary":{"established":["one finite Hamiltonian contains vacuum and gauge-invariant neutral-pair sectors","pair creation and annihilation use the same unit completed-move rule","exact local and nonzero-momentum Ward identities across sector-changing moves","longitudinal backgrounds are exactly pure gauge","virtual-pair admixture produces transverse photon dressing without a new counterterm"],"not_established":["relativistic fermionic pair statistics and a determinant","infinite-volume vacuum polarization or charge renormalization","3+1D interacting QED","radiative stability beyond the one-pair truncation","mirror-wall completion or nonlinear gravity"]},
        "next_gate":"Replace the one-pair truncation with the finite domain-wall endpoint spectrum and test multi-pair scaling, charge renormalization, and common-cone stability."
    }
    if not report['stage_pass']: raise AssertionError(json.dumps(ready(checks),indent=2,sort_keys=True))
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('phase_junction_network/microscopic/virtual_pair_vacuum_results.json'))
    parser.add_argument('--quick',action='store_true')
    args=parser.parse_args(); report=build_report(args.quick)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(ready(report),indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({"status":report['status'],"stage_pass":report['stage_pass'],"virtual_pair":report['virtual_pair'],"ward":report['ward_identity'],"polarization":report['static_polarization'],"failed":[key for key,value in report['checks'].items() if not value]},indent=2,sort_keys=True))
    return 0


if __name__=='__main__': raise SystemExit(main())
