"""Outgoing massless-mode release with explicit receiver energy bookkeeping."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
PRIOR=HERE.parent/'quantum-pulse-receiver'
spec=importlib.util.spec_from_file_location('release_prior',PRIOR/'run.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
STEP=old.DELTA

def setup(a,b,sigma,offset=0.):
    n=np.arange(200,301) if a==2 else np.arange(150,201)
    sigmaE=.14 if a==2 else .1
    Ein=a*n*STEP;Eout=b*n*STEP;shift=(a-b)*n
    amp=np.exp(-(Ein-Ein.mean()-offset)**2/(4*sigmaE**2));amp/=np.linalg.norm(amp)
    psi=amp*np.exp(1j*Ein*4.)
    r,gram=old.receiver_gram(shift,sigma,0.)
    K=len(r)+int(shift.max());R=old.ER[0]+STEP*np.arange(K)
    J=np.zeros((len(n),K),complex)
    for k,s in enumerate(shift):J[k,s:s+len(r)]=psi[k]*r
    rho=J@J.conj().T
    np.testing.assert_allclose(rho,np.outer(psi,psi.conj())*gram,atol=1e-12)
    return Ein,Eout,shift,amp,r,R,J,rho

def moments(prob,energy):
    mean=float(prob@energy)
    return mean,float(prob@((energy-mean)**2))

def release_case(a,b,sigma,offset=0.):
    Ein,Eout,shift,amp,r,R,J,rho=setup(a,b,sigma,offset)
    prob=amp**2;loss=Ein-Eout
    # Frozen release energy is calibrated to the original centered band only.
    L=2.5 if a==2 else 1.75;ell=round(L/STEP)
    no_emission=np.zeros_like(J);emission=np.zeros_like(J)
    no_emission[:,:ell]=J[:,:ell]
    emission[:,:len(R)-ell]=J[:,ell:]
    reduced=no_emission@no_emission.conj().T+emission@emission.conj().T
    rho_error=float(np.max(abs(reduced-rho)));assert rho_error<1e-11
    p0=np.sum(abs(no_emission)**2,axis=0);p1=np.sum(abs(emission)**2,axis=0)
    fail=float(p0.sum());released=L*float(p1.sum())
    R0,V0=moments(r*r,old.ER);R1,V1=moments(p0+p1,R)
    photon_before=float(prob@Ein);photon_after=float(prob@Eout);mean_loss=float(prob@loss)
    total_error=abs(photon_before+R0-photon_after-R1-released)
    total_p_error=abs(photon_before+R0-old.U-photon_after-(R1-old.U)-released)
    assert max(total_error,total_p_error)<1e-9
    assert np.min(old.U*(2*R-old.U))>0
    # Release pairs |R,vac> <-> |R-L,one L> preserve E and P exactly.
    basis_error=float(np.max(abs(R[ell:]-R[:-ell]-L)));assert basis_error<1e-10
    _,before_profile,_,_=old.profile(rho,b*STEP)
    _,after_profile,_,_=old.profile(reduced,b*STEP)
    assert np.max(abs(before_profile-after_profile))<1e-11
    success_by_input=np.sum(abs(emission)**2,axis=1)
    borrowed=float(np.sum(success_by_input[loss<L-1e-12]))
    borrowed_energy=float(np.sum(success_by_input*np.maximum(L-loss,0)))
    repaid_energy=float(np.sum(success_by_input*np.maximum(loss-L,0)))
    # Whole-receiver release to a fixed low-energy anchor moves its entire
    # spectral state to the outgoing massless mode. J is now photon x companion.
    anchor=float(R[0]);companion_energies=R-anchor
    complete_rho=J@J.conj().T
    assert np.max(abs(complete_rho-rho))<1e-12
    complete_companion=float(np.sum(abs(J)**2*companion_energies[None,:]))
    assert abs(photon_before+R0-photon_after-anchor-complete_companion)<1e-9
    assert abs(photon_before+R0-old.U-photon_after-(anchor-old.U)-complete_companion)<1e-9
    # Direct exact-loss tagging is a DIFFERENT joint conversion, not an allowed
    # replacement local operation on overlapping intermediate receiver states.
    exact_tag_rho=np.diag(prob).astype(complex)
    _,tag_profile,_,_=old.profile(exact_tag_rho,b*STEP)
    assert np.ptp(tag_profile)<1e-12
    rgram=old.receiver_gram(shift,sigma,0.)[1]
    local_reset_gram_defect=float(np.max(abs(rgram-np.eye(len(shift)))))
    exact_tag_energy_error=float(np.max(abs(Ein-Eout-loss)))
    return dict(stretch=a/b,receiver_width=sigma,input_band_center_offset=offset,
        photon_energy_before=photon_before,photon_energy_after=photon_after,photon_energy_loss=mean_loss,
        fixed_release=dict(companion_energy_per_emission=L,mean_companion_energy=released,
            no_emission_probability=fail,receiver_initial_energy=R0,receiver_final_energy=R1,
            receiver_energy_change=R1-R0,receiver_momentum_change=R1-R0,
            initial_receiver_variance=V0,final_receiver_variance=V1,
            conditional_loss_variance=float(prob@((loss-mean_loss)**2)),
            branches_borrowing_receiver_energy=borrowed,mean_borrowed_energy=borrowed_energy,
            mean_repaid_energy=repaid_energy,photon_density_matrix_change=rho_error,
            total_energy_error=total_error,total_momentum_error=total_p_error,basis_release_residual=basis_error),
        full_anchor_release=dict(anchor_energy=anchor,anchor_momentum=anchor-old.U,
            anchor_speed=(anchor-old.U)/anchor,mean_companion_energy=complete_companion,
            initial_receiver_energy_consumed=R0-anchor,photon_density_matrix_change=float(np.max(abs(complete_rho-rho)))),
        direct_exact_loss_tagging=dict(photon_profile_variation=float(np.ptp(tag_profile)),
            photon_purity=float(np.sum(prob*prob)),total_energy_error=exact_tag_energy_error,
            attempted_receiver_only_reset_gram_defect=local_reset_gram_defect))

def repeated(a,b,N=100):
    Ein,Eout,shift,amp,r,_,_,_=setup(a,b,4.)
    weights=amp**2;L=2.5 if a==2 else 1.75;ell=round(L/STEP);dmin=int(shift.min())
    assert np.all(np.diff(shift)==1)
    cut=ell-dmin;assert cut>0
    K=len(r)+N*max(int(shift.max())-ell,0)+ell+10
    energy=old.ER[0]+STEP*np.arange(K)
    p=np.zeros(K);p[:len(r)]=r*r
    R0,V0=moments(p,energy);loss=Ein-Eout;lossmean=float(weights@loss)
    variance=float(weights@((loss-lossmean)**2));cumulative=0.;rows=[]
    for step in range(1,N+1):
        captured=np.convolve(p,weights)
        # Captured index q denotes true receiver index q+dmin. Attempt constant
        # release; below the positive-energy floor retain the captured energy.
        assert captured[cut+K:].sum()<1e-20
        next_p=np.zeros(K)
        kept=captured[cut:cut+K];next_p[:len(kept)]=kept
        next_p[dmin:dmin+cut]+=captured[:cut]
        failure=float(captured[:cut].sum())
        cumulative+=L*float(captured[cut:].sum())
        p=next_p
        assert abs(p.sum()-1)<1e-10
        mean,var=moments(p,energy)
        balance=abs(step*lossmean+R0-mean-cumulative)
        assert balance<1e-8
        if step in [1,2,10,100]:
            assert abs(var-(V0+step*variance))<1e-7
            rows.append(dict(uses=step,receiver_energy=mean,receiver_momentum=mean-old.U,
                receiver_variance=var,interior_variance_prediction=V0+step*variance,
                cumulative_companion_energy=cumulative,cumulative_photon_loss=step*lossmean,
                last_step_no_emission_probability=failure,total_energy_balance_error=balance))
    return dict(stretch=a/b,receiver_width=4.,rows=rows,
        propagated='Receiver diagonal probabilities and emitted energy without resetting; not a full multiphoton joint state.')

def main():
    base=json.loads((PRIOR/'results.json').read_text())
    assert base['code_sha256']==hashlib.sha256((PRIOR/'run.py').read_bytes()).hexdigest()
    rows=[release_case(a,b,sigma) for a,b in [(2,1),(3,2)] for sigma in [.02,.1,1.,4.]]
    sensitivity=[release_case(a,b,4.,offset) for a,b in [(2,1),(3,2)] for offset in [-.2,.2]]
    histories=[repeated(a,b) for a,b in [(2,1),(3,2)]]
    result=dict(scope='Finite quantum release bookkeeping; no causal local field action or astronomical fit.',
        rows=rows,changed_input_spectrum=sensitivity,repeated_use=histories,
        exact_eventwise_loss_rule_preserved_by_fixed_release=False,
        massless_modes_have_P_equal_E=True,physical_spatial_propagation_derived=False,
        multiphoton_timing_simulated=False,graviton_identity_established=False,new_holdouts_opened=False,
        source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PRIOR/'results.json',PRIOR/'run.py']},
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    for row in rows:
        if row['receiver_width']==4:print(json.dumps(row,indent=2))
    print('Completed',len(rows),'release cases,',len(sensitivity),'spectrum changes and two 100-use histories.')

if __name__=='__main__':main()
