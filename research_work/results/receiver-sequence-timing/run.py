"""Sequential spectral-channel timing with an evolving, unreset receiver.

This remains an instantaneous mode map, not spacetime-local propagation.
"""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent;PRIOR=HERE.parent/'quantum-pulse-receiver'
spec=importlib.util.spec_from_file_location('sequence_prior',PRIOR/'run.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
STEP=old.DELTA;R=old.ER
births=[0.,30.,60.,90.]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def receiver():
    r=np.exp(-(R-100)**2/(4*4.**2));r[(R<60)|(R>140)]=0.;r/=np.linalg.norm(r)
    return r

def shifted(r,k):
    out=np.zeros(len(r),complex)
    if k>=0:out[k:]=r[:len(r)-k] if k else r
    else:out[:k]=r[-k:]
    assert abs(np.vdot(out,out).real-1)<1e-12
    return out

def gram(r,shifts):
    delta=shifts[:,None]-shifts[None,:]
    vals=np.array([np.vdot(r[k:],r[:-k]) if k else np.vdot(r,r) for k in range(int(np.ptp(shifts))+1)])
    # <T_m r | T_n r>, n-m positive, equals sum r[j-n] r[j-m]*.
    G=vals[abs(delta)];G=np.where(delta>=0,G,G.conj())
    return G

def main():
    source=json.loads((PRIOR/'results.json').read_text());assert source['code_sha256']==digest(PRIOR/'run.py')
    r=receiver();initial_variance=float((r*r)@((R-100)**2));rows=[];controls=[];envelopes=[]
    for a,b in [(2,1),(3,2)]:
        n=np.arange(200,301) if a==2 else np.arange(150,201)
        Ein=a*n*STEP;Eout=b*n*STEP;delta=(a-b)*n
        sigma=.14 if a==2 else .1
        amp=np.exp(-(Ein-Ein.mean())**2/(4*sigma*sigma));amp/=np.linalg.norm(amp)
        prob=amp*amp;S=a/b;L=2.5 if a==2 else 1.75;ell=round(L/STEP)
        increment=delta-ell;dist=np.ones(1);offset=0;G0=gram(r,delta)
        event_weights=np.array([.1,.2,.4,.3])
        input_mixture=np.zeros((len(n),len(n)),complex)
        output_mixture=np.zeros_like(input_mixture);frozen_mixture=np.zeros_like(input_mixture)
        for i,t in enumerate(births):
            # Prior photons are traced out. The exact receiver marginal is a
            # mixture of energy-translated states with their common free phase.
            # Free evolution of each translation differs only by a global phase.
            rt=r*np.exp(-1j*R*t);Gt=gram(rt,delta)
            expectedG=G0*np.exp(1j*(Ein-Eout)[:,None]*t-1j*(Ein-Eout)[None,:]*t)
            np.testing.assert_allclose(Gt,expectedG,atol=2e-12)
            # Translation invariance is checked at both support extremes and
            # the most probable carried receiver translation, without resetting.
            support=offset+np.arange(len(dist))
            for k in [int(support[0]),int(support[-1]),int(support[np.argmax(dist)])]:
                np.testing.assert_allclose(gram(shifted(r,k)*np.exp(-1j*R*t),delta),Gt,atol=2e-12)
            rho=np.outer(amp,amp)*Gt
            _,profile,local_center,width=old.profile(rho,b*STEP)
            absolute_center=t+local_center
            assert abs(absolute_center-S*t)<1e-6
            # Same state obtained at a common epoch and freely evolved to the
            # encounter: no extra time shift may be added twice.
            psi_global=amp*np.exp(1j*Ein*t)
            rho_global=np.outer(psi_global,psi_global.conj())*G0
            phase=np.exp(-1j*Eout*t)
            frame=phase[:,None]*rho_global*phase.conj()[None,:]
            frame_error=float(np.max(abs(frame-rho)));assert frame_error<2e-12
            input_mixture+=event_weights[i]*np.outer(psi_global,psi_global.conj())
            undo=np.exp(1j*Eout*t)
            output_mixture+=event_weights[i]*undo[:,None]*rho*undo.conj()[None,:]
            frozen_mixture+=event_weights[i]*undo[:,None]*(np.outer(amp,amp)*G0)*undo.conj()[None,:]
            _,_,frozen_local,_=old.profile(np.outer(amp,amp)*G0,b*STEP)
            mean_shift=float(dist@support)*STEP
            var_shift=float(dist@((support*STEP-mean_shift)**2))
            before=100+mean_shift
            dist=np.convolve(dist,prob);offset+=int(increment.min())
            after_support=offset+np.arange(len(dist));after=100+float(dist@after_support)*STEP
            np.testing.assert_allclose(dist.sum(),1,atol=1e-12)
            # All possible bank shifts remain within the positive-energy ladder.
            for k in [int(after_support[0]),int(after_support[-1])]:shifted(r,k)
            balance=abs(float(prob@Ein)+before-float(prob@Eout)-L-after)
            assert balance<1e-10
            rows.append(dict(stretch=S,sequence_index=i,encounter_time=t,local_profile_center=local_center,
                absolute_profile_center=absolute_center,output_profile_width=width,
                frozen_local_phase_control_absolute_center=t+frozen_local,frame_identity_error=frame_error,
                receiver_energy_before=before,receiver_energy_after=after,
                receiver_variance_before=initial_variance+var_shift,
                receiver_variance_after=initial_variance+float(dist@((after_support*STEP-(after-100))**2)),
                mean_photon_loss=float(prob@(Ein-Eout)),companion_energy_per_event=L,
                total_energy_and_forward_momentum_residual=balance))
        def moments_at(rho,energy,spacing,origin):
            p=np.exp(-1j*energy*origin)
            _,_,mean,width=old.profile(p[:,None]*rho*p.conj()[None,:],spacing)
            return mean+origin,width
        im,iw=moments_at(input_mixture,Ein,a*STEP,45.)
        om,ow=moments_at(output_mixture,Eout,b*STEP,S*45.)
        fm,fw=moments_at(frozen_mixture,Eout,b*STEP,45.)
        broadening=(S-1)**2/(4*4.**2)
        assert abs(om-S*im)<1e-6 and abs(fm-im)<1e-6
        assert abs(ow**2-S*S*iw**2-broadening)<1e-5
        envelopes.append(dict(stretch=S,incoherent_event_weights=event_weights.tolist(),input_mean=im,input_sd=iw,
            freely_evolving_receiver_output_mean=om,freely_evolving_receiver_output_sd=ow,
            frozen_local_phase_output_mean=fm,frozen_local_phase_output_sd=fw,
            extra_output_variance=ow**2-S*S*iw**2,gaussian_receiver_variance_prediction=broadening,
            scope='Finite Fourier-profile mixture, not a measured or causally propagated supernova light curve.'))
        # Exact reduced-state control with two photon frequencies. Enumerate
        # every distinct receiver translation generated by earlier encounters.
        shifts=np.array([200,300]) if a==2 else np.array([150,200])
        losses=shifts*STEP;energies=losses/(1-1/S);outenergies=energies/S
        pmap={0:1.}
        for t in births:
            reduced=np.zeros((2,2),complex);nextmap={}
            for k,w in pmap.items():
                rk=shifted(r,k)*np.exp(-1j*R*t)
                joint=np.zeros((2,len(R)+int(shifts.max())),complex)
                for j,s in enumerate(shifts):joint[j,s:s+len(R)]=rk/np.sqrt(2)
                reduced+=w*(joint@joint.conj().T)
                for s in shifts:
                    nk=k+int(s)-ell;nextmap[nk]=nextmap.get(nk,0)+w/2
            expected=.5*gram(r*np.exp(-1j*R*t),shifts)
            err=float(np.max(abs(reduced-expected)));assert err<2e-12
            controls.append(dict(stretch=S,encounter_time=t,receiver_mixture_terms=len(pmap),direct_trace_error=err))
            pmap=nextmap
    result=dict(scope='Sequential reduced-state spectral channel with free receiver evolution; not causal detector timing.',
        rows=rows,direct_two_frequency_controls=controls,incoherent_envelopes=envelopes,initial_receiver_energy=100.,initial_receiver_variance=initial_variance,
        receiver_reset_between_encounters=False,receiver_free_evolution='exp(-i E_R t)',
        fixed_release_is_exact_per_photon_loss=False,spacetime_local_interaction_derived=False,
        event_encounter_schedule_is_input=True,supernova_likelihood=False,holdouts_opened=False,
        source_hashes={str(p):digest(p) for p in [PRIOR/'run.py',PRIOR/'results.json']},code_sha256=digest(Path(__file__)))
    for S in [2.,1.5]:
        selected=[r for r in rows if r['stretch']==S]
        measured=np.diff([r['absolute_profile_center'] for r in selected])/np.diff(births)
        np.testing.assert_allclose(measured,S,atol=1e-7)
        print(S,'sequence profile interval factors',measured,'receiver final variance',selected[-1]['receiver_variance_after'])
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()
