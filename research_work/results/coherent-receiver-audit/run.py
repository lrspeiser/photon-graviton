"""Energy-conserving frequency-remapping diagnostic, not a cosmological model.

All energies use an arbitrary unit; hbar=1. The reservoir ladder and interaction
are stipulated. No spatial propagation, momentum, capture or gravity is derived.
"""
from pathlib import Path
import json
import numpy as np

H=Path(__file__).resolve().parent
save=lambda name,obj:(H/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def model(step=.02,maximum=80):
    photon=np.array([1.,2.,4.]);energy=np.arange(round(maximum/step)+1)*step
    quanta=np.rint(photon/step).astype(int);N=len(energy)
    # A permutation inside each total-energy eigenspace is exactly unitary and
    # energy conserving. Its inverse supplies the reverse operation.
    target=np.arange(3*N)
    for total in range(N+int(quanta.max())):
        ids={i:i*N+(total-q) for i,q in enumerate(quanta) if 0<=total-q<N}
        if len(ids)==3:
            for old,new in [(0,2),(1,0),(2,1)]:target[ids[old]]=ids[new]
        elif set(ids)=={0,1}:
            target[ids[0]]=ids[1];target[ids[1]]=ids[0]
    total_energy_quanta=(quanta[:,None]+np.arange(N)[None,:]).ravel()
    assert np.array_equal(np.sort(target),np.arange(3*N))
    assert np.array_equal(total_energy_quanta[target],total_energy_quanta)
    return photon,energy,target

def reduced(psi,env,target):
    joint=(psi[:,None]*env[None,:]).ravel()
    out=np.zeros_like(joint);out[target]=joint
    out=out.reshape(3,len(env))
    return out@out.conj().T,out

def evolve_density(rho,E,t):
    p=np.exp(-1j*E*t);return p[:,None]*rho*p.conj()[None,:]

def diagnostics(psi,env,E,Er,target):
    rho,out=reduced(psi,env,target)
    before=float(np.sum(abs(psi)**2*E)+np.sum(abs(env)**2*Er))
    finalp=float(np.sum(rho.diagonal().real*E))
    finalr=float(np.sum(abs(out)**2*Er[None,:]))
    return rho,{'photon_energy_before':float(np.sum(abs(psi)**2*E)),
       'photon_energy_after':finalp,'reservoir_energy_gain':finalr-float(np.sum(abs(env)**2*Er)),
       'total_energy_residual':finalp+finalr-before,
       'beat_visibility':float(2*abs(rho[0,1])),
       'photon_purity':float(np.trace(rho@rho).real)}

if __name__=='__main__':
    E,Er,target=model()
    psi=np.array([0,1,1],complex)/np.sqrt(2)
    ground=np.zeros(len(Er),complex);ground[0]=1
    rho0,ground_result=diagnostics(psi,ground,E,Er,target)
    assert np.allclose(np.diag(rho0),[.5,.5,0]) and abs(rho0[0,1])<1e-14
    times=[.13,.7,1.3,2.]
    ground_covariance=[]
    for t in times:
        delayed=psi*np.exp(-1j*E*t)
        left=reduced(delayed,ground,target)[0]
        ground_covariance.append(float(np.linalg.norm(left-evolve_density(rho0,E,t))))
    # Ideal coherent dilation maps input energies 2,4 to 1,2 on this subspace.
    def ideal(p):
        q=np.array([p[1],p[2],0]);return np.outer(q,q.conj())
    ideal_defects=[float(np.linalg.norm(ideal(psi*np.exp(-1j*E*t))-evolve_density(ideal(psi),E,t))) for t in times]
    widths=[.1,.25,.5,1,2,5]
    rows=[];rho_plot={}
    for sigma in widths:
        env=np.exp(-(Er-30)**2/(4*sigma*sigma)).astype(complex)
        # Keep all populated states away from truncated-ladder boundaries.
        env[(Er<8)|(Er>65)]=0;env/=np.linalg.norm(env)
        rho,result=diagnostics(psi,env,E,Er,target)
        result.update({'reservoir_energy_sd_parameter':sigma,
                       'untruncated_gaussian_overlap_prediction':float(np.exp(-1/(8*sigma*sigma)))})
        defects=[];joint=[]
        for t in times:
            delayed=psi*np.exp(-1j*E*t)
            defects.append(float(np.linalg.norm(reduced(delayed,env,target)[0]-evolve_density(rho,E,t))))
            # Evolving BOTH input and reservoir restores total time symmetry.
            all_delayed=reduced(delayed,env*np.exp(-1j*Er*t),target)[0]
            joint.append(float(np.linalg.norm(all_delayed-evolve_density(rho,E,t))))
        result['fixed_reservoir_covariance_defects']=defects
        result['joint_evolution_covariance_max_defect']=max(joint)
        # Same energy probabilities, but remove all reservoir energy coherence.
        dephased=np.zeros((3,3),complex)
        for n in np.flatnonzero(abs(env)>0):
            # Only two initial photon components are populated. Trace over the
            # final reservoir explicitly, retaining equal-reservoir cross terms.
            sources=np.array([len(Er)+n,2*len(Er)+n]);dest=target[sources]
            fph=dest//len(Er);fenv=dest%len(Er)
            for i in range(2):
                for j in range(2):
                    if fenv[i]==fenv[j]:dephased[fph[i],fph[j]]+=abs(env[n])**2/2
        result['dephased_reservoir_beat_visibility']=float(2*abs(dephased[0,1]))
        assert abs(result['total_energy_residual'])<2e-12
        assert max(joint)<1e-12
        assert result['dephased_reservoir_beat_visibility']==0
        assert abs(result['beat_visibility']-result['untruncated_gaussian_overlap_prediction'])<2e-5
        rows.append(result);rho_plot[str(sigma)]=rho
    assert max(ground_covariance)<1e-12 and max(ideal_defects)>.1
    out={'status':'Synthetic energy/coherence channel tests; no astronomical fit or propagation law',
       'photon_levels':E.tolist(),'reservoir_step':float(Er[1]-Er[0]),'reservoir_max_energy':float(Er[-1]),
       'interaction_unitary_permutation':True,'interaction_preserves_each_total_energy_eigenspace':True,
       'sharp_ground_reservoir':ground_result,'stationary_channel_covariance_defects':ground_covariance,
       'ideal_dilation_covariance_defects':ideal_defects,'coherent_reservoir_cases':rows,
       'held_out_astronomical_scores_evaluated':False}
    save('results.json',out)
    # A beat diagnostic displays coherence, not a supernova light curve.
    t=np.linspace(0,4*np.pi,501)
    curves=[]
    for name,rho in [('sharp receiver',rho0),('coherent sigma=0.5',rho_plot['0.5']),('coherent sigma=5',rho_plot['5'])]:
        values=1+2*np.real(rho[0,1]*np.exp(-1j*(E[0]-E[1])*t))
        curves.append({'model':name,'time':t.tolist(),'normalized_beat_signal':values.tolist()})
    save('beat-curves.json',curves)
    print(json.dumps(out,indent=2))
