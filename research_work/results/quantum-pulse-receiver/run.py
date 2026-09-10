"""Finite energy/momentum-conserving spectral swaps and Fourier pulse profiles.

No spacetime-local interaction or causal detector propagation is implemented.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
DELTA=.01
U=100. # Receiver light-cone energy E_R-P_R; positive and unchanged by a forward transfer.
ER=50.01+DELTA*np.arange(10001)

def profile(rho,spacing):
    period=2*np.pi/spacing
    t=np.linspace(-period/2,period/2,4096,endpoint=False)
    coeff=np.array([np.diagonal(rho,offset=k).sum() for k in range(len(rho))])
    phase=np.exp(1j*spacing*np.arange(1,len(rho))[:,None]*t)
    probability=(coeff[0].real+2*np.real(coeff[1:]@phase))/period
    dt=period/len(t)
    assert probability.min()>-1e-12 and abs(probability.sum()*dt-1)<1e-10
    mean=float(np.sum(t*probability)*dt)
    variance=float(np.sum((t-mean)**2*probability)*dt)
    # Check the fast diagonal-sum transform against a direct density-matrix sum.
    for index in [123,1101,2048,3077]:
        v=np.exp(1j*spacing*np.arange(len(rho))*t[index])
        direct=float(np.real(v.conj()@rho@v)/period)
        # This vector convention gives the same -i(E_n-E_m)t intensity.
        assert abs(direct-probability[index])<1e-10
    return t,probability,mean,np.sqrt(variance)

def receiver_gram(shifts,sigma,reference):
    r=np.exp(-(ER-100)**2/(4*sigma*sigma));r/=np.linalg.norm(r)
    maxlag=int(np.ptp(shifts));overlaps=np.empty(maxlag+1)
    overlaps[0]=1.
    for k in range(1,maxlag+1):overlaps[k]=np.dot(r[:-k],r[k:])
    difference=shifts[:,None]-shifts[None,:]
    gram=overlaps[abs(difference)]*np.exp(-1j*DELTA*difference*reference)
    return r,gram

def family(a,b):
    if (a,b)==(2,1):n=np.arange(200,301);sigmaE=.14
    else:n=np.arange(150,201);sigmaE=.1
    incoming=a*n*DELTA;outgoing=b*n*DELTA;shifts=(a-b)*n
    assert not set((a*n).tolist()) & set((b*n).tolist())
    amplitude=np.exp(-(incoming-incoming.mean())**2/(4*sigmaE**2));amplitude/=np.linalg.norm(amplitude)
    weights=amplitude**2
    receiver_size=len(ER)+int(shifts.max())
    er_all=ER[0]+DELTA*np.arange(receiver_size)
    assert np.min(U*(2*er_all-U))>0
    # Every input basis state maps to a unique disjoint output basis state.
    # Pairwise swaps plus identity elsewhere are an exactly unitary completion.
    source=(a*n)[:,None]*receiver_size+np.arange(len(ER))[None,:]
    target=(b*n)[:,None]*receiver_size+np.arange(len(ER))[None,:]+shifts[:,None]
    assert np.unique(target).size==target.size and np.intersect1d(source.ravel(),target.ravel()).size==0
    in_energy=incoming[:,None]+ER[None,:]
    out_energy=outgoing[:,None]+ER[None,:]+DELTA*shifts[:,None]
    in_momentum=incoming[:,None]+ER[None,:]-U
    out_momentum=outgoing[:,None]+ER[None,:]+DELTA*shifts[:,None]-U
    energy_error=float(np.max(abs(in_energy-out_energy)))
    momentum_error=float(np.max(abs(in_momentum-out_momentum)))
    assert max(energy_error,momentum_error)<1e-10
    rows=[];curves=[]
    for sigma in [.02,.1,1.,4.]:
        for birth in [0.,4.,8.]:
            for reference in [0.,3.]:
                psi=amplitude*np.exp(1j*incoming*birth)
                r,gram=receiver_gram(shifts,sigma,reference)
                rho=np.outer(psi,psi.conj())*gram
                assert abs(np.trace(rho)-1)<1e-12 and np.linalg.eigvalsh(rho).min()>-1e-10
                ti,pi,mi,si=profile(np.outer(psi,psi.conj()),a*DELTA)
                to,po,mo,so=profile(rho,b*DELTA)
                S=a/b;expected=S*birth-(S-1)*reference
                expected_sd=np.sqrt((S/(2*sigmaE))**2+((S-1)/(2*sigma))**2)
                assert abs(mi-birth)<1e-5 and abs(mo-expected)<1e-4
                assert abs(so-expected_sd)/expected_sd<.002
                ideal=amplitude*np.exp(1j*outgoing*expected)
                fidelity=float(np.real(ideal.conj()@rho@ideal))
                Eg0=float(weights@incoming);Eg1=float(weights@outgoing)
                R0=float((r*r)@ER);P0=R0-U
                R1=float(sum(weights[k]*np.dot(r*r,ER+DELTA*shifts[k]) for k in range(len(n))))
                P1=R1-U
                assert abs(Eg0+R0-Eg1-R1)<1e-10 and abs(Eg0+P0-Eg1-P1)<1e-10
                # Joint phase translation must translate output by the same amount.
                _,gram_both=receiver_gram(shifts,sigma,reference+2)
                psi_both=amplitude*np.exp(1j*incoming*(birth+2))
                both=np.outer(psi_both,psi_both.conj())*gram_both
                phase=np.exp(1j*outgoing*2)
                covariance=float(np.max(abs(both-phase[:,None]*rho*phase.conj()[None,:])))
                assert covariance<1e-11
                direct_error=None
                if birth==0 and reference==0 and sigma==1:
                    joint=np.zeros((len(n),receiver_size),complex)
                    for k,shift in enumerate(shifts):joint[k,shift:shift+len(r)]=psi[k]*r
                    direct=joint@joint.conj().T
                    direct_error=float(np.max(abs(direct-rho)));assert direct_error<1e-11
                    assert abs(np.sum(abs(joint)**2*er_all[None,:])-R1)<1e-10
                rows.append(dict(receiver_width=sigma,birth_phase_time=birth,receiver_reference_time=reference,
                    input_profile_center=mi,output_profile_center=mo,predicted_center=expected,
                    input_profile_sd=si,output_profile_sd=so,gaussian_predicted_sd=expected_sd,
                    ideal_pulse_fidelity=fidelity,photon_energy_before=Eg0,photon_energy_after=Eg1,
                    receiver_energy_before=R0,receiver_energy_after=R1,receiver_momentum_before=P0,receiver_momentum_after=P1,
                    joint_time_translation_residual=covariance,direct_joint_trace_residual=direct_error))
                if birth==4 and reference==0:
                    curves.append(dict(receiver_width=sigma,time=to[::8].tolist(),profile=po[::8].tolist()))
    # Dephasing the receiver, without changing its energy probabilities, removes
    # all off-diagonal output terms because different transfers leave distinct labels.
    td,pd,_,_=profile(np.diag(weights).astype(complex),b*DELTA)
    assert np.ptp(pd)<1e-12
    return dict(stretch=S,photon_levels=len(n),input_band=[float(incoming.min()),float(incoming.max())],
        output_band=[float(outgoing.min()),float(outgoing.max())],receiver_initial_rest_energy_mean=100.,
        minimum_receiver_mass_squared=float(np.min(U*(2*er_all-U))),
        basis_swaps=int(target.size),basis_energy_residual=energy_error,basis_momentum_residual=momentum_error,
        dephased_receiver_profile_variation=float(np.ptp(pd)),rows=rows,curves=curves)

def main():
    families=[family(2,1),family(3,2)]
    priors=[HERE.parent/'coherent-receiver-audit/results.json',HERE.parent/'recoil-delay-receiver/results.json']
    result=dict(scope='Finite spectral-unitary and Fourier pulse-profile test; not causal propagation or a void model.',
        units='hbar=c=1; arbitrary energies and reciprocal-energy times',families=families,
        receiver_dispersion='P_R=E_R-U, m_R^2=U*(2*E_R-U), U=100',
        receiver_contains_initial_coherence=True,energy_retained_in_receiver_not_escaping_companions=True,
        frequency_map_stipulated=True,spacetime_locality_derived=False,causal_arrival_map_derived=False,
        repeated_receiver_reuse_simulated=False,new_holdouts_opened=False,
        source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in priors},
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    for f in families:
        print('Stretch',f['stretch'],'basis swaps',f['basis_swaps'])
        for row in f['rows']:
            if row['birth_phase_time']==4 and row['receiver_reference_time']==0:
                print(row['receiver_width'],'center',row['output_profile_center'],'sd',row['output_profile_sd'],'fidelity',row['ideal_pulse_fidelity'])

if __name__=='__main__':main()
