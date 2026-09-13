"""Two labeled packets, one receiver, no resets; finite spectral channel only."""
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
step=.05
n=np.arange(96,105)
ein=2*n*step
eout=n*step
er=50+step*np.arange(201)
size=len(er)+2*int(n.max())
er_all=er[0]+step*np.arange(size)
a=np.exp(-((ein-10)/.22)**2/4)
a/=np.linalg.norm(a)
t1,t2,tref=0.,4.,1.
p1=a*np.exp(1j*ein*t1)
p2=a*np.exp(1j*ein*t2)
records=[]
for sigma in [.1,.5,1.]:
    r=np.exp(-((er-55)/sigma)**2/4)*np.exp(1j*er*tref)
    r/=np.linalg.norm(r)
    # First operation; second packet remains in its input band.
    after1=np.zeros((len(n),len(n),size),complex)
    for i,shift in enumerate(n):
        after1[i,:,shift:shift+len(er)]=p1[i]*p2[:,None]*r[None,:]
    total_mid=eout[:,None,None]+ein[None,:,None]+er_all[None,None,:]
    total_final=eout[:,None,None]+eout[None,:,None]+er_all[None,None,:]
    def second(state):
        out=np.zeros_like(state)
        for j,shift in enumerate(n):
            out[:,j,shift:]=state[:,j,:-shift]
        return out
    joint=second(after1)
    norm=float(np.vdot(joint,joint).real)
    flat=joint.reshape(len(n)**2,size)
    rho=flat@flat.conj().T
    shifts=(n[:,None]+n[None,:]).ravel()
    photons=(p1[:,None]*p2[None,:]).ravel()
    overlap=[]
    for lag in range(int(np.ptp(shifts))+1):
        overlap.append(1. if lag==0 else float(np.vdot(abs(r[:-lag]),abs(r[lag:])).real))
    diff=shifts[:,None]-shifts[None,:]
    gram=np.array(overlap)[abs(diff)]*np.exp(-1j*step*diff*tref)
    expected=photons[:,None]*photons.conj()[None,:]*gram
    trace_error=float(np.max(abs(rho-expected)))
    wait_errors=[]
    for wait in [.3,2.,11.]:
        sequential=second(after1*np.exp(-1j*total_mid*wait))
        common=joint*np.exp(-1j*total_final*wait)
        wait_errors.append(float(np.max(abs(sequential-common))))
    initial=2*float(abs(a)**2@ein)+float(abs(r)**2@er)
    final=float(np.sum(abs(joint)**2*total_final))
    receiver_gain=float(np.sum(abs(joint)**2*er_all[None,None,:])-abs(r)**2@er)
    # Partial trace over packet 2; compare with the first-use marginal.
    reduced=rho.reshape(len(n),len(n),len(n),len(n))
    marginal=np.einsum('ijkj->ik',reduced)
    firstflat=after1.reshape(len(n),-1)
    marginal_error=float(np.max(abs(marginal-firstflat@firstflat.conj().T)))
    assert abs(norm-1)<1e-12 and trace_error<1e-12
    assert max(wait_errors)<1e-11 and abs(initial-final)<1e-11 and marginal_error<1e-12
    records.append(dict(receiver_sigma=sigma,normalization=norm,
        joint_trace_formula_error=trace_error,free_evolution_errors=wait_errors,
        total_energy_error=abs(initial-final),receiver_energy_gain=receiver_gain,
        first_packet_marginal_change=marginal_error,
        continuum_gaussian_common_jitter_sd=1/(2*sigma),
        continuum_predicted_profile_centers=[2*t1-tref,2*t2-tref],
        continuum_predicted_profile_center_separation=2*(t2-t1)))
result=dict(scope='Finite spectral operations, no spatial interaction or causal arrival prediction',
            no_receiver_reset=True,waits=[.3,2.,11.],stretch_postulated=2,
            photons_labeled_by_distinct_modes=True,receiver_lightcone_energy=50,
            input_centers=[t1,t2],receiver_reference=tref,
            timing_values_are_analytic_continuum_predictions_not_detected_arrivals=True,
            rows=records)
(HERE/'sequence-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
