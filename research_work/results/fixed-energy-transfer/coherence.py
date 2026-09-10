"""Conditional quantum-channel completion of the energy ladder; not a local action."""
from pathlib import Path
import math,json
import numpy as np
from scipy.special import comb
HERE=Path(__file__).resolve().parent
N=96;m=np.arange(N+1)
psi=np.exp(-.25*((m-60)/3)**2).astype(complex);psi/=np.linalg.norm(psi)
rho=np.outer(psi,psi.conj())
def channel(r,q):
    out=np.zeros_like(r,dtype=complex)
    for k in range(N+1):
        n=np.arange(k,N+1)
        coeff=np.sqrt(comb(n,k)*q**(n-k)*(1-q)**k)
        out[:N+1-k,:N+1-k]+=r[k:,k:]*coeff[:,None]*coeff[None,:]
    return out
def phase(r,t):return r*np.exp(1j*(m[:,None]-m[None,:])*t)
def center(r):
    # Fourier time density uses exp[-i (m-n)t]; first circular moment below.
    detected=r.copy();detected[0,:]=0;detected[:,0]=0
    moment=np.trace(detected,offset=-1)
    return float(np.angle(moment)),float(abs(moment)/np.trace(detected).real)
rows=[]
for q in [.99,.8,1/1.885801]:
    out=channel(rho,q);trace=float(np.trace(out).real)
    assert abs(trace-1)<1e-12
    assert np.linalg.eigvalsh(out)[0]>-1e-12
    input_mean=float(m@np.diag(rho).real);output_mean=float(m@np.diag(out).real)
    assert abs(output_mean-q*input_mean)<1e-10
    for t in [.2,.5,1.2]:
        delayed=channel(phase(rho,t),q)
        err=float(np.max(abs(delayed-phase(out,t))))
        measured,concentration=center(delayed)
        assert err<1e-12 and abs(measured-t)<1e-10
        rows.append(dict(q=q,input_pulse_center=t,output_pulse_center=measured,
            required_dilated_center=t/q,mean_frequency_ratio=output_mean/input_mean,
            covariance_max_error=err,purity=float(np.trace(out@out).real),
            first_circular_moment_magnitude=concentration,vacuum_probability=float(out[0,0].real)))
# State-by-state energy and collinear momentum at the isometry level.
for initial in range(N+1):
    for lost in range(initial+1):
        assert initial == (initial-lost)+lost
result=dict(status='Known amplitude-damping channel mathematics reinterpreted as an energy ladder; no new local graviton action.',
    energy_units='epsilon',time_units='hbar/epsilon',frequency_levels=N+1,
    input_mean_level=float(m@np.diag(rho).real),cases=rows,
    conclusion='Mean frequency decreases but translated pulse centers retain their separation; stationary receiving bath does not generate universal temporal dilation.',
    limitations=['Phase-coordinate diagnostic, not an actual detector likelihood','Collinear momentum bookkeeping only; no angular momentum or local field stress',
        'Vacuum boundary retained','No clock coupling or astrophysical transition amplitude derived'])
(HERE/'coherence-results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
