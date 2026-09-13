"""Construct positive histories with identical current radiation but different storage."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm

OUT=Path(__file__).resolve().parent
beta=.1

def response(age):
    p=np.exp(-age)
    comp=(np.exp(-age)-np.exp(-beta*age))/(beta-1)
    return np.array([p,comp,1-p-comp])

def mean_burst(age,width=.02):
    # Finite-duration constant bursts, not idealized delta functions.
    # Integrate a unit-total source over this age interval with an augmented matrix.
    M=np.array([[-1,0,0,1],[1,-beta,0,0],[0,beta,0,0],[0,0,0,0]],float)
    initial=expm(M*width)@np.array([0.,0.,0.,1/width])
    free=M[:3,:3]
    return expm(free*(age-width/2))@initial[:3]

rows=[]
for old_age in [10.,30.,100.]:
    ages=[.1,1.,old_age]
    kernel=np.column_stack([mean_burst(t) for t in ages])
    assert np.all(kernel>=0) and np.max(abs(kernel.sum(axis=0)-1))<1e-12
    baseline=np.array([1.,1.,0.])
    # Add an old burst while exactly compensating both present radiation reservoirs.
    direction=np.r_[-np.linalg.solve(kernel[:2,:2],kernel[:2,2]),1.]
    limits=[baseline[i]/(-direction[i]) for i in [0,1] if direction[i]<0]
    added=.9*min(limits)
    altered=baseline+added*direction
    assert np.all(altered>=0)
    first=kernel@baseline;second=kernel@altered
    assert np.max(abs(first[:2]-second[:2]))<1e-12
    assert abs(sum(second)-sum(altered))<1e-10*sum(altered)
    # Independent finite-interval quadrature of analytic impulse response.
    from scipy.integrate import quad
    numeric=np.column_stack([[quad(lambda x:response(x)[j],t-.01,t+.01,epsabs=1e-12)[0]/.02 for j in range(3)] for t in ages])
    assert np.max(abs(numeric-kernel))<1e-12
    rows.append(dict(old_burst_age_alpha_c_T=old_age,burst_ages_alpha_c_T=ages,burst_width_alpha_c_T=.02,
        baseline_emitted_energies=baseline.tolist(),alternative_emitted_energies=altered.tolist(),
        baseline_current_photon_companion_deposit=first.tolist(),alternative_current_photon_companion_deposit=second.tolist(),
        deposited_energy_ratio=float(second[2]/first[2]),baseline_total_emitted=float(sum(baseline)),alternative_total_emitted=float(sum(altered)),
        present_radiation_max_absolute_difference=float(max(abs(first[:2]-second[:2])))))
result=dict(scope='Constructive identifiability test in the homogeneous constant-rate branch; histories are synthetic positive finite bursts, not inferred stars. Equal present bolometric reservoirs do not imply equal spectra.',
    beta_over_alpha=beta,energy_units='arbitrary common physical energy per volume; all emitted energy is explicitly counted',cases=rows)
(OUT/'history-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
