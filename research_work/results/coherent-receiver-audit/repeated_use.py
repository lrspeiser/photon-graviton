"""Repeat the stipulated exchange without resetting the receiving state.

Exact translated-state mixture in the populated interior of the finite ladder.
This checks energy buildup and inter-photon correlations, not spatial timing.
"""
import numpy as np
from scipy.stats import binom
from run import model,reduced,save

E,Er,target=model(maximum=300)
step=Er[1]-Er[0]
psi=np.array([0,1,1],complex)/np.sqrt(2)
results=[]
def entropy(rho):
    w=np.linalg.eigvalsh(rho);w=w[w>1e-14]
    return float(-np.sum(w*np.log2(w)))

for width in [.5,2,5]:
    env=np.exp(-(Er-30)**2/(4*width*width)).astype(complex)
    env[(Er<8)|(Er>65)]=0;env/=np.linalg.norm(env)
    initial=abs(env)**2;m0=float(initial@Er);v0=float(initial@((Er-m0)**2))
    # Execute two applications of the SAME unitary, retaining the whole state.
    first=reduced(psi,env,target)[1]
    joint=np.zeros((3,3,len(Er)),complex)
    for a in range(3):
        before=(psi[:,None]*first[a][None,:]).ravel()
        after=np.zeros_like(before);after[target]=before
        joint[a]=after.reshape(3,len(Er))
    reduced_two=joint.reshape(9,-1)@joint.reshape(9,-1).conj().T
    tensor=reduced_two.reshape(3,3,3,3)
    rho1=np.einsum('abcb->ac',tensor);rho2=np.einsum('abad->bd',tensor)
    product=np.kron(rho1,rho2)
    distance=.5*float(np.sum(abs(np.linalg.eigvalsh(reduced_two-product))))
    assert np.isclose(np.trace(reduced_two),1)
    # Compare the explicitly repeated unitary with translated-receiver overlaps.
    sums=[2,3,3,4]
    branches=np.array([np.roll(env,round(s/step)) for s in sums])/2
    predicted=branches@branches.conj().T
    active=[0,1,3,4]
    assert np.max(abs(reduced_two[np.ix_(active,active)]-predicted))<1e-12
    repeated=[]
    for n in [1,2,10,100]:
        # Each photon leaves receiver shifts +1 or +2 with equal probability.
        # After tracing photons: sum_k Binomial(n,k) T_(n+k) rho T_(n+k)^dagger.
        assert 65+2*n+4<Er[-1], 'Would hit boundary; this interior formula would fail'
        prob=np.zeros_like(initial)
        for k in range(n+1):prob+=binom.pmf(k,n,.5)*np.roll(initial,round((n+k)/step))
        mean=float(prob@Er);var=float(prob@((Er-mean)**2))
        assert abs(mean-(m0+1.5*n))<1e-10
        assert abs(var-(v0+.25*n))<1e-9
        repeated.append({'photons':n,'receiver_energy_gain':mean-m0,
                         'receiver_energy_variance_increase':var-v0,
                         'photon_energy_loss':1.5*n})
    results.append({'receiver_width':width,'two_photon_distance_from_product_state':distance,
        'two_photon_mutual_information_bits':entropy(rho1)+entropy(rho2)-entropy(reduced_two),
        'first_photon_beat_visibility':float(2*abs(rho1[0,1])),
        'second_photon_beat_visibility':float(2*abs(rho2[0,1])),
        'repeated_energy_ledger':repeated})
save('repeated-use-results.json',{'status':'No resets; exact interior translated-state model, with explicit two-step unitary verification. No free-flight or spatial dynamics specified.',
    'reservoir_max_energy':float(Er[-1]),'cases':results})
print(results)
