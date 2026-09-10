"""Canonical moving-profile completion; not a covariant graviton model."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp,simpson
from run import old,ROOT,HERE

def simulate(mu,tight=False):
    N=257;labels=np.linspace(-2,-1,N)
    F0=2*np.sin(np.pi*(labels+2))**2*np.cos(16*np.pi*(labels+2))
    E0=float(simpson(F0**2/2,x=labels));M=mu*E0
    initial=np.r_[np.array([labels,np.ones(N),F0,np.zeros(N)]).ravel(),0.,0.]
    def rhs(t,y):
        x,J,F,P=y[:-2].reshape(4,N);X,K=y[-2:]
        W,Wx=old.profile(x-X);v=1/(1+.1*t*W);vx=-.1*t*Wx*v*v;vs=-.1*W*v*v
        force=simpson(vx*F**2/2*J,x=labels)
        return np.r_[np.array([v,vx*J,-vx*F,-vx*P-vs*F**2/2]).ravel(),K/np.hypot(M,K),force]
    times=np.linspace(0,24,961)
    sol=solve_ivp(rhs,[0,24],initial,t_eval=times,method='DOP853',rtol=2e-12 if tight else 2e-10,
                  atol=2e-13 if tight else 2e-11,max_step=.025 if tight else .05)
    assert sol.success
    x,J,F,P=sol.y[:-2].T.reshape(len(times),4,N).transpose(1,0,2);X,K=sol.y[-2:]
    W,_=old.profile(x-X[:,None]);v=1/(1+.1*times[:,None]*W)
    momentum=simpson(F**2/2*J,x=labels,axis=1)
    light=simpson(v*F**2/2*J,x=labels,axis=1);reservoir=simpson(P*J,x=labels,axis=1)
    kinetic=K*K/(np.hypot(M,K)+M)
    energy_error=float(np.max(abs(light+reservoir+kinetic-E0))/E0)
    momentum_error=float(np.max(abs(momentum+K-E0))/E0)
    assert energy_error<1e-7 and momentum_error<1e-10 and J.min()>0
    assert np.all(x[-1]-X[-1]>4)
    return dict(support_rest_energy_ratio=mu,tight=tight,
        final_support_position=float(X[-1]),maximum_support_speed=float(np.max(abs(K/np.hypot(M,K)))),
        outgoing_wave_energy_fraction=float(light[-1]/E0),reservoir_energy_fraction=float(reservoir[-1]/E0),
        support_recoil_energy_fraction=float(kinetic[-1]/E0),support_momentum_fraction=float(K[-1]/E0),
        maximum_fractional_energy_error=energy_error,maximum_fractional_canonical_momentum_error=momentum_error)

def main():
    sources=[Path(__file__),HERE/'run.py',HERE.parent/'bounded-wave-reservoir/run.py']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    runs=[simulate(m) for m in [10.,100.,10000.]]+[simulate(10.,True)]
    for key in ['outgoing_wave_energy_fraction','reservoir_energy_fraction','support_recoil_energy_fraction']:
        assert abs(runs[0][key]-runs[-1][key])<1e-8
    for p,h in hashes.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    result=dict(scope='Moving-profile Hamiltonian energy and canonical momentum completion only.',
        hamiltonian_extension='W(x-X); support Hamiltonian sqrt(M^2+K^2); Xdot=K/sqrt(M^2+K^2); Kdot=integral v_x e dx.',
        runs=runs,source_hashes=hashes,profile_shape_origin_derived=False,
        symmetric_relativistic_stress_tensor_derived=False,holdouts_opened=False)
    (HERE/'moving-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(runs,indent=2))

if __name__=='__main__':main()
