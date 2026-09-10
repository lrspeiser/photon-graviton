"""Explicit impulse delivered to the prescribed conversion profile."""
from pathlib import Path
import importlib.util
import json
import hashlib
import numpy as np
from scipy.integrate import solve_ivp, simpson

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=HERE.parent/'bounded-wave-reservoir/run.py'
spec=importlib.util.spec_from_file_location('bounded_profile',SOURCE)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def simulate(a,b,N,tight=False):
    labels=np.linspace(-2,-1,N)
    F0=2*np.sin(np.pi*(labels+2))**2*np.cos(16*np.pi*(labels+2))
    E0=float(simpson(F0**2/2,x=labels))
    initial=np.r_[np.array([labels,np.ones(N),F0,np.zeros(N)]).ravel(),0.]
    def rhs(t,y):
        x,J,F,P=y[:-1].reshape(4,N);W,Wx=old.profile(x)
        v=1/(1+(b+a*t)*W);vx=-(b+a*t)*Wx*v*v;vs=-a*W*v*v
        states=np.array([v,vx*J,-vx*F,-vx*P-vs*F**2/2])
        # Opposite of the canonical wave momentum rate -integral v_x*e dx.
        impulse_rate=simpson(vx*(F**2/2)*J,x=labels)
        return np.r_[states.ravel(),impulse_rate]
    t=np.linspace(0,24,961)
    sol=solve_ivp(rhs,[0,24],initial,method='DOP853',rtol=2e-12 if tight else 2e-10,
                  atol=2e-13 if tight else 2e-11,max_step=.025 if tight else .05,t_eval=t)
    assert sol.success
    state=sol.y[:-1].T.reshape(len(t),4,N);x,J,F,P=state.transpose(1,0,2)
    wave_momentum=simpson(J*F**2/2,x=labels,axis=1)
    impulse=sol.y[-1]
    balance=float(np.max(abs(wave_momentum+impulse-E0)))
    assert balance<1e-8
    W,Wx=old.profile(x);v=1/(1+(b+a*t[:,None])*W);vx=-(b+a*t[:,None])*Wx*v*v
    sampled_force=simpson(vx*J*F**2/2,x=labels,axis=1)
    impulse_quad=float(simpson(sampled_force,x=t))
    gain=float(simpson(J[-1]*P[-1],x=labels))
    result=dict(a=a,b=b,rays=N,tight=tight,initial_energy_and_momentum=E0,
        final_wave_momentum=float(wave_momentum[-1]),
        net_profile_impulse=float(impulse[-1]),net_profile_impulse_fraction=float(impulse[-1]/E0),
        minimum_profile_impulse=float(impulse.min()),maximum_profile_impulse=float(impulse.max()),
        final_reservoir_energy=gain,canonical_momentum_balance_error=balance,
        independent_force_quadrature_error=abs(impulse_quad-impulse[-1]),
        reservoir_gain_minus_profile_impulse=gain-float(impulse[-1]))
    return result

def main():
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),SOURCE]}
    runs=[simulate(0,.4,129),simulate(.1,0,129),simulate(.1,0,257,True)]
    delta=runs[-1]['net_profile_impulse_fraction']
    assert abs(delta-(1-np.exp(-.15)))<1e-8
    assert abs(runs[-1]['reservoir_gain_minus_profile_impulse'])<1e-8
    assert abs(runs[0]['net_profile_impulse'])<1e-8
    recoils=[]
    for mu in [1.,10.,100.,10000.]:
        # Known finite-support recoil energy for rest-energy mu*E_initial.
        recoil=delta*delta/(np.hypot(mu,delta)+mu)
        final_energy=mu+recoil
        assert abs(final_energy**2-delta**2-mu**2)<max(1e-12,mu**2*1e-14)
        recoils.append(dict(support_rest_energy_over_initial_wave_energy=mu,
            recoil_energy_over_initial_wave_energy=float(recoil),
            recoil_fraction_of_claimed_transfer=float(recoil/delta),
            maximum_reservoir_fraction_if_wave_loss_unchanged=float(delta-recoil)))
    out=dict(scope='Canonical profile impulse and conditional finite-support recoil, not a completed moving-medium model.',
        runs=runs,recoil_comparisons=recoils,source_hashes=hashes,
        moving_profile_backreaction_solved=False,ordinary_companion_momentum_derived=False,
        holdouts_opened=False,total_cosmic_energy_budget_evaluated=False)
    for p,h in hashes.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    (HERE/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
