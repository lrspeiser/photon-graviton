"""Relativistic storage/recoil kinematics; no microscopic or coherent-wave claim."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent

def boost(v,b):
    g=1/np.sqrt(1-b*b)
    return g*np.array([v[0]-b*v[1],v[1]-b*v[0]])

def cycle(M,E,q,tau):
    mu=np.sqrt(M*M+2*M*E);beta=E/(M+E);gamma=(M+E)/mu
    T=gamma*tau;shift=beta*T;delay=(1-beta)*T
    initial=np.array([M+E,E]);excited=np.array([M+E,E])
    output=np.array([M,0.])+q*np.array([E,E])+(1-q)*np.array([E,E])
    assert np.max(abs(initial-excited))<1e-12 and np.max(abs(output-initial))<1e-10
    assert abs((M+E)**2-E**2-mu**2)<1e-7
    errors=[]
    for t in [-1.,0.,T/2,T,T+1.,T+7.]:
        if t<0:moment=E*t
        elif t<T:moment=(M+E)*beta*t
        else:moment=M*shift+E*(shift+t-T)
        errors.append(abs(moment/(M+E)-beta*t))
    assert max(errors)<1e-10 and abs(M*shift-E*delay)<1e-9
    boosts=[]
    for b in [-.8,-.3,0,.3,.8]:
        # Transform each physical constituent, not just the aggregate formula.
        before=boost(np.array([M,0.]),b)+boost(np.array([E,E]),b)
        during=boost(excited,b)
        after=boost(np.array([M,0.]),b)+boost(q*np.array([E,E]),b)+boost((1-q)*np.array([E,E]),b)
        err=max(np.max(abs(before-during)),np.max(abs(after-during)))
        assert err<1e-8
        boosts.append(dict(beta=b,conservation_residual=float(err)))
    return dict(rest_energy_M=M,input_photon_energy=E,photon_fraction=q,proper_storage_time=tau,
        excited_rest_energy=mu,stored_internal_energy=mu-M,kinetic_energy=M+E-mu,
        lab_storage_time=T,receiver_speed=beta,receiver_displacement=shift,arrival_delay=delay,
        final_photon_energy=q*E,final_companion_energy=(1-q)*E,final_receiver_rest_energy=M,
        final_receiver_momentum=0.,center_energy_error=max(errors),delay_displacement_residual=abs(M*shift-E*delay),boosts=boosts)

def sequence(M,kappa,q):
    # Receiver's own proper age sets an assumed residence schedule. No reset.
    E=1.;mu=np.sqrt(M*M+2*M*E);gamma=(M+E)/mu;beta=E/(M+E);a=M/mu
    position=10.;last_release=0.;proper_age=0.;rows=[]
    for emitted in [0.,12.,24.,36.]:
        absorption=emitted+position
        assert absorption>=last_release
        age_at_absorption=proper_age+absorption-last_release
        residence=.1+kappa*age_at_absorption
        duration=gamma*residence;release=absorption+duration
        position+=beta*duration
        arrival=release+100-position
        predicted=emitted+100+a*residence
        assert abs(arrival-predicted)<1e-10 and position<100
        rows.append(dict(emitted=emitted,absorbed=absorption,receiver_age_at_absorption=age_at_absorption,
            proper_residence=residence,released=release,receiver_position=position,arrival=arrival,
            assigned_spectral_factor=1/q))
        proper_age=age_at_absorption+residence;last_release=release
    intervals=[]
    for left,right in zip(rows,rows[1:]):
        dt=right['emitted']-left['emitted']
        measured=(right['arrival']-left['arrival'])/dt
        analytic=1+a*kappa*(1+(1-a)*left['proper_residence']/dt)
        assert abs(measured-analytic)<1e-10
        if not kappa:assert abs(measured-1)<1e-10
        intervals.append(dict(event_stretch=measured,analytic_event_stretch=analytic,
            assigned_spectral_factor=1/q,ratio_difference=measured-1/q))
    return dict(rest_energy_M=M,input_energy=E,kappa=kappa,photon_fraction=q,
        residence_postulate='proper_wait = 0.1 + kappa * receiver_proper_age_at_absorption',
        rows=rows,intervals=intervals)

def main():
    # Discrete controls of the analytic positive-energy null-momentum bound.
    controls=[]
    velocities=np.array([-1.,-.5,0.,.5,.9,1.])
    for maximum in [0.,.5,.9]:
        objective=-(velocities<=maximum).astype(float)
        solution=linprog(objective,A_eq=np.array([np.ones(6),velocities]),b_eq=[1.,1.],bounds=(0,None),method='highs')
        assert solution.success and abs(solution.fun)<1e-10
        controls.append(dict(maximum_forward_velocity=maximum,maximum_energy_in_those_states=-float(solution.fun),weights=solution.x.tolist()))
    # With no exactly forward null state, conserved incoming null momentum is infeasible.
    forbidden=linprog(np.zeros(5),A_eq=np.array([np.ones(5),velocities[:-1]]),b_eq=[1.,1.],bounds=(0,None),method='highs')
    assert forbidden.status==2
    cycles=[cycle(M,1.,q,tau) for M in [.1,1.,10.,10000.] for q in [.5,.8,1.] for tau in [.1,1.,10.]]
    seq=[sequence(M,k,q) for M in [1.,10.,10000.] for k in [0.,.2] for q in [.8,1.]]
    prior=HERE.parent/'null-stream-transfer/results.json'
    result=dict(scope='Conditional positive-energy kinematics and massive storage counterexample; no new fundamental time law.',
        units='c=1; M denotes receiver rest energy, not mass in kilograms',null_input_controls=controls,
        without_forward_null_state_infeasible=True,cycles=cycles,sequences=seq,
        quantum_coherence_derived=False,graviton_transition_derived=False,void_receiver_population_derived=False,
        new_holdouts_opened=False,astronomical_fit=False,
        input_hashes={str(prior):hashlib.sha256(prior.read_bytes()).hexdigest()},
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print('Cycles:',len(cycles),'boost checks:',sum(len(r['boosts']) for r in cycles))
    for r in seq:
        if r['photon_fraction']==.8:print('M',r['rest_energy_M'],'kappa',r['kappa'],'event factors',[round(x['event_stretch'],9) for x in r['intervals']],'assigned spectrum',1/.8)

if __name__=='__main__':main()
