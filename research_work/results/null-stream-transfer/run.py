"""Joined stress, transport and event-timing test of a kinetic null-stream alternative."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def window(x):
    x=np.asarray(x);v=np.zeros_like(x,dtype=float);inside=(x>1)&(x<3)
    v[inside]=np.sin(np.pi*(x[inside]-1)/2)**4
    return v

def evolve(labels,a,b,times):
    count=len(labels)
    def rhs(t,y):
        h=a*window(labels+t)*(1+b*t)
        transfer=h*y[:count]
        return np.r_[-transfer,transfer]
    sol=solve_ivp(rhs,[0,float(times[-1])],np.r_[np.ones(count),np.zeros(count)],t_eval=times,
        method='DOP853',rtol=2e-11,atol=2e-12,max_step=.03)
    assert sol.success
    return sol.y[:count].T,sol.y[count:].T

def pulse(count,a,b):
    labels=np.linspace(-2,-1,count);initial=2*np.sin(np.pi*(labels+2))**2
    times=np.linspace(0,10,201);g,c=evolve(labels,a,b,times)
    total=initial*(g+c)
    energy=np.trapezoid(total,labels,axis=1)
    centroid=np.trapezoid((labels[None,:]+times[:,None])*total,labels,axis=1)/energy
    exact=np.exp(-.75*a*(1+b*(2-labels)))
    result=dict(grid=count,rate=a,time_evolution=b,final_photon_energy=float(np.trapezoid(initial*g[-1],labels)),
        final_companion_energy=float(np.trapezoid(initial*c[-1],labels)),
        maximum_local_transfer_balance=float(np.max(abs(g+c-1))),
        maximum_global_energy_error=float(np.max(abs(energy-1))),
        maximum_centroid_error=float(np.max(abs(centroid-(-1.5+times)))),
        maximum_exit_fraction_error=float(np.max(abs(g[-1]-exact))),
        boosted=[])
    assert result['maximum_local_transfer_balance']<1e-10
    assert result['maximum_exit_fraction_error']<1e-8
    assert result['maximum_centroid_error']<1e-9
    k=np.array([1.,1.,0.,0.])
    for beta in [-.8,-.3,0,.3,.8]:
        gamma=1/np.sqrt(1-beta**2);doppler=gamma*(1-beta)
        L=np.eye(4);L[:2,:2]=[[gamma,-gamma*beta],[-gamma*beta,gamma]]
        # Check the tensor transformation and opposite source four-vectors.
        T=np.outer(k,k);transformed=L@T@L.T
        tensor_error=float(np.max(abs(transformed-doppler**2*T)))
        h=.37;Q=h*k;Qprime=L@Q
        rateprime=h/doppler
        source_error=float(np.max(abs(Qprime-rateprime*doppler**2*k)))
        # Equal boosted-time slices, not an incorrectly boosted equal-lab-time integral.
        energies=[];centroids=[]
        for tprime in [40.,41.]:
            labt=(tprime/gamma+beta*labels)/(1-beta)
            assert labt.min()>10
            labx=labels+labt;xprime=gamma*(labx-beta*labt)
            uprime=doppler**2*initial*(g[-1]+c[-1])
            E=float(np.trapezoid(uprime,xprime));X=float(np.trapezoid(xprime*uprime,xprime)/E)
            energies.append(E);centroids.append(X)
        speed=centroids[1]-centroids[0]
        assert abs(energies[0]-doppler)<1e-9 and abs(speed-1)<1e-9
        assert tensor_error<1e-12 and source_error<1e-12
        result['boosted'].append(dict(beta=beta,energy=energies[0],expected_energy=doppler,
            momentum=energies[0],centroid_speed=speed,tensor_residual=tensor_error,source_residual=source_error))
    return result

def ray(a,b,birth):
    def rhs(t,y):
        x,g,c=y;transfer=a*float(window(x))*(1+b*t)*g
        return [1.,-transfer,transfer]
    def arrived(t,y):return y[0]-6
    arrived.terminal=True;arrived.direction=1
    sol=solve_ivp(rhs,[birth,birth+7],[0,1,0],events=arrived,method='DOP853',rtol=2e-11,atol=2e-12,max_step=.03)
    assert sol.success and len(sol.t_events[0])==1
    arrival=float(sol.t_events[0][0]);x,g,c=sol.y_events[0][0]
    exact=np.exp(-.75*a*(1+b*(birth+2)))
    assert abs(arrival-birth-6)<1e-10 and abs(g-exact)<1e-8 and abs(g+c-1)<1e-10
    return dict(rate=a,time_evolution=b,birth=birth,arrival=arrival,photon_energy=float(g),
        companion_energy=float(c),conditional_spectral_factor=float(1/g))

def main():
    pulses=[pulse(n,a,b) for a,b in [(0.,0.),(.2,0.),(.2,.04)] for n in [129,257]]
    rays=[];intervals=[]
    for a in [0.,.02,.2]:
        for b in [0.,.04]:
            group=[ray(a,b,t) for t in [0.,.1,2.,5.]];rays+=group
            for left,right in zip(group,group[1:]):
                stretch=(right['arrival']-left['arrival'])/(right['birth']-left['birth'])
                assert abs(stretch-1)<1e-9
                intervals.append(dict(rate=a,time_evolution=b,first_birth=left['birth'],second_birth=right['birth'],
                    measured_event_stretch=stretch,first_conditional_spectral_factor=left['conditional_spectral_factor']))
    src=HERE.parent/'electromagnetic-audit/spectral-aging-predictions.csv'
    priorpath=HERE.parent/'companion-birth-timing/results.json'
    prior=json.loads(priorpath.read_text())
    assert digest(src)==prior['input_sha256']
    obs=pd.read_csv(src);assert len(obs)==35 and obs['object'].is_unique
    stationary=float(np.sum(((obs.observed_aging_rate-1)/obs.sigma)**2))
    matched=float(np.sum(((obs.observed_aging_rate-1/(1+obs.z))/obs.sigma)**2))
    assert abs(stationary-prior['spectral_aging']['fixed_birth_reset_chi2'])<1e-10
    assert abs(matched-prior['spectral_aging']['fixed_shared_phase_chi2'])<1e-10
    data=dict(scope='Alternative kinetic stress/transport closure; no microscopic photon interaction or graviton identification.',
        units='c=1, pulse initial energy=1, arbitrary length/time normalization',
        pulses=pulses,rays=rays,intervals=intervals,
        spectral_aging=dict(rows=35,no_stretch_chi2=stationary,matched_stretch_chi2=matched,
            data_status='Reused exposed data, exact reproduction of earlier scores, not new validation.',
            conditional_on_intrinsic_spectral_templates=True,full_covariance_included=False),
        assumptions=['straight fixed-distance light-speed propagation','fixed ordinary endpoint clocks',
            'phenomenological achromatic per-photon energy loss','photon count retained',
            'co-propagating null companion stress with no subsequent energy loss',
            'prescribed conversion rate; no derived action or capture law'],
        coherent_EM_phase_derived=False,graviton_identity_established=False,new_holdouts_opened=False,
        input_hashes={str(p):digest(p) for p in [src,priorpath]},code_sha256=digest(Path(__file__)))
    (HERE/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(pulses=[{k:r[k] for k in ['grid','rate','time_evolution','final_photon_energy','final_companion_energy','maximum_exit_fraction_error']} for r in pulses],
        aging=data['spectral_aging']),indent=2))

if __name__=='__main__':main()
