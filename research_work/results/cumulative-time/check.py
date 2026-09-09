"""Cumulative whole-signal stretch: ray, waveform and energy checks."""
from pathlib import Path
import json,hashlib,os
import numpy as np
from scipy.integrate import solve_ivp,quad

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'cumulative-time'


def main():
    protocol=json.loads((HERE/'protocol.json').read_text());tol=protocol['tolerances'];cases=protocol['ray_cases']
    rays=[]
    for length in cases['slab_lengths']:
        for alpha in cases['alpha']:
            expected=np.exp(alpha*length/2)
            for launch in cases['launch_times']:
                for frequency in cases['initial_frequencies']:
                    def rhs(x,y):
                        g=np.sin(np.pi*x/length)**2
                        gx=np.pi/length*np.sin(2*np.pi*x/length)
                        answer=[]
                        for offset in [0,3]:
                            t,k,companion=y[offset:offset+3]
                            n=1+alpha*g*t;nx=alpha*gx*t
                            # dx/dt=1/n, dk/dt=k*n_x/n^2; independent of an imposed energy-loss ODE.
                            answer.extend([n,k*nx/n,alpha*g*k/n])
                        return answer
                    delay=cases['second_flash_delay']
                    initial=[launch,frequency,0.,launch+delay,frequency,0.]
                    sol=solve_ivp(rhs,[0,length],initial,rtol=2e-12,atol=1e-13,max_step=length/80)
                    assert sol.success
                    t,k,ec,t2,k2,ec2=sol.y[:,-1]
                    # g=0 and n=1 at reception, hence omega=k at both endpoints.
                    stretch=(t2-t)/delay
                    frequency_error=max(abs(k/frequency*expected-1),abs(k2/frequency*expected-1))
                    duration_error=abs(stretch/expected-1)
                    energy_error=max(abs(k+ec-frequency),abs(k2+ec2-frequency))/frequency
                    assert max(frequency_error,duration_error)<tol['ray_relative']
                    assert energy_error<tol['ledger_absolute']
                    rays.append({'length':length,'alpha':alpha,'launch_time':launch,'initial_frequency':frequency,
                                 'predicted_stretch':float(expected),'arrival_stretch':float(stretch),
                                 'energy_remaining_fraction':float(k/frequency),'companion_gain_fraction':float(ec/frequency),
                                 'frequency_relative_error':float(frequency_error),'duration_relative_error':float(duration_error),
                                 'energy_relative_error':float(energy_error)})
    # Independent Gaussian pulse integrals: normalized photon rate and fixed per-photon energy 2.
    pulses=[]
    for stretch in [1.,1.01,1.1,2.,3.,6.]:
        rate=lambda t:np.exp(-t*t/2)/np.sqrt(2*np.pi)
        out_rate=lambda t:rate(t/stretch)/stretch
        count=quad(out_rate,-np.inf,np.inf,epsabs=1e-11)[0]
        energy=quad(lambda t:2/stretch*out_rate(t),-np.inf,np.inf,epsabs=1e-11)[0]
        variance=quad(lambda t:t*t*out_rate(t),-np.inf,np.inf,epsabs=1e-11)[0]/count
        assert abs(count-1)<tol['pulse_integral_relative']
        assert abs(energy/(2/stretch)-1)<tol['pulse_integral_relative']
        assert abs(np.sqrt(variance)/stretch-1)<tol['pulse_integral_relative']
        pulses.append({'stretch':stretch,'integrated_photon_count':count,'energy_ratio':energy/2,
                       'duration_ratio':float(np.sqrt(variance)),'peak_power_ratio':1/stretch**2})
    composition=[]
    for s1,s2,delay1,delay2 in [(1.2,1.5,3.,7.),(2.,3.,0.,5.),(1.,1.,3.,4.)]:
        t=np.array([-1.,0.,2.]);nested=s2*(s1*t+delay1)+delay2
        direct=s1*s2*t+s2*delay1+delay2
        assert np.max(abs(nested-direct))<1e-12
        composition.append({'stretch':s1*s2,'delay':s2*delay1+delay2})
    ledgers=[]
    for alpha in [.01,.1,.5]:
        for capture in [0.,alpha,2*alpha,1.]:
            def flow(x,y):
                gamma,companion,deposit=y
                return [-alpha*gamma,alpha*gamma-capture*companion,capture*companion]
            end=10.;sol=solve_ivp(flow,[0,end],[1.,0.,0.],rtol=1e-11,atol=1e-13)
            gamma,companion,deposit=sol.y[:,-1]
            ec=(alpha*end*np.exp(-alpha*end) if capture==alpha else alpha*(np.exp(-alpha*end)-np.exp(-capture*end))/(capture-alpha))
            assert max(abs(sol.y.sum(axis=0)-1))<tol['ledger_absolute']
            assert abs(companion-ec)<tol['ledger_absolute'] and min(gamma,companion,deposit)>-tol['ledger_absolute']
            ledgers.append({'alpha':alpha,'capture_rate':capture,'photon':float(gamma),'companion':float(companion),'deposit':float(deposit)})
    predictions=[]
    for z in [0.,.01,.1,1.,2.,5.]:
        s=1+z
        predictions.append({'z':z,'stretch':s,'photon_energy_ratio':1/s,'companion_energy_fraction':1-1/s,
                            'duration_ratio':s,'photon_arrival_rate_ratio':1/s,'bolometric_flux_ratio':1/s**2,'event_fluence_ratio':1/s})
    source=ROOT/'research_work/results/conversion-first/results.json'
    prior=json.loads(source.read_text());alpha=prior['fitted']['alpha_per_mly']
    illustration=[{'distance_billion_light_years':bly,'uniform_effective_stretch':float(np.exp(alpha*1000*bly))} for bly in [1,5,10]]
    result={'scope':protocol['scope'],'agreed_receiver':'Companion sector receives all energy lost by photon stretching.',
            'checks_pass':True,'ray_cases':rays,'pulse_checks':pulses,'composition_checks':composition,'capture_ledgers':ledgers,
            'predictions':predictions,'historical_rate_extrapolations':illustration,
            'maximum_ray_relative_error':max(max(r['frequency_relative_error'],r['duration_relative_error'],r['energy_relative_error']) for r in rays),
            'fundamental_time_law_derived':False,'momentum_and_capture_microphysics_derived':False,'observational_validation':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'check.py',source]}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'cumulative-time-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'rays':len(rays),'pulses':len(pulses),'ledgers':len(ledgers),'max_error':result['maximum_ray_relative_error'],'predictions':predictions,'distance_illustrations':illustration},indent=2))


if __name__=='__main__':main()
