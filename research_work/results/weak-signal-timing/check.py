"""Independent carrier and event timing in an evolving, spatial companion field."""
from pathlib import Path
import json,hashlib,os
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'weak-signal-timing'


def background(cfg,energy):
    m=cfg['modes'];mass0=cfg['inertia']*cfg['length'];mass=mass0/2
    k=2*np.pi*np.arange(1,m+1)/cfg['length'];w=np.exp(-.5*(cfg['smoothing_sigma']*k)**2)
    def fields(y,x):
        co=np.cos(k*x);si=np.sin(k*x)
        qc=y[2:2+m];qs=y[2+m:2+2*m];pc=y[2+2*m:2+3*m];ps=y[2+3*m:2+4*m]
        n=y[0]+np.sum(w*(qc*co+qs*si))
        nx=np.sum(w*k*(-qc*si+qs*co))
        nt=y[1]/mass0+np.sum(w*(pc*co+ps*si))/mass
        return n,nx,nt
    def rhs(t,y):
        n,nx,nt=fields(y,y[-2]);source=y[-1]/n**2
        co=np.cos(k*y[-2]);si=np.sin(k*y[-2]);spring=mass*cfg['wave_speed']**2*k*k
        return np.r_[y[1]/mass0,source,y[2+2*m:2+4*m]/mass,
                     -spring*y[2:2+m]+source*w*co,-spring*y[2+m:2+2*m]+source*w*si,
                     1/n,y[-1]*nx/n**2]
    initial=np.zeros(4*m+4);initial[0]=1;initial[1]=mass0*cfg['rolling_rate'];initial[-2]=cfg['packet_x'];initial[-1]=energy
    sol=solve_ivp(rhs,[0,cfg['duration']],initial,method='DOP853',dense_output=True,rtol=2e-11,atol=1e-13,max_step=.02)
    assert sol.success
    def field(t,x):
        if t < -1e-10 or t>cfg['duration']+1e-10:raise ValueError('Background extrapolation forbidden')
        return fields(sol.sol(t),x)
    return sol,field


def probe(field,source,distance,launch,frequency,end):
    n0=field(launch,source)[0]
    def hamilton(t,y):
        n,nx,nt=field(t,y[0]);return [1/n,y[1]*nx/n**2]
    def arrival(t,y):return y[0]-source-distance
    arrival.terminal=True;arrival.direction=1
    ray=solve_ivp(hamilton,[launch,end],[source,frequency*n0],events=arrival,method='DOP853',rtol=1e-10,atol=1e-12,max_step=.03)
    assert ray.success and len(ray.t_events[0])==1
    arrival_time=ray.t[-1];n1=field(arrival_time,source+distance)[0]
    frequency_stretch=frequency/(ray.y[1,-1]/n1)
    def variation(x,y):
        n,nx,nt=field(y[0],x);return [n,nt*y[1]]
    var=solve_ivp(variation,[source,source+distance],[launch,1.],method='DOP853',rtol=1e-10,atol=1e-12,max_step=.03)
    assert var.success
    return {'launch':launch,'path_length':distance,'frequency':frequency,'arrival_time':float(arrival_time),
            'frequency_stretch':float(frequency_stretch),'local_event_stretch':float(var.y[1,-1]),
            'independent_arrival_time_difference':float(abs(var.y[0,-1]-arrival_time)),
            'frequency_timing_relative_difference':float(abs(frequency_stretch/var.y[1,-1]-1)),
            'n_emitted':float(n0),'n_received':float(n1),
            'clock_scaling_controls':{str(p):float(var.y[1,-1]*(n0/n1)**p) for p in [0,1,2]}}


def main():
    protocol=json.loads((HERE/'protocol.json').read_text());cfg=protocol['background'];s=protocol['signals'];tols=protocol['checks']
    prior_path=ROOT/'research_work/results/companion-backreaction/companion-backreaction-results.json'
    prior=json.loads(prior_path.read_text());runs=[]
    for energy in cfg['packet_energies']:
        sol,field=background(cfg,energy)
        y3=sol.sol(3);ep3=y3[-1]/field(3,y3[-2])[0]
        if energy:
            old=next(r for r in prior['runs'] if r['modes']==cfg['modes'] and r['inertia']==cfg['inertia'] and r['initial_packet_energy']==energy)
            reference_error=abs(ep3-old['final_packet_energy']);assert reference_error<tols['maximum_reference_background_energy_difference']
        else:reference_error=0.
        rows=[];colors=[];events=[]
        for distance in s['path_lengths']:
            group=[]
            for launch in s['launch_times']:
                row=probe(field,s['source_x'],distance,launch,s['reference_frequency'],cfg['duration'])
                assert row['frequency_timing_relative_difference']<tols['maximum_frequency_timing_relative_difference']
                assert row['independent_arrival_time_difference']<tols['maximum_arrival_time_difference']
                rows.append(row);group.append(row)
            other=probe(field,s['source_x'],distance,0,s['color_control_frequency'],cfg['duration'])
            color_error=abs(other['frequency_stretch']/group[0]['frequency_stretch']-1)
            assert color_error<tols['maximum_color_relative_difference'];colors.append({'path_length':distance,'relative_difference':color_error})
            times=np.array([r['arrival_time'] for r in group]);launches=np.array(s['launch_times'])
            interval=np.diff(times)/np.diff(launches);local=np.array([r['local_event_stretch'] for r in group])
            events.append({'path_length':distance,'finite_interval_stretches':interval.tolist(),
                           'local_stretch_range':[float(local.min()),float(local.max())],
                           'local_stretch_fractional_range':float(local.max()/local.min()-1),
                           'affine_prediction_max_time_residual':float(np.max(abs(times-(times[0]+local[0]*(launches-launches[0])))))})
        if energy==0:
            assert max(abs(r['local_event_stretch']/np.exp(cfg['rolling_rate']*r['path_length'])-1) for r in rows)<1e-8
        result={'driving_packet_energy':energy,'background_reference_energy_difference':reference_error,'signals':rows,'color_checks':colors,'finite_events':events}
        runs.append(result);print(json.dumps({'energy':energy,'events':events}),flush=True)
    result={'scope':protocol['scope'],'checks_pass':True,'runs':runs,'finite_probe_backreaction_included':False,
            'atom_clock_law_derived':False,'astronomical_validation':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'check.py',prior_path]}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'weak-signal-timing-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('40 weak-signal timing cases and 8 color controls checked.')


if __name__=='__main__':main()
