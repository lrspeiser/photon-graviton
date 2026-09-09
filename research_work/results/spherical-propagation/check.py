"""Energy-accounted spherical radiation shells and a 3D radial companion field."""
from pathlib import Path
import gc
import hashlib
import importlib.util
import json
import os
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import CubicHermiteSpline
from scipy.special import erfc

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'spherical-propagation'
PROBE=HERE.parent/'weak-signal-timing/check.py'
spec=importlib.util.spec_from_file_location('spherical_probe',PROBE)
probe_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe_module)


def background(dr,case,cfg):
    K=cfg['inertia'];v=cfg['field_speed'];sigma=cfg['smoothing_width'];mass=4*np.pi*K
    cells=int(round(case['outer_radius']/dr))-1;r=dr*np.arange(1,cells+1)
    g=np.zeros(cells)
    for center in case['centers']:
        z=(r-center)/cfg['half_width'];g+=np.where(abs(z)<1,np.cos(np.pi*z/2)**2,0)
    active=np.flatnonzero(g>0);ra=r[active];ga=g[active]
    shell_count=int(round(case['source_duration']/case['spacing']))
    assert shell_count*case['spacing']==case['source_duration']
    shell_energy=case['power']*case['spacing'];initial_fuel=shell_count*shell_energy
    source=cfg['source_radius'];count=0;fuel=initial_fuel
    y=np.zeros(2*cells+2*shell_count);y[2*cells:2*cells+shell_count]=source
    times=[];values=[];velocities=[];ledger=[];emissions=[];max_energy_error=0.;max_boundary=0.;minimum_n=1.
    def optical(state,positions):
        offset=ra[None,:]-np.atleast_1d(positions)[:,None]
        W=ga[None,:]*np.exp(-.5*(offset/sigma)**2)/(np.sqrt(2*np.pi)*sigma*ra[None,:])
        return 1+dr*(W@state[active]),dr*((W*offset/sigma**2)@state[active]),dr*(W@state[cells+active]),W
    def interacting(X):
        tail=cfg['kernel_tail_sigmas']*sigma
        return np.flatnonzero((X>ra.min()-tail)&(X<ra.max()+tail))
    def rhs(t,state):
        w=state[:cells];u=state[cells:2*cells];X=state[2*cells:2*cells+count];P=state[2*cells+shell_count:2*cells+shell_count+count]
        idx=interacting(X);n,nx,nt,W=optical(state,X[idx]);source_weights=P[idx]/n**2
        acc=v*v*(np.r_[w[1:],0]-2*w+np.r_[0,w[:-1]])/dr**2
        acc[active]+=(source_weights@W)/mass
        speed=np.zeros(shell_count);force=np.zeros(shell_count);speed[:count]=1.;speed[idx]=1/n;force[idx]=source_weights*nx
        return np.r_[u,acc,speed,force]
    def diagnose(t,state):
        nonlocal max_energy_error,max_boundary,minimum_n
        w=state[:cells];u=state[cells:2*cells];X=state[2*cells:2*cells+count];P=state[2*cells+shell_count:2*cells+shell_count+count]
        idx=interacting(X);n,nx,nt,W=optical(state,X[idx]);E=P.copy();E[idx]/=n
        edges=np.diff(np.r_[0,w,0])/dr
        field_energy=mass*dr/2*(np.dot(u,u)+v*v*np.dot(edges,edges))
        max_energy_error=max(max_energy_error,abs(fuel+E.sum()+field_energy-initial_fuel)/initial_fuel)
        edge_energy=mass*dr/2*(np.dot(u[-8:],u[-8:])+v*v*np.dot(edges[-8:],edges[-8:]))
        max_boundary=max(max_boundary,edge_energy/initial_fuel)
        if len(n):minimum_n=min(minimum_n,float(n.min()))
        conversion=float(np.sum(P[idx]*nt/n**2))
        field_velocity=u[active]/ra
        return [t,fuel,float(E.sum()),float(field_energy),conversion,float(field_velocity.mean()),float(np.max(abs(field_velocity)))]
    endpoints=list(np.arange(shell_count)*case['spacing'])+[case['end_time']]
    for j in range(shell_count):
        t=endpoints[j]
        n0=optical(y,[source])[0][0]
        y[2*cells+shell_count+j]=shell_energy*n0
        count=j+1;fuel=initial_fuel-count*shell_energy
        emissions.append({'time':float(t),'photon_energy':shell_energy,'canonical_momentum':float(shell_energy*n0),'remaining_fuel':float(fuel)})
        sol=solve_ivp(rhs,[t,endpoints[j+1]],y,method='DOP853',rtol=2e-9,atol=2e-12,max_step=.75*dr/v)
        assert sol.success
        # Field value and velocity are continuous through each emission. Photon
        # and fuel energies jump equally and oppositely; check both sides.
        for k,(time,state) in enumerate(zip(sol.t,sol.y.T)):
            row=diagnose(float(time),state)
            if times and time==times[-1]:
                ledger[-1]=row
            else:
                times.append(float(time));values.append(state[active].copy());velocities.append(state[cells+active].copy());ledger.append(row)
        y=sol.y[:,-1].copy()
    assert max_energy_error<cfg['checks']['relative_total_energy_drift'],max_energy_error
    assert max_boundary<cfg['checks']['relative_boundary_energy'],max_boundary
    assert minimum_n>.5
    history=CubicHermiteSpline(times,np.array(values),np.array(velocities),extrapolate=False)
    def field(t,X):
        if t<0 or t>case['end_time']:raise ValueError('Background extrapolation forbidden')
        offset=ra-X;W=ga*np.exp(-.5*(offset/sigma)**2)/(np.sqrt(2*np.pi)*sigma*ra)
        return float(1+dr*np.dot(W,history(t))),float(dr*np.dot(W*offset/sigma**2,history(t))),float(dr*np.dot(W,history(t,1)))
    # Continuum static solution: source power remains Q when n_t=0.
    # The Gaussian integral from the emitter outward is retained explicitly.
    source_profile=g*.5*erfc((source-r)/(np.sqrt(2)*sigma))
    static=np.array([case['power']/(mass*v*v)*dr*np.sum(source_profile/np.maximum(radius,r)) for radius in ra])
    low,high=case['field_window'];sample_t=np.linspace(low,high,101)
    sampled_phi=history(sample_t)/ra;sampled_u=history(sample_t,1)/ra
    static_relative_error=float(np.max(abs(sampled_phi-static))/np.max(abs(static)))
    def static_optical(X):
        W=ga*np.exp(-.5*((ra-X)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        return 1+dr*np.dot(W,static)
    static_delay=quad(lambda X:static_optical(X)-1,source,case['detector_radius'],epsabs=1e-11,points=case['centers'])[0]
    rows=np.array(ledger);select=(rows[:,0]>=low)&(rows[:,0]<=high)
    late=rows[select];energy_growth=float((late[-1,3]-late[0,3])/(late[-1,0]-late[0,0]))
    summary={'case':case['name'],'power':case['power'],'packet_spacing':case['spacing'],'grid_spacing':dr,'cells':cells,
             'shell_count':shell_count,'initial_fuel_energy':initial_fuel,'final_fuel_energy':fuel,
             'final_field_energy':float(rows[-1,3]),'maximum_relative_total_energy_drift':max_energy_error,
             'maximum_relative_boundary_energy':max_boundary,'minimum_photon_optical_factor':minimum_n,
             'static_prediction':{'support_radii':ra.tolist(),'field_on_support':static.tolist(),'excess_travel_time':float(static_delay)},
             'late_field_diagnostics':{'window':[low,high],'maximum_relative_static_profile_difference':static_relative_error,
                 'maximum_absolute_field_velocity':float(np.max(abs(sampled_u))),
                 'mean_field_velocity':float(sampled_u.mean()),'mean_total_field_energy_growth':energy_growth},
             'ledger_columns':['time','remaining_emitter_fuel','photon_energy','companion_field_energy','photon_to_field_power','mean_phi_t_on_support','maximum_absolute_phi_t_on_support'],
             'sampled_energy_ledger':rows[::max(1,len(rows)//200)].tolist()+[rows[-1].tolist()],
             'emissions':emissions}
    return summary,field


def main():
    cfg=json.loads((HERE/'protocol.json').read_text());runs=[];convergence=[]
    for case in cfg['cases']:
        for dr in cfg['grid_spacings']:
            row,field=background(dr,case,cfg)
            probes=[probe_module.probe(field,cfg['source_radius'],case['detector_radius']-cfg['source_radius'],t,1.,case['end_time']) for t in case['launches']]
            assert max(p['frequency_timing_relative_difference'] for p in probes)<cfg['checks']['frequency_timing_relative_difference']
            assert max(p['independent_arrival_time_difference'] for p in probes)<cfg['checks']['independent_arrival_difference']
            late=[p for p in probes if p['launch']>=case['late_launch_minimum']]
            row.update({'probes':probes,'late_probe_summary':{
                'maximum_absolute_redshift':float(max(abs(p['frequency_stretch']-1) for p in late)),
                'mean_redshift':float(np.mean([p['frequency_stretch']-1 for p in late])),
                'mean_excess_travel_time':float(np.mean([p['arrival_time']-p['launch']-p['path_length'] for p in late]))},
                'maximum_sampled_absolute_redshift':float(max(abs(p['frequency_stretch']-1) for p in probes)),
                'finite_event_intervals':[{'launch_interval':[a['launch'],b['launch']],
                    'duration_stretch':(b['arrival_time']-a['arrival_time'])/(b['launch']-a['launch'])} for a,b in zip(probes[:-1],probes[1:])]})
            runs.append(row);del field;gc.collect()
            print(f'{case["name"]} dr={dr}: probes={row["late_probe_summary"]}, field={row["late_field_diagnostics"]}',flush=True)
        a,b=runs[-2:]
        check={'case':case['name'],'maximum_absolute_probe_stretch_change':max(abs(x['frequency_stretch']-y['frequency_stretch']) for x,y in zip(a['probes'],b['probes'])),
               'maximum_absolute_late_probe_stretch_change':max(abs(x['frequency_stretch']-y['frequency_stretch']) for x,y in zip(a['probes'],b['probes']) if x['launch']>=case['late_launch_minimum']),
               'relative_final_field_energy_change':abs(a['final_field_energy']/b['final_field_energy']-1),
               'absolute_static_delay_change':abs(a['static_prediction']['excess_travel_time']-b['static_prediction']['excess_travel_time'])}
        convergence.append(check)
        assert check['maximum_absolute_late_probe_stretch_change']<cfg['checks']['absolute_grid_change_in_late_probe_stretch'],check
        assert check['maximum_absolute_probe_stretch_change']<cfg['checks']['absolute_grid_change_in_probe_stretch'],check
        assert check['relative_final_field_energy_change']<cfg['checks']['relative_grid_change_in_final_field_energy'],check
        assert check['absolute_static_delay_change']<cfg['checks']['absolute_grid_change_in_static_delay'],check
    result={'scope':cfg['scope'],'checks_pass':True,'runs':runs,'grid_convergence':convergence,
            'full_angular_illumination_tested':False,'cosmic_stationarity_or_capture_or_gravity_proved':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'check.py',HERE/'protocol.json',PROBE]}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'spherical-propagation-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
