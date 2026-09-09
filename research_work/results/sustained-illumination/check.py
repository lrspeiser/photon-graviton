"""Longer closed photon trains compared with a derived rolling continuum branch."""
from pathlib import Path
import gc
import hashlib
import importlib.util
import json
import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from scipy.optimize import brentq
from scipy.special import erf

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'sustained-illumination'
PROBE=HERE.parent/'weak-signal-timing/check.py'
spec=importlib.util.spec_from_file_location('sustained_probe',PROBE)
probe_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe_module)


def continuum(power,cfg):
    K=cfg['inertia'];v=cfg['field_speed'];G=len(cfg['centers'])*cfg['half_width']
    # Divide the closure by u; the apparent zero root before division is spurious.
    def equation(u):
        return 4*K*v*u-power*(-np.expm1(-2*u*G))/u if u else -2*power*G
    u=brentq(equation,0,2*power*G/(2*K*v),xtol=1e-15)
    S=np.exp(u*G);out=power/S**2;field_rate=2*K*v*u*u
    assert abs(power-out-2*field_rate)<1e-13
    return {'integrated_coupling':G,'rolling_rate':u,'stretch':float(S),'redshift':float(S-1),
            'outgoing_photon_power':float(out),'photon_power_deficit':float(power-out),
            'field_energy_growth_rate':field_rate,'in_flight_photon_energy_growth_rate':field_rate,
            'closure_residual':float(equation(u))}


def background(cells,power,spacing,cfg):
    lo,hi=cfg['domain'];dx=(hi-lo)/cells;x=lo+(np.arange(cells)+.5)*dx
    K=cfg['inertia'];v=cfg['field_speed'];sigma=cfg['smoothing_width'];g=np.zeros(cells)
    for center in cfg['centers']:
        z=(x-center)/cfg['half_width'];g+=np.where(abs(z)<1,np.cos(np.pi*z/2)**2,0)
    active=np.flatnonzero(g>0);xa=x[active];ga=g[active]
    packets=int(round(cfg['source_duration']/spacing));assert packets*spacing==cfg['source_duration']
    packet_energy=power*spacing;initial_energy=packets*packet_energy
    def kernel(positions):
        offset=xa[None,:]-np.atleast_1d(positions)[:,None]
        W=ga[None,:]*np.exp(-.5*(offset/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        return W,offset
    def optical(y,positions):
        W,offset=kernel(positions);phi=y[active];u=y[cells+active]
        return 1+dx*(W@phi),dx*((W*offset/sigma**2)@phi),dx*(W@u),W
    def rhs(t,y):
        phi=y[:cells];u=y[cells:2*cells]
        X=y[2*cells:2*cells+packets];P=y[2*cells+packets:-1]
        n,nx,nt,W=optical(y,X);sources=P/n**2
        acc=v*v*(np.roll(phi,-1)-2*phi+np.roll(phi,1))/dx**2
        acc[active]+=(sources@W)/K
        forces=sources*nx
        external=forces.sum()-K*dx*np.dot(acc,(np.roll(phi,-1)-np.roll(phi,1))/(2*dx))
        return np.r_[u,acc,1/n,forces,external]
    y0=np.zeros(2*cells+2*packets+1)
    y0[2*cells:2*cells+packets]=cfg['source_x']-spacing*np.arange(packets)
    y0[2*cells+packets:-1]=packet_energy
    sol=solve_ivp(rhs,[0,cfg['end_time']],y0,method='DOP853',dense_output=False,
                  rtol=2e-9,atol=2e-12,max_step=.75*dx/v)
    assert sol.success
    # Retain only interaction-region histories, with nt the derivative of the
    # same Hermite interpolant used for n. This avoids storing full dense states.
    history=CubicHermiteSpline(sol.t,sol.y[active,:].T,sol.y[cells+active,:].T,extrapolate=False)
    time=sol.t;energy_errors=[];momentum_errors=[];boundary=[];records=[]
    a,b=cfg['control_volume'];width=cfg['control_edge_width'];rt2=np.sqrt(2)
    for t,y in zip(time,sol.y.T):
        phi=y[:cells];u=y[cells:2*cells];X=y[2*cells:2*cells+packets];P=y[2*cells+packets:-1]
        n,nx,nt,W=optical(y,X);E=P/n;gradient=(np.roll(phi,-1)-phi)/dx
        density=K/2*(u*u+v*v*gradient*gradient);field_energy=dx*density.sum()
        energy_errors.append(abs(field_energy+E.sum()-initial_energy)/initial_energy)
        momentum=P.sum()-K*dx*np.dot(u,(np.roll(phi,-1)-np.roll(phi,1))/(2*dx))
        momentum_errors.append(abs(momentum-initial_energy-y[-1])/initial_energy)
        boundary.append(dx*(density[:8].sum()+density[-8:].sum())/initial_energy)
        weight=.5*(erf((X-a)/(rt2*width))-erf((X-b)/(rt2*width)))
        left=np.exp(-.5*((X-a)/width)**2)/(np.sqrt(2*np.pi)*width)
        right=np.exp(-.5*((X-b)/width)**2)/(np.sqrt(2*np.pi)*width)
        conversion=E*nt/n
        records.append([t,field_energy,float(E.sum()),float(E@weight),float((conversion*weight).sum()),
                        float((E*left/n).sum()),float((E*right/n).sum()),float(conversion.sum()),
                        float(u[active].mean()),float(u[active].min()),float(u[active].max())])
    rows=np.array(records);low,high=cfg['energy_rate_window'];late=(time>=low)&(time<=high)
    r=rows[late];elapsed=high-low
    curve=CubicSpline(time,rows,axis=0,extrapolate=False)
    energy_curve=CubicHermiteSpline(time,rows[:,1],rows[:,7],extrapolate=False)
    flight_curve=CubicHermiteSpline(time,rows[:,3],rows[:,5]-rows[:,6]-rows[:,4],extrapolate=False)
    rates={name:float(curve.integrate(low,high)[i]/elapsed) for name,i in
           [('weighted_conversion',4),('incoming',5),('outgoing',6),('total_conversion',7)]}
    rates.update({'window':[low,high],
                  'field_energy_growth':float((energy_curve(high)-energy_curve(low))/elapsed),
                  'in_flight_photon_energy_growth':float((flight_curve(high)-flight_curve(low))/elapsed),
                  'mean_field_velocity_on_support':float(curve.integrate(low,high)[8]/elapsed),
                  'minimum_field_velocity_on_support':float(r[:,9].min()),
                  'maximum_field_velocity_on_support':float(r[:,10].max())})
    residual=rates['in_flight_photon_energy_growth']-(rates['incoming']-rates['outgoing']-rates['weighted_conversion'])
    rates['integrated_control_ledger_relative_error']=abs(residual)/max(abs(rates['total_conversion']),1e-30)
    assert rates['integrated_control_ledger_relative_error']<cfg['checks']['relative_integrated_control_ledger_error'],rates
    assert max(energy_errors)<cfg['checks']['relative_energy_drift']
    assert max(momentum_errors)<cfg['checks']['relative_energy_drift']
    assert max(boundary)<1e-12
    def field(t,X):
        if t<0 or t>cfg['end_time']:raise ValueError('Background extrapolation forbidden')
        W,offset=kernel([X]);phi=history(t);u=history(t,1)
        return float(1+dx*(W@phi)[0]),float(dx*((W*offset/sigma**2)@phi)[0]),float(dx*(W@u)[0])
    summary={'cells':cells,'power':power,'spacing':spacing,'packet_count':packets,
             'initial_total_energy':initial_energy,'final_field_energy':float(rows[-1,1]),
             'maximum_relative_energy_drift':float(max(energy_errors)),
             'maximum_relative_momentum_rate_error':float(max(momentum_errors)),
             'maximum_relative_boundary_energy':float(max(boundary)),
             'rate_diagnostics':rates,
             'history_columns':['time','field_energy','all_photon_energy','weighted_in_flight_photon_energy',
                                'weighted_conversion_rate','incoming_photon_power','outgoing_photon_power',
                                'total_conversion_rate','mean_field_velocity','minimum_field_velocity','maximum_field_velocity'],
             'sampled_energy_history':rows[::max(1,len(rows)//220)].tolist()+[rows[-1].tolist()]}
    return summary,field


def main():
    cfg=json.loads((HERE/'protocol.json').read_text());checks=cfg['checks'];runs=[];convergence=[]
    olddir=HERE.parent/'radiation-train';oldcfg=json.loads((olddir/'protocol.json').read_text())
    old=json.loads((olddir/'radiation-train-results.json').read_text())
    refcfg={**cfg,**oldcfg,'control_volume':cfg['control_volume'],'control_edge_width':cfg['control_edge_width'],
            'energy_rate_window':[18,24],'checks':{**oldcfg['checks'],**checks}}
    reference,field=background(2112,.003,.5,refcfg)
    saved=next(r for r in old['runs'] if r['cells']==2112 and r['power']==.003 and r['spacing']==.5)
    energy_difference=abs(reference['final_field_energy']/saved['final_field_energy']-1)
    probe_differences=[]
    for p in saved['probes'][::4]:
        new=probe_module.probe(field,refcfg['source_x'],refcfg['detector_x']-refcfg['source_x'],p['launch'],1.,refcfg['end_time'])
        probe_differences.append(abs(new['frequency_stretch']-p['frequency_stretch']))
    assert energy_difference<checks['relative_reference_energy_difference']
    assert max(probe_differences)<checks['absolute_reference_probe_difference']
    del field;gc.collect()
    print('Saved-background and interpolation reference passed.',flush=True)
    for case in cfg['cases']:
        theory=continuum(case['power'],cfg)
        for cells in cfg['cells']:
            row,field=background(cells,case['power'],case['spacing'],cfg)
            probes=[probe_module.probe(field,cfg['source_x'],cfg['detector_x']-cfg['source_x'],t,1.,cfg['end_time'])
                    for t in cfg['probe_launches']]
            assert max(p['frequency_timing_relative_difference'] for p in probes)<checks['frequency_timing_relative_difference']
            assert max(p['independent_arrival_time_difference'] for p in probes)<checks['independent_arrival_difference']
            late=[p for p in probes if p['launch']>=cfg['late_launch_minimum']]
            times=np.array([p['launch'] for p in late]);z=np.array([p['frequency_stretch']-1 for p in late])
            slope=np.polyfit(times-times.mean(),z,1)[0]
            row.update({'probes':probes,'continuum_prediction':theory,'late_statistics':{
                'mean_redshift':float(z.mean()),'minimum_redshift':float(z.min()),'maximum_redshift':float(z.max()),
                'linear_drift_per_time':float(slope),'relative_mean_redshift_difference_from_continuum':float(z.mean()/theory['redshift']-1)},
                'finite_event_intervals':[{'launch_interval':[a['launch'],b['launch']],
                   'duration_stretch':(b['arrival_time']-a['arrival_time'])/(b['launch']-a['launch'])} for a,b in zip(probes[:-1],probes[1:])]})
            rates=row['rate_diagnostics']
            row['continuum_rate_comparison']={key:float(rates[key]/theory['field_energy_growth_rate']-1)
                for key in ['field_energy_growth','in_flight_photon_energy_growth']}
            runs.append(row);del field;gc.collect()
            print(f'{case} cells={cells}: {row["late_statistics"]}; rates={row["continuum_rate_comparison"]}',flush=True)
        a,b=runs[-2:]
        change={'power':case['power'],'spacing':case['spacing'],
                'relative_photon_loss_change':abs(a['final_field_energy']/b['final_field_energy']-1),
                'relative_mean_redshift_change':abs(a['late_statistics']['mean_redshift']/b['late_statistics']['mean_redshift']-1),
                'maximum_relative_rate_change':max(abs(a['rate_diagnostics'][key]/b['rate_diagnostics'][key]-1) for key in ['field_energy_growth','in_flight_photon_energy_growth']),
                'maximum_absolute_probe_stretch_change':max(abs(x['frequency_stretch']-y['frequency_stretch']) for x,y in zip(a['probes'],b['probes']))}
        convergence.append(change)
        assert change['maximum_relative_rate_change']<checks['relative_grid_change_in_energy_rates'],change
        assert change['relative_photon_loss_change']<checks['relative_grid_change_in_total_photon_loss'],change
        assert change['relative_mean_redshift_change']<checks['relative_grid_change_in_mean_redshift'],change
        assert change['maximum_absolute_probe_stretch_change']<checks['absolute_grid_change_in_probe_stretch'],change
    result={'scope':cfg['scope'],'checks_pass':True,'runs':runs,'grid_convergence':convergence,
            'reference_check':{'relative_energy_difference':energy_difference,'maximum_absolute_probe_difference':max(probe_differences)},
            'continuum_stability_proved':False,'matter_clocks_capture_and_gravity_derived':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in
                             [HERE/'check.py',HERE/'protocol.json',PROBE,olddir/'protocol.json',olddir/'radiation-train-results.json']}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'sustained-illumination-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
