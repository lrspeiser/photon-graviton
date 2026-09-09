"""Coupled spatial regions and independent weak-probe carrier/event timing."""
from pathlib import Path
import gc
import hashlib
import importlib.util
import json
import os
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'multiple-regions'
PROBE_PATH=HERE.parent/'weak-signal-timing/check.py'
spec=importlib.util.spec_from_file_location('independent_probe',PROBE_PATH)
probe_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe_module)


def background(cells,count,energy,cfg):
    lo,hi=cfg['domain'];dx=(hi-lo)/cells;x=lo+(np.arange(cells)+.5)*dx
    K=cfg['inertia'];v=cfg['field_speed'];sigma=cfg['smoothing_width'];g=np.zeros(cells)
    for center in cfg['centers'][:count]:
        z=(x-center)/cfg['half_width'];g+=np.where(abs(z)<1,np.cos(np.pi*z/2)**2,0)
    def optical(y,X):
        phi=y[:cells];u=y[cells:2*cells]
        W=g*np.exp(-.5*((x-X)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        return 1+dx*np.dot(W,phi),dx*np.dot(W*(x-X)/sigma**2,phi),dx*np.dot(W,u),W
    def rhs(t,y):
        phi=y[:cells];u=y[cells:2*cells];X,P=y[-3:-1]
        n,nx,nt,W=optical(y,X)
        acceleration=v*v*(np.roll(phi,-1)-2*phi+np.roll(phi,1))/dx**2+P*W/(K*n*n)
        force=P*nx/(n*n)
        profile_grid_force=force-K*dx*np.dot(acceleration,(np.roll(phi,-1)-np.roll(phi,1))/(2*dx))
        return np.r_[u,acceleration,1/n,force,profile_grid_force]
    initial=np.zeros(2*cells+3);initial[-3:-1]=[cfg['source_x'],energy]
    sol=solve_ivp(rhs,[0,cfg['end_time']],initial,method='DOP853',dense_output=True,
                  rtol=2e-9,atol=2e-12,max_step=.75*dx/v)
    assert sol.success
    totals=[];boundary=[];mom=[];minimum_n=[]
    for t in np.linspace(0,cfg['end_time'],251):
        y=sol.sol(t);phi=y[:cells];u=y[cells:2*cells];n=optical(y,y[-3])[0]
        gradient=(np.roll(phi,-1)-phi)/dx
        density=K/2*(u*u+v*v*gradient*gradient)
        totals.append(dx*density.sum()+y[-2]/n)
        boundary.append(dx*(density[:8].sum()+density[-8:].sum())/energy)
        momentum=y[-2]-K*dx*np.dot(u,(np.roll(phi,-1)-np.roll(phi,1))/(2*dx))
        mom.append(abs(momentum-energy-y[-1])/energy);minimum_n.append(n)
    drift=max(abs(np.asarray(totals)-energy))/energy
    assert drift<cfg['checks']['relative_energy_drift']
    assert max(mom)<cfg['checks']['relative_energy_drift']
    assert max(boundary)<1e-12
    assert min(minimum_n)>.5
    y=sol.y[:,-1];final_energy=y[-2]/optical(y,y[-3])[0]
    summary={'cells':cells,'region_count':count,'initial_driver_energy':energy,
             'final_driver_energy':float(final_energy),'driver_energy_lost':float(energy-final_energy),
             'driver_energy_stretch':float(energy/final_energy),'relative_energy_drift':float(drift),
             'relative_momentum_rate_error':float(max(mom)),'relative_boundary_energy':float(max(boundary)),
             'initial_field_energy':0.0,'final_driver_x':float(y[-3]),'minimum_driver_optical_factor':float(min(minimum_n))}
    def field(t,X):
        if t<0 or t>cfg['end_time']:raise ValueError('Background extrapolation forbidden')
        return optical(sol.sol(t),X)[:3]
    return summary,field


def main():
    cfg=json.loads((HERE/'protocol.json').read_text());runs=[]
    for energy in cfg['driver_energies']:
        for count in cfg['region_counts']:
            for cells in cfg['cells']:
                summary,field=background(cells,count,energy,cfg)
                summary['probes']=[]
                for launch in cfg['probe_launches']:
                    row=probe_module.probe(field,cfg['source_x'],cfg['detector_x']-cfg['source_x'],launch,1.,cfg['end_time'])
                    assert row['frequency_timing_relative_difference']<cfg['checks']['frequency_timing_relative_difference']
                    assert row['independent_arrival_time_difference']<cfg['checks']['independent_arrival_difference']
                    summary['probes'].append(row)
                summary['finite_event_intervals']=[
                    {'launch_interval':[a['launch'],b['launch']],
                     'arrival_interval':[a['arrival_time'],b['arrival_time']],
                     'duration_stretch':(b['arrival_time']-a['arrival_time'])/(b['launch']-a['launch'])}
                    for a,b in zip(summary['probes'][:-1],summary['probes'][1:])]
                runs.append(summary)
                print(f"E={energy} regions={count} N={cells}: driver S={summary['driver_energy_stretch']:.9g}; probe S={[round(r['frequency_stretch'],9) for r in summary['probes']]}",flush=True)
                del field;gc.collect()
    convergence=[];composition=[]
    for energy in cfg['driver_energies']:
        for count in cfg['region_counts']:
            a,b=[r for r in runs if r['initial_driver_energy']==energy and r['region_count']==count]
            loss_delta=abs(a['driver_energy_lost']/b['driver_energy_lost']-1)
            probe_delta=max(abs(x['frequency_stretch']-y['frequency_stretch']) for x,y in zip(a['probes'],b['probes']))
            assert loss_delta<cfg['checks']['relative_grid_change_in_driver_loss']
            assert probe_delta<cfg['checks']['absolute_grid_change_in_probe_stretch']
            convergence.append({'energy':energy,'regions':count,'relative_loss_change':loss_delta,'max_probe_stretch_change':probe_delta})
        one,three=[next(r for r in runs if r['initial_driver_energy']==energy and r['region_count']==count and r['cells']==max(cfg['cells'])) for count in [1,3]]
        S=one['driver_energy_stretch'];actual=three['driver_energy_stretch']
        composition.append({'energy':energy,'joint_three_region_stretch':actual,'frozen_single_region_factor_cubed':S**3,
                            'reciprocal_energy_weak_loss_prediction':1+3*(S-1),
                            'scope':'Diagnostic approximations; the joint run evolves the full field and packet without resets.'})
    result={'scope':cfg['scope'],'checks_pass':True,'runs':runs,'grid_convergence':convergence,'composition_comparisons':composition,
            'finite_probe_backreaction_included':False,'matter_clock_law_derived':False,'capture_and_gravity_derived':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'check.py',HERE/'protocol.json',PROBE_PATH]}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'multiple-regions-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'grid_convergence':convergence,'composition':composition},indent=2))


if __name__=='__main__':main()
