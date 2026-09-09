"""Closed finite packet trains and weak-probe timing in their shared field."""
from pathlib import Path
import gc
import hashlib
import importlib.util
import json
import os
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'radiation-train'
PROBE=HERE.parent/'weak-signal-timing/check.py'
spec=importlib.util.spec_from_file_location('train_probe',PROBE)
probe_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe_module)


def background(cells,power,spacing,cfg):
    lo,hi=cfg['domain'];dx=(hi-lo)/cells;x=lo+(np.arange(cells)+.5)*dx
    K=cfg['inertia'];v=cfg['field_speed'];sigma=cfg['smoothing_width'];g=np.zeros(cells)
    for center in cfg['centers']:
        z=(x-center)/cfg['half_width'];g+=np.where(abs(z)<1,np.cos(np.pi*z/2)**2,0)
    active=np.flatnonzero(g>0);xa=x[active];ga=g[active]
    packets=int(round(cfg['source_duration']/spacing));assert packets*spacing==cfg['source_duration']
    packet_energy=power*spacing;initial_energy=packets*packet_energy
    def optical(y,positions):
        offset=xa[None,:]-np.atleast_1d(positions)[:,None]
        W=ga[None,:]*np.exp(-.5*(offset/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        phi=y[active];u=y[cells+active]
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
    sol=solve_ivp(rhs,[0,cfg['end_time']],y0,method='DOP853',dense_output=True,
                  rtol=2e-9,atol=2e-12,max_step=.75*dx/v)
    assert sol.success
    energy_errors=[];momentum_errors=[];boundary=[];minimum_n=[];snapshots=[]
    for t in np.linspace(0,cfg['end_time'],281):
        y=sol.sol(t);phi=y[:cells];u=y[cells:2*cells];X=y[2*cells:2*cells+packets];P=y[2*cells+packets:-1]
        n=optical(y,X)[0];gradient=(np.roll(phi,-1)-phi)/dx
        density=K/2*(u*u+v*v*gradient*gradient);field_energy=dx*density.sum();photon_energy=np.sum(P/n)
        energy_errors.append(abs(field_energy+photon_energy-initial_energy)/initial_energy)
        momentum=P.sum()-K*dx*np.dot(u,(np.roll(phi,-1)-np.roll(phi,1))/(2*dx))
        momentum_errors.append(abs(momentum-initial_energy-y[-1])/initial_energy)
        boundary.append(dx*(density[:8].sum()+density[-8:].sum())/initial_energy);minimum_n.append(n.min())
        if any(abs(t-mark)<1e-8 for mark in [0,8,16,24,28]):
            snapshots.append({'time':float(t),'field_energy':float(field_energy),'photon_energy':float(photon_energy)})
    assert max(energy_errors)<cfg['checks']['relative_energy_drift']
    assert max(momentum_errors)<cfg['checks']['relative_energy_drift']
    assert max(boundary)<1e-12
    assert min(minimum_n)>.5
    summary={'cells':cells,'power':power,'spacing':spacing,'packet_count':packets,'initial_packet_energy':packet_energy,
             'initial_total_energy':initial_energy,'final_field_energy':snapshots[-1]['field_energy'],
             'maximum_relative_energy_drift':float(max(energy_errors)),'maximum_relative_momentum_rate_error':float(max(momentum_errors)),
             'maximum_relative_boundary_energy':float(max(boundary)),'minimum_driver_optical_factor':float(min(minimum_n)),
             'energy_snapshots':snapshots}
    def field(t,X):
        if t<0 or t>cfg['end_time']:raise ValueError('Background extrapolation forbidden')
        n,nx,nt,_=optical(sol.sol(t),[X]);return float(n[0]),float(nx[0]),float(nt[0])
    return summary,field


def main():
    cfg=json.loads((HERE/'protocol.json').read_text());runs=[];convergence=[]
    for case in cfg['cases']:
        for cells in cfg['cells']:
            row,field=background(cells,case['power'],case['spacing'],cfg);probes=[]
            for launch in cfg['probe_launches']:
                p=probe_module.probe(field,cfg['source_x'],cfg['detector_x']-cfg['source_x'],launch,1.,cfg['end_time'])
                assert p['frequency_timing_relative_difference']<cfg['checks']['frequency_timing_relative_difference']
                assert p['independent_arrival_time_difference']<cfg['checks']['independent_arrival_difference']
                probes.append(p)
            t=np.array(cfg['probe_launches']);z=np.array([p['frequency_stretch']-1 for p in probes]);center=t.mean()
            design=np.column_stack([np.ones(len(t)),t-center]);fit=np.linalg.lstsq(design,z,rcond=None)[0]
            residual=z-design@fit
            row['probes']=probes
            row['sample_statistics']={'mean_redshift':float(z.mean()),'minimum_redshift':float(z.min()),'maximum_redshift':float(z.max()),
                                      'linear_drift_per_time':float(fit[1]),'detrended_rms':float(np.sqrt(np.mean(residual**2))),
                                      'detrended_rms_over_mean_redshift':float(np.sqrt(np.mean(residual**2))/abs(z.mean())),
                                      'scope':'Unweighted nonuniform launch-time sample over a finite window; drift subtraction is descriptive, not removal of physical variation.'}
            row['finite_event_intervals']=[{'launch_interval':[a['launch'],b['launch']],
                'duration_stretch':(b['arrival_time']-a['arrival_time'])/(b['launch']-a['launch'])}
                for a,b in zip(probes[:-1],probes[1:])]
            runs.append(row);del field;gc.collect()
            print(f"power={case['power']} spacing={case['spacing']} N={cells}: {row['sample_statistics']}",flush=True)
        a,b=runs[-2:];field_change=abs(a['final_field_energy']/b['final_field_energy']-1)
        mean_change=abs(a['sample_statistics']['mean_redshift']/b['sample_statistics']['mean_redshift']-1)
        rms_change=abs(a['sample_statistics']['detrended_rms']/b['sample_statistics']['detrended_rms']-1)
        probe_change=max(abs(x['frequency_stretch']-y['frequency_stretch']) for x,y in zip(a['probes'],b['probes']))
        convergence.append({**case,'relative_photon_loss_change':field_change,'relative_mean_redshift_change':mean_change,'relative_detrended_rms_change':rms_change,'maximum_absolute_probe_stretch_change':probe_change})
        assert field_change<cfg['checks']['relative_grid_change_in_total_photon_loss'],convergence[-1]
        assert mean_change<cfg['checks']['relative_grid_change_in_mean_redshift'],convergence[-1]
        assert rms_change<cfg['checks']['relative_grid_change_in_detrended_rms'],convergence[-1]
        assert probe_change<cfg['checks']['absolute_grid_change_in_probe_stretch'],convergence[-1]
    result={'scope':cfg['scope'],'checks_pass':True,'runs':runs,'grid_convergence':convergence,
            'stationary_population_or_astronomical_validity_established':False,'capture_and_gravity_derived':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'check.py',HERE/'protocol.json',PROBE]}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'radiation-train-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(convergence,indent=2))


if __name__=='__main__':main()
