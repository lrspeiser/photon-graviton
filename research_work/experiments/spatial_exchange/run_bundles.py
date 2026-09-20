"""SE-B evolving test-body bundles with same-stage background integration."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from emitter_model import ExchangeEvolution
from probes import coefficients,evaluate,derivative
from run_emitter import save

ROOT=Path(__file__).resolve().parent


def configurations():
    base=dict(n=32,length=16.,radius=.9,dt=.02,chi=200.,mix=.5,emitter_angle=np.pi/2,probe_radius=.9)
    variants=[dict(name='emission-only',chi=0.,mix=0.,emitter_angle=0.),
              dict(name='X',emitter_angle=0.),dict(name='mixed',emitter_angle=np.pi/4),dict(name='Y'),
              dict(name='time',dt=.01),dict(name='space',n=40,dt=.01),dict(name='probe-radius',probe_radius=1.2)]
    return [dict(base,**v) for v in variants]


def initial_probes():
    offsets=[(0.,0.)]
    for delta in (.05,.025):offsets += [(delta,0.),(-delta,0.),(0.,delta),(0.,-delta)]
    rays=np.array([[-1.5,1+y,z,1.,0.,0.] for y,z in offsets])
    theta=np.arange(6)*np.pi/3
    bodies=np.column_stack((2*np.cos(theta),2*np.sin(theta),np.zeros(6),-.02*np.sin(theta),.02*np.cos(theta),np.zeros(6)))
    return np.vstack((rays,bodies)),np.array([0.]*9+[1.]*6)


def bundle(events):
    matrices=[]
    for start,delta in [(1,.05),(5,.025)]:
        y=(np.array(events[start]['position'])[1:]-events[start+1]['position'][1:])/(2*delta)
        z=(np.array(events[start+2]['position'])[1:]-events[start+3]['position'][1:])/(2*delta)
        matrices.append(np.column_stack((y,z)))
    M=matrices[-1];det=float(np.linalg.det(M))
    return dict(matrices=matrices,bundle_error=float(np.max(abs(matrices[0]-M))),
                trace_distortion=float(1-np.trace(M)/2),shear_components=[float((M[1,1]-M[0,0])/2),float(-(M[0,1]+M[1,0])/2)],
                rotation=float((M[1,0]-M[0,1])/2),area_gain=1/abs(det) if det else None)


def run(cfg,out):
    wall=time.perf_counter();model=ExchangeEvolution(**{k:v for k,v in cfg.items() if k not in ['name','dt','probe_radius']})
    y=model.initial();initial=y.copy();probes,masses=initial_probes();initial_p=probes.copy()
    steps=round(4/cfg['dt']);trace=[];states=[];velocities=[];energies=[];cone=0.;extent=0.;events=[None]*9
    def stage(bg,pr):
        arrays=coefficients(model.unpack(bg)[0])
        return model.rhs(bg),derivative(pr,masses,cfg['probe_radius'],model.length,arrays)
    for step in range(steps+1):
        arrays=coefficients(model.unpack(y)[0]);vel=[];energy=[]
        for row,mass in zip(probes,masses):
            value=evaluate(row[:3],row[3:],mass,cfg['probe_radius'],model.length,*arrays)
            vel.append(value['velocity']);energy.append(value['energy'])
            cone=max(cone,float(np.linalg.norm(value['velocity']-value['bbar'])/value['cbar']-1))
        vel=np.array(vel);energy=np.array(energy);states.append(probes.copy());velocities.append(vel);energies.append(energy)
        if step:
            previous=states[-2]
            for i in range(9):
                if events[i] is None and previous[i,0]<1.5<=probes[i,0]:
                    fraction=(1.5-previous[i,0])/(probes[i,0]-previous[i,0])
                    point=previous[i]+fraction*(probes[i]-previous[i]);velocity=velocities[-2][i]+fraction*(vel[i]-velocities[-2][i])
                    events[i]=dict(time=(step-1+fraction)*cfg['dt'],position=point[:3],momentum=point[3:],velocity=velocity,
                                   coordinate_energy=float(energies[-2][i]+fraction*(energy[i]-energies[-2][i])))
        extent=max(extent,float(np.max(abs(model.unpack(y)[2]))+model.radius))
        if step%round(.1/cfg['dt'])==0 or step==steps:
            metric=model.metrics(y);metric['time']=step*cfg['dt'];trace.append(metric)
        if step==steps:break
        dt=cfg['dt'];a,ap=stage(y,probes);b,bp=stage(y+dt*a/2,probes+dt*ap/2)
        c,cp=stage(y+dt*b/2,probes+dt*bp/2);d,dp=stage(y+dt*c,probes+dt*cp)
        y=y+dt*(a+2*b+2*c+d)/6;probes=probes+dt*(ap+2*bp+2*cp+dp)/6
    drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger'])
    clearance=model.length/2-1-extent-4*model.maximum_characteristic
    background_pass=drift<1e-5 and clearance>0 and max(t['edge_amplitude'] for t in trace)<1e-5 and max(t['cone_error'] for t in trace)<1e-10
    crossed=all(e is not None for e in events);transport=bundle(events) if crossed else None
    references={'emission-only':('evidence-v1','emission-only'),'X':('evidence-v1','exchange-200'),
                'mixed':('emitter-v1','emitter-mixed-chi200'),'Y':('emitter-v1','emitter-Y-chi200'),
                'time':('emitter-v1','time-refinement'),'space':('emitter-v1','space-refinement'),
                'probe-radius':('emitter-v1','emitter-Y-chi200')}
    folder,name=references[cfg['name']];reference=ROOT/folder/(name+'.npz')
    replay=float(np.max(abs(y-np.load(reference)['final']))) if reference.exists() else None
    row=dict(config=cfg,background_passed=background_pass,energy_drift=drift,clearance=clearance,probe_cone_excess=cone,
             all_crossed=crossed,transport=transport,background_replay_error=replay,events=events,
             central_bend=float(np.arctan2(events[0]['velocity'][1],events[0]['velocity'][0])) if crossed else None,
             arrival_offset=events[0]['time']-3 if crossed else None,
             final_body_positions=probes[9:,:3],final_body_velocities=vel[9:],
             body_radial_displacement=np.linalg.norm(probes[9:,:3],axis=1)-2,
             wall_seconds=time.perf_counter()-wall,
             passed=bool(background_pass and cone<1e-10 and crossed and transport['bundle_error']<.001 and (replay is None or replay<1e-12)))
    np.savez_compressed(out/(cfg['name']+'.npz'),initial_background=initial,final_background=y,initial_probes=initial_p,
                        time=np.arange(steps+1)*cfg['dt'],probe_states=states,probe_velocities=velocities,probe_energies=energies)
    save(out/(cfg['name']+'.json'),dict(summary=row,background_trace=trace))
    print('done',cfg['name'],'passed',row['passed'],'bend',row['central_bend'],'delay',row['arrival_offset'],flush=True)
    return row


def main():
    out=ROOT/'bundle-v1';out.mkdir(exist_ok=False)
    save(out/'manifest.json',dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
         sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['run_bundles.py','probes.py','emitter_model.py','run_emitter.py','bundle-protocol.md']},configurations=configurations()))
    rows=[]
    for cfg in configurations():
        print('start',cfg['name'],flush=True);rows.append(run(cfg,out))
    base=rows[3];comparisons=[]
    for row in rows[4:6]:
        fine=row;time_case=row['config']['name']=='time'
        bend=abs(base['central_bend']-fine['central_bend'])/max(abs(fine['central_bend']),1e-8)
        delay=abs(base['arrival_offset']-fine['arrival_offset'])
        matrix=float(np.max(abs(np.array(base['transport']['matrices'][-1])-fine['transport']['matrices'][-1])))
        comparisons.append(dict(name=row['config']['name'],bend_relative_difference=bend,arrival_difference=delay,matrix_difference=matrix,
                                passed=bend<(.01 if time_case else .05) and delay<(1e-4 if time_case else .001) and matrix<(.001 if time_case else .005)))
    save(out/'summary.json',dict(runs=rows,comparisons=comparisons,passed=all(r['passed'] for r in rows+comparisons)))
    print('campaign complete',flush=True)


if __name__=='__main__':main()
