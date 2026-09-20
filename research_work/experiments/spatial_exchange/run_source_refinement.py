"""SE-R seven source-resolution runs; source comparison is a separate audit."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from emitter_model import ExchangeEvolution,rk4,COMMON
from model import ExchangeEvolution as OriginalEvolution

ROOT=Path(__file__).resolve().parent


def save(path,obj):
    def cast(value):
        if isinstance(value,np.ndarray):return value.tolist()
        if isinstance(value,np.generic):return value.item()
        raise TypeError(type(value).__name__)
    path.write_text(json.dumps(obj,indent=2,default=cast)+'\n',encoding='utf-8')


def configurations():
    common=dict(length=16.,radius=.9,dt=.01,mix=.5,chi=200.)
    variants=[dict(name=f'mixed-n{n}',n=n,emitter_angle=np.pi/4) for n in (32,40,48,56)]
    variants += [dict(name=f'Y-n{n}',n=n,emitter_angle=np.pi/2) for n in (48,56)]
    variants += [dict(name='mixed-n40-rotation',n=40,emitter_angle=np.pi/4,angle=np.pi/3)]
    return [dict(common,**v) for v in variants]


def main():
    out=ROOT/'source-refinement-v1';out.mkdir(exist_ok=False)
    manifest=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('source-refinement-protocol.md','emitter_model.py','run_source_refinement.py')},helper_source=dict(path='research_work/experiments/common_cone/model.py',sha256=hashlib.sha256((COMMON/'model.py').read_bytes()).hexdigest()),configurations=configurations(),numpy=np.__version__)
    save(out/'manifest.json',manifest)
    rows=[]
    for cfg in configurations():
        print('start',cfg['name'],flush=True);wall=time.perf_counter();cpu=time.process_time()
        model=ExchangeEvolution(**{k:v for k,v in cfg.items() if k not in ('name','dt')});y=model.initial();initial=y.copy()
        trace=[];extent=0.;steps=round(4/cfg['dt']);stride=round(.1/cfg['dt'])
        for step in range(steps+1):
            extent=max(extent,float(np.max(abs(model.unpack(y)[2]))+model.radius))
            if step%stride==0 or step==steps:
                metric=model.metrics(y);metric['time']=step*cfg['dt'];trace.append(metric)
            if step<steps:y=rk4(model,y,cfg['dt'])
        drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger'])
        cone=max(t['cone_error'] for t in trace);edge=max(t['edge_amplitude'] for t in trace)
        clearance=model.length/2-1-extent-4*model.maximum_characteristic
        null=True
        if cfg['name'] in ('no-emission','no-excitation'):null=max(max(t['radiation_amplitude'],t['companion_amplitude']) for t in trace)<1e-12
        if cfg['name']=='emission-only':null=max(t['companion_amplitude'] for t in trace)<1e-12
        row=dict(config=cfg,initial=trace[0],final=trace[-1],relative_energy_drift=drift,cone_error=cone,edge_amplitude=edge,clearance=clearance,null_control_passed=null,
                 maximum_source_extent=extent,maximum_characteristic=model.maximum_characteristic,
                 momentum_drift=max(np.linalg.norm(np.array(t['momentum'])-trace[0]['momentum']) for t in trace),angular_drift=max(np.linalg.norm(np.array(t['angular'])-trace[0]['angular']) for t in trace),
                 wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,passed=drift<1e-5 and cone<1e-10 and edge<1e-5 and clearance>0 and null)
        np.savez_compressed(out/(cfg['name']+'.npz'),initial=initial,final=y)
        save(out/(cfg['name']+'.json'),dict(summary=row,trace=trace));rows.append(row)
        print('done',cfg['name'],'passed',row['passed'],'X',row['final']['radiation_propagation'],'Y',row['final']['companion_propagation'],flush=True)
    result=dict(runs=rows,numerical_gates_passed=all(r['passed'] for r in rows),
                source_accuracy_passed=None)
    save(out/'summary.json',result);print('evolutions complete; source comparisons pending',flush=True)



if __name__=='__main__':main()
