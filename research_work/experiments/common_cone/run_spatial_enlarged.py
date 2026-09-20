"""SR-1 amendment: identical dynamics/resolution, more boundary clearance."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
import subprocess
import numpy as np
from run_spatial import ROOT, save
from evolution import rk4
from spatial_model import SpatialEvolution


def main():
    original=json.loads((ROOT/'spatial-v1/summary.json').read_text())
    checks=json.loads((ROOT/'spatial-v1/controls.json').read_text())
    if not checks['passed']:raise RuntimeError('Original controls failed')
    out=ROOT/'spatial-v2';out.mkdir(exist_ok=False)
    sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('spatial-protocol.md','spatial-amendment.md','spatial_model.py','run_spatial.py','run_spatial_enlarged.py','evolution.py','model.py')}
    save(out/'manifest.json',dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),sources=sources,numpy=np.__version__,control_sha256=hashlib.sha256((ROOT/'spatial-v1/controls.json').read_bytes()).hexdigest()))
    summaries=[]
    for entry in original['runs']:
        cfg=dict(entry['config']);fine=cfg['name'].endswith('space')
        cfg.update(n=40 if fine else 28,length=15. if fine else 14.)
        print('start',cfg['name'],flush=True)
        model=SpatialEvolution(**{k:v for k,v in cfg.items() if k not in ('name','dt')})
        y=model.initial();initial=y.copy();trace=[];steps=round(2.5/cfg['dt'])
        for step in range(steps+1):
            metrics=model.metrics(y);metrics['time']=step*cfg['dt'];trace.append(metrics)
            if step<steps:y=rk4(model,y,cfg['dt'])
        maximum_speed=max(x['max_characteristic'] for x in trace)
        clearance=model.length/2-1-max(x['source_extent'] for x in trace)-2.5*maximum_speed
        drift=max(abs(x['ledger']-trace[0]['ledger']) for x in trace)/abs(trace[0]['ledger'])
        cone=max(x['cone_residual'] for x in trace);edge=max(x['edge_amplitude'] for x in trace)
        result=dict(config=cfg,initial=trace[0],final=trace[-1],ledger_drift=drift,cone_error=cone,edge_amplitude=edge,clearance=clearance,
                    bend=trace[-1]['probe_angle']-trace[0]['probe_angle'],momentum_drift=max(np.linalg.norm(x['momentum']-trace[0]['momentum']) for x in trace),angular_drift=max(np.linalg.norm(x['angular']-trace[0]['angular']) for x in trace))
        result['passed']=drift<1e-5 and cone<1e-10 and edge<1e-5 and clearance>0
        np.savez_compressed(out/(cfg['name']+'.npz'),initial=initial,final=y)
        save(out/(cfg['name']+'.json'),dict(summary=result,trace=trace));summaries.append(result)
        print('done',cfg['name'],'passed',result['passed'],'bend',result['bend'],flush=True)
    comparisons=[]
    def compare(left,right,leftdir,rightdir,name):
        cfg=left['config'];other=right['config']
        lm=SpatialEvolution(**{k:v for k,v in cfg.items() if k not in ('name','dt')})
        rm=SpatialEvolution(**{k:v for k,v in other.items() if k not in ('name','dt')})
        ly=np.load(leftdir/(cfg['name']+'.npz'))['final'];ry=np.load(rightdir/(other['name']+'.npz'))['final']
        position=float(np.linalg.norm(lm.unpack(ly)[2][-1]-rm.unpack(ry)[2][-1]))
        relative=abs(left['final']['field']-right['final']['field'])/abs(right['final']['field'])
        comparisons.append(dict(name=name,probe_position_difference=position,field_energy_relative_difference=relative,passed=position<.01 and relative<.1))
    base=next(x for x in summaries if x['config']['name']=='s1-fast')
    for tag in ('time','space'):
        compare(base,next(x for x in summaries if x['config']['name']==f's1-fast-{tag}'),out,out,tag)
    for old,new in zip(original['runs'],summaries):
        compare(old,new,ROOT/'spatial-v1',out,'domain-'+new['config']['name'])
    save(out/'summary.json',dict(runs=summaries,comparisons=comparisons,passed=all(x['passed'] for x in summaries+comparisons)))
    print('complete',all(x['passed'] for x in summaries+comparisons),flush=True)


if __name__=='__main__':main()
