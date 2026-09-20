"""Read-only comparison of linear/cubic measurement of saved circulation."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np
from refine_local_transfer import rotation

ROOT=Path(__file__).resolve().parent


def weights(t,order):
    if order==1:return (0,1),np.stack((1-t,t))
    return (-1,0,1,2),np.stack((-t*(t-1)*(t-2)/6,(t+1)*(t-1)*(t-2)/2,-(t+1)*t*(t-2)/2,(t+1)*t*(t-1)/6))


def sample(field,position,h,order):
    grid=(position+4)/h;base=np.floor(grid).astype(int);frac=grid-base
    offsets,w=weights(frac,order);result=np.zeros((len(position),3));n=field.shape[0]
    for indices in itertools.product(range(len(offsets)),repeat=3):
        shift=np.array([offsets[k] for k in indices]);node=(base+shift)%n
        weight=w[indices[0],:,0]*w[indices[1],:,1]*w[indices[2],:,2]
        result+=weight[:,None]*field[node[:,0],node[:,1],node[:,2]]
    return result


def measure(field,h,radius,count,rot,order):
    t=np.arange(count)*2*np.pi/count
    position=np.stack((radius*np.cos(t),radius*np.sin(t),np.zeros(count)),axis=-1)@rot.T
    tangent=np.stack((-radius*np.sin(t),radius*np.cos(t),np.zeros(count)),axis=-1)@rot.T
    return float(np.sum(sample(field,position,h,order)*tangent)*2*np.pi/count)


def main():
    out=ROOT/'loop-interpolation-v1';out.mkdir(exist_ok=False)
    t=np.linspace(0,1,33);nodes,w=weights(t,3)
    polynomial_errors=[float(np.max(abs(np.sum(w*np.array(nodes)[:,None]**p,axis=0)-t**p))) for p in range(4)]
    controls=[]
    for n in (24,32,48,64):
        h=8/n;x=np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)),indexing='ij'),axis=-1)
        for rotated,rot in ((False,np.eye(3)),(True,rotation())):
            field=.02*np.exp(-np.sum(x*x,axis=-1)/1.4**2)[...,None]*np.cross(rot[:,2],x)
            for radius,order in itertools.product((1.,1.5,2.),(1,3)):
                exact=2*np.pi*.02*radius**2*np.exp(-radius**2/1.4**2)
                values=[measure(field,h,radius,count,rot,order) for count in (256,512)]
                error=abs(values[1]-exact)/abs(exact)
                controls.append(dict(n=n,rotated=rotated,radius=radius,order=order,exact=float(exact),values=values,relative_error=error,
                                     finest_gate_passed=bool(error<(.02 if order==1 else .001)) if n==64 else None))
    inputs={};observed=[]
    for name in ('n40','n48','n64','n48-rotation'):
        file=ROOT/'local-transfer-refinement-v1'/(name+'.npz');cfgfile=file.with_suffix('.json')
        inputs[file.name]=hashlib.sha256(file.read_bytes()).hexdigest();inputs[cfgfile.name]=hashlib.sha256(cfgfile.read_bytes()).hexdigest()
        cfg=json.loads(cfgfile.read_text());a=np.load(file)['final'][0];rot=np.array(cfg['rotation']);h=8/cfg['n']
        for radius,order in itertools.product((1.,1.5,2.),(1,3)):
            values=[measure(a,h,radius,count,rot,order) for count in (256,512)]
            observed.append(dict(name=name,radius=radius,order=order,values=values))
    comparisons=[]
    for order in (1,3):
        selected={r['name']:r['values'][1] for r in observed if r['order']==order and r['radius']==1.5}
        space=abs(selected['n48']-selected['n64'])/max(abs(selected['n64']),1e-10)
        turn=abs(selected['n48']-selected['n48-rotation'])/max(abs(selected['n48']),1e-10)
        comparisons.append(dict(order=order,spatial_difference=space,rotation_difference=turn,within_spatial_threshold=space<.05,within_rotation_threshold=turn<.02))
    result=dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_hashes=inputs,
                polynomial_errors=polynomial_errors,controls=controls,measurements=observed,comparisons=comparisons,
                controls_passed=max(polynomial_errors)<1e-12 and all(r['finest_gate_passed'] for r in controls if r['n']==64))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(controls_passed=result['controls_passed'],polynomial_max=max(polynomial_errors),comparisons=comparisons,
                         finest_linear_max=max(r['relative_error'] for r in controls if r['n']==64 and r['order']==1),
                         finest_cubic_max=max(r['relative_error'] for r in controls if r['n']==64 and r['order']==3))))
    assert result['controls_passed']


if __name__=='__main__':main()
