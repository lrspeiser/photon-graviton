"""SE-LR4, unchanged local transfer equations with numerical refinement."""
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.ndimage import map_coordinates

from local_transfer import rhs, ledger
from local_rotor import derivative

ROOT = Path(__file__).resolve().parent


def rotation():
    axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
    x, y, z = axis
    cross = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3)+np.sin(.573)*cross+(1-np.cos(.573))*(cross@cross)


def loop(a, h, radius, count, rot):
    theta = np.arange(count)*2*np.pi/count
    position = rot@np.array([radius*np.cos(theta), radius*np.sin(theta), np.zeros(count)])
    tangent = rot@np.array([-radius*np.sin(theta), radius*np.cos(theta), np.zeros(count)])
    values = np.array([map_coordinates(a[..., j], (position+4)/h, order=1, mode='wrap', prefilter=False) for j in range(3)])
    return float(np.sum(values*tangent)*2*np.pi/count)


def outer(state, coords, h, axis):
    a, pi, q, p = state
    ma = -np.stack([np.sum(pi*derivative(a, j, h), axis=-1) for j in range(3)], axis=-1)
    mq = -np.stack([np.sum(p*derivative(q, j, h), axis=-1) for j in range(3)], axis=-1)
    ja = (np.cross(coords, ma)+np.cross(a, pi))@axis
    total = ja+(np.cross(coords, mq)+np.cross(q, p))@axis
    squared = np.sum(coords**2, axis=-1)
    return [float(np.sum(density[squared > radius**2])*h**3) for radius in (1.2, 1.5, 2.) for density in (ja, total)]


def main():
    out = ROOT/'local-transfer-refinement-v1'; out.mkdir(exist_ok=False)
    filenames = ('local_rotor.py', 'run_local_rotor.py', 'local_transfer.py', 'refine_local_transfer.py')
    manifest = dict(commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in filenames},
                    trace_columns=['time','energy','A_energy','KQ_energy','Jx','Jy','Jz','Px','Py','Pz','edge_energy','A_Jx','A_Jy','A_Jz',
                                   'A_J_outside_1.2','total_J_outside_1.2','A_J_outside_1.5','total_J_outside_1.5','A_J_outside_2','total_J_outside_2'])
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    rows = []
    for name, n, rot in [('n40',40,np.eye(3)), ('n48',48,np.eye(3)), ('n64',64,np.eye(3)), ('n48-rotation',48,rotation())]:
        print('start',name,flush=True); start = time.monotonic()
        h=8/n; dt=.01; axis=rot[:,2]
        coords=np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)), indexing='ij'),axis=-1)
        f=.03*np.maximum(1-np.sum(coords**2,axis=-1)/1.2**2,0)**3
        state=np.zeros((4,n,n,n,3)); state[2]=f[...,None]*rot[:,0]; state[3]=f[...,None]*rot[:,1]
        initial=state.copy()
        trace=[[0.,*ledger(state,coords,h,2.,1.),*outer(state,coords,h,axis)]]
        for step in range(100):
            k1=rhs(state,h,2.,1.); k2=rhs(state+.5*dt*k1,h,2.,1.)
            k3=rhs(state+.5*dt*k2,h,2.,1.); k4=rhs(state+dt*k3,h,2.,1.)
            state += dt*(k1+2*k2+2*k3+k4)/6
            trace.append([(step+1)*dt,*ledger(state,coords,h,2.,1.),*outer(state,coords,h,axis)])
        trace=np.array(trace)
        loops=[dict(radius=r,coarse=loop(state[0],h,r,128,rot),fine=loop(state[0],h,r,256,rot)) for r in (1.,1.5,2.)]
        e=float(np.max(abs(trace[:,1]-trace[0,1]))/trace[0,1])
        j=float(np.max(np.linalg.norm(trace[:,4:7]-trace[0,4:7],axis=1))/np.linalg.norm(trace[0,4:7]))
        edge=float(np.max(trace[:,10])/trace[0,1]); j0=float(trace[0,4:7]@axis)
        ja=float(trace[-1,11:14]@axis)
        row=dict(name=name,n=n,dt=dt,rotation=rot.tolist(),energy_drift=e,angular_drift=j,edge_fraction=edge,
                 initial_energy=float(trace[0,1]),initial_J=j0,A_J=ja,A_angular_fraction=ja/j0,
                 A_energy_fraction=float(trace[-1,2]/trace[0,1]),outer_angular_fractions=(trace[-1,14:20]/j0).tolist(),loops=loops,
                 passed=bool(np.isfinite(trace).all() and np.isfinite(state).all() and e<1e-5 and j<.01 and edge<1e-5
                             and max(abs(r['coarse']-r['fine']) for r in loops)<1e-5),seconds=time.monotonic()-start)
        np.savez_compressed(out/(name+'.npz'),initial=initial,final=state,trace=trace)
        (out/(name+'.json')).write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
        print(json.dumps(row),flush=True)
    old=json.loads((ROOT/'local-transfer-v1'/'space.json').read_text())
    comparisons=[]
    for label, coarse, middle, fine in [('circulation',old['loops'][1]['fine'],rows[1]['loops'][1]['fine'],rows[2]['loops'][1]['fine']),
                                        ('angular',old['A_Jz'],rows[1]['A_J'],rows[2]['A_J'])]:
        previous=abs(coarse-middle)/max(abs(middle),1e-10); error=abs(middle-fine)/max(abs(fine),1e-10)
        comparisons.append(dict(name=label+'-space',relative_difference=error,previous_difference=previous,passed=error<.05 and error<previous))
    for label, base, rotated in [('circulation',rows[1]['loops'][1]['fine'],rows[3]['loops'][1]['fine']),('angular',rows[1]['A_J'],rows[3]['A_J'])]:
        error=abs(base-rotated)/max(abs(base),1e-10)
        comparisons.append(dict(name=label+'-rotation',relative_difference=error,passed=error<.02))
    summary=dict(runs=rows,comparisons=comparisons,passed=all(r['passed'] for r in rows+comparisons))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(dict(passed=summary['passed'],comparisons=comparisons)),flush=True)


if __name__=='__main__':main()
