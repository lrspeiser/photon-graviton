"""SE-C circulation and angular attribution; no scalar-to-swirl multiplier."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from audit_emitter import rotation

ROOT=Path(__file__).resolve().parent


def loops(beta,line,rot):
    sample=RegularGridInterpolator((line,line,line),beta.transpose(1,2,3,0),bounds_error=True)
    rows=[]
    for radius in (1,2,3,4):
        values=[]
        for n in (128,256,512):
            theta=2*np.pi*np.arange(n)/n
            point=np.column_stack((radius*np.cos(theta),radius*np.sin(theta),np.zeros(n)))@rot.T
            tangent=np.column_stack((-radius*np.sin(theta),radius*np.cos(theta),np.zeros(n)))@rot.T
            values.append(float(np.sum(sample(point)*tangent)*2*np.pi/n))
        error=abs(values[-1]-values[-2]);rows.append(dict(radius=radius,values=values,quadrature_error=error,passed=error<1e-9+.005*abs(values[-1])))
    return rows


def analyze(folder,cfg):
    path=folder/(cfg['name']+'.npz');state=np.load(path)['final'];doc=json.loads(path.with_suffix('.json').read_text())['summary']
    n=cfg['n'];h=cfg['length']/n;nf=6*n**3
    f=state[:nf].reshape(6,n,n,n);pi=state[nf:2*nf].reshape(f.shape)
    q=state[2*nf:2*nf+18].reshape(6,3);p=state[2*nf+18:2*nf+36].reshape(6,3)
    line=np.arange(n)*h-cfg['length']/2;xyz=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1)
    current=np.empty((3,n,n,n));wave=np.empty_like(current)
    for j in range(3):
        grad=(np.roll(f,-1,j+1)-np.roll(f,1,j+1))/(2*h)
        current[j]=np.sum(pi*grad,axis=0);wave[j]=np.sum(pi[4:]*grad[4:],axis=0)
    orbital=-np.cross(xyz,current.transpose(1,2,3,0)).sum(axis=(0,1,2))*h**3
    wave_angular=-np.cross(xyz,wave.transpose(1,2,3,0)).sum(axis=(0,1,2))*h**3
    spin=np.cross(f[1:4].transpose(1,2,3,0),pi[1:4].transpose(1,2,3,0)).sum(axis=(0,1,2))*h**3
    matter=np.cross(q,p).sum(axis=0);total=matter+orbital+spin
    error=float(np.max(abs(total-doc['final']['angular'])))
    rot=rotation(cfg.get('angle',0));axis=rot@np.array([0.,0.,1.]);initial=float(np.dot(doc['initial']['angular'],axis))
    beta=np.exp(.16*f[0])*.04*f[1:4]/np.sqrt(1+.08**2*np.sum(f[1:4]**2,axis=0))
    circles=loops(beta,line,rot)
    return dict(campaign=folder.name,name=cfg['name'],loops=circles,
                angular_projection={name:float(value@axis) for name,value in [('matter',matter),('field_orbital',orbital),('XY_orbital',wave_angular),('field_spin',spin),('total',total)]},
                XY_fraction_initial_angular=float(wave_angular@axis)/initial if abs(initial)>1e-15 else None,
                angular_reconstruction_error=error,passed=error<1e-10 and all(c['passed'] for c in circles),
                state_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    out=ROOT/'circulation-v1';out.mkdir(exist_ok=False)
    line=np.arange(32)*.5-8;xyz=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=0)
    beta=np.zeros_like(xyz);beta[0]=-.001*xyz[1];beta[1]=.001*xyz[0]
    controls=loops(beta,line,np.eye(3));control_error=max(abs(v-2*np.pi*.001*c['radius']**2) for c in controls for v in c['values'])
    rows=[]
    for name in ['evidence-v1','emitter-v1']:
        folder=ROOT/name;manifest=json.loads((folder/'manifest.json').read_text())
        for cfg in manifest['configurations']:
            rows.append(analyze(folder,cfg))
    result=dict(passed=bool(control_error<1e-12 and all(r['passed'] for r in rows)),count=len(rows),control_error=control_error,
                commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['circulation.py','circulation-protocol.md','audit_emitter.py']})
    (out/'results.json').write_text(json.dumps(dict(summary=result,rows=rows),indent=2)+'\n')
    print(json.dumps(result))
    for row in rows:print(row['campaign'],row['name'],'loop r2',row['loops'][1]['values'][-1],'XY angular fraction',row['XY_fraction_initial_angular'])
    if not result['passed']:raise RuntimeError('SE-C gate failed; retain archive')


if __name__=='__main__':main()
