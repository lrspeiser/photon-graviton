"""Independent raw-state SR-1 energy reconstruction and archive verification."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from spatial_model import SpatialEvolution

ROOT=Path(__file__).resolve().parent


def reconstruct(y,cfg):
    n=cfg['n'];L=cfg.get('length',12.);dx=L/n;dv=dx**3;nf=4*n**3
    f=y[:nf].reshape(4,n,n,n);pi=y[nf:2*nf].reshape(f.shape)
    q=y[2*nf:2*nf+21].reshape(7,3);p=y[2*nf+21:-1].reshape(7,3)
    U=.08*f[0];s=cfg['s'];eta=cfg['eta']
    alpha=np.exp(U);z=np.exp(2*s*U);a=np.exp((1+3*s)*U);d=np.exp((1-s)*U)
    c=np.exp((1+s)*U);beta=c*.5*eta*f[1:]/np.sqrt(1+eta**2*np.sum(f[1:]**2,axis=0))
    density=a*np.sum(pi*pi,axis=0)/2+.2**2*np.sum(f*f,axis=0)/2
    for axis in range(1,4):
        forward=(np.roll(f,-1,axis=axis)-f)/dx
        backward=(f-np.roll(f,1,axis=axis))/dx
        density+=d*np.sum(forward**2+backward**2,axis=0)/4
        density-=beta[axis-1]*np.sum(pi*(forward+backward)/2,axis=0)
    field=float(np.sum(density)*dv)
    coordinate=np.arange(n)*dx-L/2
    xyz=np.stack(np.meshgrid(coordinate,coordinate,coordinate,indexing='ij'),axis=-1)
    matter=0.;velocity=[];centers=[];speeds=[]
    for i in range(7):
        w=np.maximum(1-np.sum((xyz-q[i])**2,axis=-1)/.9**2,0)**3;w/=w.sum()
        mass=0. if i==6 else 1.
        e=np.sqrt(mass**2+z*np.dot(p[i],p[i]))
        matter+=float(np.sum(w*(alpha*e+np.einsum('jxyz,j->xyz',beta,p[i]))))
        center=np.sum(w*beta,axis=(1,2,3));speed=np.sum(w*c)
        velocity.append(p[i]*np.sum(w*alpha*z/e)+center);centers.append(center);speeds.append(speed)
    return dict(field=field,matter=matter,total=field+matter,ledger=field+matter+y[-1],q=q,p=p,
                velocity=np.array(velocity),centers=np.array(centers),speeds=np.array(speeds))


def main():
    checks=[];energy_errors=[];direction_errors=[];rng=np.random.default_rng(704)
    def check(name,value):checks.append(dict(name=name,passed=bool(value)))
    controls=json.loads((ROOT/'spatial-v1/controls.json').read_text())
    check('control counts',len(controls['local'])==240 and len(controls['grid'])==24 and len(controls['negative'])==4)
    for row in controls['local']:
        expected=row['gradient_error']<2e-6 and row['eigenvalue_error']<1e-10 and row['symmetry_error']<1e-10 and row['minimum_energy_eigenvalue']>0 and row['momentum_hessian_minimum']>-1e-10 and row['cone_error']<1e-12
        check('local '+str(row['index']),expected and row['passed'])
    for k,row in enumerate(controls['grid']):check('grid '+str(k),row['error']<2e-6 and row['passed'])
    for k,row in enumerate(controls['negative']):check('negative '+str(k),row['omitted_error']>max(100*row['correct_error'],1e-6) and row['passed'])
    check('CC-2 control equivalence',controls['cc2_equivalence']['energy']<1e-11 and controls['cc2_equivalence']['rhs']<1e-11)
    archives={};states={}
    for version in ('spatial-v1','spatial-v2'):
        directory=ROOT/version;manifest=json.loads((directory/'manifest.json').read_text());summary=json.loads((directory/'summary.json').read_text());archives[version]=summary
        for name,digest in manifest['sources'].items():
            check(version+' current source '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
            committed=subprocess.check_output(['git','show',manifest['source_commit']+':research_work/experiments/common_cone/'+name],cwd=ROOT)
            check(version+' committed source '+name,hashlib.sha256(committed).hexdigest()==digest)
        if version=='spatial-v2':check('reused controls hash',manifest['control_sha256']==hashlib.sha256((ROOT/'spatial-v1/controls.json').read_bytes()).hexdigest())
        check(version+' run count',len(summary['runs'])==8)
        for row in summary['runs']:
            cfg=row['config'];name=cfg['name'];label=version+'/'+name
            raw=np.load(directory/(name+'.npz'));doc=json.loads((directory/(name+'.json')).read_text());trace=doc['trace']
            check(label+' summary identity',row==doc['summary'])
            check(label+' time samples',len(trace)==round(2.5/cfg['dt'])+1 and trace[-1]['time']==2.5)
            for state in ('initial','final'):
                rebuilt=reconstruct(raw[state],cfg);saved=row[state]
                error=max(abs(rebuilt[k]-saved[k]) for k in ('field','matter','total','ledger'));energy_errors.append(error)
                check(label+' '+state+' independent energy',error<1e-10)
                speed=np.linalg.norm(rebuilt['velocity']-rebuilt['centers'],axis=1)/rebuilt['speeds']
                check(label+' '+state+' independent cone',speed.max()<1+1e-10 and abs(speed[-1]-1)<1e-10)
            states[(version,name)]=reconstruct(raw['final'],cfg)
            drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger'])
            cone=max(t['cone_residual'] for t in trace);edge=max(t['edge_amplitude'] for t in trace)
            clearance=cfg.get('length',12.)/2-1-max(t['source_extent'] for t in trace)-2.5*max(t['max_characteristic'] for t in trace)
            expected=drift<1e-5 and cone<1e-10 and edge<1e-5 and clearance>0
            check(label+' gate recorded honestly',expected==row['passed'] and abs(clearance-row['clearance'])<1e-12 and abs(drift-row['ledger_drift'])<1e-12)
            check(label+' bend reference',abs(row['bend']-(trace[-1]['probe_angle']-trace[0]['probe_angle']))<1e-12)
            model=SpatialEvolution(**{k:v for k,v in cfg.items() if k not in ('name','dt')});y=raw['final'];dy=model.rhs(y)
            fd,pd,qd,ppd=model.unpack(dy)
            gradient=np.concatenate(((-model.dv*(pd+model.gamma*fd)).ravel(),(model.dv*fd).ravel(),(-ppd).ravel(),qd.ravel(),[0.]))
            for k in range(3):
                direction=rng.normal(size=y.size);direction[-1]=0.;direction/=np.linalg.norm(direction)
                eps=2e-5
                numeric=(reconstruct(y+eps*direction,cfg)['total']-reconstruct(y-eps*direction,cfg)['total'])/(2*eps)
                error=abs(numeric-np.dot(gradient,direction));direction_errors.append(error)
                check(label+' final Hamiltonian derivative '+str(k),error<2e-6)
        for comparison in summary['comparisons']:
            name=comparison['name']
            if name.startswith('domain-'):
                left=states[('spatial-v1',name[7:])];right=states[('spatial-v2',name[7:])]
            else:
                left=states[(version,'s1-fast')];right=states[(version,'s1-fast-'+name)]
            position=float(np.linalg.norm(left['q'][-1]-right['q'][-1]));relative=abs(left['field']-right['field'])/abs(right['field'])
            check(version+' comparison '+name,abs(position-comparison['probe_position_difference'])<1e-12 and abs(relative-comparison['field_energy_relative_difference'])<1e-10 and comparison['passed']==(position<.01 and relative<.1))
        check(version+' overall gate',summary['passed']==all(r['passed'] for r in summary['runs']+summary['comparisons']))
    result=dict(passed=all(c['passed'] for c in checks),checks=checks,count=len(checks),max_independent_energy_error=max(energy_errors),max_independent_gradient_error=max(direction_errors),original_campaign_passed=archives['spatial-v1']['passed'],enlarged_campaign_passed=archives['spatial-v2']['passed'])
    (ROOT/'spatial-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
    if not result['passed']:raise RuntimeError('Archive audit failed')


if __name__=='__main__':main()
