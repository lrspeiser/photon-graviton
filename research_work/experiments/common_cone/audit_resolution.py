"""SR-2 raw-state energy audit; explicitly distinguishes partial evidence."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parent


def reconstruct(y,cfg):
    n=cfg['n'];L=cfg['length'];h=L/n;nf=4*n**3
    f=y[:nf].reshape(4,n,n,n);pi=y[nf:2*nf].reshape(f.shape)
    q=y[2*nf:2*nf+21].reshape(7,3);p=y[2*nf+21:-1].reshape(7,3)
    U=.08*f[0];alpha=np.exp(U);z=np.exp(2*U);a=np.exp(4*U);c=z
    beta=c*.04*f[1:]/np.sqrt(1+.08**2*np.sum(f[1:]**2,axis=0))
    density=a*np.sum(pi*pi,axis=0)/2+.2**2*np.sum(f*f,axis=0)/2
    indices=np.arange(n)
    for axis in range(1,4):
        def at(offset):return np.take(f,(indices+offset)%n,axis=axis)
        if cfg['order']==2:
            forward=(at(1)-f)/h;backward=(f-at(-1))/h
            grad=(forward+backward)/2
            density+=np.sum(forward**2+backward**2,axis=0)/4
        else:
            grad=sum(weight*at(offset) for offset,weight in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*h)
            fourth=sum(weight*at(offset) for offset,weight in ((-2,1),(-1,-4),(0,6),(1,-4),(2,1)))/h**4
            density+=np.sum(grad*grad,axis=0)/2+h**6*np.sum(fourth*fourth,axis=0)/128
        density-=beta[axis-1]*np.sum(pi*grad,axis=0)
    field=float(np.sum(density)*h**3)
    coordinate=np.arange(n)*h-L/2;xyz=np.stack(np.meshgrid(coordinate,coordinate,coordinate,indexing='ij'),axis=-1)
    matter=0.;velocities=[];centers=[];speeds=[]
    for i in range(7):
        w=np.maximum(1-np.sum((xyz-q[i])**2,axis=-1)/.9**2,0)**3;w/=w.sum()
        mass=0. if i==6 else 1.;e=np.sqrt(mass*mass+z*np.dot(p[i],p[i]))
        matter+=float(np.sum(w*(alpha*e+np.einsum('jxyz,j->xyz',beta,p[i]))))
        center=np.sum(w*beta,axis=(1,2,3));speed=np.sum(w*c)
        velocities.append(p[i]*np.sum(w*alpha*z/e)+center);centers.append(center);speeds.append(speed)
    velocity=np.array(velocities);centers=np.array(centers);speeds=np.array(speeds)
    relative=np.linalg.norm(velocity-centers,axis=1)/speeds
    return dict(field=field,matter=matter,total=field+matter,ledger=field+matter+y[-1],q=q,p=p,velocity=velocity,
                cone_error=max(0.,float(relative.max()-1),abs(float(relative[-1]-1))))


def rotation(angle):
    u=np.array([1.,2.,3.])/np.sqrt(14);x,y,z=u
    skew=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    return np.cos(angle)*np.eye(3)+(1-np.cos(angle))*np.outer(u,u)+np.sin(angle)*skew


def main():
    directory=ROOT/'resolution-v1';manifest=json.loads((directory/'manifest.json').read_text());controls=json.loads((directory/'controls.json').read_text())
    checks=[];errors=[];rows=[];states={}
    def check(name,value):checks.append(dict(name=name,passed=bool(value)))
    for name,digest in manifest['sources'].items():
        check('current '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
        blob=subprocess.check_output(['git','show',manifest['source_commit']+':research_work/experiments/common_cone/'+name],cwd=ROOT)
        check('committed '+name,hashlib.sha256(blob).hexdigest()==digest)
    check('control counts',len(controls['equivalence'])==36 and len(controls['gradients'])==48 and controls['fourier']['nonzero_modes']==4095)
    for row in controls['equivalence']:check('equivalence '+str(row['index']),max(row['rhs_error'],row['energy_error'],row['cone_error'])<1e-11 and row['passed'])
    for row in controls['gradients']:check('gradient '+str((row['order'],row['index'])),row['error']<2e-6 and row['passed'])
    check('Fourier positivity',controls['fourier']['minimum_nonzero_stiffness']>0 and controls['fourier']['zero_mode']==0 and controls['fourier']['passed'])
    for row in controls['spectral']:check('spectral '+str(row['mode']),row['relative_error']<1e-12 and row['passed'])
    for row in controls['derivative']:check('derivative order '+str(row['mode']),min(row['ratios'])>12 and row['passed'])
    for cfg in manifest['configurations']:
        name=cfg['name'];path=directory/(name+'.json')
        if not path.exists():continue
        doc=json.loads(path.read_text());row=doc['summary'];trace=doc['trace'];raw=np.load(directory/(name+'.npz'));rows.append(row)
        check(name+' configuration',row['config']==cfg)
        check(name+' horizon',trace[0]['time']==0 and trace[-1]['time']==2.5)
        for endpoint in ('initial','final'):
            value=reconstruct(raw[endpoint],cfg);err=max(abs(value[k]-row[endpoint][k]) for k in ('field','matter','total','ledger'));errors.append(err)
            check(name+' '+endpoint+' independent energy',err<1e-10)
            check(name+' '+endpoint+' independent cone',value['cone_error']<1e-10)
            check(name+' '+endpoint+' probe velocity',np.max(abs(value['velocity'][-1]-row[endpoint]['probe_velocity']))<1e-11)
        states[name]=value
        v0=np.array(row['initial']['probe_velocity']);v1=value['velocity'][-1]
        magnitude=float(np.arctan2(np.linalg.norm(np.cross(v0,v1)),np.dot(v0,v1)))
        back=v1@rotation(cfg.get('angle',0));signed=float(np.arctan2(back[1],back[0]))
        check(name+' bend definition',abs(magnitude-row['unsigned_bend'])<1e-12 and abs(signed-row['bend'])<1e-12)
        drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger'])
        cone=max(t['cone_residual'] for t in trace);edge=max(t['edge_amplitude'] for t in trace)
        clearance=cfg['length']/2-1-row['maximum_source_extent']-2.5*row['maximum_characteristic']
        check(name+' extrema bounds',row['maximum_source_extent']>=max(t['source_extent'] for t in trace)-1e-12 and row['maximum_characteristic']>=max(t['max_characteristic'] for t in trace)-1e-12)
        gate=drift<1e-5 and cone<1e-10 and edge<1e-5 and clearance>0
        check(name+' gate honesty',abs(drift-row['ledger_drift'])<1e-12 and abs(clearance-row['clearance'])<1e-12 and gate==row['passed'])
    complete=(directory/'summary.json').exists();accuracy=None
    if complete:
        summary=json.loads((directory/'summary.json').read_text());check('complete run count',len(rows)==9 and summary['runs']==rows)
        by={r['config']['name']:r for r in rows};acc=summary['accuracy']
        b40=abs(by['o4-n40']['bend']);b56=abs(by['o4-n56']['bend']);b70=abs(by['o4-n70']['bend'])
        computed=dict(finest_relative_bend_difference=abs(b70-b56)/b70,preceding_absolute_difference=abs(b56-b40),latest_absolute_difference=abs(b70-b56),finest_method_relative_difference=abs(by['o2-n70']['bend']/by['o4-n70']['bend']-1),time_relative_bend_difference=abs(by['o4-n56-time']['bend']/by['o4-n56']['bend']-1))
        vb=states['o4-n56-rotated']['velocity'][-1]@rotation(np.pi/3);vu=states['o4-n56']['velocity'][-1]
        rangle=float(np.arctan2(np.linalg.norm(np.cross(vb,vu)),np.dot(vb,vu)))
        computed.update(rotation_direction_error_radians=rangle,rotation_error_relative_to_bend=rangle/b56)
        old=json.loads((ROOT/'spatial-v2/s1-fast.json').read_text())['summary']
        computed['replay_error']=max(abs(by['o2-n28']['bend']-old['bend']),abs(by['o2-n28']['final']['total']-old['final']['total']))
        for key,value in computed.items():check('accuracy '+key,abs(value-acc[key])<1e-10)
        accuracy=computed['finest_relative_bend_difference']<.01 and computed['latest_absolute_difference']<computed['preceding_absolute_difference'] and computed['finest_method_relative_difference']<.01 and computed['time_relative_bend_difference']<.001 and computed['rotation_error_relative_to_bend']<.01 and computed['replay_error']<1e-10
        check('accuracy decision honest',accuracy==acc['passed']==summary['accuracy_passed'])
        check('numerical decision honest',summary['numerical_passed']==all(r['passed'] for r in rows))
    output=dict(passed=all(r['passed'] for r in checks),complete=complete,completed_runs=len(rows),declared_runs=9,count=len(checks),max_independent_energy_error=max(errors,default=0),accuracy_passed=accuracy,checks=checks)
    (ROOT/'resolution-audit.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in output.items() if k!='checks'}))
    if not output['passed']:raise RuntimeError('SR-2 archive audit failed')


if __name__=='__main__':main()
