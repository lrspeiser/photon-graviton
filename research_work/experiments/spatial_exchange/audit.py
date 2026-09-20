"""Independent SE-1 energy reconstruction and honest partial/full audit."""
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
    n=cfg['n'];L=cfg['length'];radius=cfg['radius'];h=L/n;nf=6*n**3
    f=y[:nf].reshape(6,n,n,n);pi=y[nf:2*nf].reshape(f.shape);offset=2*nf
    q=y[offset:offset+18].reshape(6,3);p=y[offset+18:offset+36].reshape(6,3)
    Q=y[offset+36:offset+42];P=y[offset+42:offset+48]
    U=.08*f[0];alpha=np.exp(U);z=np.exp(2*U);a=np.exp(4*U)
    beta=z*.04*f[1:4]/np.sqrt(1+.08**2*np.sum(f[1:4]**2,axis=0))
    propagation=a[None]*pi*pi/2
    for axis in range(1,4):
        forward=(np.roll(f,-1,axis)-f)/h;backward=(f-np.roll(f,1,axis))/h
        propagation+=(forward*forward+backward*backward)/4-beta[axis-1]*pi*(forward+backward)/2
    mix=cfg.get('mix',.5)*np.tanh(U/.001);mass2=.75**2*np.exp(2*cfg.get('chi',0)*U)
    potential=mass2*(f[5]-mix*f[4])**2/2
    field=float((propagation.sum()+.2**2*np.sum(f[:4]**2)/2+potential.sum())*h**3)
    line=np.arange(n)*h-L/2;xyz=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1)
    matter=0.;internal=0.;emission=cfg.get('emission',.4);vel=[];cone=[]
    for i in range(6):
        weight=np.maximum(1-np.sum((xyz-q[i])**2,axis=-1)/radius**2,0)**3;weight/=weight.sum()
        sampled_x=np.sum(weight*f[4]);mass=1+(P[i]**2+(Q[i]-emission*sampled_x)**2)/2;internal+=mass-1
        e=np.sqrt(mass*mass+z*np.dot(p[i],p[i]))
        matter+=float(np.sum(weight*(alpha*e+np.einsum('jxyz,j->xyz',beta,p[i]))))
        center=np.sum(weight*beta,axis=(1,2,3));speed=np.sum(weight*z)
        velocity=p[i]*np.sum(weight*alpha*z/e)+center;vel.append(velocity);cone.append(np.linalg.norm(velocity-center)/speed)
    return dict(total=field+matter,field=field,matter=matter,ledger=field+matter+y[-1],internal_rest=internal,
                radiation_propagation=float(propagation[4].sum()*h**3),companion_propagation=float(propagation[5].sum()*h**3),mixing_potential=float(potential.sum()*h**3),
                q=q,p=p,Q=Q,P=P,velocity=np.array(vel),cone_error=max(0.,max(cone)-1))


def rotation(angle):
    u=np.array([1.,2.,3.])/np.sqrt(14);x,y,z=u;cross=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    return np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(u,u)+np.sin(angle)*cross


def main():
    directory=ROOT/'evidence-v1';manifest=json.loads((directory/'manifest.json').read_text());controls=json.loads((directory/'controls.json').read_text())
    checks=[];errors=[];rows=[];states={}
    def check(name,value):checks.append(dict(name=name,passed=bool(value)))
    for name,digest in manifest['sources'].items():
        check('current '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
        committed=subprocess.check_output(['git','show',manifest['source_commit']+':research_work/experiments/spatial_exchange/'+name],cwd=ROOT)
        check('committed '+name,hashlib.sha256(committed).hexdigest()==digest)
    helper=manifest['helper_source'];repo=Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=ROOT,text=True).strip())
    check('helper current',hashlib.sha256((repo/helper['path']).read_bytes()).hexdigest()==helper['sha256'])
    check('helper committed',hashlib.sha256(subprocess.check_output(['git','show',manifest['source_commit']+':'+helper['path']],cwd=ROOT)).hexdigest()==helper['sha256'])
    check('control count',len(controls['rows'])==48)
    for row in controls['rows']:
        expected=row['gradient_error']<2e-6 and row['minimum_mass']>0 and row['cone_error']<1e-10 and row['principal_eigenvalue_error']<1e-10 and row['symmetrizer_error']<1e-10 and row['mixing_matrix_minimum_eigenvalue']>-1e-12
        check('control '+str(row['index']),expected and row['passed'])
    for cfg in manifest['configurations']:
        name=cfg['name'];path=directory/(name+'.json')
        if not path.exists():continue
        doc=json.loads(path.read_text());row=doc['summary'];trace=doc['trace'];raw=np.load(directory/(name+'.npz'));rows.append(row)
        check(name+' configuration',cfg==row['config'])
        check(name+' complete horizon',trace[0]['time']==0 and trace[-1]['time']==4)
        for endpoint in ('initial','final'):
            data=reconstruct(raw[endpoint],cfg)
            error=max(abs(data[k]-row[endpoint][k]) for k in ('total','field','matter','ledger','internal_rest','radiation_propagation','companion_propagation','mixing_potential'));errors.append(error)
            check(name+' '+endpoint+' raw energy',error<1e-10)
            check(name+' '+endpoint+' cone',data['cone_error']<1e-10)
        states[name]=data
        drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger']);cone=max(t['cone_error'] for t in trace);edge=max(t['edge_amplitude'] for t in trace)
        clearance=cfg['length']/2-1-row['maximum_source_extent']-4*row['maximum_characteristic']
        null=True
        if name in ('no-emission','no-excitation'):null=max(max(t['radiation_amplitude'],t['companion_amplitude']) for t in trace)<1e-12
        if name=='emission-only':null=max(t['companion_amplitude'] for t in trace)<1e-12
        expected=drift<1e-5 and cone<1e-10 and edge<1e-5 and clearance>0 and null
        check(name+' honest gates',abs(drift-row['relative_energy_drift'])<1e-12 and abs(clearance-row['clearance'])<1e-12 and null==row['null_control_passed'] and expected==row['passed'])
        check(name+' extrema bounds',row['maximum_source_extent']>=max(t['source_extent'] for t in trace)-1e-12 and row['maximum_characteristic']>=max(t['max_characteristic'] for t in trace)-1e-12)
    complete=(directory/'summary.json').exists();campaign_passed=None
    if complete:
        summary=json.loads((directory/'summary.json').read_text());check('complete count',len(rows)==11 and summary['runs']==rows);base=states['exchange-200']
        for comparison in summary['comparisons']:
            name=comparison['name'];data=states[name];cfg=next(r['config'] for r in rows if r['config']['name']==name)
            q=data['q']@rotation(cfg.get('angle',0));position=float(np.max(np.linalg.norm(q-base['q'],axis=1)))
            error=max(abs(data[k]-base[k])/max(abs(base[k]),1e-12) for k in ('radiation_propagation','companion_propagation','mixing_potential'))
            if name=='sign-mirror':
                matter=max(float(np.max(abs(data[k]-base[k]))) for k in ('q','p','Q','P'));field=abs(data['field']/base['field']-1)
                passed=matter<1e-7 and field<1e-8
                check(name+' sign symmetry',abs(matter-comparison['matter_state_difference'])<1e-12 and abs(field-comparison['total_field_relative_difference'])<1e-10)
            else:passed=position<(1e-3 if name=='time-refinement' else .01) and error<(.01 if name=='time-refinement' else .05)
            check(name+' comparison honesty',abs(position-comparison['position_difference'])<1e-12 and abs(error-comparison['channel_energy_relative_difference'])<1e-8 and passed==comparison['passed'])
        campaign_passed=all(r['passed'] for r in summary['runs']+summary['comparisons']);check('aggregate honesty',campaign_passed==summary['passed'])
    output=dict(passed=all(c['passed'] for c in checks),complete=complete,completed_runs=len(rows),declared_runs=11,campaign_passed=campaign_passed,count=len(checks),max_independent_energy_error=max(errors,default=0),checks=checks)
    (ROOT/'audit.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in output.items() if k!='checks'}))
    if not output['passed']:raise RuntimeError('SE-1 archive audit failed')


if __name__=='__main__':main()
