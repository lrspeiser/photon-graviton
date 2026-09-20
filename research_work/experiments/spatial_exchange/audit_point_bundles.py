"""Independent point-bundle detector and cardinal-basis reconstruction."""
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_emitter import reconstruct

ROOT=Path(__file__).resolve().parent


def main():
    directory=ROOT/'point-bundle-v1';manifest=json.loads((directory/'manifest.json').read_text())
    checks=[];rows=[];derived={}
    for name,digest in manifest['sources'].items():checks.append(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    for cfg in manifest['configurations']:
        path=directory/(cfg['name']+'.json')
        if not path.exists():continue
        doc=json.loads(path.read_text());row=doc['summary'];rows.append(row)
        raw=np.load(directory/(cfg['name']+'.npz'));t=raw['time'];states=raw['probe_states'];vel=raw['probe_velocities'];energies=raw['probe_energies']
        checks.append(states.shape==(len(t),15,6) and len(t)==round(4/cfg['dt'])+1 and t[-1]==4)
        for endpoint in ['initial','final']:
            energy=reconstruct(raw[endpoint+'_background'],cfg)['ledger']
            reference=doc['background_trace'][0 if endpoint=='initial' else -1]['ledger']
            checks.append(abs(energy-reference)<1e-10)
            n=cfg['n'];h=cfg['length']/n;nf=6*n**3
            fields=raw[endpoint+'_background'][:nf].reshape(6,n,n,n)
            u=.08*fields[0];alpha=np.exp(u);z=np.exp(2*u)
            beta=.04*z*fields[1:4]/np.sqrt(1+.08**2*np.sum(fields[1:4]**2,axis=0))
            line=np.arange(n)*h-cfg['length']/2
            xyz=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1)
            index=0 if endpoint=='initial' else -1
            for i in range(15):
                q,p=states[index,i,:3],states[index,i,3:];mass=0 if i<9 else 1
                distance=abs((xyz-q)/h)
                basis=np.where(distance<1,(4-6*distance**2+3*distance**3)/6,np.where(distance<2,(2-distance)**3/6,0.))
                weights=np.prod(basis,axis=-1)
                checks.append(abs(weights.sum()-1)<1e-12)
                E=np.sqrt(mass*mass+z*np.dot(p,p))
                velocity=p*np.sum(weights*alpha*z/E)+np.sum(weights*beta,axis=(1,2,3))
                H=np.sum(weights*(alpha*E+np.einsum('jxyz,j->xyz',beta,p)))
                checks.append(np.max(abs(velocity-vel[index,i]))<1e-12 and abs(H-energies[index,i])<1e-12)
        event_positions=[];event_velocities=[];event_times=[]
        for i in range(9):
            hits=np.flatnonzero((states[:-1,i,0]<1.5)&(states[1:,i,0]>=1.5))
            if len(hits)==0:
                checks.append(row['events'][i] is None);continue
            j=hits[0];x=states[j:j+2,i,0]
            event_time=float(np.interp(1.5,x,t[j:j+2]))
            position=np.array([np.interp(event_time,t[j:j+2],states[j:j+2,i,k]) for k in range(3)])
            velocity=np.array([np.interp(event_time,t[j:j+2],vel[j:j+2,i,k]) for k in range(3)])
            energy=float(np.interp(event_time,t[j:j+2],energies[j:j+2,i]))
            recorded=row['events'][i]
            checks.append(abs(event_time-recorded['time'])<1e-12 and np.max(abs(position-recorded['position']))<1e-12 and np.max(abs(velocity-recorded['velocity']))<1e-12 and abs(energy-recorded['coordinate_energy'])<1e-12)
            event_positions.append(position);event_velocities.append(velocity);event_times.append(event_time)
        if len(event_positions)==9:
            matrices=[]
            for start,step in [(1,.05),(5,.025)]:
                matrices.append(np.array([(event_positions[start]-event_positions[start+1])[1:],(event_positions[start+2]-event_positions[start+3])[1:]]).T/(2*step))
            checks.append(np.max(abs(np.array(matrices)-row['transport']['matrices']))<1e-11)
            bend=float(np.arctan2(event_velocities[0][1],event_velocities[0][0]));delay=event_times[0]-3
            checks.append(abs(bend-row['central_bend'])<1e-12 and abs(delay-row['arrival_offset'])<1e-12)
            M=matrices[-1]
            checks.append(abs(1-np.trace(M)/2-row['transport']['trace_distortion'])<1e-11)
            checks.append(abs(1/abs(np.linalg.det(M))-row['transport']['area_gain'])<1e-10)
            derived[cfg['name']]=(bend,delay,M)
        checks.append(np.max(abs(states[-1,9:,:3]-row['final_body_positions']))<1e-12)
        trace=doc['background_trace']
        drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger'])
        background_pass=drift<1e-5 and row['clearance']>0 and max(t['edge_amplitude'] for t in trace)<1e-5 and max(t['cone_error'] for t in trace)<1e-10
        checks.append(abs(drift-row['energy_drift'])<1e-12 and background_pass==row['background_passed'])
        checks.append((len(event_positions)==9)==row['all_crossed'])
        if len(event_positions)==9:
            bundle_error=float(np.max(abs(matrices[0]-matrices[1])))
            checks.append(abs(bundle_error-row['transport']['bundle_error'])<1e-11)
            expected=background_pass and row['probe_cone_excess']<1e-10 and bundle_error<.001 and (row['background_replay_error'] is None or row['background_replay_error']<1e-12)
            checks.append(expected==row['passed'])
        else:checks.append(not row['passed'])
    available_comparisons=[]
    for emitter in ('control','Y'):
        names=[f'{emitter}-n{n}' for n in (32,48,64)]
        if not all(n in derived for n in names):continue
        coarse,middle,fine=[derived[n] for n in names]
        difference=abs(middle[0]-fine[0])/max(abs(fine[0]),1e-8)
        previous=abs(coarse[0]-middle[0])/max(abs(middle[0]),1e-8)
        matrix=float(np.max(abs(middle[2]-fine[2])))
        available_comparisons.append(dict(name=emitter+'-space',bend_relative_difference=difference,previous_difference=previous,matrix_difference=matrix,passed=difference<.01 and difference<previous and matrix<.001))
    if all(k in derived for k in ('Y-n48','control-n48','Y-n64','control-n64')):
        gain48=derived['Y-n48'][0]-derived['control-n48'][0];gain64=derived['Y-n64'][0]-derived['control-n64'][0]
        available_comparisons.append(dict(name='paired-effect',coarse_difference=gain48,fine_difference=gain64,passed=abs(gain64-gain48)<=max(.2*abs(gain64),1e-7)))
    if all(k in derived for k in ('Y-n32','Y-time')):
        base,fine=derived['Y-n32'],derived['Y-time']
        bend=abs(base[0]-fine[0])/max(abs(fine[0]),1e-8);delay=abs(base[1]-fine[1]);matrix=float(np.max(abs(base[2]-fine[2])))
        available_comparisons.append(dict(name='time',bend_relative_difference=bend,arrival_difference=delay,matrix_difference=matrix,passed=bend<.001 and delay<1e-5 and matrix<1e-4))
    complete=(directory/'summary.json').exists();campaign=None
    if complete:
        summary=json.loads((directory/'summary.json').read_text());checks.append(len(rows)==7 and rows==summary['runs'])
        checks.append(len(available_comparisons)==4)
        for expected,recorded in zip(available_comparisons,summary['comparisons']):
            checks.append(expected['name']==recorded['name'] and expected['passed']==recorded['passed'])
            checks.append(all(abs(expected[k]-recorded[k])<1e-9 for k in expected if k not in ('name','passed')))
        campaign=summary['passed'];checks.append(campaign==all(r['passed'] for r in rows+available_comparisons))
    result=dict(passed=bool(all(checks)),checks=len(checks),completed_runs=len(rows),declared_runs=7,complete=complete,campaign_passed=campaign,available_comparisons=available_comparisons)
    (ROOT/'point-bundle-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result));assert result['passed']


if __name__=='__main__':main()
