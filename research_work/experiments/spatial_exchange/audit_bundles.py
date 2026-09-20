"""Independent detector reconstruction from archived every-step probe traces."""
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_emitter import reconstruct

ROOT=Path(__file__).resolve().parent


def main():
    directory=ROOT/'bundle-v1';manifest=json.loads((directory/'manifest.json').read_text())
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
    complete=(directory/'summary.json').exists();campaign=None
    if complete:
        summary=json.loads((directory/'summary.json').read_text());checks.append(len(rows)==7 and rows==summary['runs'])
        base=derived['Y']
        for c in summary['comparisons']:
            fine=derived[c['name']];time_case=c['name']=='time'
            bend=abs(base[0]-fine[0])/max(abs(fine[0]),1e-8);delay=abs(base[1]-fine[1]);matrix=float(np.max(abs(base[2]-fine[2])))
            expected=bend<(.01 if time_case else .05) and delay<(1e-4 if time_case else .001) and matrix<(.001 if time_case else .005)
            checks.append(abs(bend-c['bend_relative_difference'])<1e-10 and abs(delay-c['arrival_difference'])<1e-12 and abs(matrix-c['matrix_difference'])<1e-11 and expected==c['passed'])
        campaign=summary['passed'];checks.append(campaign==all(r['passed'] for r in rows+summary['comparisons']))
    result=dict(passed=bool(all(checks)),checks=len(checks),completed_runs=len(rows),declared_runs=7,complete=complete,campaign_passed=campaign)
    (ROOT/'bundle-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result));assert result['passed']


if __name__=='__main__':main()
