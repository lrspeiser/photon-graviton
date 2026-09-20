"""SE-R independent energy/source checks and explicitly partial comparisons."""
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_emitter import reconstruct
from total_source import analyze

ROOT=Path(__file__).resolve().parent


def main():
    directory=ROOT/'source-refinement-v1'
    manifest=json.loads((directory/'manifest.json').read_text());checks=[];states={};completed=0
    for name,digest in manifest['sources'].items():
        checks.append(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    for cfg in manifest['configurations']:
        path=directory/(cfg['name']+'.json')
        if not path.exists():continue
        completed+=1;doc=json.loads(path.read_text());row=doc['summary'];trace=doc['trace']
        raw=np.load(directory/(cfg['name']+'.npz'))
        checks.append(row['config']==cfg and trace[0]['time']==0 and trace[-1]['time']==4)
        for endpoint in ['initial','final']:
            independent=reconstruct(raw[endpoint],cfg)
            checks.append(max(abs(independent[k]-row[endpoint][k]) for k in ['field','matter','ledger','internal_rest','radiation_propagation','companion_propagation','mixing_potential'])<1e-10)
        drift=max(abs(t['ledger']-trace[0]['ledger']) for t in trace)/abs(trace[0]['ledger'])
        clearance=cfg['length']/2-1-row['maximum_source_extent']-4*row['maximum_characteristic']
        passed=drift<1e-5 and clearance>0 and max(t['cone_error'] for t in trace)<1e-10 and max(t['edge_amplitude'] for t in trace)<1e-5
        checks.append(passed==row['passed']);checks.append(passed)
        states[cfg['name']]=analyze(directory,cfg);checks.append(states[cfg['name']]['passed'])
    reused={'emitter-mixed-chi200':'mixed-coarse','emitter-Y-chi200':'Y-coarse',
            'time-refinement':'Y-time','space-refinement':'Y-n40','rotation':'Y-rotation'}
    original=ROOT/'emitter-v1'
    for cfg in json.loads((original/'manifest.json').read_text())['configurations']:
        if cfg['name'] in reused and (original/(cfg['name']+'.json')).exists():
            key=reused[cfg['name']];states[key]=analyze(original,cfg);checks.append(states[key]['passed'])
    comparisons=[]
    def difference(left,right,key):
        return abs(states[left][key]-states[right][key])/max(abs(states[right][key]),1e-12)
    for label,a,b,wlimit,tlimit,preceding in [
        ('mixed-space','mixed-n48','mixed-n56',.05,.001,('mixed-n40','mixed-n48')),
        ('Y-space','Y-n48','Y-n56',.05,.001,('Y-n40','Y-n48')),
        ('mixed-time','mixed-coarse','mixed-n32',.001,.0001,None),
        ('Y-time','Y-coarse','Y-time',.001,.0001,None),
        ('mixed-rotation','mixed-n40-rotation','mixed-n40',.02,.001,None),
        ('Y-rotation','Y-rotation','Y-coarse',.02,.001,None)]:
        needed=[a,b]+(list(preceding) if preceding else [])
        if not all(k in states for k in needed):continue
        wave=difference(a,b,'wave_source');total=difference(a,b,'total_coupling_source')
        previous=difference(*preceding,'wave_source') if preceding else None
        passed=wave<wlimit and total<tlimit and (previous is None or wave<previous)
        comparisons.append(dict(name=label,wave_relative_difference=wave,total_relative_difference=total,previous_wave_difference=previous,passed=passed))
    complete=completed==7 and len(comparisons)==6
    result=dict(audit_passed=bool(all(checks)),checks=len(checks),completed_new_runs=completed,declared_new_runs=7,
                complete=complete,source_accuracy_passed=all(c['passed'] for c in comparisons) if complete else None,
                comparisons=comparisons,states=states,
                analysis_sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['audit_source_refinement.py','total_source.py','audit_emitter.py']})
    (ROOT/'source-refinement-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['states','analysis_sources']}))
    assert result['audit_passed']


if __name__=='__main__':main()
