"""Read-only weak-field evidence and physical-target audit."""
import json
import hashlib
import subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];OUT=HERE/'weak-field-v1'


def main():
    read=lambda name:json.loads((OUT/name).read_text())
    sha=lambda raw:hashlib.sha256(raw).hexdigest()
    checks=[]
    def check(name,passed):checks.append(dict(name=name,passed=bool(passed)))
    for name,h in read('hashes.json').items():check('evidence '+name,sha((OUT/name).read_bytes())==h)
    manifest=read('manifest.json')
    for name,h in manifest['hashes'].items():
        check('working '+name,sha((ROOT/name).read_bytes())==h)
        check('pinned '+name,sha(subprocess.check_output(['git','show',manifest['git_head']+':'+name],cwd=ROOT))==h)
    orbits=read('orbits.json');lenses=read('lensing.json');summary=read('summary.json')
    check('orbit count',len(orbits)==404);check('lens count',len(lenses)==48)
    check('zero flat targets',sum(abs(r['log_speed_slope'])<=.1 for r in orbits)==summary['flat_target_passes']==0)
    check('numerical orbit gates',all(max(r['numerical_errors'])<=1e-7 for r in orbits))
    check('numerical lens gates',all(r['relative_error']<=1e-7 for r in lenses))
    check('unscreened light factor',all(abs(r['ratio_to_unscreened_spatial_curvature_benchmark']-.5)<1e-12 for r in lenses if r['omega']==0))
    check('stability threshold',all(r['stable_circular_orbit']==(r['omega']*r['r']<(1+np.sqrt(5))/2) for r in orbits))
    result=dict(passed=all(c['passed'] for c in checks),checks=checks,scope='Evidence integrity and calculation verification; physical flat-curve target remains failed.')
    (HERE/'weak-field-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(dict(passed=result['passed'],checks=len(checks))))
    if not result['passed']:raise AssertionError('Weak-field archive audit failed')


if __name__=='__main__':main()
