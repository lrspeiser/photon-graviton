"""Read-only provenance and independent interior-error reconstruction."""
import json
import hashlib
import subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]


def main():
    checks=[]
    def check(name,passed):checks.append(dict(name=name,passed=bool(passed)))
    for folder,matched in (('boundary-v1',False),('boundary-matched-v1',True)):
        out=HERE/folder;read=lambda name:json.loads((out/name).read_text())
        sha=lambda data:hashlib.sha256(data).hexdigest()
        for name,h in read('hashes.json').items():check(folder+' hash '+name,sha((out/name).read_bytes())==h)
        manifest=read('manifest.json')
        for name,h in manifest['hashes'].items():
            check(folder+' working '+name,sha((ROOT/name).read_bytes())==h)
            check(folder+' pinned '+name,sha(subprocess.check_output(['git','show',manifest['git_head']+':'+name],cwd=ROOT))==h)
        rows=read('runs.json');summary=read('summary.json');arrays=np.load(out/'profiles.npz')
        check(folder+' count',len(rows)==24)
        for row in rows:
            n=row['n'];dx=12/n;key=f'n{n}-w{row["width"]}-g{row["strength"]}'
            profile=arrays[key];x=profile[:,0];mask=abs(x)<3.5
            z=np.maximum(1-(x/.6)**2,0);initial_p=8*x/.6**2*z**3
            if matched:
                initial=dx*np.sum(initial_p**2)
                difference=.5*dx*np.sum((profile[:,3]**2+profile[:,4]**2)[mask])
            else:
                initial=.5*dx*np.sum(initial_p**2+((np.roll(z**4,-1)-z**4)/dx)**2)
                df=profile[:,3];dp=profile[:,4]
                difference=.5*dx*np.sum((dp**2+((np.roll(df,-1)-df)/dx)**2)[mask])
            error=np.sqrt(difference/initial)
            check(folder+' error '+key,abs(error-row['relative_interior_state_error'])<1e-12)
            passed=error<=.01 and row['ledger_drift']<=1e-5 and (not matched or row['opposite_characteristic_energy']<=1e-12)
            check(folder+' gate '+key,passed==row['passed'])
        check(folder+' summary count',sum(r['passed'] for r in rows)==summary['passes'])
    result=dict(passed=all(c['passed'] for c in checks),checks=checks,
                scope='Archive and metric audit passes are distinct from failed OB-1 absorption gates.')
    (HERE/'boundary-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(dict(passed=result['passed'],checks=len(checks))))
    if not result['passed']:raise AssertionError('Boundary archive audit failed')


if __name__=='__main__':main()
