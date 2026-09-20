"""Read-only source, archive and scaling audit for CC-2S."""
import json
import hashlib
import subprocess
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];OUT=HERE/'slow-motion-v1'


def main():
    read=lambda name:json.loads((OUT/name).read_text())
    sha=lambda data:hashlib.sha256(data).hexdigest()
    checks=[]
    def check(name,passed):checks.append(dict(name=name,passed=bool(passed)))
    for name,h in read('hashes.json').items():check('archive '+name,sha((OUT/name).read_bytes())==h)
    manifest=read('manifest.json')
    for name,h in manifest['hashes'].items():
        check('working '+name,sha((ROOT/name).read_bytes())==h)
        check('pinned '+name,sha(subprocess.check_output(['git','show',manifest['git_head']+':'+name],cwd=ROOT))==h)
    rows=read('cases.json');ref=read('refinements.json');quad=read('quadratures.json')
    check('96 cases and refinements',len(rows)==len(ref)==96)
    check('refinement gates',all(r['relative_error']<=1e-6 for r in ref))
    check('four quadratures',len(quad)==4 and all(q['relative_error']<=1e-7 for q in quad))
    for row in rows:
        base=next(r for r in rows if r['speed']==.2 and r['sense']==1 and r['coupling_ratio']==1 and r['radius_ratio']==row['radius_ratio'])
        star=base['star_fraction']*(row['speed']/.2)**2*row['sense']*row['coupling_ratio']**2
        light=base['photon_fraction']*(row['speed']/.2)*row['sense']*row['coupling_ratio']**2
        check('stellar scaling '+str(len(checks)),abs(star-row['star_fraction'])/max(abs(star),1e-30)<1e-12)
        check('photon scaling '+str(len(checks)),abs(light-row['photon_fraction'])/max(abs(light),1e-30)<1e-12)
    result=dict(passed=all(c['passed'] for c in checks),checks=checks)
    (HERE/'slow-motion-audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(dict(passed=result['passed'],checks=len(checks))))
    if not result['passed']:raise AssertionError('CC-2S archive/scaling audit failed')


if __name__=='__main__':main()
