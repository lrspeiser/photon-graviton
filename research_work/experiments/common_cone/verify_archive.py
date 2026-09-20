"""Read-only CC-1 provenance and independently reconstructed summary checks."""
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=HERE/'evidence-v1'


def main():
    read=lambda name:json.loads((OUT/name).read_text(encoding='utf8'))
    sha=lambda data:hashlib.sha256(data).hexdigest()
    checks=[]
    def check(name, passed):
        checks.append(dict(name=name,passed=bool(passed)))
    for name,expected in read('hashes.json').items():
        check('hash '+name,sha((OUT/name).read_bytes())==expected)
    manifest=read('manifest.json')
    for name,expected in manifest['hashes'].items():
        committed=subprocess.check_output(['git','show',manifest['git_head']+':'+name],cwd=ROOT)
        check('pinned '+name,sha(committed)==expected)
        current=(ROOT/name).read_bytes()
        if name.endswith('/run_audit.py'):
            # Only post-save stdout formatting was fixed after the first run.
            old=b'print(json.dumps(summary,indent=2),flush=True)'
            new=b"print((out/'summary.json').read_text(encoding='utf8'),flush=True)"
            check('runner output-only correction',current==committed.replace(old,new))
        else: check('unchanged '+name,current==committed)
    rows=read('local.json');sources=read('sources.json');summary=read('summary.json')
    check('600 local samples',len(rows)==600)
    check('108 source samples',len(sources['cases'])==108)
    for row in rows:
        check(f'local gates {row["index"]}',row['passed']==all(row['residuals'][k]<=row['limits'][k] for k in row['residuals']))
        f=np.array(row['F']);p=np.array(row['p']);g=row['g'];eta=row['eta'];kappa=row['kappa']
        alpha=np.exp(g*f[0]);beta=alpha*kappa*eta*f[1:]/np.sqrt(1+eta**2*np.dot(f[1:],f[1:]))
        e=np.sqrt(row['mass']**2+p@p)
        check(f'independent energy {row["index"]}',abs(alpha*e+beta@p-row['energy'])<1e-12)
        check(f'independent kinetic bound {row["index"]}',abs(alpha-np.linalg.norm(beta)-row['minimum_field_kinetic_eigenvalue'])<1e-12)
    check('local pass count',sum(r['passed'] for r in rows)==summary['local_passes']==600)
    check('source pass count',sum(r['passed'] for r in sources['cases'])==summary['source_passes']==108)
    for r in sources['rotation']:
        values=[v['sample'] for v in sources['cases'] if v['n']==r['n'] and v['radius']==r['radius']]
        spread=max(abs(v-values[0]) for v in values)/abs(values[0])
        check(f'rotation spread {r["n"]}/{r["radius"]}',abs(spread-r['rotation_spread'])<1e-14)
    result=dict(passed=all(c['passed'] for c in checks),checks=checks,
                first_run_exit='Nonzero after all evidence was saved: NumPy int64 in stdout JSON. Scientific artifacts intact; printing fixed separately.')
    print(json.dumps(dict(passed=result['passed'],checks=len(checks),first_run_exit=result['first_run_exit']),indent=2))
    if not result['passed']: raise AssertionError([c for c in checks if not c['passed']])


if __name__=='__main__': main()
