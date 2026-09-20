"""Independent reduced phase-mode ledger reconstruction and figure."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import hashlib
import subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];OUT=HERE/'phase-budget-v1'


def main():
    read=lambda name:json.loads((OUT/name).read_text())
    sha=lambda value:hashlib.sha256(value).hexdigest()
    checks=[]
    def check(name,passed):checks.append(dict(name=name,passed=bool(passed)))
    for name,h in read('hashes.json').items():check('archive '+name,sha((OUT/name).read_bytes())==h)
    manifest=read('manifest.json')
    for name,h in manifest['hashes'].items():
        check('working '+name,sha((ROOT/name).read_bytes())==h)
        check('pinned '+name,sha(subprocess.check_output(['git','show',manifest['git_head']+':'+name],cwd=ROOT))==h)
    rows=read('runs.json');branches=read('branches.json')
    check('16 unique trajectories',len(rows)==len({r['id'] for r in rows})==16)
    for row in rows:
        a=np.load(OUT/(row['id']+'.npz'));b,p,q=a['state'].T;t=row['ratio']
        terms=np.column_stack((.5*p*p*np.exp(-2*b*b),.5*b*b,t*np.exp(-.5*b*b)))
        error=np.max(abs(terms-a['energy']))
        check('independent terms '+row['id'],error<1e-12)
        drift=np.max(abs(terms.sum(axis=1)+q-row['initial_energy']))/max(1,abs(row['initial_energy']))
        check('ledger '+row['id'],abs(drift-row['energy_ledger_drift'])<1e-12 and drift<=1e-6)
        check('outlet monotonic '+row['id'],np.min(np.diff(q))>=-1e-12)
    for branch in branches:
        check('minimum energy split '+str(branch['ratio']),abs(branch['matter_fraction']+branch['field_potential_fraction']+branch['releasable_fraction']-1)<1e-12)
    summary=dict(passed=all(c['passed'] for c in checks),checks=checks)
    (HERE/'phase-budget-audit.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf8',newline='\n')
    if not summary['passed']:raise AssertionError('PF-0 independent audit failed')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for row in rows:
        if row['damping']!=.05 or row['ratio']<=1:continue
        a=np.load(OUT/(row['id']+'.npz'))
        axes[0].loglog(a['times'][1:],np.maximum(abs(a['state'][1:,0]),1e-12),label=f't={row["ratio"]:.9g}')
    axes[0].set(title='Seeded local amplitude with explicit outlet',xlabel='Dimensionless time',ylabel='|B|')
    axes[0].legend(fontsize=8)
    t=np.geomspace(1+1e-8,10,500)
    axes[1].semilogx(t,1/t,label='Matter rest-energy fraction')
    axes[1].semilogx(t,np.log(t)/t,label='Field potential fraction')
    axes[1].semilogx(t,1-(1+np.log(t))/t,label='Available kinetic/outlet fraction')
    axes[1].set(title='Stationary-branch energy budget',xlabel='Source-to-stiffness parameter t',ylabel='Fraction of zero-field matter energy',ylim=(-.02,1.02))
    axes[1].legend(fontsize=8)
    fig.suptitle('PF-0: reduced feedback diagnostic — not a 3D swirl or astronomical fit')
    fig.savefig(HERE/'phase-budget.png',dpi=170);plt.close(fig)
    print(json.dumps(dict(passed=summary['passed'],checks=len(checks))))


if __name__=='__main__':main()
