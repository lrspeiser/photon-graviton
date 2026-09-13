"""Audit one atomic refinement snapshot without assuming its process is live."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.special import ndtr

H=Path(__file__).resolve().parent;R=H.parents[2]
sys.path.insert(0,str(H.parent/'timing-continuous-scatter'))
from integration import averaged_likelihood
cache=R/'research_work/generated/timing-fresh-refinement'
raw=(cache/'checkpoint.json').read_bytes();state=json.loads(raw)
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
for path,digest in state['hashes'].items():assert sha(R/path)==digest
base=json.loads((H/'results.json').read_text());lookup={r['label']:r for r in base['cases']}
expected=[]
for family in ['split_gaussian','split_gaussian_shoulder']:
    for snr in [5,20]:
        for truth in [0,1]:
            eligible=[r for r in base['cases'] if (r['family'],r['snr'],r['truth_b'])==(family,snr,truth)]
            expected.append(max(eligible,key=lambda r:(abs(r['fit']['best']['parameters'][1]-truth),r['label']))['label'])
assert state['selection']==expected
assert len({r['label'] for r in state['cases']})==len(state['cases'])
checked=[]
for row in state['cases']:
    label=row['label'];assert label in expected
    file=cache/(label+'.npz');assert sha(file)==row['array_sha256']
    oldfile=R/'research_work/generated/timing-fresh-calibration'/(label+'.npz')
    assert sha(oldfile)==lookup[label]['array_sha256']
    with np.load(file) as data,np.load(oldfile) as old:
        x=data['log_width'];z=data['redshift'];logs=data['event_log_likelihoods']
        assert np.array_equal(x,old['log_width']) and np.array_equal(z,old['redshift'])
        assert logs.shape==(98,321) and np.all(np.isfinite(logs))
        L=np.exp(logs-logs.max(axis=1)[:,None])
        a,b,s=row['fit']['best']['parameters'];mu=a+b*np.log1p(z)
        assert np.log(5)<=a<=np.log(100) and -2<=b<=3 and 0<=s<=.6
        outside=((mu<x[0])|(mu>x[-1])).astype(float) if s==0 else ndtr((x[0]-mu)/s)+ndtr((mu-x[-1])/s)
        assert max(outside)<=1.00001e-6 and row['fit']['best']['success']
        value=averaged_likelihood(x,L,mu,s);assert np.all(value>0)
        objective=float(-np.log(value).sum())
        objective_error=abs(objective-row['fit']['best']['objective']);assert objective_error<1e-8
        D=abs((logs-logs.max(axis=1)[:,None])-(old['event_log_likelihoods']-old['event_log_likelihoods'].max(axis=1)[:,None]))
        metric=float(np.mean(averaged_likelihood(x,L*D,mu,s)/value))
        delta=abs(b-lookup[label]['fit']['best']['parameters'][1])
        assert abs(metric-row['posterior_weighted_mean_event_curve_change'])<1e-10
        assert abs(delta-row['absolute_b_change'])<1e-10
        assert row['numerical_gate_pass']==bool(delta<=.05 and metric<=.1)
        checked.append(dict(label=label,objective_recomputation_error=objective_error,
            delta_b=delta,curve_change=metric,gate_pass=row['numerical_gate_pass']))
result=dict(snapshot_sha256=hashlib.sha256(raw).hexdigest(),completed=len(checked),expected=8,
    complete=len(checked)==8,checked=checked,
    limitations='Uses the shared population integrator; does not prove its independent correctness, global optimization, coverage, width-grid convergence or physical validity. Snapshot is not process-liveness evidence.')
destination=H/'refinement-audit.json' if result['complete'] else cache/'audit-snapshot.json'
destination.write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
