"""Verify one immutable checkpoint snapshot; summarize only a complete batch."""
from pathlib import Path
import hashlib, itertools, json, sys
import numpy as np
from scipy.stats import beta, chi2
from scipy.special import ndtr
H=Path(__file__).resolve().parent;R=H.parents[2]
sys.path.insert(0,str(H.parent/'timing-continuous-scatter'))
from integration import averaged_likelihood

path=R/'research_work/generated/timing-fresh-calibration/checkpoint.json'
raw=path.read_bytes();state=json.loads(raw)
for name,digest in state['hashes'].items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
p=json.loads((H/'protocol.json').read_text())
generation=json.loads((H.parent/'timing-coverage-calibration/protocol.json').read_text())['calibration']
expected={(seed,family,snr,b) for seed,family,snr,b in itertools.product(p['seeds'],generation['families'],generation['nominal_snr'],generation['truth_b'])}
found=set();max_error=0.
for row in state['cases']:
    key=tuple(row[k] for k in ['seed','family','snr','truth_b'])
    assert key in expected and key not in found
    found.add(key)
    a_path=path.parent/(row['label']+'.npz')
    assert hashlib.sha256(a_path.read_bytes()).hexdigest()==row['array_sha256']
    with np.load(a_path) as d:
        x=d['log_width'];logs=d['event_log_likelihoods'];z=d['redshift']
    assert logs.shape==(98,321) and z.shape==(98,) and np.all(np.diff(x)>0)
    assert np.isfinite(logs).all() and np.isfinite(z).all() and np.all(z>-1)
    L=np.exp(logs-logs.max(axis=1)[:,None]);lz=np.log1p(z)
    best=row['fit']['best'];a,b,s=best['parameters']
    assert best['success'] and np.log(5)<=a<=np.log(100) and -2<=b<=3 and 0<=s<=.6
    def support(mu,sigma):
        return float(np.max((mu<x[0])|(mu>x[-1]))) if sigma==0 else float(np.max(ndtr((x[0]-mu)/sigma)+ndtr((mu-x[-1])/sigma)))
    assert abs(support(a+b*lz,s)-best['maximum_outside'])<1e-12
    assert support(a+b*lz,s)<=1.00001e-6
    recomputed=float(-np.log(averaged_likelihood(x,L,a+b*lz,s)).sum())
    delta=abs(recomputed-best['objective']);max_error=max(max_error,delta)
    assert delta<1e-8
    good=[]
    for t in row['truth_profiles']:
        aa,ss=t['parameters']
        assert abs(support(aa+row['truth_b']*lz,ss)-t['outside'])<1e-12
        if t['success'] and t['outside']<=1.00001e-6 and t['objective']<1e11:
            aa,ss=t['parameters']
            value=float(-np.log(averaged_likelihood(x,L,aa+row['truth_b']*lz,ss)).sum())
            assert abs(value-t['objective'])<1e-8
            good.append(value)
    lr=2*(min(good)-recomputed) if good else None
    assert (lr is None)==(row['raw_lr'] is None)
    if lr is not None:assert abs(lr-row['raw_lr'])<4e-8
    valid=bool(good) and lr>=-2e-6 and min(a-np.log(5),np.log(100)-a,b+2,3-b,.6-s)>1e-5
    assert bool(valid)==row['valid']
    assert row['accepted_95']==bool(valid and max(0,lr)<=chi2.ppf(.95,1))
complete=found==expected
audit=dict(checkpoint_sha256=hashlib.sha256(raw).hexdigest(),verified_cases=len(found),expected_cases=len(expected),complete=complete,
           source_and_array_hashes_verified=True,objectives_recomputed=True,maximum_objective_discrepancy=max_error,
           physical_validation=False)
print(json.dumps(audit))
if not complete:
    # A bounded snapshot is evidence of integrity, not proof the process is live.
    (H/'audit-snapshot.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8',newline='\n')
    sys.exit(0)
cells=[]
for family,snr,b in itertools.product(generation['families'],generation['nominal_snr'],generation['truth_b']):
    rows=[v for v in state['cases'] if (v['family'],v['snr'],v['truth_b'])==(family,snr,b)]
    n=len(rows);assert n==20
    errors=np.array([v['fit']['best']['parameters'][1]-b for v in rows]);k=sum(v['accepted_95'] for v in rows);invalid=sum(not v['valid'] for v in rows)
    cells.append(dict(family=family,snr=snr,truth_b=b,n=n,mean_bias=float(errors.mean()),monte_carlo_se=float(errors.std(ddof=1)/np.sqrt(n)),accepted=k,invalid=invalid,
       binomial_95=[float(beta.ppf(.025,k,n-k+1)) if k else 0.,float(beta.ppf(.975,k+1,n-k)) if k<n else 1.],
       screen_pass=bool(abs(errors.mean())<=p['screens']['maximum_absolute_mean_bias'] and k/n>=p['screens']['minimum_95_acceptance_fraction'] and invalid<=p['screens']['maximum_invalid_per_cell'])))
(H/'summary.json').write_text(json.dumps(dict(audit=audit,cells=cells),indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(cells,indent=2))
