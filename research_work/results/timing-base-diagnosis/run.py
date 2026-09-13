from pathlib import Path
import json,hashlib
import numpy as np
from scipy.stats import chi2
H=Path(__file__).resolve().parent;R=H.parents[2];p=R/'research_work/generated/timing-coverage-calibration/checkpoint.json'
r=json.loads(p.read_text(encoding='utf-8'));cases=list(r['cases'].values());assert len(cases)==160
for name,expected in r['hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==expected
# Immutable base-only content, independent of refinements still being written.
base=json.dumps(r['cases'],sort_keys=True,separators=(',',':')).encode();sha=hashlib.sha256(base).hexdigest()
if (H/'results.json').exists():assert json.loads((H/'results.json').read_text(encoding='utf-8'))['base_cases_sha256']==sha
protocol=json.loads((H.parent/'timing-coverage-calibration/protocol.json').read_text(encoding='utf-8'));gate=protocol['calibration']['diagnostic_thresholds'];out=[]
for key in sorted(set((v['family'],v['snr'],v['truth_b']) for v in cases)):
 vs=[v for v in cases if (v['family'],v['snr'],v['truth_b'])==key];bad=[v for v in vs if not v['valid_fit']];bias=float(np.mean([v['fit']['b']-v['truth_b'] for v in vs]));accepted=sum(v['coverage']['0.95'] for v in vs);lr=sum(v['profile_lr_at_truth'] is not None and v['profile_lr_at_truth']<=chi2.ppf(.95,1) for v in vs)
 out.append(dict(family=key[0],snr=key[1],truth_b=key[2],n=len(vs),mean_bias=bias,invalid=len(bad),accepted_truth=accepted,lr_only_truth=int(lr),screen_pass=abs(bias)<=gate['maximum_absolute_mean_bias'] and accepted/len(vs)>=gate['minimum_nominal_95_coverage'] and len(bad)<=gate['maximum_optimizer_or_boundary_failures']))
invalid=[dict(label=v['label'],b=v['fit']['b'],sigma=v['fit']['sigma'],optimizer_success=v['fit']['optimizer_success'],interior_solution=v['fit']['interior_solution'],outside_mass=v['fit']['maximum_population_mass_outside_width_bounds'],profile_lr_at_truth=v['profile_lr_at_truth']) for v in cases if not v['valid_fit']]
assert all(abs(x['sigma']-.03)<1e-8 and x['optimizer_success'] and not x['interior_solution'] for x in invalid)
result=dict(scope='Completed frozen base-case diagnosis only; final numerical refinements excluded',base_cases_sha256=sha,source_hashes=r['hashes'],cells=out,invalid_cases=invalid,invalid_count=len(invalid),all_invalid_at_scatter_floor=True,all_cells_pass=all(x['screen_pass'] for x in out))
(H/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(out))
