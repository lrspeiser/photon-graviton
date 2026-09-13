"""Free-scatter refits on independently evaluated finer event grids."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[2]
sys.path.insert(0,str(H.parent/'timing-boundary-capable'))
from run import fit
source=H.parent/'timing-scatter-grid-refinement/results.json';ref=json.loads(source.read_text(encoding='utf-8'))
oldpath=H.parent/'timing-boundary-capable/results.json';old={v['label']:v for v in json.loads(oldpath.read_text(encoding='utf-8'))['cases']}
rows=[]
for c in ref['cases']:
 coarse=R/'research_work/generated/timing-coverage-calibration'/(c['label']+'.npz')
 fine=R/'research_work/generated/timing-scatter-grid-refinement'/(c['label']+'.npz')
 assert hashlib.sha256(coarse.read_bytes()).hexdigest()==c['input_sha256']
 assert hashlib.sha256(fine.read_bytes()).hexdigest()==c['refined_array_sha256']
 with np.load(coarse) as a,np.load(fine) as b:
  assert np.array_equal(a['log_width'],b['log_width'][::2])
  assert np.array_equal(a['event_log_likelihoods'],b['event_log_likelihoods'][:,::2])
  assert np.array_equal(a['redshift'],b['redshift']) and np.array_equal(a['ids'],b['ids'])
  result=fit(b['log_width'],b['event_log_likelihoods'],b['redshift'],old[c['label']]['original'])
 q0=old[c['label']]['revised']['best']['parameters'];q1=result['best']['parameters']
 rows.append(dict(label=c['label'],coarse_sha256=c['input_sha256'],fine_sha256=c['refined_array_sha256'],coarse_fit=old[c['label']]['revised'],fine_fit=result,b_change=q1[1]-q0[1],sigma_change=q1[2]-q0[2]))
 print(json.dumps(dict(label=c['label'],b_change=q1[1]-q0[1],sigma_change=q1[2]-q0[2],fine_best=result['best'])),flush=True)
paths=[Path(__file__),H.parent/'timing-boundary-capable/run.py',H.parent/'timing-continuous-scatter/integration.py',source,oldpath]
out=dict(scope='Two exposed artificial samples; free-scatter width-grid check, not population calibration',hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},cases=rows)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
