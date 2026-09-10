"""Test extra-force resolution on all saved candidate-library positions."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CODE=HERE.parent/'full-bar-completion/run.py'
spec=importlib.util.spec_from_file_location('full_field',CODE)
full=importlib.util.module_from_spec(spec);spec.loader.exec_module(full)
path=ROOT/'research_work/data-cache/training-orbit-launches/full-field-orbits.npz'
f=np.load(path);states=f['trajectories'];ids=f['source_id']
axispath=full.CACHE/'axisymmetric-reference-80-L128-N512.npz'
coarsepath=full.CACHE/'fine.npz';finepath=full.CACHE/'finer.npz'
axis=full.load_field(axispath);coarse=full.load_field(coarsepath);fine=full.load_field(finepath)
rows=[]
for j in range(states.shape[1]):
    points=states[:,j,:3]
    valid=np.isfinite(points).all(axis=1)
    points=points[valid]
    errors=[]
    for start in range(0,len(points),32):
        x=points[start:start+32]
        a=axis.evaluate(x)[1];c=coarse.evaluate(x)[1];b=fine.evaluate(x)[1]
        err=np.linalg.norm(b-c,axis=1)/np.linalg.norm(a+b,axis=1)
        assert np.isfinite(err).all()
        errors.extend(err.tolist())
    if errors:
        index=int(np.argmax(errors))
        rows.append(dict(seed_index=j,source_id=str(ids[j]),positions=len(points),
                         maximum_fraction=float(max(errors)),median_fraction=float(np.median(errors)),
                         worst_position_kpc=points[index].tolist(),field_gate_pass=bool(max(errors)<.01)))
    else:rows.append(dict(seed_index=j,source_id=str(ids[j]),positions=0,field_gate_pass=False,
                         reason='No complete trajectory available; not evaluated'))
out=dict(scope='Saved-path extra-force refinement only; not continuum or stellar-likelihood validation.',
         rows=rows,total_positions=sum(r['positions'] for r in rows),
         passing_seeds=sum(r['field_gate_pass'] for r in rows),
         threshold_fraction=.01,
         input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       [Path(__file__),CODE,path,axispath,coarsepath,finepath]},
         holdouts_opened=False)
(HERE/'field-check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print('Field check passed',out['passing_seeds'],'of',len(rows),'at',out['total_positions'],'positions',flush=True)
