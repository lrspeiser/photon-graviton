"""Check current result hashes and preserve the earlier failure as history."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
checked=[]
for name in ['results.json','verification.json','verification-unsplit-results.json']:
    data=json.loads((HERE/name).read_text())
    for rel,expected in data['input_hashes'].items():
        actual=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
        assert actual==expected,(name,rel)
        checked.append(rel)
old=json.loads((HERE/'initial-unsplit-results.json').read_text())
new=json.loads((HERE/'verification-unsplit-results.json').read_text())
error=float(np.max(abs(np.array([r['fine_extra_acceleration'] for r in old['rows']])-
                       np.array([r['fine_extra_acceleration'] for r in new['rows']]))))
assert error<1e-10
# Additional imported dependencies and empirical-fit source are fingerprinted
# explicitly, beyond the per-run generated-cache hashes.
dependencies=[HERE/'run.py',HERE/'check.py',HERE/'report.py',Path(__file__),
    HERE.parent/'rotating-bar-orbits/fast_multipole.py',
    HERE.parent/'bar-field-foundation/field.py',
    HERE.parent/'joint-galaxy-audit/results.json',
    HERE.parent/'conservative-field-completion/run.py']
out=dict(current_input_hashes_verified=len(checked),
         original_unsplit_force_reproduction_max_absolute=error,
         original_unsplit_provenance='Historical implementation hashes retained; current --unsplit reproduces its forces.',
         dependency_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies})
(HERE/'integrity.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))
