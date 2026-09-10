"""Verify the published audit artifacts and their declared scope."""
from pathlib import Path
import csv
import hashlib
import json
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
manifest=json.loads((HERE/'input-manifest.json').read_text())
for rel,expected in manifest['inputs'].items():
    path=HERE/rel if rel in ['protocol.json','observations.json'] else ROOT/rel
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expected, rel
assert hashlib.sha256((HERE/'run.py').read_bytes()).hexdigest()==manifest['script_sha256']
result=json.loads((HERE/'results.json').read_text())
assert result['exact_time_candidate']['amendment_sha256']==hashlib.sha256((HERE/'time-candidate-amendment.md').read_bytes()).hexdigest()
counts={}
for name,expected in [('redshift-predictions.csv',170),('spectral-aging-predictions.csv',35),
                       ('radio-chromaticity.csv',5),('firas-predictions.csv',430),
                       ('seven-band-predictions.csv',70),('prescribed-time-field.csv',48),
                       ('exact-time-candidate.csv',12)]:
    with (HERE/name).open() as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==expected, name
    counts[name]=len(rows)
coverage=(HERE/'coverage.md').read_text(encoding='utf-8')
ids=re.findall(r'^\| (R\d\d) \|',coverage,re.M)
assert ids==[f'R{i:02d}' for i in range(1,33)]
docs=[HERE/'report.md',HERE/'coverage.md',ROOT/'research_plan/electromagnetic-transfer-specification.md']
links=[]
for path in docs:
    for target in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
        if target.startswith(('https://','http://','#')):
            continue
        resolved=(path.parent/target.split('#')[0]).resolve()
        assert resolved.is_relative_to(ROOT) and resolved.exists(), (path,target)
        links.append(str(resolved.relative_to(ROOT)))
contract=json.loads((ROOT/'research_plan/universe-contract.json').read_text())
assert contract['active_branch']['specification']=='research_plan/electromagnetic-transfer-specification.md'
assert contract['active_branch']['additional_time_field_required'] is None
verification=dict(status='Numerical/artifact checks passed; physical theory and observational coverage incomplete',
                  source_hash_checks=len(manifest['inputs']),row_counts=counts,requirements=len(ids),
                  verified_local_links=len(links),figure_visually_inspected=True,
                  preserved_paper_pdf_sha256=hashlib.sha256((ROOT/'papers/cumulative-time-companions/manuscript.pdf').read_bytes()).hexdigest(),
                  unresolved=result['exact_time_candidate']['failures'])
(HERE/'verification.json').write_text(json.dumps(verification,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(verification,indent=2))
