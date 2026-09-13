"""Audit original candidate lists without treating ranked lists as parent pools."""
from pathlib import Path
import hashlib
import json
import re
import requests
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache/selection-fields'
previous = ROOT / 'research_work/results/selection-plate-history/results.json'
plates = json.loads(previous.read_text())['plates']
manifest, records, counts = [], [], []
for plate in plates:
    holes = CACHE / f"plateHolesSorted-{plate['plate_id']:06d}.par"
    text = holes.read_text()
    inputs = re.findall(r'^plateinput\d+ (\S+_SCI_\d+\.par)',text,re.M)
    assert len(inputs) == 1
    url = 'https://svn.sdss.org/public/data/sdss/platelist/trunk/inputs/' + inputs[0]
    path = CACHE / inputs[0].split('/')[-1]
    if not path.exists():
        response = requests.get(url,timeout=(15,45))
        response.raise_for_status()
        assert response.text.startswith('targettype SCIENCE')
        path.write_bytes(response.content)
    text = path.read_text()
    assert int(re.search(r'^designid (\d+)',text,re.M)[1]) == plate['design_id']
    body, name = re.search(r'typedef struct \{(.*?)\}\s*(\w+);',text,re.S).groups()
    columns = re.findall(r'\b(?:char|double|int|float|long)\s+(\w+)(?:\[\d+\])?;',body)
    n = 0
    for line in text.splitlines():
        if not line.startswith(name+' '):
            continue
        values = re.findall(r'"[^"]*"|\{[^}]*\}|[^\s]+',line)[1:]
        assert len(values) == len(columns)
        row = dict(zip(columns,[v.strip('"').strip() for v in values]))
        keep = ['targetids','priority','sourcetype','apogee2_target1','apogee2_target2','apogee2_target3','qa','selbits']
        records.append(dict(design_id=plate['design_id'], **{k:row[k] for k in keep}))
        n += 1
    counts.append(dict(design_id=plate['design_id'],input_rows=n))
    manifest.append(dict(url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),path=str(path.relative_to(ROOT))))
    print(counts[-1],flush=True)
t = pd.DataFrame(records)
assert t.targetids.str.fullmatch(r'2MASS-J\d{8}[+-]\d{7}').all()
t['APOGEE_ID'] = t.targetids.str.replace('2MASS-J','2M',regex=False)
t['normal'] = t.apogee2_target1.astype('int64').map(lambda b:bool(b & (1<<14)))
failure_path = ROOT / 'research_work/results/selection-parent-eligibility/results.json'
failures = json.loads(failure_path.read_text())['failures']
ids = [r['APOGEE_ID'] for r in failures]
assert all(r['failed'] == ['selected_midIR_error_0_to_0p1'] for r in failures)
exception = t.loc[t.APOGEE_ID.isin(ids)]
assert set(exception.APOGEE_ID) == set(ids)
drilled_path = ROOT / 'research_work/data-cache/selection-plate-history/drilled-science-metadata.parquet'
d = pd.read_parquet(drilled_path,columns=['design_id','APOGEE_ID'])
duplicates = t.loc[t.duplicated(['design_id','APOGEE_ID'],keep=False)]
joined = d.merge(t,on=['design_id','APOGEE_ID'],how='left',indicator=True,validate='one_to_many')
assert joined._merge.eq('both').all()
summary = []
for sid,part in exception.groupby('APOGEE_ID'):
    summary.append(dict(APOGEE_ID=sid,designs=part.design_id.tolist(),normal_all=bool(part.normal.all()),
                        qa=sorted(part.qa.unique().tolist()),selbits=sorted(part.selbits.unique().tolist())))
out = ROOT/'research_work/data-cache/selection-science-inputs'
out.mkdir(exist_ok=True)
t.to_parquet(out/'ranked-candidates.parquet',index=False)
result = dict(scope='Ranked science-input membership; not complete pre-selection parent or inclusion probabilities',
    manifest=manifest,designs=counts,total_input_rows=len(t),unique_input_targets=int(t.APOGEE_ID.nunique()),
    unique_normal_targets=int(t.loc[t.normal,'APOGEE_ID'].nunique()),drilled_assignments_matched=len(d),
    repeated_design_target_groups=duplicates[['design_id','APOGEE_ID']].drop_duplicates().to_dict('records'),
    infrared_exception_count=len(ids),infrared_exceptions=summary,
    input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [previous,failure_path,drilled_path]},
    heldout_kinematics_read=False)
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print({k:result[k] for k in ['total_input_rows','unique_input_targets','unique_normal_targets','drilled_assignments_matched','infrared_exception_count']})
