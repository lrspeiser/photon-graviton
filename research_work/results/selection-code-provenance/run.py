"""Archive public legacy targeting code and test its distinct missing-error guards."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET
import requests
import numpy as np
from astropy.io import fits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache/selection-code-provenance'
CACHE.mkdir(exist_ok=True)
BASE='https://svn.sdss.org/public/repo/apogee/apogeetarget/'
files=['tags/','trunk/pro/check_data_quality.pro','trunk/pro/select_science.pro','trunk/pro/select_science2.pro','trunk/pro/compute_second_quant.pro']
manifest=[]
texts={}
for name in files:
    path=CACHE/('tags.xml' if name=='tags/' else name.split('/')[-1])
    if not path.exists():
        r=requests.get(BASE+name,timeout=(15,30));r.raise_for_status();path.write_bytes(r.content)
    texts[name]=path.read_text()
    manifest.append(dict(url=BASE+name,path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
tags=ET.fromstring(texts['tags/']).find('index')
assert 'finite(mean(data.e_4_5_m,/nan))' in texts['trunk/pro/select_science.pro']
assert 'finite(mean(data.e_4_5_m))' in texts['trunk/pro/select_science2.pro']
# Numerical illustration of the two guards, not execution of IDL or a survey rerun.
example=np.array([.05,np.nan,.15])
ordinary_guard=bool(np.isfinite(np.mean(example)))
nan_aware_guard=bool(np.isfinite(np.nanmean(example)))
assert not ordinary_guard and nan_aware_guard
parent=ROOT/'research_work/data-cache/selection-fields/apogee2Object_300+00.fits'
with fits.open(parent) as hd:
    t=hd[1].data
    eligible_ext=np.isfinite(t['AK_TARG']) & (t['AK_TARG'] > -50)
    error=np.array(t['TARG_4_5_ERR'][eligible_ext])
    summary=dict(rows=len(error),nonfinite_errors=int((~np.isfinite(error)).sum()),
                 finite_errors_above_0p1=int((np.isfinite(error)&(error>.1)).sum()),
                 finite_error_min=float(np.nanmin(error)),finite_error_max=float(np.nanmax(error)))
result=dict(scope='Legacy code provenance and guard illustration; 2017 targeting implementation not established',
    manifest=manifest,svn_index_revision=tags.get('rev'),listed_tags=[x.get('name') for x in tags if x.tag=='dir'],
    illustrative_guard_test=dict(errors=[.05,None,.15],ordinary_mean_guard=ordinary_guard,nan_aware_mean_guard=nan_aware_guard,
      ordinary_branch_applies_error_threshold=False,nan_aware_branch_applies_error_threshold=True),
    parent_extinction_available_error_summary=summary,
    parent_sha256=hashlib.sha256(parent.read_bytes()).hexdigest(),
    same_input_as_historical_routine_proven=False,explains_42_exceptions=False,heldout_kinematics_read=False)
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k!='manifest'},indent=2))
