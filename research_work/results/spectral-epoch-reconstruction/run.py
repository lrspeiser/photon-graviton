"""Recover published observer dates; compare equivalent age-slope likelihoods."""
from pathlib import Path
import csv
import hashlib
import json
import re
import urllib.request
import numpy as np
from bs4 import BeautifulSoup

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
source=ROOT/'research_work/generated/snid-template-audit/blondin2008.html'
if not source.exists():
    source.parent.mkdir(parents=True,exist_ok=True)
    with urllib.request.urlopen('https://arxiv.org/html/0804.3595v1',timeout=40) as response:
        source.write_bytes(response.read())
if (HERE/'results.json').exists():
    previous=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    assert hashlib.sha256(source.read_bytes()).hexdigest()==previous['hashes'][str(source.relative_to(ROOT))]
catalog=HERE.parent/'electromagnetic-audit/spectral-aging-predictions.csv'
with catalog.open(newline='',encoding='utf-8') as f:
    summary={r['object']:r for r in csv.DictReader(f)}
table=BeautifulSoup(source.read_bytes(),'html.parser').find(id='S3.T2')
epochs=[]
name=None
for tr in table.select('tr')[2:]:
    cells=[td.get_text(' ',strip=True) for td in tr.find_all(['td','th'],recursive=False)]
    if len(cells)!=5: continue
    if cells[0]: name=cells[0]
    assert name in summary
    nums=[float(v) for v in re.findall(r'\d+(?:\.\d+)?',cells[2])]
    assert len(nums)==2
    phase=(-1 if '-' in cells[2] or '\u2212' in cells[2] else 1)*nums[0]
    epochs.append(dict(object=name,JD=2450000+float(cells[1]),phase_days=phase,
                       phase_error_days=nums[1],published_elapsed_days=float(cells[3])))
assert len(epochs)==35
results=[]
for name in dict.fromkeys(e['object'] for e in epochs):
    rows=[e for e in epochs if e['object']==name]
    t=np.array([r['JD'] for r in rows]);t-=t[0]
    age=np.array([r['phase_days'] for r in rows]);sigma=np.array([r['phase_error_days'] for r in rows])
    assert np.max(abs(t-np.array([r['published_elapsed_days'] for r in rows])))<1e-8
    X=np.column_stack([np.ones(len(t)),t]);w=sigma**-2
    covariance=np.linalg.inv(X.T@(w[:,None]*X))
    fit=covariance@(X.T@(w*age))
    dt,dy=t[1:],age[1:]-age[0]
    C=np.diag(sigma[1:]**2)+sigma[0]**2*np.ones((len(dt),len(dt)))
    inverse_dt=np.linalg.solve(C,dt)
    slope=float(inverse_dt@dy/(dt@inverse_dt))
    variance=float(1/(dt@inverse_dt))
    assert abs(slope-fit[1])<1e-10 and abs(variance-covariance[1,1])<1e-10
    # Deliberately incomplete comparison, retained to expose origin-error effects.
    naive=float(np.sum(dt*dy/sigma[1:]**2)/np.sum(dt**2/sigma[1:]**2))
    naive_error=float(1/np.sqrt(np.sum(dt**2/sigma[1:]**2)))
    published=summary[name]
    results.append(dict(object=name,n_epochs=len(rows),observer_span_days=float(t[-1]),
                        z=float(published['z']),slope=slope,error=float(np.sqrt(variance)),
                        naive_fixed_origin_slope=naive,naive_fixed_origin_error=naive_error,
                        published_slope=float(published['observed_aging_rate']),
                        published_error=float(published['sigma'])))
out=dict(scope='Reconstruction from rounded published epoch ages; not raw-spectral replication',
         source='https://arxiv.org/html/0804.3595v1#S3.T2',epochs=epochs,results=results,
         assumptions='Independent absolute spectral-age errors before differencing; unknown template covariance omitted',
         hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,catalog,Path(__file__)]})
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(results,indent=2))
