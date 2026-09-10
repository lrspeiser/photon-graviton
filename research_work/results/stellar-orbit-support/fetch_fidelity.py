"""Retrieve public astrometric-quality estimates for the already matched IDs."""
import json,time
import pandas as pd
import pyvo
from astropy.table import Table
from run import CACHE,OUT,HERE,save,digest

endpoint='https://dc.g-vo.org/tap'
query='SELECT t.source_id, f.fidelity_v1, f.fidelity_v2 FROM TAP_UPLOAD.targets AS t LEFT OUTER JOIN gedr3spur.main AS f ON t.source_id=f.source_id'
path=OUT/'gaia-astrometric-fidelity.fits';state=HERE/'fidelity-job.json'
ids=pd.read_parquet(CACHE/'matched-with-gaia-covariance.parquet',columns=['source_id'])
assert ids.source_id.is_unique and ids.source_id.dtype.kind in 'iu'
if not path.exists():
    if state.exists():
        old=json.loads(state.read_text())
        assert old['query']==query and old['requested_sources']==len(ids)
        job=pyvo.dal.AsyncTAPJob(old['job_url'])
    else:
        service=pyvo.dal.TAPService(endpoint)
        job=service.submit_job(query,uploads={'targets':Table.from_pandas(ids)},maxrec=200000)
        save('fidelity-job.json',dict(endpoint=endpoint,job_url=job.url,query=query,requested_sources=len(ids)))
    if job.phase=='PENDING':job.run()
    last=None
    while True:
        phase=job.phase
        if phase!=last:print('Fidelity TAP phase:',phase,flush=True);last=phase
        if phase in ['COMPLETED','ERROR','ABORTED']:break
        time.sleep(5)
    job.raise_if_error()
    assert phase=='COMPLETED'
    table=job.fetch_result().to_table()
    assert len(table)==len(ids)
    table.write(path,overwrite=False)
table=Table.read(path).to_pandas()
assert table.source_id.is_unique and set(table.source_id)==set(ids.source_id)
for key in ['fidelity_v1','fidelity_v2']:
    assert table.loc[table[key].notna(),key].between(0,1).all()
save('fidelity-download.json',dict(endpoint=endpoint,table='gedr3spur.main',query=query,
    requested_sources=len(ids),returned_sources=len(table),missing_fidelity_v1=int(table.fidelity_v1.isna().sum()),
    missing_fidelity_v2=int(table.fidelity_v2.isna().sum()),sha256=digest(path),
    source_documentation='https://dc.g-vo.org/tableinfo/gedr3spur.main',
    status='Public source-ID join; no velocity or gravity selection used'))
print('Fidelity file verified:',len(table),'sources',flush=True)
