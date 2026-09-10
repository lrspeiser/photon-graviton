"""Retrieve exact 2MASS epoch/position records for the training review list."""
from pathlib import Path
import io,json,hashlib,re
import requests,pandas as pd
HERE=Path(__file__).resolve().parent
OUT=HERE.parents[2]/'research_work/data-cache/stellar-association-window'
candidates=pd.read_csv(HERE/'review-candidates.csv',dtype={'source_id':str,'two_mass_designation':str})
ids=sorted(candidates.two_mass_designation.unique())
assert all(re.fullmatch(r'\d{8}[+-]\d{7}',x) for x in ids)
frames=[]; records=[]
for start in range(0,len(ids),50):
    values=','.join("'"+x+"'" for x in ids[start:start+50])
    query='SELECT designation,ra,dec,jdate,err_maj,err_min,err_ang,j_m,h_m,k_m FROM fp_psc WHERE designation IN ('+values+')'
    path=OUT/f'2mass-epochs-{start//50}.csv'
    archived=HERE/'inputs'/path.name
    if not path.exists() and archived.exists(): path.write_bytes(archived.read_bytes())
    if not path.exists():
        response=requests.post('https://irsa.ipac.caltech.edu/TAP/sync',data={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query},timeout=45)
        response.raise_for_status()
        parsed=pd.read_csv(io.StringIO(response.text),dtype={'designation':str})
        assert 'jdate' in parsed
        path.write_text(response.text,encoding='utf-8',newline='\n')
    frame=pd.read_csv(path,dtype={'designation':str});frames.append(frame)
    records.append({'query':query,'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':len(frame)})
    print('Retrieved batch',start//50, len(frame),flush=True)
allrows=pd.concat(frames,ignore_index=True)
assert allrows.designation.is_unique and set(allrows.designation)==set(ids)
allrows.to_csv(OUT/'2mass-review-epochs.csv',index=False)
(HERE/'download.json').write_text(json.dumps({'endpoint':'https://irsa.ipac.caltech.edu/TAP/sync','requested':len(ids),'returned':len(allrows),'batches':records},indent=2)+'\n',encoding='utf-8',newline='\n')
