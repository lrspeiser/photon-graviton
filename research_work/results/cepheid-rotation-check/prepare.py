"""Archive and parse the published Cepheid-only rotation table, not a halo fit."""
from pathlib import Path
import hashlib
import json
import requests
from bs4 import BeautifulSoup

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache/cepheid-rotation-check'
CACHE.mkdir(parents=True,exist_ok=True)
URL='https://oup.silverchair-cdn.com/article-minimal/8416425'
path=CACHE/'article.html'
if not path.exists():
    response=requests.get(URL,timeout=45);response.raise_for_status()
    path.write_bytes(response.content)
soup=BeautifulSoup(path.read_text(encoding='utf-8'),'html.parser')
copies=[]
for table in soup.find_all('table'):
    header=table.get_text(' ',strip=True)
    if not header.startswith('R (kpc)'):continue
    rows=[]
    for tr in table.find_all('tr'):
        cells=[td.get_text(' ',strip=True) for td in tr.find_all('td')]
        if len(cells)!=3:continue
        try:values=list(map(float,cells))
        except ValueError:continue
        rows.append(dict(R_kpc=values[0],vc_kms=values[1],error_kms=values[2]))
    if rows:copies.append(rows)
assert len(copies)==4 and all(rows==copies[0] for rows in copies)
rows=copies[0]
assert len(rows)==12 and len({r['R_kpc'] for r in rows})==12
assert all(6<=r['R_kpc']<=18 and r['vc_kms']>0 and r['error_kms']>0 for r in rows)
out=dict(title='Feng et al. 2026, Cepheid-only Table 1',doi='10.1093/mnras/stag011',
    official_article='https://academic.oup.com/mnras/article/546/2/stag011/8416425',
    archive_url=URL,article_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    identical_html_table_copies=4,rows=rows,eligible_rows=12,excluded_rows=[],
    error_scope='Published bootstrap errors from 100 resamples; not all systematic uncertainty or bin covariance',
    frame=dict(R_sun_kpc=8.275,z_sun_kpc=.025,V_R_sun_kms=11.1,V_phi_sun_kms=250.2,V_z_sun_kms=7.9),
    independence='Cepheid population newly evaluated here, shared Gaia and Galaxy; individual-star independence not verified',
    scope='Final drift-corrected Cepheid-only rotation table; no combined-tracer curve or fitted halo imported')
(HERE/'observations.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
print(json.dumps(out,indent=2))
