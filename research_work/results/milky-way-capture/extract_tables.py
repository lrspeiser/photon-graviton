"""Extract inputs from locally saved primary-paper HTML (optional beautifulsoup4).

python extract_tables.py --eilers-html eilers.html --bovy-html bovy.html --retrieved YYYY-MM-DD
Download the URLs below separately; extraction validates table size and sample cells.
Running this overwrites inputs.json, so review a refreshed source before reanalysis.
"""
import argparse
import hashlib
import json
from pathlib import Path
from bs4 import BeautifulSoup

p=argparse.ArgumentParser()
p.add_argument('--eilers-html',type=Path,required=True)
p.add_argument('--bovy-html',type=Path,required=True)
p.add_argument('--retrieved',required=True)
a=p.parse_args()
sources={}
for name,index,width,columns,url in [
 ('eilers',8,4,['R_kpc','vc_kms','err_minus_kms','err_plus_kms'],'https://arxiv.org/html/1810.09466'),
 ('bovy',19,8,['Fe_H','alpha_Fe','R_kpc','surface_density_Msun_pc2','surface_density_error','R0_minus_R_kpc','Kz_over_2piG_Msun_pc2','Kz_error'],'https://arxiv.org/html/1309.0809')]:
    raw=getattr(a,name+'_html').read_bytes()
    table=BeautifulSoup(raw,'html.parser').find_all('table')[index]
    rows=[]
    for tr in table.find_all('tr'):
        cells=tr.find_all(['td','th'],recursive=False)
        if len(cells)!=width: continue
        try: values=[float(c.get_text(' ',strip=True)) for c in cells]
        except ValueError: continue
        rows.append(dict(zip(columns,values)))
    assert len(rows)==(38 if name=='eilers' else 43),(name,len(rows))
    assert rows[0]['R_kpc']==(5.27 if name=='eilers' else 4.63)
    sources[name]={'url':url,'retrieved':a.retrieved,'source_html_sha256':hashlib.sha256(raw).hexdigest(),'table':1 if name=='eilers' else 3,'rows':rows}
sources['user_table']={'provenance':'Six rows pasted by user; stellar baseline and remaining three claimed rows unavailable. Rounded percentages are not published observational uncertainties.','rows':[dict(R_kpc=r,speed_excess_percent=p) for r,p in [(5.27,16),(8.19,25),(12.25,36),(15.22,44),(20.27,50),(24.82,63)]]}
Path(__file__).with_name('inputs.json').write_text(json.dumps(sources,indent=2)+'\n',encoding='utf-8',newline='\n')
