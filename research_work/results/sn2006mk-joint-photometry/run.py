from pathlib import Path
import requests,json,hashlib
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[2]
url='https://cdsarc.cds.unistra.fr/ftp/J/ApJS/224/3/lcs/q106.W6yr.clean.nn2.Wstd.dat'
p=R/'research_work/generated/essence-dr3/q106.W6yr.clean.nn2.Wstd.dat'
if not p.exists():
 response=requests.get(url,timeout=60);response.raise_for_status();p.write_bytes(response.content)
sha=hashlib.sha256(p.read_bytes()).hexdigest()
if (H/'results.json').exists():assert json.loads((H/'results.json').read_text(encoding='utf-8'))['source']['sha256']==sha
text=p.read_text(encoding='utf-8');assert 'zeropoint of 25' in text
rows=[]
for line in text.splitlines():
 if not line.strip() or line.startswith('#'):continue
 f=line.split();assert len(f)==6
 row=dict(observation=f[0],mjd=float(f[1]),band=f[2],flux=float(f[3]),error_lower=float(f[4]),error_upper=float(f[5]));assert row['error_lower']>0 and row['error_upper']>0;rows.append(row)
assert len(set(x['observation'] for x in rows))==len(rows)
auditpath=H.parent/'essence-spectrum-acquisition/results.json';audit=json.loads(auditpath.read_text(encoding='utf-8'));epochs=sorted([x for x in audit['spectra'] if x['target']=='q106'],key=lambda x:x['mid_mjd'])
summary=[];links=[]
for band in sorted(set(x['band'] for x in rows)):
 rr=sorted([x for x in rows if x['band']==band],key=lambda x:x['mjd']);summary.append(dict(band=band,n=len(rr),negative_flux_rows=sum(x['flux']<0 for x in rr),mjd_range=[rr[0]['mjd'],rr[-1]['mjd']]))
 near=[]
 for ep in epochs:
  x=min(rr,key=lambda x:abs(x['mjd']-ep['mid_mjd']));assert abs(x['mjd']-ep['mid_mjd'])<.2;assert x['flux']>0 and x['error_lower']==x['error_upper'];near.append(x)
 a,b=near;ratio=b['flux']/a['flux'];sigma=ratio*np.sqrt((a['error_upper']/a['flux'])**2+(b['error_upper']/b['flux'])**2)
 delta=-2.5*np.log10(ratio);sigma_delta=2.5/np.log(10)*sigma/ratio
 links.append(dict(band=band,nearby_photometry=near,photometry_minus_spectral_midpoint_days=[x['mjd']-ep['mid_mjd'] for x,ep in zip(near,epochs)],later_to_earlier_flux_ratio=float(ratio),diagonal_linearized_ratio_sigma=float(sigma),magnitude_change=float(delta),diagonal_linearized_magnitude_sigma=float(sigma_delta)))
# These formulas use within-band relative flux and do not require an absolute zeropoint.
rband=next(x for x in links if x['band']=='R4m');iband=next(x for x in links if x['band']=='I4m')
color=dict(change_R_minus_I_magnitudes=rband['magnitude_change']-iband['magnitude_change'],diagonal_linearized_sigma=float(np.hypot(rband['diagonal_linearized_magnitude_sigma'],iband['diagonal_linearized_magnitude_sigma'])))
out=dict(scope='Same-object measured photometry linked to spectral epochs; no source-luminosity or distance fit',source=dict(url=url,sha256=sha,bytes=p.stat().st_size),header=text.split('#Observation')[0],photometry=rows,bands=summary,spectral_links=links,color_change=color,redshift_versions=dict(lightcurve_header_heliocentric=0.4750,updated_catalog_host_heliocentric=0.4754),hashes={str(q.relative_to(R)):hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__),auditpath]})
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:out[k] for k in ['bands','spectral_links','color_change']}))
