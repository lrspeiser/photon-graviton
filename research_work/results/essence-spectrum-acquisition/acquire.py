"""Acquire identified development spectra; never infer identity from date alone."""
from pathlib import Path
import csv, json, hashlib, requests
import numpy as np
from astropy.io import fits
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/generated/essence-dr3'
CACHE.mkdir(parents=True,exist_ok=True)
OLD=json.loads((HERE/'results.json').read_text()) if (HERE/'results.json').exists() else None
oldhash={x['filename']:x['sha256'] for x in OLD['files']} if OLD else {}
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def acquire(name,url,params=None):
 p=CACHE/name
 if not p.exists():
  r=requests.get(url,params=params,timeout=60);r.raise_for_status();p.write_bytes(r.content)
 if name in oldhash: assert digest(p)==oldhash[name],f'Changed cached input: {name}'
 return p
params={'REQUEST':'doQuery','FORMAT':'csv','LANG':'ADQL','MAXREC':1000,'QUERY':"SELECT * FROM ivoa.ObsCore WHERE obs_collection='ESSENCE'"}
meta=acquire('obscore.csv','https://archive.eso.org/tap_obs/sync',params)
cat=acquire('table6.dat','https://cdsarc.cds.unistra.fr/ftp/J/ApJS/224/3/table6.dat')
# Served catalog is pipe-delimited despite its fixed-width ReadMe description.
catalog={}
for line in cat.read_text().splitlines():
 f=[v.strip() for v in line.split('|')]
 assert len(f)==13
 catalog[f[0]]={'iau_name':f[1],'classification':f[2],'coordinates':f[6],'z_snid':f[7],'z_snid_error':f[8],'z_host':f[9],'z_host_error':f[10]}
# Inspect actual field layout below: coordinates are a single column.
rows=list(csv.DictReader(meta.open()))
names={'b013','b027','d093','q106','r190','s347'}
selected=[x for x in rows if 'v3.0' in x['obs_creator_did'] and x['target_name'].split('-')[0] in names]
files=[{'filename':meta.name,'sha256':digest(meta),'url':'https://archive.eso.org/tap_obs/sync','params':params},{'filename':cat.name,'sha256':digest(cat),'url':'https://cdsarc.cds.unistra.fr/ftp/J/ApJS/224/3/table6.dat'}]
records=[]
epochs=json.loads((ROOT/'research_work/results/spectral-epoch-reconstruction/results.json').read_text())['epochs']
for x in selected:
 name=x['obs_creator_did'].split('?')[-1];url='https://dataportal.eso.org/dataPortal/file/'+x['dp_id'];p=acquire(name,url)
 files.append({'filename':name,'sha256':digest(p),'bytes':p.stat().st_size,'url':url})
 with fits.open(p,checksum=True) as h:
  h.verify('exception');checks=[{'checksum':a.verify_checksum(),'datasum':a.verify_datasum()} for a in h];assert all(a['checksum']==1 and a['datasum']==1 for a in checks)
  ph,sh=h[0].header,h[1].header;d=h[1].data[0];w,f,e=[np.asarray(d[k],float) for k in ['WAVE','FLUX','ERR']]
  assert np.all(np.isfinite(w)) and np.all(np.diff(w)>0)
  assert sh['TUNIT1'].lower()=='angstrom'
  assert abs(w[0]/10-ph['WAVELMIN'])<.1 and abs(w[-1]/10-ph['WAVELMAX'])<.1
  good=np.isfinite(f)&np.isfinite(e)&(e>0)
  base=x['target_name'].split('-')[0];identity=catalog[base]
  obj='b027' if base=='b027' else identity['iau_name'];pe=[a for a in epochs if a['object']==obj]
  nearest=min(pe,key=lambda a:abs(a['JD']-(sh['TMID']+2400000.5))) if pe else None
  rec={'filename':name,'target':x['target_name'],'catalog_identity':identity,'metadata':x,'header':{k:ph.get(k) for k in ['OBJECT','DATE-OBS','MJD-OBS','MJD-END','SPECSYS','TIMESYS','FLUXCAL','TOT_FLUX','FLUXERR','CONTNORM','SPEC_RES','SNR','BARYVEL','REFERENC','PROG_ID']},'mid_mjd':sh['TMID'],'wavelength_unit':sh['TUNIT1'],'wavelength_comment':sh.get('TCOMM1'),'flux_unit':sh['TUNIT2'],'n_pixels':len(w),'valid_error_pixels':int(good.sum()),'wavelength_range_angstrom':[float(w[0]),float(w[-1])],'checksum_checks':checks,'nearest_published_epoch':nearest,'midpoint_minus_published_days':float(sh['TMID']+2400000.5-nearest['JD']) if nearest else None}
  records.append(rec)
 print(name,records[-1]['valid_error_pixels'],flush=True)
result={'scope':'Acquisition and metadata audit only; no distant spectral-age fit or causal validation','catalog_rows':len(catalog),'archive_rows':len(rows),'dr3_rows':sum('v3.0' in x['obs_creator_did'] for x in rows),'files':files,'spectra':records,'excluded_alias_note':'E142 is described as D093 in DR3 documentation, but contradictory E132 text also occurs. Not included until separate position/provenance audit. Date agreement alone never establishes object identity.'}
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

