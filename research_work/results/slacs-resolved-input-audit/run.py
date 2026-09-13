"""Audit training radial arrays and pinned observational metadata only."""
from pathlib import Path
import urllib.request,json,hashlib,csv,io
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
prior=json.loads((HERE.parent/'slacs-resolved-data/manifest.json').read_text(encoding='utf-8'))
CACHE=ROOT/'research_work/generated/slacs-resolved-data';commit=prior['commit'];repo=prior['repository']
tree=json.load(urllib.request.urlopen('https://api.github.com/repos/'+repo+'/git/trees/'+commit+'?recursive=1'));assert not tree['truncated']
lookup={r['path']:r for r in tree['tree'] if r['type']=='blob'}
extra=[]
for name in ['ExternalLenses/SLACS/kinematic_sample_slacs_preprocessing.ipynb','ExternalLenses/SLACS/slacs_all_params.csv']:
 url='https://raw.githubusercontent.com/'+repo+'/'+commit+'/'+name
 raw=urllib.request.urlopen(url).read();entry=lookup[name]
 assert len(raw)==entry['size'] and hashlib.sha1(('blob '+str(len(raw))+'\0').encode()+raw).hexdigest()==entry['sha']
 path=CACHE/name;path.write_bytes(raw)
 extra.append(dict(path=name,url=url,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),git_blob_sha1=entry['sha']))
# Verify every original file; numerical inspection remains restricted to training.
for f in prior['files']:assert hashlib.sha256((ROOT/f['cache_path']).read_bytes()).hexdigest()==f['sha256']
metadata={r['name']:r for r in csv.DictReader(io.StringIO((CACHE/'ExternalLenses/SLACS/slacs_all_params.csv').read_text(encoding='utf-8')))}
rows=[]
for system in prior['systems']:
 if system['role']!='training':continue
 name='SDSS'+system['Name'];folder=CACHE/'ExternalLenses/SLACS/slacs_kcwi_data'/name
 with (folder/(name+'_kinmap.csv')).open(encoding='utf-8') as f:
  reader=csv.DictReader(f);assert reader.fieldnames==['bin_inner_edge','bin_outer_edge','binned_Vrms'];data=list(reader)
 inner=np.array([float(r['bin_inner_edge']) for r in data]);outer=np.array([float(r['bin_outer_edge']) for r in data]);v=np.array([float(r['binned_Vrms']) for r in data])
 cov=np.loadtxt(folder/(name+'_kinmap_cov.csv'),delimiter=',',comments='#')
 assert cov.shape==(len(v),len(v)) and np.all(np.isfinite(cov)) and np.all(np.isfinite(v)) and np.all(v>0)
 assert inner[0]==0 and np.all(outer>inner) and np.all(np.diff(inner)>0) and np.allclose(inner[1:],outer[:-1],rtol=0,atol=1e-12)
 asym=float(np.max(abs(cov-cov.T)));assert asym<1e-10
 eig=np.linalg.eigvalsh(cov);assert eig.min()>0;np.linalg.cholesky(cov)
 sd=np.sqrt(np.diag(cov));cor=cov/sd[:,None]/sd[None,:];off=cor[np.triu_indices(len(v),1)]
 m=metadata[name];psf=float(m['psf_ifu_kcwi']);assert np.isfinite(psf) and psf>0
 rows.append(dict(Name=system['Name'],role='training',n_bins=len(v),inner_arcsec=inner.tolist(),outer_arcsec=outer.tolist(),vrms_kms=v.tolist(),covariance_kms_squared=cov.tolist(),psf_fwhm_arcsec=psf,release_use_flag=m['flag_ifu_kcwi'],maximum_covariance_asymmetry=asym,minimum_covariance_eigenvalue=float(eig.min()),covariance_condition_number=float(eig.max()/eig.min()),off_diagonal_correlation_range=[float(off.min()),float(off.max())],fractional_diagonal_error_range=[float((sd/v).min()),float((sd/v).max())]))
out=dict(scope='Training input validity only; no model fit or reserved radial score',source_commit=commit,additional_sources=extra,systems=rows,total_training_bins=sum(r['n_bins'] for r in rows),notebook_interpretation='Loads V_rms and full covariance directly; IFU_shells centered at zero, Gaussian FWHM from psf_ifu_kcwi',unresolved='Original Voronoi-to-radial aggregation and covariance Monte Carlo regeneration not provided by inspected preprocessing notebook; README formula typo not used')
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(systems=len(rows),bins=out['total_training_bins'],psf_range=[min(r['psf_fwhm_arcsec'] for r in rows),max(r['psf_fwhm_arcsec'] for r in rows)],maximum_correlation=max(r['off_diagonal_correlation_range'][1] for r in rows),maximum_condition=max(r['covariance_condition_number'] for r in rows),flags={r['Name']:r['release_use_flag'] for r in rows})))
