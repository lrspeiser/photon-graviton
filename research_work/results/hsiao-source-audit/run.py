from pathlib import Path
import json,hashlib,tarfile,io,requests
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[2]
a=R/'research_work/generated/essence-passbands/hsiao_template.tar.gz';url='https://astrophysics.physics.fsu.edu/~hsiao/data/hsiao_template.tar.gz'
if not a.exists():
 response=requests.get(url,timeout=60);response.raise_for_status();a.parent.mkdir(parents=True,exist_ok=True);a.write_bytes(response.content)
sha=hashlib.sha256(a.read_bytes()).hexdigest()
if (H/'results.json').exists():assert json.loads((H/'results.json').read_text(encoding='utf-8'))['author_archive_sha256']==sha
with tarfile.open(a) as t:
 raw=t.extractfile('snflux_1a.dat').read();lcraw=t.extractfile('lc_template.dat').read()
x=np.loadtxt(io.BytesIO(raw));lc=np.loadtxt(io.BytesIO(lcraw),comments=';')
p=R/'research_work/generated/des-photometric-calibration/selected/snsed/Hsiao07.dat';y=np.loadtxt(p)
manifest=json.loads((H.parent/'des-photometric-calibration/dependencies.json').read_text(encoding='utf-8'))
record=next(v for v in manifest['selected'] if v['archive_path']=='snsed/Hsiao07.dat');assert hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256']
assert np.array_equal(x,y)
phase=np.unique(x[:,0]);wave=np.unique(x[:,1]);assert len(x)==len(phase)*len(wave);assert np.array_equal(x[:,:2],np.array([(t,w) for t in phase for w in wave]))
assert np.all(np.isfinite(x)) and np.all(x[:,2]>=0);assert np.all(np.diff(lc[:,0])>0)
out=dict(scope='Author-release identity and conditional source-model readiness; not independent temporal calibration',author_archive_url=url,author_archive_sha256=sha,author_members={'snflux_1a.dat':hashlib.sha256(raw).hexdigest(),'lc_template.dat':hashlib.sha256(lcraw).hexdigest()},cached_source=str(p.relative_to(R)),cached_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),numerically_identical=True,rows=len(x),phase_count=len(phase),phase_range_days=phase[[0,-1]].tolist(),wavelength_count=len(wave),wavelength_range_angstrom=wave[[0,-1]].tolist(),phase_step_days=np.unique(np.diff(phase)).tolist(),wavelength_step_angstrom=np.unique(np.diff(wave)).tolist(),nonfinite_values=int(np.sum(~np.isfinite(x))),negative_flux_values=int(np.sum(x[:,2]<0)),first_phase_flux_range=[float(x[x[:,0]==phase[0],2].min()),float(x[x[:,0]==phase[0],2].max())],lightcurve_columns=lcraw.decode().splitlines()[0],lightcurve_shape=list(lc.shape),lightcurve_phase_range_days=lc[[0,-1],0].tolist(),lightcurve_phase_grid_matches_spectra=bool(np.array_equal(lc[:,0],phase)),source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(out))
