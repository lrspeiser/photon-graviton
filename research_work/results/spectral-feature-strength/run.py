"""Measure residual spectral structure; do not interpret it as supernova-only signal."""
from pathlib import Path
import json,hashlib
import numpy as np
from astropy.io import fits
H=Path(__file__).resolve().parent;R=H.parents[2]
ap=H.parent/'essence-spectrum-acquisition/results.json';audit=json.loads(ap.read_text())
records=sorted([x for x in audit['spectra'] if x['target']=='q106'],key=lambda x:x['mid_mjd'])
hashes={str(ap.relative_to(R)):hashlib.sha256(ap.read_bytes()).hexdigest()};rows=[]
def read(rec):
 p=R/'research_work/generated/essence-dr3'/rec['filename'];sha=hashlib.sha256(p.read_bytes()).hexdigest();assert sha==next(x['sha256'] for x in audit['files'] if x['filename']==p.name);hashes[str(p.relative_to(R))]=sha
 with fits.open(p) as hd:
  d=hd[1].data[0];return [np.asarray(d[k],float).copy() for k in ['WAVE','FLUX','ERR']]
for rec in records:
 w,f,e=read(rec);host=next(x for x in audit['spectra'] if x['target']=='q106-gal' and abs(x['mid_mjd']-rec['mid_mjd'])<1e-8)
 hw,hf,he=read(host);assert np.array_equal(hw,w)
 rest=w/(1+float(rec['catalog_identity']['z_host']));mask=(rest>=4000)&(rest<=5500)
 x=rest[mask];x=(x-x.mean())/np.ptp(x);f=f[mask];e=e[mask];hf=hf[mask];he=he[mask];n=len(x)
 assert np.all(np.isfinite(f)) and np.all(e>0);scale=float(np.median(f));assert scale>0
 for degree in [2,3,5]:
  polynomial=np.column_stack([x**j for j in range(degree+1)])
  for host_basis in [False,True]:
   design=np.column_stack([polynomial,hf/np.median(hf)]) if host_basis else polynomial
   q=np.linalg.qr(design/e[:,None],mode='reduced')[0];rank=design.shape[1];assert np.linalg.matrix_rank(design/e[:,None])==rank
   residual=e*(f/e-q@(q.T@(f/e)))
   noise_trace=float(np.sum(e**2*(1-np.sum(q*q,axis=1))))
   energy=float(residual@residual);excess=(energy-noise_trace)/n
   # Independent direct least-squares implementation checks the projection.
   coeff=np.linalg.lstsq(design/e[:,None],f/e,rcond=None)[0]
   assert np.max(np.abs(residual-(f-design@coeff)))<1e-9*max(1.,np.max(abs(f)))
   weighted_excess=(float(np.sum((residual/e)**2))-(n-rank))/float(np.sum(1/e**2))
   rows.append(dict(weighted_positive_excess_rms_fraction=float(np.sqrt(max(0.,weighted_excess))/scale),filename=rec['filename'],mid_mjd=rec['mid_mjd'],degree=degree,host_basis=host_basis,n_pixels=n,rank=rank,median_flux=scale,median_error_fraction=float(np.median(e)/scale),residual_rms_fraction=float(np.sqrt(energy/n)/scale),expected_noise_rms_fraction=float(np.sqrt(noise_trace/n)/scale),signed_excess_variance_fraction=float(excess/scale**2),positive_excess_rms_fraction=float(np.sqrt(max(0.,excess))/scale),chi_square_per_dof=float(np.sum((residual/e)**2)/(n-rank)),host_basis_coefficient=float(coeff[-1]) if host_basis else None))
# Verify the known trace formula with independent diagonal-noise Monte Carlo.
rng=np.random.default_rng(4101);draw=rng.normal(size=(4000,n));projected=(draw-(draw@q)@q.T)*e
mc=float(np.mean(np.sum(projected**2,axis=1)));relative=abs(mc/noise_trace-1);assert relative<.02
out=dict(scope='Observed residual structure after polynomial/nuisance-host projection; not a spectral age or uncontaminated feature estimate',window_rest_angstrom=[4000,5500],rows=rows,noise_trace_check=dict(draws=4000,analytic=noise_trace,monte_carlo=mc,relative_error=relative),hashes=hashes,limitations=['Diagonal pixel covariance assumed','Host basis treated fixed despite its own extraction uncertainty and shared errors','Host coefficient is a nuisance projection, not measured contamination or physical mixing fraction','Polynomial may remove real broad supernova features','Residual structure can contain sky/telluric artifacts and host residuals','Barycentric/heliocentric frame approximation retained','Pre-projection injected feature RMS is not identical to measured post-projection residual RMS'])
out['hashes'][str(Path(__file__).relative_to(R))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
for row in rows:print(json.dumps({k:row[k] for k in ['mid_mjd','degree','host_basis','positive_excess_rms_fraction','chi_square_per_dof']}))
print(out['noise_trace_check'])
