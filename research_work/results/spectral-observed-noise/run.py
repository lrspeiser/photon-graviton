from pathlib import Path
import hashlib,json
import numpy as np
from scipy.ndimage import gaussian_filter1d
from astropy.io import fits
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
p=json.loads((HERE/'protocol.json').read_text(encoding='utf-8-sig'))
loader=HERE.parent/'spectral-object-exclusion/run.py';s=loader.read_text(encoding='utf-8');assert s.count('rows=[]')==1
ns={'__file__':str(loader),'__name__':'library_loader'};exec(compile(s.split('rows=[]')[0],str(loader),'exec'),ns)
library,names,phase,targets,wave=[ns[k] for k in ['library','names','phase','targets','wave']]
pairs=[]
for obj in sorted(set(names[targets])):
 ix=[i for i in targets if names[i]==obj];a=min(ix,key=lambda i:phase[i]);b=max(ix,key=lambda i:phase[i])
 if phase[b]-phase[a]>=p['pair_minimum_phase_span_days']:pairs.append((a,b))
assert pairs
objects=np.unique(names);groups=[np.flatnonzero(names==o) for o in objects]
raw=np.array([r['values'] for r in library]);lnstep=float(np.median(np.diff(np.log(wave))))
blur_sigma=1/(p['additional_gaussian_resolving_power']*2.354820045*lnstep)
blurred=gaussian_filter1d(raw,blur_sigma,axis=1)
auditpath=HERE.parent/'essence-spectrum-acquisition/results.json';audit=json.loads(auditpath.read_text())
records=sorted([x for x in audit['spectra'] if x['target']=='q106'],key=lambda x:x['mid_mjd']);assert len(records)==2
hashes={str(v.relative_to(ROOT)):hashlib.sha256(v.read_bytes()).hexdigest() for v in [Path(__file__),HERE/'protocol.json',loader,auditpath]}
operators=[]
for record in records:
 path=ROOT/'research_work/generated/essence-dr3'/record['filename'];expected=next(v['sha256'] for v in audit['files'] if v['filename']==path.name);assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
 hashes[str(path.relative_to(ROOT))]=expected
 with fits.open(path) as h:
  d=h[1].data[0];z=float(record['catalog_identity']['z_host']);rw=np.asarray(d['WAVE'],float)/(1+z)
  m=(rw>=p['rest_window_angstrom'][0])&(rw<=p['rest_window_angstrom'][1]);rw=rw[m];err=np.asarray(d['ERR'],float)[m];f=np.asarray(d['FLUX'],float)[m]
  assert np.all(np.isfinite(err)) and np.all(err>0)
  scale=float(np.median(f));assert scale>0
  # Resampling templates onto observed pixels preserves the supplied data errors.
  native=np.array([np.interp(rw,wave,a) for a in blurred]);xx=(rw-rw.mean())/np.ptp(rw);design=np.column_stack([xx**j for j in range(p['continuum_polynomial_degree']+1)])
  q=np.linalg.qr(design/err[:,None],mode='reduced')[0]
  def project(a):
   whitened=a/err
   return whitened-(whitened@q)@q.T
  centered=native-native.mean(axis=1,keepdims=True);rms=np.sqrt(np.mean(centered**2,axis=1));assert np.all(rms>0)
  signal=centered/rms[:,None]*scale
  basis=project(signal);basis/=np.linalg.norm(basis,axis=1)[:,None]
  # An arbitrary quadratic continuum must be removed by this projection.
  residual=project(np.broadcast_to(design@np.array([2.,-.7,.4]),(1,len(rw))))
  assert np.max(np.abs(residual))<1e-8
  operators.append(dict(signal=signal,basis=basis,err=err,q=q,pixels=len(rw),scale=scale,median_noise=float(np.median(err)/scale)))

def estimate(op,values,indices):
 y=values/op['err'];y-=(y@op['q'])@op['q'].T;y/=np.linalg.norm(y,axis=1)[:,None];score=y@op['basis'].T
 objscore=np.empty((len(indices),len(objects)));objbest=np.empty(objscore.shape,int)
 for j,ix in enumerate(groups):
  best=ix[np.argmax(score[:,ix],axis=1)];objbest[:,j]=best;objscore[:,j]=score[np.arange(len(indices)),best]
 objscore[names[indices,None]==objects[None,:]]=-np.inf
 top=np.argsort(objscore,axis=1)[:,-p['top_distinct_objects']:];best=np.take_along_axis(objbest,top,axis=1);assert np.all(names[best]!=names[indices,None])
 return np.median(phase[best],axis=1)
rows=[]
for seed in [None]+p['seeds']:
 draws=[np.random.default_rng(np.random.SeedSequence([seed,e])).normal(size=(len(pairs),op['pixels'])) if seed is not None else np.zeros((len(pairs),op['pixels'])) for e,op in enumerate(operators)]
 for amp in ([1.] if seed is None else p['feature_rms_fraction_of_median_flux']):
  estimates=[]
  for e,op in enumerate(operators):
   indices=np.array([pair[e] for pair in pairs]);values=amp*op['signal'][indices]+draws[e]*op['err'];estimates.append(estimate(op,values,indices))
  per=[]
  for j,(a,b) in enumerate(pairs):
   dt=float(phase[b]-phase[a]);got=float(estimates[1][j]-estimates[0][j]);per.append(dict(object=str(names[a]),early_phase=float(phase[a]),late_phase=float(phase[b]),recovered_early=float(estimates[0][j]),recovered_late=float(estimates[1][j]),true_span=dt,recovered_span=got,span_error=got-dt,ratio=got/dt))
  ratios=np.array([x['ratio'] for x in per]);errors=np.array([x['span_error'] for x in per])
  row=dict(seed=seed,feature_rms_fraction=amp,median_span_ratio=float(np.median(ratios)),span_error_rms_days=float(np.sqrt(np.mean(errors**2))),nonpositive_spans=int(sum(x['recovered_span']<=0 for x in per)),objects=per);rows.append(row);print(json.dumps({k:v for k,v in row.items() if k!='objects'}),flush=True)
out=dict(scope=p['scope'],protocol=p,n_objects=len(pairs),library_spectra=len(library),additional_blur_sigma_template_pixels=blur_sigma,observations=[{k:v for k,v in op.items() if k in ['pixels','scale','median_noise']} for op in operators],rows=rows,hashes=hashes)
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
