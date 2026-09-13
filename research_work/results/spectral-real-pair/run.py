from pathlib import Path
import json,hashlib
import numpy as np
from astropy.io import fits
H=Path(__file__).resolve().parent;ROOT=H.parents[2]
proto=json.loads((H/'protocol.json').read_text(encoding='utf-8-sig'))
loader=H.parent/'spectral-observed-noise/run.py';s=loader.read_text(encoding='utf-8');assert s.count(chr(10)+'rows=[]'+chr(10))==1
ns={'__file__':str(loader),'__name__':'operator_loader'};exec(compile(s.split(chr(10)+'rows=[]'+chr(10))[0],str(loader),'exec'),ns)
operators,records,pairs,names,phase,groups,objects,estimate,audit=[ns[k] for k in ['operators','records','pairs','names','phase','groups','objects','estimate','audit']]
assert not any('2006mk' in n or n=='q106' for n in names)
actual=[];powers=[];fitsrows=[]
for epoch,(op,rec) in enumerate(zip(operators,records)):
 path=ROOT/'research_work/generated/essence-dr3'/rec['filename']
 with fits.open(path) as h:
  d=h[1].data[0];rw=np.asarray(d['WAVE'],float)/(1+float(rec['catalog_identity']['z_host']));m=(rw>=4000)&(rw<=5500);rw=rw[m];f=np.asarray(d['FLUX'],float)[m]
 host=next(x for x in audit['spectra'] if x['target']=='q106-gal' and x['mid_mjd']==rec['mid_mjd'])
 hostpath=ROOT/'research_work/generated/essence-dr3'/host['filename'];assert hashlib.sha256(hostpath.read_bytes()).hexdigest()==next(x['sha256'] for x in audit['files'] if x['filename']==host['filename'])
 with fits.open(hostpath) as h:hostflux=np.asarray(h[1].data[0]['FLUX'],float)[m]
 actual.append(f);white=f/op['err'];res=white-op['q']@(op['q'].T@white);power=float(res@res-(op['pixels']-op['q'].shape[1]));assert power>0;powers.append(power)
 for degree in proto['sensitivity_degrees']:
  for host_basis in proto['sensitivity_host_basis']:
   x=(rw-rw.mean())/np.ptp(rw);design=np.column_stack([x**j for j in range(degree+1)])
   if host_basis:design=np.column_stack([design,hostflux/np.median(hostflux)])
   q=np.linalg.qr(design/op['err'][:,None],mode='reduced')[0]
   basis=op['signal']/op['err'];basis-=(basis@q)@q.T;basis/=np.linalg.norm(basis,axis=1)[:,None]
   y=white-q@(q.T@white);y/=np.linalg.norm(y);score=basis@y
   best=[int(ix[np.argmax(score[ix])]) for ix in groups];best=sorted(best,key=lambda j:score[j],reverse=True)[:5]
   fitsrows.append(dict(epoch=epoch,degree=degree,host_basis=host_basis,estimated_phase=float(np.median(phase[best])),matches=[dict(object=str(names[j]),phase=float(phase[j]),correlation=float(score[j])) for j in best]))
calibration=[]
for seed in proto['seeds']:
 for rho in proto['adjacent_pixel_correlations']:
  estimates=[];amp_ranges=[]
  for e,op in enumerate(operators):
   ix=np.array([pair[e] for pair in pairs]);signal=op['signal'][ix];w=signal/op['err'];proj=w-(w@op['q'])@op['q'].T;amps=np.sqrt(powers[e]/np.sum(proj*proj,axis=1));assert np.all(np.isfinite(amps))
   # AR(1) draws have unit marginal variance, retaining each supplied error.
   draw=np.random.default_rng(np.random.SeedSequence([seed,e])).normal(size=signal.shape)
   if rho:
    for k in range(1,draw.shape[1]):draw[:,k]=rho*draw[:,k-1]+np.sqrt(1-rho*rho)*draw[:,k]
   values=amps[:,None]*signal+draw*op['err'];estimates.append(estimate(op,values,ix));amp_ranges.append([float(amps.min()),float(amps.max())])
  per=[dict(object=str(names[a]),true_span=float(phase[b]-phase[a]),estimated_span=float(estimates[1][j]-estimates[0][j])) for j,(a,b) in enumerate(pairs)]
  error=np.array([v['estimated_span']-v['true_span'] for v in per]);calibration.append(dict(seed=seed,rho=rho,objects=per,rms_error_days=float(np.sqrt(np.mean(error**2))),mean_error_days=float(error.mean()),amplitude_ranges=amp_ranges))
spans=[]
for degree in proto['sensitivity_degrees']:
 for host in proto['sensitivity_host_basis']:
  ff=[x for x in fitsrows if x['degree']==degree and x['host_basis']==host];spans.append(dict(degree=degree,host_basis=host,estimated_span=ff[1]['estimated_phase']-ff[0]['estimated_phase']))
summary=[]
for rho in proto['adjacent_pixel_correlations']:
 err=np.array([v['estimated_span']-v['true_span'] for c in calibration if c['rho']==rho for v in c['objects']]);summary.append(dict(rho=rho,n_reused_pairs=len(err),mean_error_days=float(err.mean()),rms_error_days=float(np.sqrt(np.mean(err**2))),error_quantiles_05_50_95=np.quantile(err,[.05,.5,.95]).tolist()))
paths=[Path(__file__),H/'protocol.json',loader,H.parent/'spectral-feature-strength/results.json',H.parent/'essence-spectrum-acquisition/results.json']
out=dict(scope=proto['scope'],protocol=proto,observed_projected_excess_power=powers,real_epoch_fits=fitsrows,real_spans=spans,synthetic_summary=summary,synthetic_cases=calibration,observer_span_days=records[1]['mid_mjd']-records[0]['mid_mjd'],hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'real_spans':spans,'synthetic_summary':summary}))
