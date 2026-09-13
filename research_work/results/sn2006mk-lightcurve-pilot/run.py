from pathlib import Path
import numpy as np,json,hashlib
from scipy.interpolate import RegularGridInterpolator
H=Path(__file__).resolve().parent;R=H.parents[2]
paths=[H.parent/'spectral-real-pair/results.json',H.parent/'sn2006mk-joint-photometry/results.json',H.parent/'essence-spectrum-acquisition/results.json',R/'research_work/generated/des-photometric-calibration/selected/snsed/Hsiao07.dat']
spec,phot,audit=[json.loads(p.read_text(encoding='utf-8')) for p in paths[:3]]
d=np.loadtxt(paths[3]);t=np.unique(d[:,0]);w=np.unique(d[:,1]);interp=RegularGridInterpolator((t,w),d[:,2].reshape(len(t),len(w)),bounds_error=True)
assert hashlib.sha256(paths[3].read_bytes()).hexdigest()==json.loads((H.parent/'hsiao-source-audit/results.json').read_text())['cached_sha256']
epoch=min(x['mid_mjd'] for x in audit['spectra'] if x['target']=='q106');fits=sorted([x for x in spec['real_epoch_fits'] if x['degree']==2 and not x['host_basis']],key=lambda x:x['epoch']);start=fits[0]['estimated_phase']
stretches={'spectral_pair':spec['observer_span_days']/(fits[1]['estimated_phase']-start),'event_stretch_1_plus_z':1.4754,'unchanged_arrival_intervals':1.}
def phase(date,A):return start+(date-epoch)/A
# Same eligibility for all timing hypotheses, derived only from source support.
rows=phot['photometry'];eligible=[x for x in rows if all(-19<=phase(x['mjd'],A)<=85 for A in stretches.values())]
excluded=[dict(x,reason='Outside common source-phase support [-19,85]',phases={k:phase(x['mjd'],A) for k,A in stretches.items()}) for x in rows if x not in eligible]
results=[]
for band in ['R','I']:
 path=R/'research_work/generated/essence-passbands'/f'CTIO4m_{band}.dat';paths.append(path);data=np.loadtxt(path);wo,energy=data.T
 curve=sorted([x for x in eligible if x['band']==band+'4m'],key=lambda x:x['mjd']);anchor=min(curve,key=lambda x:abs(x['mjd']-epoch));others=[x for x in curve if x!=anchor]
 assert all(x['error_upper']==x['error_lower'] for x in curve)
 for name,A in stretches.items():
  def bandflux(date):
   p=phase(date,A);f=interp(np.column_stack([np.full(len(wo),p),wo/1.4754]));return float(np.trapezoid(f*energy,wo))
  norm=bandflux(anchor['mjd']);ratio=np.array([bandflux(x['mjd'])/norm for x in others]);y=np.array([x['flux'] for x in others]);errors=np.array([x['error_upper'] for x in others]);pred=anchor['flux']*ratio;res=y-pred
  C=np.diag(errors**2)+anchor['error_upper']**2*np.outer(ratio,ratio)
  chi=float(res@np.linalg.solve(C,res));u=anchor['error_upper']*ratio;diag=errors**2
  alternate=float(np.sum(res**2/diag)-np.sum(res*u/diag)**2/(1+np.sum(u*u/diag)));assert abs(chi-alternate)<1e-8
  results.append(dict(band=band,case=name,event_stretch=A,anchor=anchor,n_scored=len(others),chi_square_conditional=chi,chi_square_independent_check=alternate,predictions=[dict(observation=x,predicted_flux=float(p),source_phase=phase(x['mjd'],A),residual=float(e)) for x,p,e in zip(others,pred,res)],residual_covariance=C.tolist()))
out=dict(scope='Conditional full in-support light-curve shape prediction anchored once per band; not absolute brightness or calibrated model selection',common_phase_support=[-19,85],input_rows=len(rows),eligible_rows=len(eligible),excluded_rows=excluded,cases=stretches,results=results,hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[Path(__file__)]})
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({'eligible':len(eligible),'excluded':len(excluded),'results':[{k:v for k,v in x.items() if k in ['band','case','n_scored','chi_square_conditional']} for x in results]}))
