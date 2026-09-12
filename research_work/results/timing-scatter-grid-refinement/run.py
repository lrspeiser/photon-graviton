"""Recompute interleaved width nodes for two exposed artificial samples."""
from pathlib import Path
import importlib.util,json,hashlib,sys
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'timing-population'))
from estimator import event_log_likelihood
sys.path.insert(0,str(HERE.parent/'timing-continuous-scatter'))
from integration import log_likelihood
OUT=ROOT/'research_work/generated/timing-scatter-grid-refinement';OUT.mkdir(parents=True,exist_ok=True)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
p=json.loads((HERE.parent/'timing-coverage-calibration/protocol.json').read_text(encoding='utf-8'))
c=p['calibration'];cases=json.loads((HERE.parent/'timing-boundary-diagnosis/protocol.json').read_text(encoding='utf-8'))['cases']
spec=importlib.util.spec_from_file_location('injection',HERE.parent/'timing-injection/run_injections.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
audit=json.loads((ROOT/p['input_audit']).read_text(encoding='utf-8'));paths={Path(v['path']).name:ROOT/v['cache_relative_path'] for v in audit['input_manifest']}
for v in audit['input_manifest']:assert sha(ROOT/v['cache_relative_path'])==v['sha256']
hd=fits.getdata(paths['DES-SN5YR_DES_HEAD.FITS.gz'],1);pd=fits.getdata(paths['DES-SN5YR_DES_PHOT.FITS.gz'],1)
heads={str(v['SNID']).strip():v for v in hd};inventory=json.loads((ROOT/p['source_inventory']).read_text(encoding='utf-8'))['runs'][0]['events'];slots=[]
for item in inventory:
 h=heads[item['id']];obs=pd[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])]
 t=np.asarray(obs['MJD'])-float(h['PEAKMJD']);e=np.asarray(obs['FLUXCALERR'],float)
 keep=(np.char.strip(obs['BAND'].astype(str))==item['band'])&(abs(t)<=80)&np.isfinite(t)&np.isfinite(e)&(e>0)
 t=t[keep];e=e[keep];order=np.argsort(t)
 slots.append(dict(id=item['id'],z=item['z'],t=t[order],e=e[order]/np.median(e)))
assert len(slots)==98
results=[]
for case in cases:
 label=case['label'];seed=int(label.split('seed')[-1]);b=case['truth_b'];base=ROOT/'research_work/generated/timing-coverage-calibration'/(label+'.npz')
 assert sha(base)==case['array_sha256']
 with np.load(base) as data:x=data['log_width'];base_logs=data['event_log_likelihoods'];z=data['redshift'];ids=data['ids']
 assert np.array_equal(ids,np.array([s['id'] for s in slots]))
 fine=np.empty(2*len(x)-1);fine[::2]=x;fine[1::2]=(x[:-1]+x[1:])/2
 rng=np.random.default_rng(seed);newlogs=[];reconstruction_error=None
 for i,slot in enumerate(slots):
  scale=np.exp(rng.normal(0,c['intrinsic_log_scatter']))*(1+slot['z'])**b;offset=rng.uniform(-5,5);e=slot['e']/5
  flux=old.injected(slot['t'],scale,offset,'split_gaussian')+rng.normal(0,e)
  if i==0:
   check=event_log_likelihood(slot['t'],flux,e,x,c['quadrature'],p);reconstruction_error=float(np.max(abs(check-base_logs[0])))
   assert reconstruction_error<1e-10
  newlogs.append(event_log_likelihood(slot['t'],flux,e,fine[1::2],c['quadrature'],p))
  if (i+1)%25==0:print(json.dumps(dict(case=label,events=i+1,total=98)),flush=True)
 combined=np.empty((98,len(fine)));combined[:,::2]=base_logs;combined[:,1::2]=np.array(newlogs)
 array=OUT/(label+'.npz');np.savez_compressed(array,log_width=fine,event_log_likelihoods=combined,redshift=z,ids=ids)
 profiles=[]
 for step in [2,1]:
  logs=combined[:,::step];logs=logs-logs.max(axis=1)[:,None]
  for sigma in [0.,.01,.03]:
   objective=lambda q:-log_likelihood(fine[::step],logs,z,q[0],q[1],sigma)
   starts=[[case['fit']['a'],case['fit']['b']],[np.log(30),0],[np.log(30),1]]
   fits=[minimize(objective,q,method=method,bounds=[np.log([5,80]),[-.5,1.5]],options={'maxiter':1000}) for method in ['L-BFGS-B','Powell'] for q in starts]
   valid=[f for f in fits if f.success and np.isfinite(f.fun)];assert valid
   best=min(valid,key=lambda f:f.fun)
   profiles.append(dict(nodes=len(fine[::step]),sigma=sigma,a=float(best.x[0]),b=float(best.x[1]),successful_starts=len(valid),trials=[dict(success=bool(f.success),objective=float(f.fun),a=float(f.x[0]),b=float(f.x[1])) for f in fits]))
 results.append(dict(label=label,input_sha256=sha(base),refined_array_sha256=sha(array),first_event_reconstruction_error=reconstruction_error,profiles=profiles))
 print(json.dumps(dict(case=label,profiles=[{k:v for k,v in r.items() if k!='trials'} for r in profiles])),flush=True)
(HERE/'results.json').write_text(json.dumps(dict(scope='Exposed-case width discretization check, not coverage or observational validation',cases=results,source_hashes={str(f.relative_to(ROOT)):sha(f) for f in [HERE/'run.py',HERE.parent/'timing-population/estimator.py',HERE.parent/'timing-continuous-scatter/integration.py',HERE.parent/'timing-coverage-calibration/protocol.json',HERE.parent/'timing-injection/run_injections.py']}),indent=2)+'\n',encoding='utf-8',newline='\n')
