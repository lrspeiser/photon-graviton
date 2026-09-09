import json,hashlib,math,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar,brentq
from scipy.special import lambertw
from scipy.integrate import cumulative_trapezoid,quad
P=Path(__file__).resolve().parent
C=299792.458
MPC_PER_MLY=0.3066013938
MODELS=['static_flux','flrw','uniform_proxy','pilot_static_fixed']
ZG=np.linspace(0,.5,100001)
IG=cumulative_trapezoid(1/np.sqrt(.3*(1+ZG)**3+.7),ZG,initial=0)
FG=(1+ZG)*IG

def save(n,x):P.joinpath(n).write_text(json.dumps(x,indent=2,allow_nan=False))
def read(n):return json.loads(P.joinpath(n).read_text())
def sha(n):return hashlib.sha256(P.joinpath(n).read_bytes()).hexdigest()
def prepare():
 old={r['group'] for r in read('pilot_features.json')};rows=[];labels={}
 for l in P.joinpath('table2.dat').read_text().splitlines():
  if not l[77:83].strip() or not l[84:89].strip():continue
  dm=float(l[77:83]);err=float(l[84:89]);d=10**((dm-25)/5);group=int(l[8:15]);pgc=int(l[:7])
  if not (10<=d<=150 and 0<err<=.30 and group not in old):continue
  rows.append(dict(pgc=pgc,group=group,dm=dm,dm_err=err,distance_mpc=d,sd_mpc=d*math.log(10)/5*err,ra=float(l[137:145]),dec=float(l[146:154])))
  labels[str(pgc)]=float(l[22:27])
 groups={}
 for r in sorted(rows,key=lambda r:(r['dm_err'],r['pgc'])):groups.setdefault(r['group'],r)
 rows=sorted(groups.values(),key=lambda r:r['pgc'])
 for r in rows:
  tile=int(r['ra']//15)*4+min(3,int((math.sin(math.radians(r['dec']))+1)*2))
  bucket=int(hashlib.sha256(f'time-redshift-sbf-v2:{tile}'.encode()).hexdigest()[:8],16)%10
  r.update(tile=tile,split='train' if bucket<6 else 'validation' if bucket<8 else 'test')
 save('features.json',rows)
 for s in ['train','validation','test']:save(s+'_labels.json',{str(r['pgc']):labels[str(r['pgc'])] for r in rows if r['split']==s})
 save('provenance.json',dict(source_sha256=sha('table2.dat'),protocol_sha256=sha('protocol.json'),source_rows=len(P.joinpath('table2.dat').read_text().splitlines()),excluded_pilot_groups=sorted(old)))
 print(json.dumps({s:dict(n=sum(r['split']==s for r in rows),tiles=len({r['tile'] for r in rows if r['split']==s})) for s in ['train','validation','test']},indent=2))
def pred(m,k,rows):
 d=np.array([r['distance_mpc'] for r in rows]);sd=np.array([r['sd_mpc'] for r in rows]);x=k*d/C
 if m in ['static_flux','pilot_static_fixed']:
  w=lambertw(x).real;v=C*np.expm1(w);der=k/(1+w)
 elif m=='uniform_proxy':v=C*np.expm1(x);der=k*np.exp(x)
 else:
  if np.max(x)>FG[-1]:raise ValueError('FLRW interpolation range exceeded')
  z=np.interp(x,FG,ZG);integ=np.interp(z,ZG,IG)
  v=C*z;der=k/(integ+(1+z)/np.sqrt(.3*(1+z)**3+.7))
 return v,np.sqrt(300**2+(der*sd)**2)
def objective(m,k,rr,y):
 v,s=pred(m,k,rr);return float(np.sum(np.log(s)+.5*((y-v)/s)**2+.5*np.log(2*np.pi)))
def fit(m,rr,y):return minimize_scalar(lambda k:objective(m,k,rr,y),bounds=(10,150),method='bounded')
def train():
 rr=[r for r in read('features.json') if r['split']=='train'];labels=read('train_labels.json');y=np.array([labels[str(r['pgc'])] for r in rr]);out={}
 rng=np.random.default_rng(9102026);tiles=sorted({r['tile'] for r in rr})
 bootids=[[i for t in rng.choice(tiles,len(tiles),replace=True) for i,r in enumerate(rr) if r['tile']==t] for _ in range(1000)]
 for m in MODELS:
  if m=='pilot_static_fixed':
   k=read('pilot_frozen_parameters.json')['static_flux']['k_kms_mpc'];out[m]=dict(k_kms_mpc=k,nll=objective(m,k,rr,y),fitted_on_this_sample=False);continue
  o=fit(m,rr,y);k=float(o.x)
  ci=[brentq(lambda q:objective(m,q,rr,y)-o.fun-.5,10,k),brentq(lambda q:objective(m,q,rr,y)-o.fun-.5,k,150)]
  bk=[float(fit(m,[rr[i] for i in ix],y[ix]).x) for ix in bootids]
  out[m]=dict(k_kms_mpc=k,profile_68_k=ci,bootstrap_95_k=np.percentile(bk,[2.5,97.5]).tolist(),nll=float(o.fun),fitted_on_this_sample=True)
  if m!='flrw':out[m].update(alpha_per_mpc=k/C,alpha_per_million_ly=k/C*MPC_PER_MLY,percent_stretch_per_million_ly=100*math.expm1(k/C*MPC_PER_MLY))
 save('frozen_parameters.json',out);print(json.dumps(out,indent=2))
def predictions():
 out=[]
 for r in read('features.json'):
  item=dict(r)
  for m,f in read('frozen_parameters.json').items():
   v,s=pred(m,f['k_kms_mpc'],[r]);item[m]=dict(predicted_z=float(v[0]/C),sigma_z=float(s[0]/C))
  out.append(item)
 save('frozen_predictions.json',out)
 save('prediction_seal.json',dict(protocol_sha256=sha('protocol.json'),parameters_sha256=sha('frozen_parameters.json'),predictions_sha256=sha('frozen_predictions.json')))
 print('Predictions and hashes saved before held-out scoring.')
def score():
 out={};scored=[]
 for s in ['validation','test']:
  rr=[r for r in read('frozen_predictions.json') if r['split']==s];lab=read(s+'_labels.json');y=np.array([lab[str(r['pgc'])] for r in rr]);stats={};errs={};nlls={}
  for m in MODELS:
   v=C*np.array([r[m]['predicted_z'] for r in rr]);sig=C*np.array([r[m]['sigma_z'] for r in rr]);e=v-y;nll=np.log(sig)+.5*(e/sig)**2+.5*np.log(2*np.pi)
   errs[m]=e;nlls[m]=nll
   stats[m]=dict(rmse_kms=float(np.sqrt(np.mean(e**2))),rmse_z=float(np.sqrt(np.mean(e**2))/C),mae_kms=float(np.mean(abs(e))),bias_kms=float(np.mean(e)),mean_nll=float(np.mean(nll)),coverage_95=float(np.mean(abs(e)<=1.96*sig)))
  out[s]=dict(n=len(rr),tiles=len({r['tile'] for r in rr}),models=stats)
  for r in rr:scored.append(dict(r,observed_z=lab[str(r['pgc'])]/C))
  if s=='test':
   rng=np.random.default_rng(9202026);tiles=sorted({r['tile'] for r in rr});deltas=[];dnll=[]
   for _ in range(5000):
    ix=[i for t in rng.choice(tiles,len(tiles),replace=True) for i,r in enumerate(rr) if r['tile']==t]
    deltas.append(np.sqrt(np.mean(errs['static_flux'][ix]**2))-np.sqrt(np.mean(errs['flrw'][ix]**2)))
    dnll.append(np.mean(nlls['static_flux'][ix]-nlls['flrw'][ix]))
   out[s]['static_minus_flrw_rmse_95']=np.percentile(deltas,[2.5,97.5]).tolist();out[s]['static_minus_flrw_mean_nll_95']=np.percentile(dnll,[2.5,97.5]).tolist()
 save('scored_predictions.json',scored);save('scores.json',out);print(json.dumps(out,indent=2))
def verify():
 rr=read('features.json');assert len(rr)==len({r['group'] for r in rr})
 assert not {r['group'] for r in rr}&set(read('provenance.json')['excluded_pilot_groups'])
 for t in {r['tile'] for r in rr}:assert len({r['split'] for r in rr if r['tile']==t})==1
 for d in [10,30,100,150]:
  row=dict(distance_mpc=d,sd_mpc=1)
  z=pred('static_flux',70,[row])[0][0]/C
  assert abs((1+z)*math.log1p(z)-70*d/C)<1e-12
  z=pred('flrw',70,[row])[0][0]/C
  exact=(1+z)*quad(lambda u:1/math.sqrt(.3*(1+u)**3+.7),0,z)[0]
  assert abs(exact-70*d/C)<1e-10
 seal=read('prediction_seal.json');assert seal['predictions_sha256']==sha('frozen_predictions.json');assert seal['protocol_sha256']==sha('protocol.json');assert seal['parameters_sha256']==sha('frozen_parameters.json')
 save('verification.json',dict(group_and_tile_separation=True,prior_group_exclusion=True,static_flux_identity=True,flrw_quadrature_agreement=True,frozen_hashes_unchanged=True));print('Verified model identities, split independence by group/tile, exclusions, and frozen file hashes.')
if __name__=='__main__':globals()[sys.argv[1]]()
