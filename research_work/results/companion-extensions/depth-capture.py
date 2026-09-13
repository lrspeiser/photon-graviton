"""Depth-dependent interception, stationary deposited-field feedback, fixed source."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

P=Path(__file__).resolve().parent;ROOT=P.parents[2];OLD=P.parent/'isotropic-galaxy-transfer'
parser=argparse.ArgumentParser();parser.add_argument('--refine',action='store_true');parser.add_argument('--limit',type=int);args=parser.parse_args()
paths=[OLD/'third-radiation-retention-results.json',OLD/'model-comparison-predictions.json',ROOT/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt',P/'depth-capture-protocol.md']
cp=json.loads(paths[0].read_text())['models']['attenuated'];saved=json.loads(paths[1].read_text());meta={}
for line in paths[2].read_text().splitlines():
 f=line.split()
 if len(f)==19:
  try:meta[f[0]]=(float(f[11]),float(f[7]))
  except ValueError:pass
data=[]
for row in saved:
 if row['model']!='companion_third':continue
 b=next(v for v in saved if v['model']=='baryons' and v['galaxy']==row['galaxy']);rd,L=meta[row['galaxy']];X=L/rd**2;eta=X**(1/3)/(1+X**(1/3))
 data.append(dict(name=row['galaxy'],split=row['split'],r=np.array(row['R_kpc']),y=np.array(row['observed_kms']),vb2=np.array(b['predicted_kms'])**2,ref=np.array(row['predicted_kms']),a=rd*cp['scale_to_disk'],rd=rd,C=2*cp['C_Msun_kpc3']*eta))
if args.limit:data=data[:args.limit]
nr,na,ns=(640,32,48) if args.refine else (320,16,24)
x=np.geomspace(1e-6,3000,nr);lx=np.log(x);mu,mw=leggauss(na);nu,nw=leggauss(ns)
impact=x[:,None]*np.sqrt(1-mu*mu);Q=np.sqrt(1+impact*impact);theta0=np.arctan(x[:,None]*mu/Q);span=np.pi/2-theta0
theta=theta0[:,:,None]+span[:,:,None]*(nu+1)/2
path=np.sqrt(impact[:,:,None]**2+(Q[:,:,None]*np.tan(theta))**2)
weights=span[:,:,None]*nw*np.cos(theta)**2/(2*Q[:,:,None]**3)
analytic=(np.pi/2-theta0-.5*np.sin(2*theta0))/(2*Q**3)
ray_error=float(np.max(abs(weights.sum(axis=2)-analytic)/np.maximum(analytic,1e-20)));assert ray_error<1e-8,ray_error
lp=np.log(np.maximum(path,1e-300));idx=np.clip(np.searchsorted(lx,lp)-1,0,nr-2);frac=np.clip((lp-lx[idx])/(lx[idx+1]-lx[idx]),0,1)
f0=(1+x*x)**-2;G=4.30091727003628e-6

def baryonic_depth(d,rr,tail=1.):
 r=d['r'];v=d['vb2'];lr=np.log(r);segments=.5*(v[:-1]+v[1:])*np.diff(lr)
 edge=np.r_[np.cumsum(segments[::-1])[::-1],0.]+tail*v[-1]
 z=np.log(np.maximum(rr,1e-300));k=np.clip(np.searchsorted(lr,z)-1,0,len(r)-2)
 h=np.clip(z-lr[k],0,np.diff(lr)[k]);slope=(v[k+1]-v[k])/np.diff(lr)[k]
 answer=edge[k]-v[k]*h-.5*slope*h*h
 answer=np.where(rr<r[0],edge[0]+.5*v[0]*(1-(rr/r[0])**2),answer)
 answer=np.where(rr>r[-1],tail*v[-1]*r[-1]/rr,answer)
 return answer

def grav(rho):
 mass=cumulative_trapezoid(x*x*rho,x,initial=0)
 outer=-cumulative_trapezoid((x*rho)[::-1],x[::-1],initial=0)[::-1]
 return mass,mass/x+outer
def mod(B,family,beta):
 S=B/(B+150**2)
 return 1+beta*S if family=='boost' else 1-beta+beta*S if family=='gate' else 1+beta*(2*S-1)

def setup(d,tail=1.):
 T=cp['k0_per_kpc']*d['a'];reference=f0*(np.exp(-T*analytic)@mw/2)
 return dict(T=T,reference=reference,Bb=baryonic_depth(d,d['a']*x,tail),Bbpath=baryonic_depth(d,d['a']*path,tail),factor=4*np.pi*G*d['C']*d['a']**2)

def solve(d,inp,family,beta,seed=1.,feedback=True):
 rho=inp['reference']*seed;residual=np.inf
 for it in range(100):
  mass,B=grav(rho if feedback else inp['reference']);Bp=B[idx]*(1-frac)+B[idx+1]*frac
  Bp=np.where(path>x[-1],mass[-1]/path,Bp)
  modifier=mod(inp['Bb']+inp['factor']*B,family,beta)
  raymod=mod(inp['Bbpath']+inp['factor']*Bp,family,beta)
  tau=inp['T']*np.sum(weights*raymod,axis=2)
  target=f0*modifier*(np.exp(-tau)@mw/2)
  residual=float(np.max(abs(target-rho))/max(np.max(target),1e-30))
  if residual<2e-7:rho=target;break
  rho=.5*rho+.5*target
 mass,B=grav(rho);base_mass,_=grav(inp['reference']);r=d['r']/d['a']
 v=np.sqrt(d['vb2']+inp['factor']*PchipInterpolator(lx,mass/x)(np.log(r)))
 return dict(predicted_kms=v.tolist(),converged=bool(residual<2e-7),iterations=it+1,residual=residual,mass_ratio=float(mass[-1]/base_mass[-1]),central_depth_kms2=float(inp['Bb'][0]+inp['factor']*B[0]))

coarse=json.loads((P/'depth-capture-results.json').read_text()) if args.refine else None
cases=[(family,beta) for family in ['boost','gate','tilt'] for beta in [.25,.5,1.]] if not args.refine else [(s['family'],s['beta']) for s in coarse['selected'] if s['beta']>0]
# If the identity wins every family, verify the weak gate's metric-dependent gain;
# this is a numerical diagnostic, not a replacement selection.
if args.refine and not cases:cases=[('gate',.25)]
results={f'{family}_{beta}':[] for family,beta in cases};baseline=[]
for i,d in enumerate(data):
 inp=setup(d);r0=solve(d,inp,'boost',0);baseline.append(dict(galaxy=d['name'],split=d['split'],**r0))
 for family,beta in cases:
  row=solve(d,inp,family,beta)
  if args.refine:
   alt=solve(d,inp,family,beta,seed=.1)
   row['alternate_seed_max_speed_difference']=float(np.max(abs(np.array(row['predicted_kms'])-alt['predicted_kms'])))
   row['alternate_seed_converged']=alt['converged']
  results[f'{family}_{beta}'].append(dict(galaxy=d['name'],split=d['split'],**row))
 if (i+1)%10==0 or i+1==len(data):print('Depth galaxies',i+1,'/',len(data),flush=True)

def scores(rows):
 return {split:dict(n=sum(d['split']==split for d in data),RMSE_kms=float(np.sqrt(np.mean([np.mean((np.array(row['predicted_kms'])-d['y'])**2) for d,row in zip(data,rows) if d['split']==split]))),log_RMS=float(np.sqrt(np.mean([np.mean(np.log10(np.array(row['predicted_kms'])/d['y'])**2) for d,row in zip(data,rows) if d['split']==split])))) for split in ['train','validation','test'] if any(d['split']==split for d in data)}
branches=[]
for family,beta in cases:
 rows=results[f'{family}_{beta}'];branches.append(dict(family=family,beta=beta,rows=rows,scores=scores(rows),all_converged=all(r['converged'] for r in rows)))
selected=[]
if not args.refine:
 for family in ['boost','gate','tilt']:
  candidates=[dict(family=family,beta=0.,scores=scores(baseline))]+[b for b in branches if b['family']==family and b['all_converged']]
  best=min(candidates,key=lambda x:x['scores']['train']['log_RMS']);selected.append(dict(family=family,beta=best['beta'],scores=best['scores']))
out=dict(scope='Spherical capture with sphericalized midplane baryonic depth, fixed 150 km/s saturation scale and stationary deposited-field feedback; no full 3D self-consistency or stability proof.',resolution=dict(nr=nr,na=na,ns=ns),ray_integral_relative_error=ray_error,baseline=baseline,baseline_scores=scores(baseline),baseline_max_archived_speed_difference=float(max(np.max(abs(np.array(r['predicted_kms'])-d['ref'])) for r,d in zip(baseline,data))),branches=branches,selected=selected,source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if args.limit:print('Pilot',out['baseline_max_archived_speed_difference'],[(b['family'],b['beta'],b['all_converged'],max(r['iterations'] for r in b['rows'])) for b in branches]);raise SystemExit
target=P/('depth-capture-refined.json' if args.refine else 'depth-capture-results.json');target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print('Selected',selected);print('reference',out['baseline_scores']);print('all cases',[(b['family'],b['beta'],b['scores'],b['all_converged']) for b in branches])
