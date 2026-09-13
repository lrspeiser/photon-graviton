"""Axisymmetric attenuated capture, multipole force, and shared mixture tests."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_legendre
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar

P=Path(__file__).resolve().parent;ROOT=P.parents[2];OLD=P.parent/'isotropic-galaxy-transfer'
args=argparse.ArgumentParser();args.add_argument('--refine',action='store_true');args.add_argument('--limit',type=int);aargs=args.parse_args()
files=[OLD/'third-radiation-retention-results.json',OLD/'model-comparison-predictions.json',OLD/'model-comparison-results.json',ROOT/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt',P/'directional-capture-protocol.md']
cp=json.loads(files[0].read_text())['models']['attenuated'];pred=json.loads(files[1].read_text());control=json.loads(files[2].read_text())
meta={}
for line in files[3].read_text().splitlines():
 f=line.split()
 if len(f)==19:
  try:meta[f[0]]=(float(f[11]),float(f[7]))
  except ValueError:pass
data=[]
for v in pred:
 if v['model']!='companion_third':continue
 b=next(x for x in pred if x['model']=='baryons' and x['galaxy']==v['galaxy']);rd,L=meta[v['galaxy']];eta=(L/rd**2)**(1/3)/(1+(L/rd**2)**(1/3))
 data.append(dict(name=v['galaxy'],split=v['split'],r=np.array(v['R_kpc']),y=np.array(v['observed_kms']),vb2=np.array(b['predicted_kms'])**2,ref=np.array(v['predicted_kms']),rd=rd,a=rd*cp['scale_to_disk'],C=2*cp['C_Msun_kpc3']*eta))
if aargs.limit:data=data[:aargs.limit]
nr,nt,nv,nphi,lmax=(768,48,32,64,24) if aargs.refine else (384,24,16,32,12)
x=np.geomspace(1e-6,3000,nr);u,uw=leggauss(nt);v,vw=leggauss(nv)
# Exact reflection symmetry halves azimuth and position evaluations.
nh=nphi//2;phi=np.pi*(np.arange(nh)+.5)/nh
nz=np.repeat(v,nh);nx=(np.sqrt(1-v*v)[:,None]*np.cos(phi)).ravel();iw=np.repeat(vw/(2*nh),nh)
W=np.array([np.ones(len(nz)),3*nz*nz,1.5*(1-nz*nz)])
assert np.allclose(W@iw,1)
ls=np.arange(0,lmax+1,2);basis=np.array([eval_legendre(l,u)*(2*l+1)*uw/2 for l in ls]);p0=np.array([eval_legendre(l,0) for l in ls])
G=4.30091727003628e-6

def density(T,q,which):
 up=u[nt//2:]
 A=1+nz*nz*(q**-2-1);dot=np.sqrt(1-up*up)[:,None]*nx+up[:,None]*nz/q**2
 angular_shape=1+up*up*(q**-2-1);result=np.empty((len(which),nr,nt//2))
 for j in range(0,nr,16):
  xx=x[j:j+16,None,None];b=xx*dot;c=1+xx*xx*angular_shape[None,:,None];h=np.maximum(c-b*b/A,1.)
  y=b/np.sqrt(A*h);integral=(np.arctan2(np.ones_like(y),y)-y/(1+y*y))/(2*np.sqrt(A)*h**1.5)
  tau=T*np.maximum(integral,0)
  for z,k in enumerate(which):result[z,j:j+16]=(np.exp(-tau*W[k])*(W[k]*iw)).sum(axis=2)/(1+x[j:j+16,None]**2*angular_shape)**2
 return np.concatenate([result[:,:,::-1],result],axis=2)

def field(rho):
 coeff=rho@basis.T;inside=np.zeros_like(coeff);outside=np.zeros_like(coeff)
 for i in range(1,nr):
  ratio=(x[i-1]/x[i])**(ls+1)
  inside[i]=ratio*inside[i-1]+.5*(x[i]-x[i-1])*(coeff[i-1]*x[i-1]*ratio+coeff[i]*x[i])
 for i in range(nr-2,-1,-1):
  ratio=(x[i]/x[i+1])**ls
  outside[i]=ratio*outside[i+1]+.5*(x[i+1]-x[i])*(coeff[i]*x[i]+coeff[i+1]*x[i+1]*ratio)
 force=np.sum(((ls+1)*inside-ls*outside)*p0/(2*ls+1),axis=1)
 return PchipInterpolator(np.log(x),force),float(inside[-1,0]*x[-1])

# Known spherical Plummer source supplies an independent field-shape control.
test=(1+x*x)**(-2.5);ff,mm=field(np.broadcast_to(test[:,None],(nr,nt)))
use=(x>.001)&(x<100);exact=x*x/(3*(1+x*x)**1.5)
plummer=float(np.max(abs(ff(np.log(x[use]))/exact[use]-1)));assert plummer<.004,plummer
names=['flat_isotropic','flat_vertical','flat_sides','thin_isotropic'];cache=[]
mu,mw=leggauss(192 if aargs.refine else 96)
for i,d in enumerate(data):
 T=cp['k0_per_kpc']*d['a'];xx=x[:,None];tt=xx*mu;bb=1+xx*xx*(1-mu*mu)
 tau=T*(tt/(2*bb*(bb+tt*tt))+(np.arctan(tt/np.sqrt(bb))+np.pi/2)/(2*bb**1.5))
 rho_ref=(np.exp(-np.maximum(tau,0))@mw/2)/(1+x*x)**2
 refmass=float(np.trapezoid(x*x*rho_ref,x))
 reflut=PchipInterpolator(np.log(x),cumulative_trapezoid(x*x*rho_ref,x,initial=0)/x)
 factor=4*np.pi*G*d['C']*d['a']**2;lr=np.log(d['r']/d['a'])
 calc_ref=np.sqrt(d['vb2']+factor*reflut(lr));ans=dict(galaxy=d['name'],split=d['split'],reference_max_difference_kms=float(np.max(abs(calc_ref-d['ref']))),variants={})
 shapes=list(density(T,.5,[0,1,2]))+list(density(T,.25,[0]))
 for name,rho,q in zip(names,shapes,[.5,.5,.5,.25]):
  lut,mass=field(rho);vv=factor*lut(lr);ratio=mass/refmass
  # Positive density <= ellipsoidal f because <W>=1; upper exterior mass bound.
  shape_avg=float(np.sum(uw/(1+u*u*(q**-2-1))**2)/2)
  tail=shape_avg/x[-1]/mass
  ans['variants'][name]=dict(raw_vcomp2=vv.tolist(),fixed_mass_vcomp2=(vv/ratio).tolist(),mass_ratio=ratio,tail_upper_fraction_of_grid_mass=tail,
      equator_to_pole_density_at_a=float(np.interp(1,x,rho[:,nt//2])/np.interp(1,x,rho[:,-1])))
 cache.append(ans)
 if (i+1)%10==0 or i+1==len(data):print('Geometry',i+1,'/',len(data),flush=True)
out=dict(resolution=dict(nr=nr,nt=nt,nv=nv,nphi=nphi,lmax=lmax),plummer_max_relative_force_error=plummer,geometry=cache,source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
if aargs.limit:
 print('Pilot',plummer,max(z['reference_max_difference_kms'] for z in cache));print({n:cache[0]['variants'][n]['mass_ratio'] for n in names});raise SystemExit

def scores(curves):
 return {s:dict(n=sum(d['split']==s for d in data),RMSE_kms=float(np.sqrt(np.mean([np.mean((v-d['y'])**2) for d,v in zip(data,curves) if d['split']==s]))),log_RMS=float(np.sqrt(np.mean([np.mean(np.log10(v/d['y'])**2) for d,v in zip(data,curves) if d['split']==s])))) for s in ['train','validation','test']}
def makecurves(f,extra):
 vv=[d['ref']**2+f*(e-(d['ref']**2-d['vb2'])) for d,e in zip(data,extra)]
 if any(np.any(v<=0) for v in vv):return None
 return [np.sqrt(v) for v in vv]
def loss(f,extra):
 curves=makecurves(f,extra)
 if curves is None:return 1e30
 return float(np.mean([np.mean(np.log10(v/d['y'])**2) for d,v in zip(data,curves) if d['split']=='train']))
coarse_file=P/'directional-capture-results.json'
frozen=json.loads(coarse_file.read_text()) if aargs.refine else None
out['reference_scores']=scores([d['ref'] for d in data]);out['branches']=[]
for name in names:
 for mode,key in [('raw','raw_vcomp2'),('fixed_mass','fixed_mass_vcomp2')]:
  extra=[np.array(g['variants'][name][key]) for g in cache];label=name+'_'+mode
  if aargs.refine:f=next(b['mixture'] for b in frozen['branches'] if b['label']==label);scan=None
  else:
   grid=np.linspace(0,1,101);ys=[loss(f,extra) for f in grid];candidates=[(0.,ys[0]),(1.,ys[-1])]
   intervals=[(0,.01),(.99,1)]+[(grid[i-1],grid[i+1]) for i in range(1,100) if ys[i]<=ys[i-1] and ys[i]<=ys[i+1]]
   for lo,hi in intervals:
    opt=minimize_scalar(lambda f:loss(f,extra),bounds=(lo,hi),method='bounded',options={'xatol':1e-8});assert opt.success;candidates.append((float(opt.x),float(opt.fun)))
   f=min(candidates,key=lambda x:x[1])[0];scan=[dict(mixture=float(g),loss=y) for g,y in zip(grid,ys)]
  curves=makecurves(f,extra);assert curves is not None
  radial=[]
  for k in range(3):
   ee=[v[np.digitize(d['r']/d['rd'],[1,3])==k]-d['y'][np.digitize(d['r']/d['rd'],[1,3])==k] for d,v in zip(data,curves) if np.any(np.digitize(d['r']/d['rd'],[1,3])==k)]
   radial.append(dict(bin=k,galaxies=len(ee),bias_kms=float(np.mean([e.mean() for e in ee])),RMSE_kms=float(np.sqrt(np.mean([np.mean(e*e) for e in ee])))))
  branch=dict(label=label,mixture=f,scores=scores(curves),radial=radial,scan=scan,predictions=[dict(galaxy=d['name'],split=d['split'],R_kpc=d['r'].tolist(),observed_kms=d['y'].tolist(),predicted_kms=v.tolist()) for d,v in zip(data,curves)])
  if aargs.refine:
   previous=next(b for b in frozen['branches'] if b['label']==label)
   branch['max_prediction_refinement_change_kms']=max(float(np.max(abs(v-np.array(p['predicted_kms'])))) for v,p in zip(curves,previous['predictions']))
  out['branches'].append(branch);print(label,round(f,5),branch['scores'],flush=True)
out['scope']='Axisymmetric opacity hypotheses; one shared trained mixture per branch; exposed data; no absolute supply, support, or lensing validation.'
target=P/('directional-capture-refined.json' if aargs.refine else 'directional-capture-results.json')
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
