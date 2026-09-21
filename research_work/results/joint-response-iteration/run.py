#!/usr/bin/env python3
"""JR-1: exploratory, source-linked, finite companion-reservoir fitting.

Run inside the photon-graviton checkout. Fresh output directories preserve old
results. --resume may add predeclared stages; it must not be used to choose new
models based on already exposed transfer results. No claim of microscopic
production or observational validation is made by optimizer convergence.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, io, json, time, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares, brentq
from scipy.special import expit
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
G=4.30091727003628e-6; CLIGHT=299792.458; KPC_M=3.0856775814913673e19
ARCSEC=206264.80624709636; MREF=1e10
AREF=np.sqrt(G*MREF*6.54e-11*KPC_M/1e6)
TRAIN_LENSES=('J0037-0942','J1112+0826','J1204+0358','J1402+6321')
TRANSFER_LENSES=('J1621+3931','J1630+4520'); LENSES=TRAIN_LENSES+TRANSFER_LENSES
INPUT_NAMES=[
 'temporal_candidate_audit/data/Rotmod_LTG.zip',
 'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt',
 'research_work/results/companion-extensions/mond-inventory-results.json',
 'research_work/results/lensing-data-readiness/lens-observations-and-image-models.json',
 'research_work/results/lensing-data-readiness/conditional-geometry.json',
 'research_work/results/slacs-light-profile-audit/results.json',
 'research_work/results/slacs-resolved-input-audit/results.json',
 'research_work/results/lens-photometric-audit/normalization-sensitivity.json',
 'research_work/results/slacs-component-refit/model.py',
 'research_work/results/slacs-resolved-fit/model.py']
DEFAULTS=dict(logA=float(np.log10(AREF)),logc=0.,p=.5,dA=0.,dc=0.,q=2.,
 logt=float(np.log10(20.)),logu=0.,beta=0.,logP0=-2.,betainf=0.,logra=0.,logsph=0.,core_sph=0.)
BOUNDS=dict(logA=(2.5,5.5),logc=(-1.5,1.5),p=(.15,1.),dA=(-1.,1.),dc=(-.8,.8),
 q=(1.2,4.),logt=(.7,2.5),logu=(-.3,.4),beta=(-1.,.45),logP0=(-5.,5.),
 betainf=(-1.,.65),logra=(-1.,1.),logsph=(-1.,2.),core_sph=(-2.,2.))
def nuisance_key(prefix,name): return prefix+'_'+name.replace('-','_').replace('+','p')
DEFAULTS['q_sph']=2.
BOUNDS['q_sph']=(.3,4.)
for _name in TRAIN_LENSES:
 DEFAULTS[nuisance_key('dm',_name)]=0.;BOUNDS[nuisance_key('dm',_name)]=(-.2,.2)
 DEFAULTS[nuisance_key('b',_name)]=0.;BOUNDS[nuisance_key('b',_name)]=(-1.,.35)
SPECS=[
 dict(name='R0_initial_root_reservoir',free=[],saturation=False,radial_beta=False),
 dict(name='R1_amplitude_core',free=['logA','logc'],saturation=False,radial_beta=False),
 dict(name='R2_source_compactness',free=['logA','logc','p','dA'],saturation=False,radial_beta=False),
 dict(name='R3_spatial_shape',free=['logA','logc','p','dA','dc','q','logt'],saturation=False,radial_beta=False),
 dict(name='R4_shared_stellar_orbits',free=['logA','logc','p','dA','dc','q','logt','logu','beta'],saturation=False,radial_beta=False),
 dict(name='R5_production_release',free=['logA','logc','p','dA','dc','q','logt','logu','beta','logP0'],saturation=True,radial_beta=False),
 dict(name='R6_radial_orbits',free=['logA','logc','p','dA','dc','q','logt','logu','beta','logP0','betainf','logra'],saturation=True,radial_beta=True),
 dict(name='R7_geometry_production',free=['logA','logc','p','dA','dc','q','logt','logu','beta','logsph','core_sph'],saturation=False,radial_beta=False),
 dict(name='R8_geometry_rate_orbits',free=['logA','logc','p','dA','dc','q','logt','logu','beta','logP0','betainf','logra','logsph','core_sph'],saturation=True,radial_beta=True)]

SPECS.append(dict(name='R9_spheroid_shape_stellar_priors',
 free=['logA','logc','p','dA','dc','q','logt','logu','logsph','core_sph','q_sph']+
 [nuisance_key(kind,name) for kind in ('dm','b') for name in TRAIN_LENSES],
 saturation=False,radial_beta=False,separate_q=True,stellar_nuisance=True))

def read_json(relative): return json.loads((ROOT/relative).read_text())

def stellar_disk_mass(rows,rd):
 """Independent surface-density integral; no observed velocity or force mass."""
 use=(rows[:,0]>0)&(rows[:,6]>0)&np.isfinite(rows[:,6]); r=rows[use,0]; s=.5e6*rows[use,6]
 if len(r)==0: raise ValueError('Missing positive stellar surface profile')
 x,w=np.polynomial.legendre.leggauss(24); m=np.pi*r[0]**2*s[0]
 if len(r)>1:
  v=r[:-1,None]+np.diff(r)[:,None]*(x+1)/2
  ln=np.log(s[:-1,None])+(np.log(s[1:])-np.log(s[:-1]))[:,None]*(x+1)/2
  m+=float(np.sum(np.pi*np.diff(r)[:,None]*w*v*np.exp(ln)))
 m+=2*np.pi*s[-1]*rd*(r[-1]+rd)
 return float(m)

def load_sparc():
 guide=read_json(INPUT_NAMES[2])['rows']; cat={}; out=[]
 for line in (ROOT/INPUT_NAMES[1]).read_text().splitlines():
  f=line.split()
  if len(f)==19:
   try:cat[f[0]]=dict(rd=float(f[11]),Re=float(f[9]),L=float(f[7])*1e9,Mgas=1.33*float(f[13])*1e9)
   except ValueError:pass
 with zipfile.ZipFile(ROOT/INPUT_NAMES[0]) as z:
  for item in guide:
   name=item['galaxy']; c=cat[name]
   d=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(name+'_rotmod.dat'))));d=d[d[:,0]>0]
   if np.any(d[:,2]<=0) or c['Re']<=0:raise ValueError('Invalid observational errors/size: '+name)
   md=stellar_disk_mass(d,c['rd']);mb=float(np.max(.7*d[:,0]*d[:,5]**2/G))
   out.append(dict(name=name,split=item['split'],r=d[:,0],observed=d[:,1],error=d[:,2],
    star_v2=.5*d[:,4]**2+.7*d[:,5]**2,gas_v2=d[:,3]*abs(d[:,3]),Mstar=md+mb,Mgas=c['Mgas'],
    Re=c['Re'],L36=c['L'],Mdisk=md,Mbulge=mb))
 assert len(out)==149
 return out

def source_parameters(Mstar,Mgas,Re,p,saturation=False,spheroid_fraction=0.):
 true_mass=10**p['logu']*np.asarray(Mstar)+np.asarray(Mgas)
 S=10**p['logu']*np.asarray(Mstar)*(1+(10**p.get('logsph',0.)-1)*np.asarray(spheroid_fraction))+np.asarray(Mgas)
 s=S/MREF;compactness=(true_mass/MREF)/np.asarray(Re)**2
 if saturation:
  x=s**(2*p['p']);t=10**p['logP0'];occupation=x*(np.sqrt(1+4/t)+1)/(np.sqrt(1+4*x/t)+1)
 else:occupation=s**p['p']
 A=10**p['logA']*occupation*compactness**p['dA']
 rc=10**p['logc']*np.asarray(Re)*compactness**p['dc']*10**(-p.get('core_sph',0.)*np.asarray(spheroid_fraction));rt=10**p['logt']*np.asarray(Re)
 return A,rc,rt

def companion_force(r,A,rc,rt,q):
 r=np.asarray(r)
 if np.any(r<=0):raise ValueError('Radius must be positive')
 return A/r*expit(q*np.log(r/rc))/np.hypot(1,r/rt)

class Lens:
 def __init__(self,name,n=1601,order=64,deproj_order=128):
  self.name=name;self.training=name in TRAIN_LENSES
  obs={x['Name']:x for x in read_json(INPUT_NAMES[3])}[name]
  geo={x['Name']:x for x in read_json(INPUT_NAMES[4])}[name]
  light={x['Name']:x for x in read_json(INPUT_NAMES[5])['rows']}[name]
  kin={x['Name']:x for x in read_json(INPUT_NAMES[6])['systems']}[name]
  masses=[x for x in read_json(INPUT_NAMES[7]) if x['Name']==name and x['imf']=='Chabrier' and x['propagation_branch']=='energy_loss_and_event_stretch']
  assert masses and float(kin['release_use_flag'])==1.
  assert np.ptp([x['conditional_log10_stellar_mass'] for x in masses])<1e-12
  self.Mstar=10**masses[0]['conditional_log10_stellar_mass'];self.Mgas=0.
  self.mass_log_error=masses[0]['published_log10_mass_error']
  self.Dl=1000*geo['conditional_Dl_Mpc'];self.ratio=geo['conditional_Dls_over_Ds'];self.theta=obs['bSIE']
  self.b=self.Dl*self.theta/ARCSEC;self.Re=light['computed_equal_area_half_light_arcsec']*self.Dl/ARCSEC
  self.y=np.array(kin['vrms_kms']);self.cov=np.array(kin['covariance_kms_squared']);self.chol=np.linalg.cholesky(self.cov)
  self.inner_arcsec=np.array(kin['inner_arcsec']);self.outer_arcsec=np.array(kin['outer_arcsec'])
  edges=np.r_[self.inner_arcsec,self.outer_arcsec[-1]]*self.Dl/ARCSEC
  assert np.allclose(self.inner_arcsec[1:],self.outer_arcsec[:-1])
  psf=kin['psf_fwhm_arcsec']*self.Dl/ARCSEC/np.sqrt(8*np.log(2))
  comps=[dict(R=c['R_arcsec']*self.Dl/ARCSEC,n=c['n'],amp=c['amp_at_R'],bn=c['bn']) for c in light['components']]
  s=importlib.util.spec_from_file_location('jr1_component',ROOT/INPUT_NAMES[8]);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
  self.base=m.ComponentModel(self.Re/1.8153,edges,psf,0.,.5,1.,20,comps,n=n,order=order,deproj_order=deproj_order)
  self.r=self.base.r;self.nu=self.base.nu;self.lr=np.log(self.r);self.frac=np.asarray(self.base.mass_interp(self.lr))
  self.gb=G*self.Mstar*self.frac/self.r**2
  x,w=np.polynomial.legendre.leggauss(order*2);self.lens_t=(x+1)*np.pi/4;self.lens_w=w*np.pi/4
  self.lens_r=self.b/np.cos(self.lens_t);self.lens_gb=self.baryon_force(self.lens_r)
  self.total_light_check=self.base.total_mass_check
 def baryon_force(self,r):
  r=np.asarray(r);f=self.base.mass_interp(np.clip(np.log(r),self.lr[0],self.lr[-1]))
  f=np.where(r<self.r[0],self.frac[0]*(r/self.r[0])**self.base.inner_power,f);f=np.where(r>self.r[-1],1.,f)
  return G*self.Mstar*f/r**2
 def orbit_factors(self,p,radial):
  beta=np.full_like(self.r,p['beta']);factor=(self.r/self.Re)**(2*p['beta'])
  if radial:
   ra=10**p['logra']*self.Re;beta=p['beta']+(p['betainf']-p['beta'])*self.r**2/(self.r**2+ra**2)
   factor*=(1+(self.r/ra)**2)**(p['betainf']-p['beta'])
  return factor,beta
 def moments(self,g,p,radial):
  factor,beta=self.orbit_factors(p,radial)
  pressure=-cumulative_trapezoid((self.nu*g*factor)[::-1],self.r[::-1],initial=0)[::-1]/factor
  second=np.trapezoid(self.r**2*pressure*(self.base.W-beta*self.base.T),self.r,axis=1)/self.base.den
  if np.any(second<=0) or np.any(~np.isfinite(second)):raise ValueError('Invalid second moment')
  return np.sqrt(second)
 def local_parameters(self,p,spec):
  if not spec.get('stellar_nuisance',False):return p
  local=p.copy()
  local['logu']+=p[nuisance_key('dm',self.name)] if self.training else 0.
  local['beta']=p[nuisance_key('b',self.name)] if self.training else float(np.mean([p[nuisance_key('b',n)] for n in TRAIN_LENSES]))
  return local
 def prediction(self,p,spec,baryons_only=False):
  p=self.local_parameters(p,spec);q=p['q_sph'] if spec.get('separate_q',False) else p['q']
  if baryons_only:A,rc,rt=0.,self.Re,20*self.Re
  else:A,rc,rt=source_parameters(self.Mstar,0.,self.Re,p,spec['saturation'],1.)
  g=10**p['logu']*self.gb+companion_force(self.r,A,rc,rt,q);v=self.moments(g,p,spec['radial_beta'])
  ray=10**p['logu']*self.lens_gb+companion_force(self.lens_r,A,rc,rt,q)
  alpha=4/CLIGHT**2*np.dot(self.lens_w,ray*self.lens_r);frac=self.ratio*alpha/(self.theta/ARCSEC)-1
  return v,float(frac),(float(A),float(rc),float(rt))
 def angle(self,p,spec,baryons_only=False):
  p=self.local_parameters(p,spec);q=p['q_sph'] if spec.get('separate_q',False) else p['q']
  if baryons_only:A,rc,rt=0.,self.Re,20*self.Re
  else:A,rc,rt=source_parameters(self.Mstar,0.,self.Re,p,spec['saturation'],1.)
  def f(b):
   r=b/np.cos(self.lens_t);g=10**p['logu']*self.baryon_force(r)+companion_force(r,A,rc,rt,q)
   return self.ratio*4/CLIGHT**2*np.dot(self.lens_w,g*r)-b/self.Dl
  grid=self.Re*np.geomspace(1e-6,1e4,200);vals=np.array([f(b) for b in grid]);ii=np.flatnonzero((vals[:-1]>0)&(vals[1:]<=0))
  if not len(ii):return 0.
  i=ii[-1];return float(brentq(f,grid[i],grid[i+1],xtol=1e-10)/self.Dl*ARCSEC)

class Experiment:
 def __init__(self,n=1601,order=64,deproj_order=128):
  self.galaxies=load_sparc();self.lenses=[Lens(name,n,order,deproj_order) for name in LENSES]
  self.r=np.concatenate([d['r'] for d in self.galaxies]);self.y=np.concatenate([d['observed'] for d in self.galaxies])
  self.error=np.concatenate([d['error'] for d in self.galaxies]);self.vs2=np.concatenate([d['star_v2'] for d in self.galaxies]);self.vg2=np.concatenate([d['gas_v2'] for d in self.galaxies])
  self.ids=np.concatenate([np.full(len(d['r']),i) for i,d in enumerate(self.galaxies)])
  self.Msph=np.array([d['Mbulge'] for d in self.galaxies]);self.Ms=np.array([d['Mstar'] for d in self.galaxies]);self.Mg=np.array([d['Mgas'] for d in self.galaxies]);self.Re=np.array([d['Re'] for d in self.galaxies])
  self.split=np.array([d['split'] for d in self.galaxies]);self.train=self.split[self.ids]=='train';self.calls=0
  size=np.array([len(d['r']) for d in self.galaxies]);self.sweight=1/self.error/np.sqrt(size[self.ids]*np.sum(self.split=='train'))
 def sparc_predictions(self,p,spec,baryons_only=False):
  v2=np.maximum(self.vg2+10**p['logu']*self.vs2,0.)
  if not baryons_only:
   A,rc,rt=source_parameters(self.Ms,self.Mg,self.Re,p,spec['saturation'],self.Msph/self.Ms)
   q=p['q']+(p['q_sph']-p['q'])*self.Msph/self.Ms if spec.get('separate_q',False) else np.full_like(self.Ms,p['q'])
   v2+=self.r*companion_force(self.r,A[self.ids],rc[self.ids],rt[self.ids],q[self.ids])
  return np.sqrt(v2)
 def residual(self,x,spec):
  self.calls+=1;p=DEFAULTS|dict(zip(spec['free'],x));y=self.sparc_predictions(p,spec);parts=[((y-self.y)*self.sweight)[self.train]]
  for l in self.lenses:
   if not l.training:continue
   v,f,_=l.prediction(p,spec);parts.extend([solve_triangular(l.chol,v-l.y,lower=True)/np.sqrt(len(l.y)*len(TRAIN_LENSES)),np.array([f/.05/np.sqrt(len(TRAIN_LENSES))])])
  if spec.get('stellar_nuisance',False):
   for l in self.lenses:
    if l.training:
     parts.append(np.array([p[nuisance_key('dm',l.name)]/l.mass_log_error/np.sqrt(len(TRAIN_LENSES)),p[nuisance_key('b',l.name)]/.3/np.sqrt(len(TRAIN_LENSES))]))
  return np.concatenate(parts)
 def metrics(self,p,spec,include_transfer=True,baryons_only=False):
  pred=self.sparc_predictions(p,spec,baryons_only);rows=[]
  for i,d in enumerate(self.galaxies):
   if not include_transfer and d['split']=='test':continue
   ix=self.ids==i;v=pred[ix];e=v-d['observed'];A,rc,rt=source_parameters(d['Mstar'],d['Mgas'],d['Re'],p,spec['saturation'],d['Mbulge']/d['Mstar'])
   rows.append(dict(name=d['name'],split=d['split'],points=int(ix.sum()),RMSE_kms=float(np.sqrt(np.mean(e**2))),fractional_RMS=float(np.sqrt(np.mean((e/d['observed'])**2))),chi2=float(np.sum((e/d['error'])**2)),mean_standardized_square=float(np.mean((e/d['error'])**2)),Mstar_input=d['Mstar'],Mgas_input=d['Mgas'],Re_kpc=d['Re'],A_km2_s2=0. if baryons_only else float(A),rc_kpc=float(rc),rt_kpc=float(rt),equivalent_reservoir_mass_Msun=0. if baryons_only else float(A*rt/G),r_kpc=d['r'].tolist(),observed_kms=d['observed'].tolist(),error_kms=d['error'].tolist(),predicted_kms=v.tolist()))
  groups={}
  for split in sorted(set(r['split'] for r in rows)):
   g=[r for r in rows if r['split']==split];groups[split]=dict(galaxies=len(g),points=sum(r['points'] for r in g),mean_galaxy_RMSE_kms=float(np.mean([r['RMSE_kms'] for r in g])),median_galaxy_fractional_RMS=float(np.median([r['fractional_RMS'] for r in g])),mean_standardized_square=float(np.mean([r['mean_standardized_square'] for r in g])),raw_chi2=sum(r['chi2'] for r in g),within10pct_galaxies=sum(r['fractional_RMS']<.10 for r in g),within20pct_galaxies=sum(r['fractional_RMS']<.20 for r in g))
  lrows=[]
  for l in self.lenses:
   if not include_transfer and not l.training:continue
   v,f,(A,rc,rt)=l.prediction(p,spec,baryons_only);e=v-l.y;wh=solve_triangular(l.chol,e,lower=True);angle=l.angle(p,spec,baryons_only);local=l.local_parameters(p,spec)
   lrows.append(dict(name=l.name,role='fit' if l.training else 'out_of_fit_transfer',n_bins=len(v),Mstar_input=l.Mstar,Mstar_used=l.Mstar*10**local['logu'],orbital_beta0_used=local['beta'],stellar_log_mass_error=l.mass_log_error,stellar_offset_dex=local['logu']-p['logu'],Re_kpc=l.Re,Dl_kpc=l.Dl,Dls_over_Ds=l.ratio,observed_vrms_kms=l.y.tolist(),predicted_vrms_kms=v.tolist(),covariance=l.cov.tolist(),inner_arcsec=l.inner_arcsec.tolist(),outer_arcsec=l.outer_arcsec.tolist(),RMSE_kms=float(np.sqrt(np.mean(e**2))),fractional_RMS=float(np.sqrt(np.mean((e/l.y)**2))),chi2=float(wh@wh),mean_standardized_square=float(wh@wh/len(v)),theta_observed_arcsec=l.theta,theta_predicted_arcsec=angle,theta_fractional_error=angle/l.theta-1,deflection_fractional_error_at_observed=f,A_km2_s2=A,rc_kpc=rc,rt_kpc=rt,equivalent_reservoir_mass_Msun=A*rt/G))
  lg={}
  for role in sorted(set(r['role'] for r in lrows)):
   g=[r for r in lrows if r['role']==role];lg[role]=dict(galaxies=len(g),bins=sum(r['n_bins'] for r in g),mean_RMSE_kms=float(np.mean([r['RMSE_kms'] for r in g])),mean_fractional_RMS=float(np.mean([r['fractional_RMS'] for r in g])),raw_chi2=sum(r['chi2'] for r in g),mean_standardized_square=float(np.mean([r['mean_standardized_square'] for r in g])),theta_fractional_RMS=float(np.sqrt(np.mean([r['theta_fractional_error']**2 for r in g]))),deflection_working_mean_square=float(np.mean([(r['deflection_fractional_error_at_observed']/.05)**2 for r in g])))
  return dict(sparc=groups,lenses=lg,sparc_rows=rows,lens_rows=lrows)

def simple_checks():
 r=np.geomspace(1e-4,1e5,20001);f=companion_force(r,100.,2.,50.,2.);m=r*r*f/G;dm=np.gradient(m,r)
 analytic=m/r*(1+2/(1+(r/2)**2)-r*r/(r*r+50**2));mask=(r>1e-3)&(r<1e4)
 P=np.geomspace(1e-6,1e6,31);a=.3;b=.4;N=2*P/(a+np.sqrt(a*a+8*b*P));k=1e9
 N2=2*k*P/(k*a+np.sqrt((k*a)**2+8*k*b*k*P))
 return dict(minimum_dM_dr=float(dm.min()),finite_mass_relative_error=float(abs(m[-1]/(100*50/G)-1)),density_derivative_error=float(np.max(abs(dm[mask]/analytic[mask]-1))),minimum_spherical_epicyclic_frequency_squared=float((3*f/r+np.gradient(f,r)).min()),clock_rate_rescaling_population_error=float(np.max(abs(N2/N-1))),stationary_rate_balance_relative_error=float(np.max(abs(P-a*N-2*b*N*N)/P)))

def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--max-stage',type=int,default=9);ap.add_argument('--starts',type=int,default=3);ap.add_argument('--max-nfev',type=int,default=250);ap.add_argument('--resume',action='store_true');ap.add_argument('--no-final',action='store_true');args=ap.parse_args();out=args.output_dir
 if out.exists() and not args.resume:raise FileExistsError('Use fresh output directory or explicit --resume')
 out.mkdir(parents=True,exist_ok=True)
 manifest=dict(experiment='JR-1',baseline='c22d188949874ea0ef05b7e782d341973d421385',protocol_commit='7267c709ace23af6afa56cadd04688c85285cf12',defaults=DEFAULTS,bounds=BOUNDS,specifications=SPECS,hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUT_NAMES},source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),lens_fit_names=TRAIN_LENSES,lens_transfer_names=TRANSFER_LENSES,data_roles='Historically exposed samples; transfer objects excluded from fitting, not genuinely new blind data',objective='Equal mean standardized-square blocks for rotation, lens kinematics and deflection; not a calibrated joint likelihood',lens_working_fractional_scale=.05,measured_lens_angle_errors_available=False)
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');start=time.monotonic();E=Experiment();print('PREPARED',len(E.r),'rotation rows',sum(len(l.y) for l in E.lenses),'lens stellar bins',flush=True)
 checks=simple_checks();(out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n');previous=DEFAULTS.copy();records=[]
 for j,spec in enumerate(SPECS[:args.max_stage+1]):
  path=out/(spec['name']+'.json')
  if path.exists() and args.resume:rec=json.loads(path.read_text());previous=DEFAULTS|rec['parameters'];records.append(rec);continue
  rng=np.random.default_rng(20260921+j);free=spec['free'];runs=[]
  if not free:p=DEFAULTS.copy();r=E.residual([],spec);cost=float(r@r)
  else:
   lo=np.array([BOUNDS[k][0] for k in free]);hi=np.array([BOUNDS[k][1] for k in free]);x=np.array([previous[k] for k in free])
   if spec['radial_beta']:x[free.index('betainf')]=previous['beta']
   starts=[np.clip(x,lo+1e-8,hi-1e-8)]
   for k in range(args.starts-1):starts.append(np.clip(starts[0]+rng.normal(0,.1 if k<2 else .25,len(free))*(hi-lo),lo+1e-6,hi-1e-6))
   fits=[]
   for k,x in enumerate(starts):
    t0=time.monotonic();fit=least_squares(E.residual,x,args=(spec,),bounds=(lo,hi),max_nfev=args.max_nfev,ftol=2e-8,xtol=2e-8,gtol=2e-7,x_scale='jac');fits.append(fit)
    runs.append(dict(start=x.tolist(),success=bool(fit.success),status=int(fit.status),objective=float(fit.fun@fit.fun),nfev=fit.nfev,optimality=float(fit.optimality),parameters=fit.x.tolist(),seconds=time.monotonic()-t0));print('START',spec['name'],k,'objective',round(fit.fun@fit.fun,5),'success',fit.success,flush=True)
   fit=min(fits,key=lambda f:float(f.fun@f.fun));p=DEFAULTS|dict(zip(free,fit.x.tolist()));cost=float(fit.fun@fit.fun)
  metrics=E.metrics(p,spec,include_transfer=False);rec=dict(specification=spec,parameters=p,free_parameter_count=len(free),objective=cost,optimizer_starts=runs,bounds_touched=[k for k in free if min(abs(p[k]-np.array(BOUNDS[k])))<1e-3],metrics=metrics)
  path.write_text(json.dumps(rec,indent=2,allow_nan=False)+'\n');records.append(rec);previous=p
  print('RESULT',spec['name'],json.dumps(dict(objective=cost,SPARC=metrics['sparc'],LENSES=metrics['lenses'],parameters=p)),flush=True)
 if args.no_final:print('Transfer data not scored: development stage only',flush=True);return
 criteria=[]
 for rec in records:
  m=rec['metrics'];l=m['lenses']['fit'];score=m['sparc']['validation']['mean_standardized_square']+l['mean_standardized_square']+l['deflection_working_mean_square'];criteria.append(dict(model=rec['specification']['name'],selection_score=score))
 chosen=min(range(len(criteria)),key=lambda i:criteria[i]['selection_score']);selected=records[chosen]
 decision=dict(criteria=criteria,selected_model=selected['specification']['name'],selection_rule='SPARC validation plus fitted lens kinematic/deflection blocks; transfer objects excluded',created_before_transfer_evaluation=True)
 (out/'selection.json').write_text(json.dumps(decision,indent=2)+'\n');full=E.metrics(selected['parameters'],selected['specification'],True)
 (out/'selected-full-results.json').write_text(json.dumps(full,indent=2,allow_nan=False)+'\n');baselines={}
 for label,only in [('baryons_only',True),('initial_root_reservoir',False)]:
  m=E.metrics(DEFAULTS,SPECS[0],True,only);(out/(label+'-full-results.json')).write_text(json.dumps(m,indent=2,allow_nan=False)+'\n');baselines[label]={k:m[k] for k in ('sparc','lenses')}
 summary=dict(experiment='JR-1',selection=decision,selected_parameters=selected['parameters'],selected_results={k:full[k] for k in ('sparc','lenses')},baselines=baselines,model_summaries=[dict(name=r['specification']['name'],nparams=r['free_parameter_count'],objective=r['objective'],parameters=r['parameters'],bounds_touched=r['bounds_touched'],sparc=r['metrics']['sparc'],lenses=r['metrics']['lenses']) for r in records],checks=checks,total_forward_evaluations=E.calls,seconds=time.monotonic()-start,limitations=['Source-linked effective spherical companion envelope, not a microscopic production solution','Conditional static lens distances and population-mass rescalings; spherical lens geometry','SPARC tabulated baryonic forces and integrated source mass are not one exact density-PDE reconstruction','No full lens image likelihood; 5 percent is a working optimization scale, not a measured uncertainty','No Solar-system, cluster, propagation-time or source-energy validation','Static observables do not identify a common exchange-clock rate'])
 (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n');print('FINAL',json.dumps(summary),flush=True)
if __name__=='__main__':main()
