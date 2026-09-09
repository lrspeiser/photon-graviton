"""Derived-data comparisons. Previously examined splits are explicitly exploratory."""
from pathlib import Path
import numpy as np,pandas as pd,json,zipfile,io,hashlib
from scipy.optimize import minimize,minimize_scalar,curve_fit,least_squares
from scipy.special import lambertw
P=Path(__file__).resolve().parent;D=P/'data';O=P/'results'
c=299792.458; gamma=7.731496595524618e-11;KPC=3.085677581491367e19;AK=299792458**2*7.7315e-5/(1e6*9.4607304725808e15)
def dump(n,v):(O/n).write_text(json.dumps(v,indent=2,allow_nan=False))
def csv(n,rows):pd.DataFrame(rows).to_csv(O/n,index=False)

# A. Redshift: old split, all parameters trained on 104 groups.
rows=json.loads((D/'redshift_features.json').read_text()); labels={}
for split in ['train','validation','test']:labels.update(json.loads((D/f'{split}_labels.json').read_text()))
# Catalog membership is a selection-dependent endpoint proxy, not a 3D density map.
counts={}
for line in (D/'cf4_table2.dat').read_text().splitlines():
 g=int(line[8:15]);counts[g]=counts.get(g,0)+1
for r in rows:
 ra,dec=np.deg2rad([r['ra'],r['dec']]);r['nhat']=np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)]);r['logN']=np.log(counts[r['group']])
tr=[r for r in rows if r['split']=='train'];emean=np.mean([r['logN'] for r in tr]);esd=np.std([r['logN'] for r in tr])
for r in rows:r['env']=(r['logN']-emean)/esd

def redpred(model,p,rr):
 d=np.array([r['distance_mpc'] for r in rr]);sd=np.array([r['sd_mpc'] for r in rr]);nh=np.array([r['nhat'] for r in rr]);k=np.full(len(rr),p[0])
 if model=='directional_rate':k*=np.exp(nh@p[1:4])
 if model=='endpoint_proxy':k*=np.exp(p[1]*np.array([r['env'] for r in rr]))
 x=k*d/c
 if model=='rolling':z=2*x/(1+np.sqrt(1+4*x));der=k/np.sqrt(1+4*x)
 else:
  a=.5 if model=='energy_only' else 1.;w=lambertw(a*x).real;z=np.expm1(w/a);der=k*np.exp((1-a)*w/a)/(1+w)
 v=c*z
 if model=='bulk_velocity':v+=nh@p[1:4]
 return v,np.sqrt(300**2+(der*sd)**2)
red={};redrows=[];preds={}
for model in ['exponential','rolling','energy_only','directional_rate','endpoint_proxy','bulk_velocity']:
 bounds=[(10,150)]+([(-.5,.5)]*3 if model=='directional_rate' else [(-1500,1500)]*3 if model=='bulk_velocity' else [(-.5,.5)] if model=='endpoint_proxy' else [])
 def obj(p):
  v,s=redpred(model,p,tr);y=np.array([labels[str(r['pgc'])] for r in tr]);return float(np.sum(np.log(s)+.5*((v-y)/s)**2))
 # Scaling makes bulk velocity numerical optimization well conditioned.
 scale=np.array([75]+([300]*3 if model=='bulk_velocity' else [1]*(len(bounds)-1)))
 opt=minimize(lambda q:obj(q*scale),np.array([75]+[0]*(len(bounds)-1))/scale,bounds=[(a/s,b/s) for (a,b),s in zip(bounds,scale)],method='L-BFGS-B',options={'ftol':1e-13,'gtol':1e-8,'maxiter':2000})
 p=opt.x*scale;res={'parameters':p.tolist(),'optimizer_success':bool(opt.success),'scores':{}}
 for split in ['train','validation','test']:
  rr=[r for r in rows if r['split']==split];v,s=redpred(model,p,rr);y=np.array([labels[str(r['pgc'])] for r in rr]);err=v-y
  res['scores'][split]={'N':len(rr),'RMSE_kms':float(np.sqrt(np.mean(err**2))),'mean_NLL_without_constant':float(np.mean(np.log(s)+.5*(err/s)**2))}
  if split=='test':preds[model]=err
  for r,yy,vv,ss in zip(rr,y,v,s):redrows.append(dict(model=model,split=split,pgc=r['pgc'],tile=r['tile'],group=r['group'],catalog_group_entries=counts[r['group']],distance_mpc=r['distance_mpc'],observed_cz=yy,predicted_cz=vv,sigma_cz=ss))
 red[model]=res
# sky-tile bootstrap: descriptive uncertainty conditional on frozen training fit.
test=[r for r in rows if r['split']=='test'];tiles=sorted({r['tile'] for r in test});rng=np.random.default_rng(124);boots=[]
for b in range(3000):boots.append(np.array([i for t in rng.choice(tiles,len(tiles),replace=True) for i,r in enumerate(test) if r['tile']==t]))
for model in red:
 deltas=[np.sqrt(np.mean(preds[model][ix]**2))-np.sqrt(np.mean(preds['exponential'][ix]**2)) for ix in boots]
 red[model]['test_delta_RMSE_vs_exponential_95_descriptive']=np.quantile(deltas,[.025,.975]).tolist()
red['metadata']={'status':'Exploratory reused split, not a new blind test','env_definition':'Standardized ln(number of entries in CF4 table2 sharing group ID)','env_train_mean':float(emean),'env_train_sd':float(esd),'flux_assumption':'D_L=R*(1+z) except energy_only D_L=R*sqrt(1+z); these are phenomenological branches, not derived screened-field flux laws.'}
dump('redshift.json',red);csv('redshift_predictions.csv',redrows)

# B. Timing: same observable exponent, two published-derived methods.
timing={}
for name in ['DES','spectral']:
 if name=='DES':
  d=pd.read_csv(D/'DES_event_averages.csv');z=d.z.values;y=d.width.values;s=d.err.values;sign=1
 else:
  d=np.array([r[1:] for r in json.loads((D/'spectral_aging.json').read_text())]);z,y,s=d.T;sign=-1
 f=lambda z,b:(1+z)**(sign*b)
 b,cov=curve_fit(f,z,y,sigma=s,absolute_sigma=True,p0=[1])
 timing[name]={'N':len(z),'b':float(b[0]),'formal_sigma_b':float(np.sqrt(cov[0,0])),'chi2_b1':float(np.sum(((y-f(z,1))/s)**2)),'chi2_b0':float(np.sum(((y-f(z,0))/s)**2)),'caveat':'Derived templates and shared calibration; diagonal errors are not a full likelihood.'}
dump('timing.json',timing)

# C. GW and clocks: published-summary comparisons, no raw detector fitting.
ext=json.loads((D/'published_measurements.json').read_text());gw=[]
secyr=365.25*86400;mpcly=3.085677581491367e22/9.4607304725808e15
for mpc in ext['GW170817']['distance_sensitivity_Mpc']:
 T=mpc*mpcly;x=gamma*T;tlight=np.log1p(x)/gamma;gap=T-tlight
 # derivative of GW arrival time wrt p at p=1 for n_GW=exp(p*gamma*t).
 derivative=-((1+x)*np.log1p(x)-x)*secyr/gamma
 gw.append({'distance_Mpc':mpc,'photon_first_gap_years_p0':float(gap),'shared_p1_propagation_gap_seconds':0.,'dp_for_10second_gap_linearized':float(10/abs(derivative))})
clock={'Galileo_epsilonG_95_approx':[ext['Galileo']['epsilonG']-1.96*ext['Galileo']['sigma'],ext['Galileo']['epsilonG']+1.96*ext['Galileo']['sigma']],'alpha_model':'alpha_dot/alpha=q_eff*gamma','q_eff_abs_bound_conservative_95':(abs(ext['alpha_clock']['alpha_dot_per_year'])+1.96*ext['alpha_clock']['sigma_per_year'])/gamma,'cases':[]}
for q,scr in [(0,1),(1,1),(1,1e-6),(1,1e-8)]:
 predicted=q*scr*gamma;clock['cases'].append({'q_intrinsic':q,'screen_fraction':scr,'predicted_alpha_dot_per_year':predicted,'standardized_residual':(predicted-ext['alpha_clock']['alpha_dot_per_year'])/ext['alpha_clock']['sigma_per_year']})
dump('GW_clocks.json',{'GW':gw,'clocks':clock,'caveat':'GW homogeneous law and constant-c alternative, illustrative distance and emission assumptions. q is an explicit alpha coupling, not universal time rescaling. Gaussian standardized residuals are diagnostics, not calibrated extreme-tail probabilities.'})

# D. CMB spectrum and angular tables.
f=np.loadtxt(D/'firas.txt');nu=f[:,0]*100*299792458.;res=f[:,2]/1000;sig=f[:,3]/1000;gal=f[:,4]/1000
h=6.62607015e-34;kb=1.380649e-23;cc=299792458.;T0=2.725
B=lambda T:2*h*nu**3/cc**2/np.expm1(h*nu/(kb*T))/1e-20
base=B(T0);mix=[]
for delta in [0,.00001,.001,.003,.01,.05]:
 def fun(p):return (.5*(B(T0*(1+p[0])*(1-delta))+B(T0*(1+p[0])*(1+delta)))-base+p[1]*gal-res)/sig
 opt=least_squares(fun,[0,0],xtol=1e-12,ftol=1e-12,gtol=1e-10)
 mix.append({'fractional_temperature_split':delta,'chi2_diagonal':float(np.sum(opt.fun**2)),'fitted_T_K':float(T0*(1+opt.x[0])),'galaxy_nuisance':float(opt.x[1])})
angular={}
for kind in ['TT','EE']:
 if not (D/f'planck_{kind}.txt').exists():continue
 a=np.loadtxt(D/f'planck_{kind}.txt');ell,y,sl,sh,best=a.T;s=.5*(sl+sh)
 angular[kind]={'N_bins':len(a),'ell_range':[float(ell.min()),float(ell.max())],'zero_signal_chi2_diagonal':float(np.sum((y/s)**2)),'published_bestfit_chi2_diagonal':float(np.sum(((y-best)/s)**2)),'bins_above_5_quoted_errors':int(np.sum(y/s>5)),'caveat':'Binned released bandpowers and diagonal approximation, not Planck likelihood. BestFit is a fitted standard-model reference. A uniform scalar acting on an initially uniform unpolarized sky predicts zero TT and EE; general field perturbations remain unspecified.'}
 csv(f'planck_{kind}_comparison.csv',[dict(ell=l,observed=d,sigma=sig,published_bestfit=b) for l,d,sig,b in zip(ell,y,s,best)])
cmb={'mixture':mix,'HFLS3':{'z':6.34,'observed_1sigma_range_K':[16.4,30.2],'constant_temperature_prediction_K':2.72548,'coherent_scaling_prediction_K':2.72548*7.34},'angular':angular}
dump('CMB.json',cmb)

# E. Rotation curves: same trained split; new surface-density proxy branch.
meta={}
for line in (D/'SPARC_Lelli2016c.mrt').read_text().splitlines():
 f=line.split()
 if len(f)!=19:continue
 try:meta[f[0]]={'method':int(f[4]),'inc':float(f[5]),'rd':float(f[11]),'q':int(f[17])}
 except ValueError:pass
sparc={}
with zipfile.ZipFile(D/'Rotmod_LTG.zip') as z:
 for fn in sorted(z.namelist()):
  name=fn.replace('_rotmod.dat','');m=meta[name];a=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(fn))))
  if m['q']>2 or m['inc']<30 or m['rd']<=0:continue
  vb2=a[:,3]*abs(a[:,3])+.5*a[:,4]*abs(a[:,4])+.7*a[:,5]*abs(a[:,5]);good=np.isfinite(a).all(axis=1)&(a[:,0]>0)&(a[:,1]>0)&(a[:,2]>0)&(vb2>0)
  if good.sum()<5:continue
  a=a[good];r=a[:,0]*KPC;vb2=vb2[good];sigma=.5*a[:,6]+.7*a[:,7]
  sparc[name]={'name':name,**m,'r':r,'v':a[:,1],'vb2':vb2,'gb':vb2*1e6/r,'sigma':np.maximum(sigma,0)}
frozen=json.loads((D/'sparc_frozen.json').read_text());splits=frozen['split'];sptr=[sparc[n] for n in splits['train']]
def sppred(model,p,d):
 x=d['gb']/AK
 if model=='ordinary':return np.sqrt(d['vb2'])
 if model=='RAR':return np.sqrt(d['r']*d['gb']/(-np.expm1(-np.sqrt(x/(10**p[0])))))/1000
 if model=='power':return np.sqrt(d['r']*(d['gb']+10**p[0]*AK*x**p[1]))/1000
 # Surface stellar mass density only. Excludes gas density; not total volume density.
 A,S=10**p[0],10**p[1];power=1 if model=='density_p1' else p[2]
 return np.sqrt(d['vb2']*(1+A/(1+(d['sigma']/S)**power)))
def sploss(model,p,ds):return np.mean([np.mean(np.log10(sppred(model,p,d)/d['v'])**2) for d in ds])
spfits={m:frozen['params'][m] for m in ['ordinary','power','RAR']};diags={}
for model in ['density_p1','density_free_p']:
 bounds=[(-3,3),(-4,5)]+([(.1,4)] if model=='density_free_p' else [])
 opts=[minimize(lambda p:sploss(model,p,sptr),[a,s]+([1] if model=='density_free_p' else []),bounds=bounds,method='L-BFGS-B',options={'ftol':1e-13,'gtol':1e-8}) for a,s in [(0,0),(1,1),(2,2)]]
 opt=min(opts,key=lambda o:o.fun);spfits[model]=opt.x.tolist();diags[model]=bool(opt.success)
sp={};sprows=[]
for model,p in spfits.items():
 result={'parameters':p,'optimizer_success':diags.get(model,True),'scores':{}}
 for split,names in splits.items():
  ds=[sparc[n] for n in names]
  for suffix,sub in [('',ds),('_independent_distance',[d for d in ds if d['method'] in [2,3,5]])]:
   lm=[];km=[]
   for d in sub:
    v=sppred(model,p,d);lm.append(np.mean(np.log10(v/d['v'])**2));km.append(np.mean((v-d['v'])**2))
   result['scores'][split+suffix]={'N_galaxies':len(sub),'N_points':sum(len(d['v']) for d in sub),'galaxy_weighted_log_RMSE_dex':float(np.sqrt(np.mean(lm))),'galaxy_weighted_RMSE_kms':float(np.sqrt(np.mean(km)))}
  if split=='test':
   for d in ds:
    for i,v in enumerate(sppred(model,p,d)):sprows.append(dict(model=model,galaxy=d['name'],r_kpc=d['r'][i]/KPC,stellar_surface_density_Msun_pc2=d['sigma'][i],observed_kms=d['v'][i],predicted_kms=v))
 sp[model]=result
sp['metadata']={'status':'New proxy fits, exploratory reuse of prior split, no new blind test','density_definition':'0.5*SBdisk + 0.7*SBbulge, Msun/pc2; not total volume density','physical_status':'Dynamical enhancement templates, not derived light-propagation or pure apparent-clock effects','mass_to_light':'Fixed disk 0.5 and bulge 0.7; inclinations and distances fixed; errors correlated within galaxies.'}
dump('SPARC.json',sp);csv('SPARC_predictions.csv',sprows)
summary={'redshift_test_RMSE':{m:r['scores']['test']['RMSE_kms'] for m,r in red.items() if m!='metadata'},'timing':timing,'GW':gw,'CMB':cmb,'SPARC_test_RMSE':{m:r['scores']['test']['galaxy_weighted_RMSE_kms'] for m,r in sp.items() if m!='metadata'}}
dump('summary.json',summary);print(json.dumps(summary,indent=2))
