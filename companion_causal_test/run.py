"""Exploratory causal tests; python with numpy/scipy/pandas. No network needed."""
from pathlib import Path
import numpy as np, json,csv,io,zipfile,re,hashlib
from scipy.optimize import minimize
P=Path(__file__).resolve().parent;D=P.parent/'companion_wave_test/data'
G=4.30091727003628e-6;C=299792458.;KPC=3.085677581491367e19;MSUN=1.98847e30;LSUN=3.828e26;YR=31557600.;T=1e10*YR
KAPPA=7.7315e-5*3.261563777167433 # per Mpc
canon=lambda s:re.sub(r'(?<=\D)0+(?=\d)','',re.sub(r'[^A-Z0-9]','',s.upper()))
def read_dust(fn):
 out={}
 for l in (P/'data'/fn).read_text().splitlines():
  try:
   d=dict(name=l[:23].strip(),ra=float(l[24:33]),dec=float(l[34:43]),dist=float(l[49:61]),lum=float(l[100:109]),elum=float(l[110:119]),mass=float(l[80:89]),sfr=float(l[62:70]))
   if d['lum']>0 and d['dist']>0:
    a,b=np.deg2rad([d['ra'],d['dec']]);d['xyz']=d['dist']*np.array([np.cos(b)*np.cos(a),np.cos(b)*np.sin(a),np.sin(b)]);out[canon(d['name'])]=d
  except ValueError:pass
 return out
dust=read_dust('themis.dat');dust_alt=read_dust('dl14.dat');meta={}
for l in (D/'SPARC_Lelli2016c.mrt').read_text().splitlines():
 f=l.split()
 if len(f)!=19:continue
 try:meta[f[0]]=dict(inc=float(f[5]),rd=float(f[11]),q=int(f[17]),lum=float(f[7])*1e9,dist=float(f[2]),mhi=float(f[13])*1e9)
 except ValueError:pass
z=zipfile.ZipFile(D/'Rotmod_LTG.zip');data={}
for fn in sorted(z.namelist()):
 name=fn.replace('_rotmod.dat','');m=meta[name];a=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(fn))))
 vb2=a[:,3]*abs(a[:,3])+.5*a[:,4]*abs(a[:,4])+.7*a[:,5]*abs(a[:,5]);good=np.isfinite(a).all(axis=1)&(a[:,0]>0)&(a[:,1]>0)&(a[:,2]>0)&(vb2>0)
 if m['q']>2 or m['inc']<30 or m['rd']<=0 or good.sum()<5:continue
 a=a[good];vb2=vb2[good];data[name]=dict(**m,r=a[:,0],v=a[:,1],vb2=vb2,mb=.5*m['lum']+1.33*m['mhi'])
# Baryonic depth proxy is calculated without observed rotation. Bulge M/L not resolved in total mass proxy.
for n,d in data.items():d['vdepth']=np.sqrt(2*G*d['mb']/d['rd'])
matched=sorted(n for n in data if canon(n) in dust)
split=json.loads((D/'sparc_frozen.json').read_text())['split'];msplit={s:[n for n in ns if n in matched] for s,ns in split.items()}
# Commit candidate definitions before fits/score evaluation.
protocol={'status':'Exploratory; reused galaxies and old split. Newly downloaded independent SED observable, not blind discovery.', 'sources':len(dust),'matched':len(matched),'split':msplit,'duration_yr':1e10,'halo_radius_Rd':10,'energy_response':'ordinary eta=1 baseline and fitted universal gain','retention':'f=x^m/(1+x^m), x=vdepth/200 km/s; m>=0','candidates':['local','external','mix','local_retention','external_retention','mix_retention'],'geometry':'local softened inverse-square production; uniform external deposition; truncate at 10 Rd','cluster':'same acceleration parameters from 149 galaxy training objects, not disk-template transfer','cv':'5 deterministic name-hash folds; all candidate scores reported; exploratory'}
(P/'protocol.json').write_text(json.dumps(protocol,indent=2))
# The near-source companion flux: F_c = L * log(S)/(4 pi d^2 S^2), S=exp(k*d).
# Includes energy loss and interval stretching. Exact exponential history is only a local approximation to revised history.
energy=[]
for n in matched:
 d=data[n];s=dust[canon(n)];L=s['lum']*(d['dist']/s['dist'])**2*LSUN;d['lbol']=L
 xyz=s['xyz'];fc=0.;fc_no_stretch=0.;nearest=1e9
 for key,t in dust.items():
  if key==canon(n):continue
  sep=float(np.linalg.norm(t['xyz']-xyz));sep=max(sep,.03);nearest=min(nearest,sep)
  x=KAPPA*sep;fc+=t['lum']*LSUN*x*np.exp(-2*x)/(4*np.pi*(sep*1000*KPC)**2)
  fc_no_stretch+=t['lum']*LSUN*x/(4*np.pi*(sep*1000*KPC)**2)
 d['uc']=fc/C
 R=d['r'][-1];mx=R*(d['v'][-1]**2-d['vb2'][-1])/G;E=mx*MSUN*C*C
 f=-np.expm1(-KAPPA*R/1000)
 full=L*T;local=full*f;ext=np.pi*(R*KPC)**2*fc*T
 alt=dust_alt.get(canon(n),s);Lalt=alt['lum']*(d['dist']/alt['dist'])**2*LSUN
 energy.append(dict(galaxy=n,Lbol_Lsun=L/LSUN,Lbol_fractional_SED_error=s['elum']/s['lum'],SPARC_distance_Mpc=d['dist'],DustPedia_distance_Mpc=s['dist'],Lbol_DL14_over_THEMIS=Lalt/L,Rmax_kpc=R,extra_mass_Msun=mx,full_conversion_shortfall=E/full,local_conversion_shortfall=E/local,catalog_external_shortfall=E/ext,combined_local_external_shortfall=E/(local+ext),external_companion_u_J_m3=fc/C,nearest_catalog_source_Mpc=nearest,external_no_stretch_over_stretched=fc_no_stretch/fc))
def templates(d,rhfac=10):
 rd=d['rd'];rh=rhfac*rd;rr=np.minimum(d['r'],rh);r=d['r'];h=rr-rd*np.arctan(rr/rd)
 # mass corresponding to perfect retention for T. Local softened q integrates k L h.
 ml=d['lbol']*T/C**2/MSUN*(KAPPA/1000)*h
 me=np.pi*(rh*KPC)**2*C*d['uc']*T/C**2/MSUN*(rr/rh)**3
 return G*ml/r,G*me/r
models={'local':([8.],[(0,15)]),'external':([6.],[(0,15)]),'mix':([8.,6.],[(0,15),(0,15)]),'local_retention':([8.,1.],[(0,15),(0,6)]),'external_retention':([6.,1.],[(0,15),(0,6)]),'mix_retention':([8.,6.,1.],[(0,15),(0,15),(0,6)])}
def predict(model,p,d,rhfac=10):
 l,e=templates(d,rhfac)
 f=1.
 if model.endswith('retention'):
  xx=(d['vdepth']/200)**p[-1];f=xx/(1+xx)
 if model.startswith('local'):extra=10**p[0]*l
 elif model.startswith('external'):extra=10**p[0]*e
 else:extra=10**p[0]*l+10**p[1]*e
 return np.sqrt(d['vb2']+f*extra)
def scores(pred,ns):
 return dict(n=len(ns),RMSE_kms=float(np.sqrt(np.mean([np.mean((pred(n)-data[n]['v'])**2) for n in ns]))),log_RMSE_dex=float(np.sqrt(np.mean([np.mean(np.log10(pred(n)/data[n]['v'])**2) for n in ns]))))
def fit(model,ns,rhfac=10):
 start,bounds=models[model]
 def loss(p):return np.mean([np.mean(np.log10(predict(model,p,data[n],rhfac)/data[n]['v'])**2) for n in ns])
 opts=[]
 for shift in [-1.,0.,1.]:
  x=np.array(start);x[0]+=shift;o=minimize(loss,x,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-12,'maxiter':800});opts.append(o)
 o=min(opts,key=lambda o:o.fun);return o.x.tolist(),bool(o.success)
results={};predrows=[]
fold={n:int(hashlib.sha256(n.encode()).hexdigest()[:8],16)%5 for n in matched}
for m in models:
 p,ok=fit(m,msplit['train']);sc={s:scores(lambda n:predict(m,p,data[n]),ns) for s,ns in msplit.items()}
 cv={};pars=[]
 for k in range(5):
  tr=[n for n in matched if fold[n]!=k];te=[n for n in matched if fold[n]==k];pp,success=fit(m,tr);pars.append(dict(fold=k,parameters=pp,success=success))
  for n in te:cv[n]=predict(m,pp,data[n])
 results[m]=dict(parameters=p,success=ok,scores=sc,cv_scores=scores(lambda n:cv[n],matched),cv_parameters=pars)
 for n in matched:
  for r,v,vp in zip(data[n]['r'],data[n]['v'],cv[n]):predrows.append(dict(model=m,galaxy=n,fold=fold[n],r_kpc=r,observed_kms=v,predicted_kms=vp))
print('MATCHED MODELS',json.dumps(results,indent=2),flush=True)
# Reference on identical matched sample; frozen previously fitted global model.
prior=json.loads((P.parent/'companion_deposition_fit/results.json').read_text())['models']
p=prior['lum_size_core']['parameters']
def reference(n):
 d=data[n];r=d['r'];rc=10**p[3]*d['rd'];amp=1e4*10**p[0]*(d['lum']/1e10)**p[1]*(d['rd']/3)**p[2];return np.sqrt(d['vb2']+amp*(1-rc/r*np.arctan(r/rc)))
refs={'baryons':scores(lambda n:np.sqrt(data[n]['vb2']),matched),'previous_empirical_core':scores(reference,matched),'ordinary_local_plus_catalog':scores(lambda n:predict('mix',[0,0],data[n]),matched)}
# Robustness of finite boundary and SED calibration: refit local candidate, same train, no model selection.
sensitivity={}
for rhfac in [5,20]:
 pp,ok=fit('local',msplit['train'],rhfac);sensitivity[str(rhfac)]=dict(parameters=pp,test=scores(lambda n:predict('local',pp,data[n],rhfac),msplit['test']))
# Modified gravity response: an extra depth factor shared across galaxies/clusters.
# gx = a0 * 10^a * (gb/a0)^q * (vdepth/200)^m; this is a force-law surrogate, not companion interaction.
a0=3085.677581491367
am={'power':([-.2,.5],[(-3,3),(0,1)]),'power_depth_positive':([-.2,.5,.2],[(-3,3),(0,1),(0,3)]),'power_depth_signed':([-.2,.5,0],[(-3,3),(0,1),(-3,3)])}
def apred(p,d):
 gb=d['vb2']/d['r'];gx=10**p[0]*a0*(gb/a0)**p[1]
 if len(p)>2:gx=gx*(d['vdepth']/200)**p[2]
 return np.sqrt(d['vb2']+d['r']*gx)
acc={};cr=list(csv.DictReader((P.parent/'companion_wave_test/cluster_comparison.csv').open()));crout=[]
for m,(st,bo) in am.items():
 o=minimize(lambda p:np.mean([np.mean(np.log10(apred(p,data[n])/data[n]['v'])**2) for n in split['train']]),st,bounds=bo,method='L-BFGS-B',options={'ftol':1e-13,'maxiter':1500});p=o.x
 ratios=[]
 for c in cr:
  r=float(c['R500_Mpc'])*1000;mt=float(c['M500_1e14Msun'])*1e14;mb=mt*(float(c['fgas500'])+.02);gb=G*mb/r**2;vd=np.sqrt(2*G*mb/r);gx=10**p[0]*a0*(gb/a0)**p[1]*(vd/200)**(p[2] if len(p)>2 else 0);mx=gx*r*r/G;ratio=(mt-mb)/mx;ratios.append(ratio);crout.append(dict(model=m,cluster=c['cluster'],required_over_predicted=ratio,vdepth_proxy_kms=vd))
 acc[m]=dict(parameters=p.tolist(),success=bool(o.success),scores={s:scores(lambda n:apred(p,data[n]),ns) for s,ns in split.items()},cluster_ratio_median=float(np.median(ratios)),cluster_ratio_range=[min(ratios),max(ratios)])
# Paired galaxy bootstrap on out-of-fold residuals for local vs retention/mix.
by={}
for row in predrows:by.setdefault((row['model'],row['galaxy']),[]).append(np.log10(row['predicted_kms']/row['observed_kms'])**2)
rng=np.random.default_rng(20260909);ix=rng.integers(0,len(matched),size=(10000,len(matched)));boot={}
for m in ['local_retention','external','mix','mix_retention']:
 diff=np.array([np.mean(by[(m,n)])-np.mean(by[('local',n)]) for n in matched]);bs=diff[ix].mean(axis=1);boot[m]=dict(mean_delta_log_MSE=float(diff.mean()),percentile95=np.quantile(bs,[.025,.975]).tolist(),note='paired held-out galaxy residual resampling; ignores overlapping-fold training uncertainty')
summary={'protocol':protocol,'energy_medians':{k:float(np.median([r[k] for r in energy])) for k in energy[0] if k!='galaxy'},'energy_shortfall_ranges_10_90':{k:np.quantile([r[k] for r in energy],[.1,.9]).tolist() for k in ['full_conversion_shortfall','local_conversion_shortfall','catalog_external_shortfall','combined_local_external_shortfall']},'models':results,'references_matched_all':refs,'boundary_sensitivity':sensitivity,'acceleration_models':acc,'paired_cv_bootstrap':boot}
(P/'results.json').write_text(json.dumps(summary,indent=2))
for fn,rows in [('energy.csv',energy),('cv_predictions.csv',predrows),('cluster_predictions.csv',crout)]:
 with (P/fn).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(P/'hashes.json').write_text(json.dumps({str(p.relative_to(P.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'data/themis.dat',P/'data/dl14.dat',D/'SPARC_Lelli2016c.mrt',D/'Rotmod_LTG.zip',D/'sparc_frozen.json']},indent=2))
print('SUMMARY',json.dumps({k:v for k,v in summary.items() if k not in ['models','protocol']},indent=2))
