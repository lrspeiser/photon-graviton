from pathlib import Path
import numpy as np,json,csv,io,zipfile
from scipy.optimize import minimize
P=Path(__file__).resolve().parent; D=P.parent/'companion_wave_test/data'; G=4.30091727003628e-6
meta={}
for line in (D/'SPARC_Lelli2016c.mrt').read_text().splitlines():
 f=line.split()
 if len(f)!=19:continue
 try:meta[f[0]]=dict(inc=float(f[5]),rd=float(f[11]),q=int(f[17]),lum=float(f[7])*1e9)
 except ValueError:pass
z=zipfile.ZipFile(D/'Rotmod_LTG.zip');data={}
for fn in sorted(z.namelist()):
 name=fn.replace('_rotmod.dat','');m=meta[name];a=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(fn))))
 vb2=a[:,3]*abs(a[:,3])+.5*a[:,4]*abs(a[:,4])+.7*a[:,5]*abs(a[:,5]);good=np.isfinite(a).all(axis=1)&(a[:,0]>0)&(a[:,1]>0)&(a[:,2]>0)&(vb2>0)
 if m['q']>2 or m['inc']<30 or m['rd']<=0 or good.sum()<5:continue
 a=a[good];vb2=vb2[good];data[name]=dict(**m,r=a[:,0],v=a[:,1],vb2=vb2)
split=json.loads((D/'sparc_frozen.json').read_text())['split']
models={'area_core':([0,0],[(-3,3),(-2,2)]),'lum_core':([0,.5,0],[(-3,3),(-.5,1.5),(-2,2)]),'lum_size_core':([0,.5,0,0],[(-3,3),(-.5,1.5),(-2,2),(-2,2)]),'acceleration':([0,.5],[(-3,3),(0,1)]),'baryons':([],[])}
(P/'protocol.json').write_text(json.dumps({'status':'Exploratory reuse of existing split, not blind','split_counts':{s:len(ns) for s,ns in split.items()},'models':models,'loss':'mean per-galaxy mean squared log10 velocity residual','selection':'minimum validation loss among 3 positive-density cored models; acceleration separate comparison','no_per_object_fitted_parameters':True},indent=2))
def predict(model,p,d):
 r=d['r']; vb=d['vb2'];L=d['lum']/1e10;rd=d['rd'];r0=3
 if model=='baryons':return np.sqrt(vb)
 if model=='acceleration':
  gb=vb/r; gx=10**p[0]*3085.677581491367*(gb/3085.677581491367)**p[1];return np.sqrt(vb+r*gx)
 if model=='area_core':amp=1e4*10**p[0]*(rd/r0);rc=10**p[1]*rd
 if model=='lum_core':amp=1e4*10**p[0]*L**p[1];rc=10**p[2]*rd
 if model=='lum_size_core':amp=1e4*10**p[0]*L**p[1]*(rd/r0)**p[2];rc=10**p[3]*rd
 return np.sqrt(vb+amp*(1-rc/r*np.arctan(r/rc)))
def loss(model,p,names):return float(np.mean([np.mean(np.log10(predict(model,p,data[n])/data[n]['v'])**2) for n in names]))
results={};predrows=[]
for m,(start,bounds) in models.items():
 if start:
  opts=[]
  for shift in [0,-.5,.5]:
   initial=np.array(start,float);initial[0]+=shift
   o=minimize(lambda p:loss(m,p,split['train']),initial,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-13,'gtol':1e-8,'maxiter':1500});opts.append(o)
  o=min(opts,key=lambda o:o.fun);params=o.x.tolist();success=bool(o.success)
 else:params=[];success=True
 scores={}
 for s,names in split.items():
  mse=[np.mean((predict(m,params,data[n])-data[n]['v'])**2) for n in names]
  scores[s]={'n':len(names),'RMSE_kms':float(np.sqrt(np.mean(mse))),'log_RMSE_dex':float(np.sqrt(loss(m,params,names)))}
  if s=='test':
   for n in names:
    d=data[n]
    for r,v,pred in zip(d['r'],d['v'],predict(m,params,d)):predrows.append(dict(model=m,galaxy=n,r_kpc=r,observed_kms=v,predicted_kms=pred))
 results[m]={'parameters':params,'optimizer_success':success,'scores':scores}
winner=min(['area_core','lum_core','lum_size_core'],key=lambda m:results[m]['scores']['validation']['log_RMSE_dex'])
# Transfer the acceleration formula, frozen at its galaxy training fit, to the 11-cluster HSE table.
cp=list(csv.DictReader((P.parent/'companion_wave_test/cluster_comparison.csv').open()));ap=results['acceleration']['parameters']; cr=[]
for d in cp:
 r=float(d['R500_Mpc'])*1000;M=float(d['M500_1e14Msun'])*1e14;Mb=M*(float(d['fgas500'])+.02);gb=G*Mb/r**2;gx=10**ap[0]*3085.677581491367*(gb/3085.677581491367)**ap[1];mx=gx*r*r/G
 cr.append(dict(cluster=d['cluster'],required_extra_Msun=M-Mb,predicted_extra_Msun=mx,required_over_predicted=(M-Mb)/mx))
# Effective density diagnostic for acceleration branch: derivative of r^2 gx on observed radii.
neg=0;total=0
for d in data.values():
 r=d['r'];v=predict('acceleration',ap,d);mh=r*(v*v-d['vb2'])/G;gradient=np.diff(mh)/np.diff(r);neg+=int((gradient<0).sum());total+=len(gradient)
summary={'models':results,'validation_selected_cored_model':winner,'cluster_acceleration_median_required_over_predicted':float(np.median([x['required_over_predicted'] for x in cr])),'cluster_acceleration_range_required_over_predicted':[min(x['required_over_predicted'] for x in cr),max(x['required_over_predicted'] for x in cr)],'acceleration_negative_mass_slope_intervals':neg,'acceleration_total_intervals':total,'note':'Negative slopes are finite-difference diagnostics at fixed noisy baryon inputs, not statistical exclusions.'}
(P/'results.json').write_text(json.dumps(summary,indent=2))
for fname,rows in [('predictions.csv',predrows),('cluster_predictions.csv',cr)]:
 with (P/fname).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(summary,indent=2))
