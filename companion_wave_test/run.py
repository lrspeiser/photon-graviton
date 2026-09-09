from pathlib import Path
import numpy as np, json, zipfile, io, csv, hashlib
from scipy.optimize import minimize_scalar
P=Path(__file__).resolve().parent; SRC=P.parent/'temporal_candidate_audit/data'; D=P/'data';D.mkdir(exist_ok=True)
for f in ['SPARC_Lelli2016c.mrt','Rotmod_LTG.zip','sparc_frozen.json']:(D/f).write_bytes((SRC/f).read_bytes())
G=4.30091727003628e-6 # kpc (km/s)^2 / solar mass
C=299792458.; MS=1.98847e30; LS=3.828e26; YEAR=31557600.; KPC=3.085677581491367e19
T=1e10; KAPPA=7.7315e-5; p=.26390611655880836
meta={}
for line in (D/'SPARC_Lelli2016c.mrt').read_text().splitlines():
 f=line.split()
 if len(f)!=19:continue
 try:meta[f[0]]=dict(inc=float(f[5]),rd=float(f[11]),q=int(f[17]),lum=float(f[7])*1e9,method=int(f[4]))
 except ValueError:pass
zs=zipfile.ZipFile(D/'Rotmod_LTG.zip'); data={}; excluded=[]
for fn in sorted(zs.namelist()):
 name=fn.replace('_rotmod.dat','');m=meta[name];a=np.atleast_2d(np.loadtxt(io.BytesIO(zs.read(fn))))
 vb2=a[:,3]*abs(a[:,3])+.5*a[:,4]*abs(a[:,4])+.7*a[:,5]*abs(a[:,5])
 good=np.isfinite(a).all(axis=1)&(a[:,0]>0)&(a[:,1]>0)&(a[:,2]>0)&(vb2>0)
 if m['q']>2 or m['inc']<30 or m['rd']<=0 or good.sum()<5:
  excluded.append(name);continue
 a=a[good];vb2=vb2[good];order=np.argsort(a[:,0]);a=a[order];vb2=vb2[order]
 data[name]=dict(**m,r=a[:,0],v=a[:,1],err=a[:,2],vg=a[:,3],vd=a[:,4],vbul=a[:,5],vb2=vb2)
rows=[]; shape_rows=[]
for name,d in data.items():
 r=d['r'];v=d['v'];vb2=d['vb2'];R=r[-1];dex=v[-1]**2-vb2[-1];M=R*dex/G
 f=-np.expm1(-KAPPA*R*.003261563777167433) # original fixed path-law candidate
 rec=dict(galaxy=name,rmax_kpc=R,L36_solar=d['lum'],v_outer_kms=v[-1],extra_v2_kms2=dex,extra_mass_signed_Msun=M,energy_J=M*MS*C*C,full_luminosity_years_B1=M*MS*C*C/(d['lum']*LS*YEAR),local_conversion_fraction=f)
 # E_required / energy supplied over 10 Gyr; B=1, capture=1; negative excess is kept, not silently discarded.
 rec['required_gain_full_conversion_B1']=rec['full_luminosity_years_B1']/T
 rec['required_gain_local_conversion_B1']=rec['required_gain_full_conversion_B1']/f
 rec['external_u_required_J_m3_eta1_tau10Gyr']=M*MS*C*C/(np.pi*(R*KPC)**2*C*T*YEAR)
 for ups in [.3,.8]:
  vb_alt=d['vg']*abs(d['vg'])+ups*d['vd']*abs(d['vd'])+(ups+.2)*d['vbul']*abs(d['vbul'])
  rec[f'extra_mass_diskML{ups}_Msun']=R*(v[-1]**2-vb_alt[-1])/G
 rows.append(rec)
 x=r/R; rc=d['rd']; shapes={'central':1/x,'uniform':x*x,'rminus2':np.ones_like(x),'cored_rminus2':(1-rc/r*np.arctan(r/rc))/(1-rc/R*np.arctan(R/rc)),'outer_shell':np.maximum(x**3-.9**3,0)/(1-.9**3)/x}
 for model,h in shapes.items():
  fun=lambda A:float(np.mean((np.sqrt(vb2+A*h)-v)**2))
  opt=minimize_scalar(fun,bounds=(0,4*max(v*v)),method='bounded',options={'xatol':1e-8})
  A=opt.x if opt.fun<fun(0) else 0
  shape_rows.append(dict(galaxy=name,model=model,amplitude_kms2=A,mse_kms2=fun(A)))
 shape_rows.append(dict(galaxy=name,model='baryons_only',amplitude_kms2=0,mse_kms2=float(np.mean((np.sqrt(vb2)-v)**2))))

def write_csv(name,items):
 with (P/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
write_csv('galaxy_energy.csv',rows);write_csv('shape_fits.csv',shape_rows)
pos=[r for r in rows if r['extra_mass_signed_Msun']>0];split=json.loads((D/'sparc_frozen.json').read_text())['split']; byname={r['galaxy']:r for r in pos}
scalings={}
for model in ['luminosity','area','volume']:
 def feature(r):return {'luminosity':r['L36_solar'],'area':r['rmax_kpc']**2,'volume':r['rmax_kpc']**3}[model]
 train=[byname[n] for n in split['train'] if n in byname]
 intercept=float(np.mean([np.log10(r['extra_mass_signed_Msun']/feature(r)) for r in train]))
 scalings[model]={'log10_normalization':intercept,'scores':{}}
 for s,names in split.items():
  rr=[byname[n] for n in names if n in byname];res=np.array([np.log10(10**intercept*feature(r)/r['extra_mass_signed_Msun']) for r in rr])
  scalings[model]['scores'][s]={'n':len(rr),'logmass_RMSE_dex':float(np.sqrt(np.mean(res*res))),'mean_residual_dex':float(res.mean())}
summary={'selection':{'n':len(data),'points':sum(len(d['v']) for d in data.values()),'excluded':excluded,'positive_outer_excess':len(pos),'nonpositive_outer_excess':[r['galaxy'] for r in rows if r['extra_mass_signed_Msun']<=0]},'assumptions':{'time_years':T,'bolometric_proxy_B':1,'disk_ML':.5,'bulge_ML':.7,'ordinary_gravity_eta':1,'p':p,'kappa_per_Mly':KAPPA},'positive_excess_percentiles_10_50_90':{key:np.percentile([r[key] for r in pos],[10,50,90]).tolist() for key in ['extra_mass_signed_Msun','full_luminosity_years_B1','required_gain_full_conversion_B1','required_gain_local_conversion_B1','external_u_required_J_m3_eta1_tau10Gyr','local_conversion_fraction']},'shape_descriptive_galaxy_equal_RMSE_kms':{m:float(np.sqrt(np.mean([r['mse_kms2'] for r in shape_rows if r['model']==m]))) for m in sorted(set(r['model'] for r in shape_rows))},'cross_galaxy_outer_mass_scalings':scalings,'examples':[r for r in rows if r['galaxy'] in ['DDO154','NGC2403','NGC3198','NGC7331']],'mass_to_light_sensitivity_median_gain':{str(ups):float(np.median([max(r[f'extra_mass_diskML{ups}_Msun'],0)*MS*C*C/(r['L36_solar']*LS*YEAR*T) for r in rows])) for ups in [.3,.8]},'data_hashes':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in D.iterdir() if f.is_file()}}
(P/'results.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
