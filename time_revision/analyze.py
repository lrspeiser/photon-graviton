from pathlib import Path
import pickle,importlib,json,hashlib
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
P=Path(__file__).resolve().parent; D=P/'data'
allowed={('pandas.core.frame','DataFrame'),('pandas.core.internals.managers','BlockManager'),('pandas._libs.internals','_unpickle_block'),('numpy.core.multiarray','_reconstruct'),('numpy','ndarray'),('numpy','dtype'),('builtins','slice'),('pandas.core.indexes.base','_new_Index'),('pandas.core.indexes.base','Index'),('pandas.core.indexes.range','RangeIndex')}
class Restricted(pickle.Unpickler):
 def find_class(self,module,name):
  if (module,name) not in allowed: raise ValueError((module,name))
  return getattr(importlib.import_module(module),name)
def fit(z,w,e):
 b,cov=curve_fit(lambda z,b:(1+z)**b,z,w,p0=[1.],sigma=e,absolute_sigma=True)
 return {'n':len(z),'b':float(b[0]),'formal_sigma_b':float(np.sqrt(cov[0,0])),'chi2_dof':float(np.sum(((w-(1+z)**b[0])/e)**2)/(len(z)-1)),'rmse_b1':float(np.sqrt(np.mean((w-(1+z))**2))),'rmse_b0':float(np.sqrt(np.mean((w-1)**2)))}
out={'DES':{},'status':'Derived-data reanalysis; published reference curves already assume dilation. No independent raw-light-curve or blind validation claim.'}; allrows=[]
for band in 'griz':
 with (D/(band+'BAND_pickled_data')).open('rb') as f: df=Restricted(f).load()
 df=df[['CID','z','Width','Width_err']].copy(); df['band']=band; allrows.append(df)
 ok=np.isfinite(df[['z','Width','Width_err']]).all(axis=1)&(df.Width_err>0)&(df.Width_err<df.Width)
 d=df[ok];out['DES'][band]=fit(d.z.to_numpy(),d.Width.to_numpy(),d.Width_err.to_numpy())
df=pd.concat(allrows,ignore_index=True);df.to_csv(D/'DES_widths.csv',index=False)
d=df[np.isfinite(df[['z','Width','Width_err']]).all(axis=1)&(df.Width_err>0)&(df.Width_err<df.Width)].copy()
g=d.groupby('CID'); avg=g.agg(z=('z','first'),width=('Width','mean'),bands=('band','count'));avg['err']=g.Width_err.apply(lambda x:np.sqrt(np.sum(x*x))/len(x));avg.to_csv(D/'DES_event_averages.csv')
out['DES']['all_event_average']=fit(avg.z.to_numpy(),avg.width.to_numpy(),avg.err.to_numpy())
# Object bootstrap retains all bands of an event together; shared-reference covariance remains unavailable.
rng=np.random.default_rng(20260908); boot=[]
for _ in range(400):
 ix=rng.integers(0,len(avg),len(avg));a=avg.iloc[ix];boot.append(fit(a.z.to_numpy(),a.width.to_numpy(),a.err.to_numpy())['b'])
out['DES']['event_bootstrap_95_b']=np.percentile(boot,[2.5,97.5]).tolist()
# Released FIRAS residuals, not the rounded reconstructed total-intensity column.
a=np.loadtxt(D/'firas.txt');nu=a[:,0]*100*299792458.; h=6.62607015e-34;kb=1.380649e-23;c=299792458.;T=2.725
x=h*nu/(kb*T);B=2*h*nu**3/c**2/np.expm1(x)/1e-20 # MJy/sr
G=B*x/(1-np.exp(-x));Y=G*(x/np.tanh(x/2)-4)
r=a[:,2]/1000;err=a[:,3]/1000
# Temperature shift + y distortion + residual Galactic-template coefficient.
X=np.array([G,Y,a[:,4]/1000]).T;A=X/err[:,None];v=r/err
cov=np.linalg.inv(A.T@A);theta=cov@A.T@v;res=r-X@theta
out['FIRAS']={'channels':len(r),'theta_names':['deltaT_over_T','y','galactic_template_coefficient'],'theta':theta.tolist(),'diagonal_only_sigma':np.sqrt(np.diag(cov)).tolist(),'chi2':float(np.sum((res/err)**2)),'dof':len(r)-3,'published_residual_rms_over_peak':float(np.sqrt(np.mean(r*r))/max(B)),'note':'Released residual data and diagonal errors only; foreground and calibration processing already applied. Not a replacement for full-covariance FIRAS limits. Coherent Planck-preserving model is exactly degenerate with blackbody once present T is chosen; no kappa measured.'}
# Screened path: n(u,y)=1+a*y*sin(pi*u)^2, y=c*t/R, source/detector n=1.
# dy/du=n, dlnomega/du=-a*sin(pi*u)^2.
a=2*7.731496595524618e-5*100
sol=[]
for y0 in [0.,1e-4]:
 z=solve_ivp(lambda u,y:[1+a*y[0]*np.sin(np.pi*u)**2,-a*np.sin(np.pi*u)**2],[0,1],[y0,0],rtol=1e-11,atol=1e-13)
 sol.append(z.y[:,-1])
S=np.exp(a/2);pulse=(sol[1][0]-sol[0][0])/1e-4
out['screened_ray']={'R_Mly':100,'active_path_fraction':.5,'local_rate_twice_mean':True,'predicted_S':S,'numerical_pulse_stretch':float(pulse),'numerical_frequency_stretch':float(np.exp(-sol[0][1])),'relative_pulse_error':float(abs(pulse/S-1)),'status':'Synthetic ray-equation consistency test, not observational evidence or microscopic matter theory.'}
out['sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in D.iterdir() if f.is_file()}
(P/'results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
