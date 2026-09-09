from pathlib import Path
import numpy as np,json,re,hashlib
from scipy.optimize import curve_fit,least_squares,brentq,minimize_scalar
from scipy.integrate import solve_ivp
P=Path(__file__).resolve().parent;D=P/'data'
text=(D/'blondin2008.txt').read_text();section=text[text.index('Aging rate measurements'):];section=section[:section.index('4. TESTING') ] if '4. TESTING' in section else section
rows=[]
for line in section.splitlines():
 m=re.match(r'\s+(\S+)\s+(0\.\d+)\s+(0\.\d+)\s+([01]\.\d+)\s+\((0\.\d+)\)',line)
 if m: rows.append([m[1],*map(float,[m[2],m[4],m[5]])])
 if len(rows)==35:break
assert len(rows)==35,len(rows)
(D/'spectral_aging.json').write_text(json.dumps(rows,indent=2))
out={}
for name,rr in [('all',rows),('high_z',[r for r in rows if r[1]>.2])]:
 z,a,s=np.array([r[1:] for r in rr]).T
 b,cov=curve_fit(lambda z,b:(1+z)**(-b),z,a,sigma=s,absolute_sigma=True,p0=[1.])
 out[name]={'N':len(z),'b':float(b[0]),'sigma_b':float(np.sqrt(cov[0,0])),'chi2_b1':float(np.sum(((a-1/(1+z))/s)**2)),'chi2_b0':float(np.sum(((a-1)/s)**2))}
# FIRAS same released residual convention, nonlinear equal-weight thermal mixture.
f=np.loadtxt(D/'firas.txt');nu=f[:,0]*100*299792458.;r=f[:,2]/1000;s=f[:,3]/1000;gal=f[:,4]/1000
h=6.62607015e-34;k=1.380649e-23;c=299792458.;T0=2.725
B=lambda T:2*h*nu**3/c**2/np.expm1(h*nu/(k*T))/1e-20
base=B(T0)
def mixture(delta):
 def fun(p):
  T=T0*(1+p[0]); mod=(B(T*(1-delta))+B(T*(1+delta)))/2-base+p[1]*gal
  return (r-mod)/s
 o=least_squares(fun,[0,0],xtol=1e-12,ftol=1e-12,gtol=1e-10)
 return float(np.sum(o.fun**2)),o.x.tolist()
best=minimize_scalar(lambda d:mixture(d)[0],bounds=(0,.02),method='bounded',options={'xatol':1e-10})
chi0=mixture(0)[0]
out['mixture']={'best_delta':float(best.x),'best_chi2':float(best.fun),'zero_delta_chi2':chi0,'N':43,'assumptions':'equal-weight two-temperature mixture, temperature and Galaxy nuisance; diagonal errors; illustrative profile threshold, not calibrated confidence bound.','grid':{str(d):{'chi2':mixture(d)[0],'delta_chi2_vs_best':mixture(d)[0]-best.fun} for d in [0,1e-5,.001,.003,.01,.05]},'delta_at_Dchi2_3_84':float(brentq(lambda d:mixture(d)[0]-best.fun-3.84,.001,.02))}
# Canonical positive-energy field with pinned endpoints: lowest normal mode.
# n=1+epsilon/(pi*v_ratio)*sin(pi*v_ratio*y+phase)*sin(pi*u).
# v_ratio=field propagation speed / photon speed, y=c*t/R, u=x/R.
eps=7.731496595524618e-5*100
modes=[]
for v in [1,.1,.01]:
 for phase in [0,np.pi/2,np.pi]:
  def fun(u,state):
   y,l=state; f=np.sin(np.pi*u)
   n=1+eps/(np.pi*v)*np.sin(np.pi*v*y+phase)*f
   return [n,-eps*np.cos(np.pi*v*y+phase)*f]
  sol=solve_ivp(fun,[0,1],[0,0],rtol=1e-11,atol=1e-13)
  modes.append({'vfield_over_c':v,'initial_phase':float(phase),'z':float(np.expm1(-sol.y[1,-1]))})
out['source_free_pinned_field']=modes
out['gravitational_clock_constraint']={'epsilon_G_mean':.19e-5,'epsilon_G_sigma':2.48e-5,'approx_95_interval':[.19e-5-1.96*2.48e-5,.19e-5+1.96*2.48e-5],'source':'Delva et al 2018; published summary, no orbit/clock reprocessing'}
out['interpretation']='Published-table reanalysis and conditional model calculations. These are not new blind data. Models sharing b=1 are not distinguished by spectral aging. Pinned-field normal-mode calculations are synthetic, not astronomical observations.'
out['hashes']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in D.iterdir() if f.is_file()}
(P/'results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
