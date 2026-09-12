"""Sensitivity of polar trajectories to a shrinking numerical core continuation."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole
from axis import evaluate
bar=FastMultipole.load(ROOT/'research_work/data-cache/bar-field/bar-L40.npz')
zero=bar.m==0;Y=np.sqrt((2*bar.l[zero]+1)/(4*np.pi))
def ordinary(r):
    c=bar.spline(np.log(r))[zero];dc=bar.spline(np.log(r),1)[zero]/r
    return float(c@Y),float(dc@Y)
rows=[];times=np.linspace(0,.25,301)
for R in (1.,3.):
    d=json.loads((HERE/f'prepared-R{R:g}.json').read_text(encoding='utf8'))
    r=next(r for r in d['records'] if np.hypot(*r['direction'][:2])==0)
    path=ROOT/r['cache'];assert hashlib.sha256(path.read_bytes()).hexdigest()==r['cache_sha256']
    c=np.load(path);old=CubicHermiteSpline(c['ages'],c['position'][:,2],c['velocity'][:,2])
    for cut in (1e-4,1e-5,1e-6):
        boundary,derivative=ordinary(cut);kappa=derivative/cut
        def field(z):
            radius=abs(z)
            if radius<cut:return boundary+.5*kappa*(z*z-cut*cut),-kappa*z
            phi,g=ordinary(radius);return phi,-np.sign(z)*g
        def rhs(t,y):return [y[1],field(y[0])[1]]
        sol=solve_ivp(rhs,(0,.25),[R,-50.],method='DOP853',rtol=2e-13,atol=2e-15,max_step=.0002,dense_output=True)
        assert sol.success
        y=sol.sol(times).T;E=.5*y[:,1]**2+np.array([field(z)[0] for z in y[:,0]])
        dx=float(max(abs(y[:,0]-old(times))));dv=float(max(abs(y[:,1]-old(times,1))));de=float(max(abs(E-E[0]))/220**2)
        rows.append(dict(R=R,cut_kpc=cut,matched_curvature=kappa,position_change=dx,velocity_change=dv,
                         energy_drift=de,passes=dx<1e-5 and dv<.05 and de<1e-5,z=y[:,0].tolist(),v=y[:,1].tolist()))
        print(R,cut,dx,dv,de,rows[-1]['passes'],flush=True)
comparisons=[]
for R in (1.,3.):
    r=[r for r in rows if r['R']==R]
    for a,b in zip(r[:-1],r[1:]):
        dx=float(max(abs(np.array(a['z'])-b['z'])));dv=float(max(abs(np.array(a['v'])-b['v'])))
        comparisons.append(dict(R=R,cut_a=a['cut_kpc'],cut_b=b['cut_kpc'],position_change=dx,velocity_change=dv,passes=dx<1e-5 and dv<.05))
(HERE/'axis-core-checks.json').write_text(json.dumps(dict(records=rows,comparisons=comparisons),indent=2)+'\n',encoding='utf8',newline='\n')
