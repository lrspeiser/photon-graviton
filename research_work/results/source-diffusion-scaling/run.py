from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.optimize import minimize_scalar
H=Path(__file__).resolve().parent
width=json.loads((H.parent/'sn2006mk-intrinsic-width/results.json').read_text(encoding='utf-8'))['effective_stretch']
# Illustrative radioactive-heating/diffusion model, not a fit to SN 2006mk.
ni,co=8.8,111.3;ratio=.2;td=13.5;T=150.
def heat(t,scale=1.):return (1-ratio)*np.exp(-t/(ni*scale))+ratio*np.exp(-t/(co*scale))
def model(diffusion,clock=1.):
 sol=solve_ivp(lambda t,y:[2*t/diffusion**2*(heat(t,clock)-y[0])],(0,T),[0.],rtol=1e-10,atol=1e-12,dense_output=True,max_step=.5);assert sol.success;return sol
base=model(td);physical=model(td*width);rescaled=model(td*width,width)
grid=np.linspace(5,120,461);target=base.sol(grid/width)[0]
# A complete time rescaling requires heating clocks to be rescaled too.
maxerr=float(np.max(abs(rescaled.sol(grid)[0]-target)));assert maxerr<1e-8
norm=float(np.max(target));change=physical.sol(grid)[0]
normalized_rms=float(np.sqrt(np.mean((change/change.max()-target/norm)**2)))
# Allow an amplitude and diffusion-time fit while retaining fixed decay times.
def score(x):
 y=model(x).sol(grid)[0];a=float(y@target/(y@y));return float(np.mean((a*y-target)**2)/norm**2)
opt=minimize_scalar(score,bounds=(8,40),method='bounded');assert opt.success;fitted=model(opt.x).sol(grid)[0];amp=float(fitted@target/(fitted@fitted))
# Integrated energy balance for E=td^2 L/(2t), zero initial stored radiation.
checks=[]
for diffusion,sol in [(td,base),(td*width,physical)]:
 end=120.;emitted=quad(lambda t:float(sol.sol(t)[0]),0,end,epsabs=1e-8)[0]
 work=quad(lambda t:diffusion**2*float(sol.sol(t)[0])/(2*t*t) if t>0 else .5,0,end,epsabs=1e-8)[0]
 stored=diffusion**2*float(sol.sol(end)[0])/(2*end);supplied=quad(heat,0,end,epsabs=1e-8)[0];residual=abs((emitted+work+stored)/supplied-1);assert residual<1e-7
 checks.append(dict(diffusion_days=diffusion,emitted=emitted,expansion_work=work,stored=stored,supplied=supplied,relative_balance_error=residual))
out=dict(scope='Conditional diffusion-model scaling test; not observed bolometric or spectral fit',parameters=dict(ni_mean_life_days=ni,co_mean_life_days=co,adopted_slow_heating_fraction=ratio,reference_diffusion_days=td,target_time_stretch=width),naive_diffusion_days=td*width,required_kappa_mass_over_velocity_ratio=width**2,complete_clock_rescaling_max_error=maxerr,naive_peak_normalized_rms=normalized_rms,fixed_clock_best_fit=dict(diffusion_days=float(opt.x),amplitude=amp,rms_as_fraction_of_target_peak=float(np.sqrt(opt.fun))),energy_checks=checks,samples=[dict(time=float(t),target=float(y),diffusion_only=float(v),best_fixed_clock=float(a)) for t,y,v,a in zip(grid,target,change,amp*fitted)])
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({k:v for k,v in out.items() if k!='samples'}))
