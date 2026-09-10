"""Inverse stationary wave-support diagnostic, not a field fit."""
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sp
import mpmath as mp

here=Path(__file__).resolve().parent
files=[here.parent/s for s in ['deposit-boundary-admissibility/configurations.json',
    'joint-galaxy-audit/results.json','lens-photometric-audit/training-photometry.json',
    'lensing-data-readiness/conditional-geometry.json']]
configs,fit,phot,geo=[json.loads(f.read_text()) for f in files]
pars=fit['sparc']['parameters']; p=pars['p']; A=pars['A']
phot={r['SDSS']:r for r in phot}; geo={r['Name']:r for r in geo if r['role']=='training'}
assert len(configs)==99 and all(r['role']=='training' for r in configs)
x,q=sp.symbols('x q',positive=True)
density=(1+x)**(-2*q)*(2/x-2*q/(1+x))
ell=sp.diff(sp.log(density),x)/2
lap=sp.factor(sp.diff(ell,x)+ell**2+2*ell/x)
derivative=sp.factor(sp.diff(lap,x))
fn=sp.lambdify((x,q),derivative,'numpy')
# Independent high precision differentiation of the amplitude, not log derivative.
mp.mp.dps=45
pm=mp.mpf(str(p))
def amp(t): return mp.sqrt((1+t)**(-2*pm)*(2/t-2*pm/(1+t)))
def lap_direct(t): return (mp.diff(amp,t,2)+2*mp.diff(amp,t)/t)/amp(t)
checks=[]
for t in [.18153,.5,1.,1.8153,4.,7.2612]:
    reference=float(mp.diff(lap_direct,mp.mpf(str(t))))
    error=abs(fn(t,p)/reference-1)
    assert error<1e-10
    checks.append({'r_over_a':t,'relative_error':error})
kpc=3.085677581491367e19; hbar=1.054571817e-34; ev_mass=1.7826619216279e-36
astar=pars['a_star_m_s2']
fractions=np.geomspace(.1,4.,80)
rows=[]
for cfg in configs:
    name=cfg['Name']; c=cfg['c_dimensionless']
    a=geo[name]['conditional_Dl_Mpc']*1000*phot[name]['Re(I)']*np.pi/(180*3600)/1.8153*kpc
    # c = A*(a_star/B)^(1-p), where B is the baryonic acceleration scale.
    B=astar*(A/c)**(1/(1-p))
    xx=fractions*1.8153
    gg=B*(1/(1+xx)**2+c/(1+xx)**(2*p))
    numerator=fn(xx,p)/a**3
    valid=numerator>0
    masses=np.where(valid,hbar*np.sqrt(np.maximum(numerator,0)/(2*gg))/ev_mass,np.nan)
    assert np.all(valid)
    # Best constant support amplitude in unweighted log acceleration, diagnostic only.
    # g_wave/g = (m_required/m_chosen)^2.
    chosen=float(np.exp(np.mean(np.log(masses))))
    ratio=(masses/chosen)**2
    rows.append({'Name':name,'role':'training','cutoff_over_Re':cfg['cutoff_over_Re'],
        'required_mass_eV_c2_at_0p1_Re':float(masses[0]),
        'required_mass_eV_c2_at_4_Re':float(masses[-1]),
        'max_over_min_required_mass':float(np.max(masses)/np.min(masses)),
        'log_optimal_constant_mass_eV_c2':chosen,
        'min_support_over_required_acceleration':float(np.min(ratio)),
        'max_support_over_required_acceleration':float(np.max(ratio))})
summary={'classification':'Inverse single stationary no-flow wave test; no observational likelihood',
    'p':p,'training_configurations':len(rows),'radial_probes_per_configuration':len(fractions),
    'range_over_Re':[.1,4.],'positive_wave_support_at_all_probes':True,
    'max_over_min_required_mass':{key:float(fun([r['max_over_min_required_mass'] for r in rows]))
         for key,fun in [('min',np.min),('median',np.median),('max',np.max)]},
    'best_constant_acceleration_ratio_range':[min(r['min_support_over_required_acceleration'] for r in rows),
         max(r['max_support_over_required_acceleration'] for r in rows)],
    'derivative_checks':checks,
    'input_sha256':{str(f.relative_to(here.parent)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}
(here/'results.json').write_text(json.dumps(summary,indent=2)+'\n',newline='\n')
(here/'configurations.json').write_text(json.dumps(rows,indent=2)+'\n',newline='\n')
(here/'symbolic-derivative.txt').write_text('d/dx [(laplacian_x sqrt(rho))/sqrt(rho)] =\n'+str(derivative)+'\n',newline='\n')
print(json.dumps(summary,indent=2))
