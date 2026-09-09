from pathlib import Path
import os
import sys,json,csv,shutil
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',Path(__file__).resolve().parents[3]/'research_work/generated'))/'radial-capture'
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/'research_work/energy'))
from radial_capture import fractions,deposited_velocity_km_s
checks=[]
r=np.linspace(0,20,201)
for a,b in [(1,.2),(1,1),(.2,1),(0,1),(1,0),(1e-7,1e-3)]:
    actual=np.array(fractions(r,a,b)).T
    ode=solve_ivp(lambda x,y:[-a*y[0],a*y[0]-b*y[1],b*y[1]],(0,20),[1,0,0],t_eval=r,rtol=1e-11,atol=1e-14)
    error=float(abs(actual-ode.y.T).max())
    residual=float(abs(actual.sum(axis=1)-1).max())
    assert error<1e-9 and residual<1e-13 and actual.min()>-1e-14
    checks.append(dict(check='independent_radial_ODE_and_energy',alpha=a,beta=b,error=error,residual=residual))
small=fractions([1e-5],.2,.3)[2][0]
assert abs(small/(.2*.3*1e-10/2)-1)<2e-6
checks.append(dict(check='short_path_double_transfer_limit',relative_error=float(abs(small/(.2*.3*1e-10/2)-1))))
alpha=7.7315e-5*3.261563777167433/1000
data=list(csv.DictReader((ROOT/'companion_causal_test/energy.csv').open(encoding='utf-8-sig')))
galaxies=[]
for row in data:
    radius=float(row['Rmax_kpc']); lum=float(row['Lbol_Lsun'])
    out=dict(galaxy=row['galaxy'],radius_kpc=radius,luminosity_Lsun=lum,alpha_R=alpha*radius)
    branches=[]
    for capture_length in [.01,1.,10.,100.]:
        beta=1/capture_length
        d=float(fractions([radius],alpha,beta)[2][0])
        branches.append(dict(capture_length_kpc=capture_length,deposited_fraction=d,required_energy_or_response_multiplier=float(row['full_conversion_shortfall'])/d,deposited_velocity_km_s=float(deposited_velocity_km_s([radius],lum,1e10,alpha,beta)[0])))
    out['scenarios']=branches
    # Fast-capture limit has exactly the archived all-local-conversion budget.
    bound=float(row['full_conversion_shortfall'])/(-np.expm1(-alpha*radius))
    assert np.isclose(bound,float(row['local_conversion_shortfall']),rtol=1e-12)
    out['instant_capture_shortfall']=bound
    galaxies.append(out)
checks.append(dict(check='instant_capture_limit_matches_archived_local_supply',galaxies=len(data),rtol=1e-12))

# Dimensionless morphology test, using fixed illustrative alpha R.
x=np.geomspace(.001,1,300); aR=1e-5
curves=[]
fig,axes=plt.subplots(1,2,figsize=(11,4.5))
for bR in [.1,1,10,1000]:
    p,c,d=fractions(x,aR,bR)
    norm=d/(aR*x)
    shape=np.sqrt(d/x)/np.sqrt(d[-1])
    axes[0].plot(x,norm,label=f'capture length = {1/bR:g} R')
    axes[1].plot(x,shape)
    curves.append(dict(beta_R=bR,outer_deposited_fraction=float(d[-1]),velocity_ratio_at_half_radius=float(np.sqrt(fractions([.5],aR,bR)[2][0]/.5/d[-1]))))
axes[0].set(xscale='log',xlabel='Radius / outer radius R',ylabel='Deposited / locally converted energy',ylim=(0,1.05))
axes[1].set(xscale='log',xlabel='Radius / outer radius R',ylabel='Deposited circular-speed contribution / value at R',ylim=(0,1.1))
axes[0].legend(fontsize=8,loc='upper left')
for ax in axes: ax.grid(alpha=.2)
fig.suptitle('Conditional point-source model: fast capture can give a flat outer contribution')
fig.text(.5,.01,'Shape comparison only. No halo fit, binding mechanism or adequate energy supply established.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.04,1,.95)); fig.savefig(OUT/'radial-capture.png',dpi=170); plt.close(fig)
result=dict(scope='Point source; radial straight transport; steady luminosity; constant coefficients; permanent fixed deposits; no secondary loss; deposited Newtonian contribution only',alpha_per_kpc=alpha,duration_years=1e10,checks=checks,shape_scenarios=curves,galaxies=galaxies)
(OUT/'radial-capture-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
for source in [ROOT/'research_work/energy/radial_capture.py',Path(__file__)]:shutil.copy2(source,OUT/source.name)
print(json.dumps(dict(checks_passed=len(checks),galaxies=len(galaxies),alpha_per_kpc=alpha,median_instant_capture_shortfall=float(np.median([g['instant_capture_shortfall'] for g in galaxies])),shape_scenarios=curves),indent=2))
