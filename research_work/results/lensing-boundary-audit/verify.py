"""Independent projected-source calculation and interpolation eligibility."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad
H=Path(__file__).resolve().parent;J=H.parent/'joint-galaxy-audit'
load=lambda p:json.loads(p.read_text())
pred=load(H/'predictions.json');src=load(J/'galaxy-rotation-predictions.json');par=load(J/'results.json')['sparc']['parameters']
G=4.30091727003628e-6;C=299792.458;KPC=3.085677581491367e19;ARCSEC=180/np.pi*3600
A,p,astar=par['A'],par['p'],par['a_star_m_s2']*KPC/1e6;q=2*p
errors=[];eligibility=[]
for name in sorted(set(x['galaxy'] for x in src)):
    d=sorted([x for x in src if x['galaxy']==name and x['model']=='baryons'],key=lambda x:x['R_kpc'])
    r=np.array([x['R_kpc'] for x in d]);gb=np.array([x['predicted_kms']**2/x['R_kpc'] for x in d]);gc=A*astar*(gb/astar)**p
    f=PchipInterpolator(r,gc,extrapolate=False);df=f.derivative();mass=r*r*gc/G
    sampled=np.unique(np.r_[r,np.geomspace(r[0],r[-1],1000)])
    min_derivative=float(np.min(2*sampled*f(sampled)+sampled**2*df(sampled)))
    eligibility.append(dict(galaxy=name,nonmonotonic_enclosed_mass_at_original_nodes=bool((np.diff(mass)<0).any()),
        negative_density_with_refined_probe_grid=bool(min_derivative<0),
        classification='node-level spherical positive-source conflict' if (np.diff(mass)<0).any() else ('interpolant-level conflict; monotone mass interpolation is available' if min_derivative<0 else 'positive on sampled grid only')))
    for x in [y for y in pred if y['galaxy']==name]:
        b=x['b_kpc'];rt=x['outer_radius_over_Rmax']*r[-1]
        # A spherical shell at s contributes fraction 1-sqrt(1-b²/s²)
        # of its mass to the cylinder of radius b. Stable form below.
        def fraction(s):
            u=(b/s)**2
            return u/(1+np.sqrt(max(0.,1-u)))
        inside=quad(lambda s:float((2*s*f(s)+s*s*df(s))/G)*fraction(s),b,r[-1],
            points=r[(r>b)&(r<r[-1])],epsabs=.001,epsrel=1e-9,limit=1000)[0] if b<r[-1] else 0.
        outside=quad(lambda s:(2-q)*gc[-1]*r[-1]**q*s**(1-q)/G*fraction(s),r[-1],rt,
            epsabs=.001,epsrel=1e-9)[0] if rt>r[-1] else 0.
        projected=float(b*b*f(b)/G)+inside+outside
        independent=4*G*projected/(C*C*b)*ARCSEC
        errors.append(abs(independent-x['extra_deflection_arcsec'])/max(abs(independent),1e-12))
assert max(errors)<1e-7
assert len(errors)==1788
summary={'independent_projected_mass_comparisons':len(errors),'maximum_relative_deflection_difference':float(max(errors)),
         'node_level_mass_monotonicity_failures':[x['galaxy'] for x in eligibility if x['nonmonotonic_enclosed_mass_at_original_nodes']],
         'interpolant_density_failures':[x['galaxy'] for x in eligibility if x['negative_density_with_refined_probe_grid']],
         'scope':'Numerical identity and source eligibility, not a lensing observation test. No failed template repaired or silently dropped.'}
positive_names={x['galaxy'] for x in eligibility if not x['negative_density_with_refined_probe_grid']}
summary['positive_sampled_templates']=len(positive_names)
summary['boundary_changes_positive_sampled_templates']={}
for frac in [.25,.5,1.]:
    vals=[]
    for n in positive_names:
        x=[v for v in pred if v['galaxy']==n and v['b_over_Rmax']==frac]
        vals.append(x[-1]['extra_deflection_arcsec']/x[0]['extra_deflection_arcsec'])
    summary['boundary_changes_positive_sampled_templates'][str(frac)]=float(np.median(vals))
for name,obj in [('verification.json',summary),('source-eligibility-refined.json',eligibility)]:
    (H/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
