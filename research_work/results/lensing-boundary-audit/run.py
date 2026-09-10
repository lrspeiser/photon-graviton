"""Conditional companion deflections from frozen SPARC radial templates.

No observed lensing data or lensing fit. Psi_c=Phi_c is an explicit additional
metric assumption, not a derived photon-companion coupling.
"""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import quad
from scipy.special import gamma
H=Path(__file__).resolve().parent;J=H.parent/'joint-galaxy-audit'
datafile=J/'galaxy-rotation-predictions.json';fitfile=J/'results.json'
data=json.loads(datafile.read_text());par=json.loads(fitfile.read_text())['sparc']['parameters']
G=4.30091727003628e-6;C=299792.458;KPC=3.085677581491367e19;ARCSEC=180/np.pi*3600
A,p,astar=par['A'],par['p'],par['a_star_m_s2']*KPC/1e6;q=2*p
def extra(g):return A*astar*(g/astar)**p
rows=[];checks=[];skipped=0
for galaxy in sorted(set(r['galaxy'] for r in data)):
    d=sorted([r for r in data if r['galaxy']==galaxy and r['model']=='baryons'],key=lambda r:r['R_kpc'])
    radius=np.array([r['R_kpc'] for r in d]);gb=np.array([r['predicted_kms']**2/r['R_kpc'] for r in d]);gc=extra(gb)
    assert np.all(np.diff(radius)>0) and np.all(gb>0)
    interp=PchipInterpolator(radius,gc,extrapolate=False)
    rmax=radius[-1];gmax=gc[-1]
    # Positive effective spherical density is a separate eligibility check.
    probes=np.unique(np.r_[radius,np.geomspace(radius[0],rmax,400)])
    dens=(2*interp(probes)/probes+interp.derivative()(probes))/(4*np.pi*G)
    checks.append(dict(galaxy=galaxy,split=d[0]['split'],minimum_sampled_effective_density_Msun_kpc3=float(dens.min()),negative_density_probe_count=int((dens<0).sum())))
    for fraction in [.25,.5,1.]:
        b=fraction*rmax
        if b<radius[0]:skipped+=1;continue
        for cut in [1.,2.,5.,10.]:
            rt=cut*rmax;gt=gmax*cut**(-q);tmax=np.arccos(b/rmax);tt=np.arccos(b/rt)
            # Split at every interpolation knot to resolve rapid radial changes.
            knots=np.arccos(b/radius[radius>b]);inside=quad(lambda t:float(interp(min(rmax,b/np.cos(t))))*b/np.cos(t),0,tmax,
                 points=knots[(knots>0)&(knots<tmax)],epsabs=1e-7,epsrel=1e-9,limit=1000)[0] if tmax>0 else 0.
            extension=quad(lambda t:gmax*(b/(rmax*np.cos(t)))**(-q)*b/np.cos(t),tmax,tt,epsabs=1e-7,epsrel=1e-10)[0] if tt>tmax else 0.
            tail=gt*rt*rt/b*(b/rt)**2/(1+np.sin(tt))
            angle=4*(inside+extension+tail)/C**2
            rows.append(dict(galaxy=galaxy,original_split=d[0]['split'],b_kpc=float(b),b_over_Rmax=fraction,Rmax_kpc=float(rmax),
                outer_radius_over_Rmax=cut,extra_deflection_arcsec=float(angle*ARCSEC),
                interior_radial_template_has_negative_density_probe=bool((dens<0).any())))
# Independent exact and quadrature comparisons for a dimensionless power law.
exact=np.sqrt(np.pi)*gamma(q/2)/(2*gamma((q+1)/2))
numeric=quad(lambda t:np.cos(t)**(q-1),0,np.pi/2,epsabs=1e-11)[0]
assert abs(numeric-exact)<1e-9
assert abs(quad(lambda t:np.cos(t),0,np.pi/2)[0]-1)<1e-12
# At b=Rmax and cutoff=Rmax the complete ray is exterior to the source.
for r in rows:
    if r['b_over_Rmax']==1 and r['outer_radius_over_Rmax']==1:
        g=next(x for x in data if x['model']=='baryons' and x['galaxy']==r['galaxy'] and x['R_kpc']==r['Rmax_kpc'])
        expect=4*extra(g['predicted_kms']**2/r['Rmax_kpc'])*r['Rmax_kpc']/C**2*ARCSEC
        assert abs(r['extra_deflection_arcsec']-expect)<1e-12
changes={}
for fraction in [.25,.5,1.]:
    ratios=[]
    for galaxy in sorted(set(r['galaxy'] for r in rows)):
        x=[r for r in rows if r['galaxy']==galaxy and r['b_over_Rmax']==fraction]
        if x:ratios.append(x[-1]['extra_deflection_arcsec']/x[0]['extra_deflection_arcsec'])
    changes[str(fraction)]={'galaxies':len(ratios),'median_ratio_cut10_to_cut1':float(np.median(ratios)),
        'minimum_ratio':float(min(ratios)),'maximum_ratio':float(max(ratios))}
summary=dict(galaxies=len(checks),predictions=len(rows),skipped_inner_impact_parameters=skipped,parameters=par,
    lensing_metric_assumption='Psi_c=Phi_c, static weak field, straight unperturbed rays and asymptotically distant endpoints',
    geometry='Spherical effective extra source; radial template over measured interval, point-baryon power-law continuation to trial cutoff, fixed enclosed mass beyond cutoff.',
    negative_density_template_galaxies=sum(x['negative_density_probe_count']>0 for x in checks),
    conditional_boundary_changes=changes,analytic_powerlaw_integral=exact,numerical_powerlaw_integral=numeric,
    new_parameters_fitted=False,observed_lensing_rows=0,heldout_likelihood_scores_opened=False,
    limitations=['No capture-predicted outer radius; trial cutoffs are sensitivity cases.',
    'Spherical extension of disk rotation is an added approximation; some templates require negative effective source density.',
    'Extra deflection only, not total lensing, shear, image separation or a lens equation.',
    'Neither the lens metric nor the companion origin is derived by this calculation.'],
    input_sha256={x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in [datafile,fitfile]})
for name,obj in [('results.json',summary),('predictions.json',rows),('source-eligibility.json',checks)]:
    (H/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
