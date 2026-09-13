"""Boundary-conditioned straight-ray capture; fixed third-law parameters."""
from pathlib import Path
import json, hashlib
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import linprog
P=Path(__file__).resolve().parent
alpha=.0002488993286382367/1000 # per kpc
inputs=['halo-deposition-map-results.json','third-radiation-retention-results.json']
halos=json.loads((P/inputs[0]).read_text())['rows']
cp=json.loads((P/inputs[1]).read_text())['models']['attenuated']
k=cp['k0_per_kpc']; assert cp['q']==1/3
cases=[1.1,2.,5.,10.,100.,None]

def primitive(t,b2,a):
    B2=a*a+b2; B=np.sqrt(B2)
    return k*a**4*(t/(2*B2*(B2+t*t))+np.arctan(t/B)/(2*B**3))

def profile(R,a,eta,ratio,nr,na):
    x=np.r_[0.,np.geomspace(1e-5,1.,nr-1)]; r=R*x
    mu,w=leggauss(na); rr=r[:,None]
    if ratio is None:
        t=rr*mu; b2=rr**2-t*t; factor=np.ones_like(t)
    else:
        D=ratio*R; d=np.sqrt(D*D+rr*rr-2*D*rr*mu)
        t=(rr*rr-D*rr*mu)/d; b2=np.maximum(rr*rr-t*t,0)
    te=-np.sqrt(np.maximum(R*R-b2,0))
    if ratio is not None:
        de=d-(t-te)
        assert np.all(de>0)
        factor=(D/d)**2*(-np.expm1(-alpha*de))/(-np.expm1(-alpha*D))
    tau=np.maximum(primitive(t,b2,a)-primitive(te,b2,a),0)
    q=eta*k/(1+(rr/a)**2)**2*factor*np.exp(-tau)
    avg=q@w/2; cumulative=cumulative_trapezoid(4*np.pi*r*r*avg,r,initial=0)
    assert cumulative[-1]>0 and np.all(np.diff(cumulative)>=0)
    cdf=cumulative/cumulative[-1]
    # Absolute angular first moment divided by scalar retained input.
    angular=cumulative_trapezoid(4*np.pi*r*r*(q@(w*mu)/2),r,initial=0)[-1]
    return x,cdf,float(cumulative[-1]),float(abs(angular)/cumulative[-1])

# Independent path integral verification, including nearly tangential rays.
errors=[]
for a,R,b,t in [(5,30,0,12),(17,50,20,10),(2,30,29,2)]:
    te=-np.sqrt(R*R-b*b)
    exact=primitive(t,b*b,a)-primitive(te,b*b,a)
    numeric=quad(lambda z:k/(1+(b*b+z*z)/(a*a))**2,te,t,epsabs=1e-11)[0]
    errors.append(abs(exact-numeric));assert abs(exact-numeric)<1e-9
rows=[]
for h in halos:
    R=h['finite_target_radius_kpc'];a=h['original_capture_scale_kpc'];eta=h['eta']
    curves=[]; diagnostics=[];maxdrift=0
    for ratio in cases:
        x,cdf,power,dipole=profile(R,a,eta,ratio,401,96)
        xx,fine,fp,fd=profile(R,a,eta,ratio,801,192)
        drift=float(np.max(np.abs(cdf-np.interp(x,xx,fine))));maxdrift=max(maxdrift,drift)
        assert drift<.002
        curves.append(fine)
        diagnostics.append(dict(D_over_R=ratio,half_deposit_radius_over_R=float(np.interp(.5,fine,xx)),retained_volume_integral_per_unattenuated_center_flux_kpc2=fp,single_source_deposit_dipole=fd if ratio is not None else None,coarse_fine_max_cdf_difference=drift))
    u=R*xx/h['rs_kpc'];target=np.log1p(u)-u/(1+u);target/=target[-1]
    curves=np.array(curves)
    for diag,c in zip(diagnostics,curves):diag['max_cumulative_fraction_error']=float(np.max(np.abs(c-target)))
    # Mixture weights describe fractions of deposited power, not source luminosities.
    A=np.vstack([np.c_[curves.T,-np.ones(len(xx))],np.c_[-curves.T,-np.ones(len(xx))]])
    b=np.r_[target,-target]
    fit=linprog(np.r_[np.zeros(len(cases)),1.],A_ub=A,b_ub=b,A_eq=[np.r_[np.ones(len(cases)),0]],b_eq=[1],bounds=[(0,None)]*(len(cases)+1),method='highs')
    assert fit.success
    mix=fit.x[:-1]@curves
    # Far source must approach a distant incident field as its distance grows.
    _,far,_,_=profile(R,a,eta,1e5,801,192)
    far_drift=float(np.max(np.abs(far-curves[-1])));assert far_drift<1e-4
    row=dict(Name=h['Name'],geometry=h['geometry'],population=h['population'],R_kpc=R,a_kpc=a,eta=eta,rs_kpc=h['rs_kpc'],target_half_mass_radius_over_R=float(np.interp(.5,target,xx)),cases=diagnostics,mixture_max_cumulative_fraction_error=float(np.max(np.abs(mix-target))),mixture_deposit_weights=fit.x[:-1].tolist(),max_resolution_drift=maxdrift,far_limit_cdf_difference=far_drift,r_over_R=xx.tolist(),target_cdf=target.tolist(),case_cdfs=curves.tolist(),best_mixture_cdf=mix.tolist())
    rows.append(row)
out=dict(scope='Boundary-conditioned external companion capture shape, not an absolute supply or new lens fit',k0_per_kpc=k,alpha_per_kpc=alpha,q=cp['q'],cases_D_over_R=cases,quadrature_optical_depth_max_absolute_error=max(errors),input_sha256={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in inputs},rows=rows)
(P/'stream-capture-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
lines=['# Do companion rivers deposit the required halo shape?','',
'Fixed capture parameters and the exact one-third law were used throughout. Straight incoming rays were attenuated along their paths inside each finite 5 Re diagnostic region. This is a conditional external-input shape test, not an absolute energy supply test or a new fit to observed lens images or stellar velocities.','',
'## Result','',
'Under the standard comparison geometry, source mixing reduces cumulative-profile differences from 14.46 to 2.44 percentage points for J1112, from 18.83 to 5.52 for J1621, and from 18.91 to 5.83 for J1630 (Chabrier proxy). It does not improve J0037, J1204 or J1402: their differences remain 17.97, 49.94 and 65.11 points. Their required profiles are more centrally concentrated than these straight-stream deposits. Under our retained geometry none of the six source mixtures improves on distant illumination in this basis. These are conditional shape comparisons, not observed lens residuals.','',
'The successful shape adjustments require the nearest tested shell: J1112 assigns about 93.1% of deposited power to D/R=1.1, while J1621 and J1630 put all power there. Thus they hit the source-distance boundary and have no independently measured source justification. Sources closer than 1.1 R, internal sources, different capture regions or subsequent migration were not tested; no general exclusion follows.','',
'## Executed transport','',
'Known transport: dF_c/ds=-kappa F_c after geometric dilution is factored out, with tau=integral kappa ds and q_ret=eta kappa F_c. Proposed capture profile: kappa=k0/[1+(r/a)^2]^2, retaining the previously calibrated k0 and a. The Hill expression eta=X^(1/3)/(1+X^(1/3)) is known mathematics with our proposed retention interpretation. Unretained absorbed energy belongs to another reservoir; no energy creation is assumed.','',
'For a point source, F_c(x)=L[1-exp(-alpha*d_entry)] exp(-tau)/(4*pi*d(x)^2). The entry distance is source to boundary; d(x) is source to deposition point. These rays include finite-distance divergence. The external companion component is supplied at the boundary. Conversion inside the region and internal stars are excluded from this particular test. Exterior opacity is unspecified; the boundary field is conditional. R=5 Re is a diagnostic boundary, not a measured capture cutoff.','',
'For each of 24 geometry/population combinations, sources were placed at 1.1, 2, 5, 10 and 100 boundary radii, plus a distant isotropic field. Spherical averages also represent uniformly populated source shells. A single beam can make an asymmetric deposit; a shell averages that asymmetry away. The one-third factor multiplies all rates and cancels only when comparing normalized radial shapes. An isotropic distant field has zero dipole.','',
'## Halo shape comparison','',
'Both curves are normalized to their mass or deposited energy inside R. The error below is the largest difference between cumulative fractions, in percentage points; it is not a velocity error, lens error, confidence level or statistical significance. Constant illumination history with no migration is required to identify deposited-power shape with stored-energy shape.','',
'| Galaxy | Geometry | Distant field error, pp | Best distance-mixture error, pp | Target half-mass radius / R | Distant half-deposit radius / R |',
'|---|---|---:|---:|---:|---:|']
for h in rows:
    if h['population']!='Chabrier':continue
    lines.append(f"| {h['Name']} | {h['geometry']} | {100*h['cases'][-1]['max_cumulative_fraction_error']:.2f} | {100*h['mixture_max_cumulative_fraction_error']:.2f} | {h['target_half_mass_radius_over_R']:.3f} | {h['cases'][-1]['half_deposit_radius_over_R']:.3f} |")
lines += ['', 'The mixture is deliberately generous: six nonnegative deposited-power weights are fitted separately to each target. It is an inverse diagnostic, not a measured source population or shared prediction. No source luminosities, durations, or energy budgets follow from these weights. Both population cases, every distance, angular dipoles, normalized curves and weights are retained in the JSON.','',
'These NFW profiles are our restricted target fits, not unique measured halos. Some have parameter boundaries and poor stellar fits. Standard geometry remains a comparison only; original capture scales are retained when evaluating that target, so it is not a recalibration under FLRW. Matching this enclosed spherical profile would still not prove projected lens agreement, support, or a viable energy history.','',
'## Interpretation and next decision','',
'Changing source distance changes which paths illuminate the receiver, but leaves its local capture coefficient fixed. A finite smooth capture law with nonsingular external sources gives finite central deposited density, whereas NFW has a central cusp. This does not alone exclude agreement over a finite measured range, but source mixing cannot arbitrarily prescribe the inner profile. The table quantifies the discrepancy over this declared diagnostic region.','',
'No bending, scattering, binding, or subsequent inward migration was added. If distance mixtures leave a deficit of central deposition, the next physical variable is an explicit transport/capture change that delivers energy inward. It must account for momentum and retained-state support; simply renaming rays rivers does not supply it. A target-fit mixture is not a reason to change the one-third exponent.','',
'Actual nearby-source positions for these six receivers are not available in this calculation. The cached M87 catalog was not reassigned to them. Thus this test is a controlled distance-family test against existing halo targets, not an observed source-to-halo prediction. No new galaxy, cluster, redshift, or lensing success is claimed.','',
'## Verification','',f"Independent optical-depth quadrature agrees within {max(errors):.3g}. Doubling radial/angular resolution changes cumulative fractions by at most {max(h['max_resolution_drift'] for h in rows):.3g}; the D/R=100000 limit agrees with distant illumination within {max(h['far_limit_cdf_difference'] for h in rows):.3g}. All cumulative profiles are monotonic, and mixture weights are nonnegative and sum to one. See stream-capture-protocol.md, executable stream-capture.py and input hashes in the results. All six project goals remain open."]
(P/'stream-capture-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
for ax,h in zip(axes.flat,[h for h in rows if h['population']=='Chabrier' and h['geometry']=='standard_flat_FLRW']):
    ax.plot(h['r_over_R'],h['target_cdf'],color='black',label='Fitted halo target')
    ax.plot(h['r_over_R'],h['case_cdfs'][-1],label='Distant input')
    ax.plot(h['r_over_R'],h['best_mixture_cdf'],'--',label='Target-fitted source mixture')
    ax.set(title=h['Name'],xlabel='Radius / diagnostic boundary',ylabel='Enclosed fraction',xlim=(0,1),ylim=(0,1));ax.grid(alpha=.2)
axes.flat[0].legend(fontsize=8)
fig.suptitle('Fixed capture law: normalized deposition vs standard-geometry halo targets')
fig.savefig(P/'stream-capture.png',dpi=160)
print('\n'.join(lines[14:33]))
