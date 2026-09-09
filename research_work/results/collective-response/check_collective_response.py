"""Resolved Gaussian target form factor for the existing scalar comparison.

No physical target, special clock law or graviton identity is adopted.
"""
from pathlib import Path
import json,os
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'collective-response'

def integrate(b,order=256):
    nodes,weights=leggauss(order)
    x=(nodes+1)/2; wx=weights/2; y=1-x
    if b==0:
        mu=nodes
        w=weights*(1+mu**2)/2
        theta2=np.sum(w*np.arccos(mu)**2)/np.sum(w)
        return dict(rate_relative_point=1.,mean_loss=1/3,second_loss=1/7,
                    photon_angle_rms=float(np.sqrt(theta2)),mean_one_minus_cosine=1.)
    # Integrate companion angles analytically, then use t=b(v-y) with
    # v=|n_initial-x*n_final|. This resolves the forward peak at large b.
    tn,tw=leggauss(96)
    limit=np.minimum(12.,2*b*x)
    t=limit[:,None]*(tn[None,:]+1)/2
    wt=limit[:,None]*tw[None,:]/2
    mu=1-(2*y[:,None]*t/b+(t/b)**2)/(2*x[:,None])
    mu=np.clip(mu,-1,1)
    a=(1+mu**2)/2
    difference=np.exp(-t*t)*(-np.expm1(-4*b*y[:,None]*t-4*b*b*y[:,None]**2))
    # Overall 2pi^2/b^3 cancels from normalized moments. The x*y factors
    # in the angular integral cancel two factors in x^3*y phase space.
    base=np.sum(wt*a*difference,axis=1)
    density=x*x*base
    total=np.sum(wx*density)
    # Point contact integral: (32pi^2/3)*integral x^3(1-x) dx = 8pi^2/15.
    ratio=(15/(4*b**3))*total
    moment_angle=np.sum(wx*x*x*np.sum(wt*a*difference*np.arccos(mu)**2,axis=1))/total
    moment_cos=np.sum(wx*x*x*np.sum(wt*a*difference*(1-mu),axis=1))/total
    return dict(rate_relative_point=float(ratio),mean_loss=float(np.sum(wx*density*y)/total),
                second_loss=float(np.sum(wx*density*y*y)/total),
                photon_angle_rms=float(np.sqrt(moment_angle)),mean_one_minus_cosine=float(moment_cos))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    fourier=[]
    for q in [0.,.3,1.,3.]:
        value=2*quad(lambda z:np.exp(-z*z/2)*np.cos(q*z)/np.sqrt(2*np.pi),0,12,epsabs=1e-12)[0]
        expected=np.exp(-q*q/2)
        assert abs(value-expected)<1e-12
        fourier.append(dict(dimensionless_momentum=q,amplitude=value,expected=float(expected)))
    results=[]; convergence=[]
    for b in [0.,.01,.1,1.,10.,100.,1000.]:
        low=integrate(b,256); high=integrate(b,512)
        errors={key:abs(low[key]-high[key])/max(abs(high[key]),1e-30) for key in high}
        assert max(errors.values())<2e-5,(b,errors)
        high['energy_times_size']=b
        high['large_b_rate_asymptote']=None if b==0 else float(5*np.sqrt(np.pi)/(8*b**3))
        high['large_b_angle_rms_asymptote']=None if b==0 else float(1/np.sqrt(b*np.sqrt(np.pi)))
        results.append(high); convergence.append(dict(b=b,relative_differences=errors))
    assert abs(results[1]['rate_relative_point']-1)<1e-3
    assert abs(results[-1]['mean_loss']-.25)<.003
    assert abs(results[-1]['second_loss']-.1)<.003
    assert abs(results[-1]['rate_relative_point']/results[-1]['large_b_rate_asymptote']-1)<.02
    assert abs(results[-1]['photon_angle_rms']/results[-1]['large_b_angle_rms_asymptote']-1)<.02
    colors=[]
    for reference_b in [0.,1.,10.,100.]:
        ref=integrate(reference_b)
        rows=[]
        for energy in [.5,1.,2.,4.]:
            actual=integrate(reference_b*energy)
            alpha_ratio=energy**6*actual['rate_relative_point']*actual['mean_loss']/(ref['rate_relative_point']*ref['mean_loss'])
            rows.append(dict(energy_ratio=energy,initial_fractional_loss_ratio=float(alpha_ratio)))
        colors.append(dict(reference_energy_times_size=reference_b,predictions=rows))
    # Collinear outgoing momenta have Q=0 at every energy split, even for a
    # very large target. Spatial coherence by itself does not suppress them.
    collinear=[]
    for retained in [.1,.5,.9,.999999999]:
        q=1-retained-(1-retained)
        assert q==0
        collinear.append(dict(retained_fraction=retained,recoil_momentum=0.,form_factor_squared=1.))
    # Illustrative large-b size for a one-arcsecond RMS angle PER EVENT.
    # This does not establish an observed bound, a viable target, or many-event blur.
    angle=np.pi/(180*3600)
    b_required=1/(np.sqrt(np.pi)*angle**2)
    size_m=b_required*1.973269804e-7  # hbar*c in eV m at E=1 eV
    result=dict(scope='Conditional spatial coherence correction to the scalar polarizability candidate; no adopted physical target or complete theory',
                gaussian_fourier_checks=fourier,form_factor_results=results,quadrature_convergence=convergence,
                color_dependence=colors,collinear_energy_splits=collinear,
                illustrative_1ev_1arcsec_per_event=dict(required_b=float(b_required),gaussian_size_m=float(size_m)),
                checks=dict(gaussian_transform=True,quadrature_refinement=True,point_limit=True,
                            large_target_spectrum_and_rate=True,large_target_angular_scaling=True,
                            spatial_coherence_does_not_force_soft_energy=True),
                limitations=['Scalar companion comparison; no spin-two amplitude.',
                             'Coherent elastic internal target channel, negligible center-of-mass recoil energy, fixed total polarizability coefficient.',
                             'Gaussian density profile stipulated; binding, constituent density, causal response and population not derived.',
                             'No temporal susceptibility, inverse reactions, transient dilation, capture or gravity response calculated.'])
    (OUT/'collective-response-results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],form_factor_results=results,color_dependence=colors,
                         illustrative_size=result['illustrative_1ev_1arcsec_per_event']),indent=2))

if __name__=='__main__': main()
