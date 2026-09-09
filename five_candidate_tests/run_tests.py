"""Five candidate families: explicit restricted proxies, not five completed theories.
Real released/derived data; all exploratory and previously inspected.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit,least_squares,brentq
P=Path(__file__).resolve().parent
def fit(z,w,e):
    bb,cv=curve_fit(lambda z,b:(1+z)**b,z,w,p0=[1.],sigma=e,absolute_sigma=True)
    return dict(b=float(bb[0]),sigma_formal=float(np.sqrt(cv[0,0])),
        rmse_stretch=float(np.sqrt(np.mean((w-1-z)**2))),
        rmse_no_stretch=float(np.sqrt(np.mean((w-1)**2))),
        chi_no_stretch=float(np.sum(((w-1)/e)**2)),
        chi_unit_stretch=float(np.sum(((w-1-z)/e)**2)),N=len(w))
d=pd.read_csv(P/'data/DES_event_averages.csv')
des=fit(d.z.to_numpy(),d.width.to_numpy(),d.err.to_numpy())
b=pd.read_csv(P/'data/DES_widths.csv')
b=b[np.isfinite(b[['z','Width','Width_err']]).all(axis=1)&(b.Width_err>0)&(b.Width_err<b.Width)].copy()
bands={}
for band in 'griz':
    s=b[b.band==band];bands[band]=fit(s.z.to_numpy(),s.Width.to_numpy(),s.Width_err.to_numpy())
# Paired event bootstrap for g-i contrast; duplicate references still unmodeled.
wide=b.pivot(index='CID',columns='band',values=['Width','Width_err','z'])
wide=wide.dropna(subset=[('Width','g'),('Width','i')]);rng=np.random.default_rng(20260909)
def contrast(w):
    vals=[]
    for k in ['g','i']:
        vals.append(fit(w['z'][k].to_numpy(),w['Width'][k].to_numpy(),w['Width_err'][k].to_numpy())['b'])
    return vals[0]-vals[1]
paired=contrast(wide);boots=[]
for _ in range(500):
    sample=wide.iloc[rng.integers(0,len(wide),len(wide))]
    boots.append(contrast(sample))

# FIRAS reconstructed from exact reference Planck function + unrounded residual.
a=np.loadtxt(P/'data/firas.txt');nu=a[:,0]*100*299792458.
h=6.62607015e-34;kb=1.380649e-23;c=299792458.
def BB(T):return 2*h*nu**3/c**2/np.expm1(h*nu/(kb*T))/1e-20
ref=BB(2.725);obs=ref+a[:,2]/1000;err=a[:,3]/1000;gal=a[:,4]/1000
def mix(delta):
    def res(t):
        T=2.725+t[0]*.001
        pred=.5*(BB(T*np.exp(delta))+BB(T*np.exp(-delta)))+t[1]*gal
        return (obs-pred)/err
    sol=least_squares(res,[0.,0.],bounds=([-725.,-100.],[275.,100.]),xtol=1e-12,ftol=1e-12,gtol=1e-10)
    return dict(log_temperature_half_separation=delta,T_center_K=float(2.725+sol.x[0]*.001),
        galaxy_coefficient=float(sol.x[1]),chi_squared=float(sol.fun@sol.fun))
base=mix(0.);mixtures=[mix(x) for x in [.001,.01,.1]]
limit=brentq(lambda x:mix(x)['chi_squared']-base['chi_squared']-3.841458820694124,1e-5,.05)
# Representative detailed-balance birth/death rates for one boson mode.
# Gamma_minus/Gamma_plus=exp(x) gives f_eq=1/(exp(x)-1).
xx=np.array([.1,1.,5.,10.]);f=1/np.expm1(xx)
collision=1.*(1+f)-np.exp(xx)*f
out=dict(status='Restricted candidates; data diagnostics, analytical tests and synthetic checks distinguished',
    candidate1=dict(model='Stationary fixed-speed path memory: dlnnu/dell=-kappa; arrival map t_o=t_e+R/c',
        prediction='Frequency redshift but no envelope dilation',DES=des,
        verdict='Fails the width diagnostic under its published source/reference assumptions; richer time-dependent memory remains open'),
    candidate2=dict(model='Coherent evolving-index propagation n=1+gamma*t; H_wave=c|p|/n',
        prediction='Spectral and envelope factors both exp(kappa R)',
        DES=des,other_checks='Prior fixed-atom electromagnetic completion fails leading hyperfine invariance; exact shared tensor action absent',
        alternative='A kinetic normalization alone changes leading amplitude, not the null cone or leading geometric-optics frequency on a static metric'),
    candidate3=dict(model='Band-dependent envelope law width=(1+z)^b_band; a proxy for chromatic temporal transport, not a spectral-redshift model',
        bands=bands,paired_g_i_N=len(wide),paired_g_minus_i_b=paired,
        paired_event_bootstrap_95_g_minus_i=np.percentile(boots,[2.5,97.5]).tolist(),
        limits='No interband/reference covariance provided; bootstrap addresses event sampling only. No inference across atomic, radio or gravitational-wave frequencies.'),
    candidate4=dict(model='Separate material and propagation metrics',
        conformal_subfamily='Prior derivation: stable local matter ratios and shared characteristics, but fixed material scale implies zero homogeneous redshift',
        genuinely_nonconformal='No full action satisfying matter, tensor and stable background requirements; no new fit score is assignable',
        previous_brightness={'flat_chi_squared':898.4745863939072,'negative_chi_squared':873.041277470174,'N':960},
        interpretation='The earlier geometry fit is kinematic and does not validate this unspecified action'),
    candidate5=dict(model='Thermal reservoir with detailed balance; compare blackbody and equal-weight mixtures of redshifted thermal components',
        FIRAS_channels=len(nu),single_temperature=base,equal_weight_mixtures=mixtures,
        conditional_profile_threshold_delta=limit,
        threshold_convention='Delta chi2=3.841 relative to delta=0; approximate diagnostic, not exact boundary-calibrated confidence limit',
        equilibrium_collision_max_abs=float(max(abs(collision))),
        identifiability='Equilibrium fit does not identify reservoir origin, rate, temperature derivation, kappa or energy budget',
        caveat='FIRAS released residuals and diagonal errors only. No full covariance, absolute temperature recalibration or angular CMB test.'),
    common_caveats=['DES derived widths inherit reference-curve dilation assumptions; no raw-light-curve fit or independent confirmation.',
        'All data reused; no newly blinded observations.',
        'No candidate has a completed stable nonexpanding matter-light-gravity action.',
        'Stationary transfer with fixed path maps delayed copies, not a stretch of the entire time axis; arbitrary broadening is not the same observable.'])
(P/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
