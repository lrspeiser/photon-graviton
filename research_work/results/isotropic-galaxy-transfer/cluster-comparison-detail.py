"""Coma full-profile and leave-one-bin-out diagnostic; same shape assumptions."""
import ast
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent
# Load only the already tested projection function, without rerunning its driver.
tree=ast.parse((HERE/'cluster-model-comparison.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='shape'],type_ignores=[]),'<existing shape>','exec'))
raw=json.loads((HERE.parent/'cluster-observation-readiness/kubo-figure-data.json').read_text())
rows=raw['rows'];R=np.array([r['published_radius_h_inverse_Mpc'] for r in rows]);R/=R[0]
y=np.array([r['shear_t'] for r in rows]);err=np.array([r['plotted_sigma_t'] for r in rows])
control=json.loads((HERE.parent/'cluster-observation-readiness/kubo-controls-results.json').read_text())['blank_rows']
yc=y-np.array([r['shear_t'] for r in control]);ec=np.hypot(err,[r['sigma_t'] for r in control])

def nfw_exact(x):
    # Known projected NFW formula: DeltaSigma/(rho_s*r_s).
    x=np.asarray(x);A=np.empty_like(x)
    lo=x<1;hi=x>1;eq=x==1
    A[lo]=np.arccosh(1/x[lo])/np.sqrt(1-x[lo]**2)
    A[hi]=np.arccos(1/x[hi])/np.sqrt(x[hi]**2-1)
    A[eq]=1
    f=np.empty_like(x);f[~eq]=(1-A[~eq])/(x[~eq]**2-1);f[eq]=1/3
    return 4*(np.log(x/2)+A)/x**2-2*f

checkx=np.geomspace(.02,500,201)
analytic_error=float(np.max(abs(nfw_exact(checkx)/shape(checkx,'nfw')-1)))
assert analytic_error<1e-6

def companion(T):
    # Identical fixed-T density and projection to coma.py.
    x=np.geomspace(1e-7,1e8,16001);mu,w=leggauss(96)
    if T:
        u=x[:,None]*mu;B2=1+x[:,None]**2*(1-mu**2);B=np.sqrt(B2)
        tau=(2*T/np.pi)*(u/(2*B2*(B2+u*u))+(np.arctan(u/B)+np.pi/2)/(2*B**3))
        J=.5*np.sum(np.exp(-np.maximum(tau,0))*w,axis=1)
    else:J=np.ones_like(x)
    logrho=PchipInterpolator(np.log(x),np.log(np.maximum(J/(1+x*x)**2,1e-300)))
    t,wt=leggauss(160);theta=(t+1)*np.pi/4;wt*=np.pi/4
    b=np.geomspace(1e-5,1e5,6001);s=np.sqrt(1+b*b)[:,None]
    rad=np.sqrt(b[:,None]**2+(s*np.tan(theta))**2)
    sigma=2*np.sum(np.exp(logrho(np.log(rad)))*s/np.cos(theta)**2*wt,axis=1)
    integ=sigma[0]*b[0]**2/2+cumulative_trapezoid(sigma*b,b,initial=0)
    F=PchipInterpolator(np.log(b),2*integ/b**2-sigma)
    return lambda b:F(np.log(b))

functions={'transparent':companion(0),'strong_interception':companion(100),
           'nfw':lambda b:shape(b,'nfw'),'mond_point_baryons':lambda b:shape(b,'mond_point_baryons'),
           'point_mass_control':lambda b:b**-2}
grid=np.linspace(-2,2,401)
# Precompute a fine scale grid once; optimize the actual function locally.
grids={k:np.array([f(R/10**a) for a in grid]) for k,f in functions.items()}

def fit(kind,indices,obs,errors):
    idx=np.asarray(indices);fn=functions[kind]
    def evaluate(a):
        f=fn(R/10**a);amp=max(0.,float(np.sum(f[idx]*obs[idx]/errors[idx]**2)/np.sum((f[idx]/errors[idx])**2)))
        pred=amp*f;return float(np.sum(((pred[idx]-obs[idx])/errors[idx])**2)),amp,pred
    if kind=='point_mass_control':a=0.
    else:
        f=grids[kind][:,idx];amp=np.maximum(0,np.sum(f*obs[idx]/errors[idx]**2,axis=1)/np.sum((f/errors[idx])**2,axis=1))
        losses=np.sum(((amp[:,None]*f-obs[idx])/errors[idx])**2,axis=1)
        i=int(np.argmin(losses));opt=minimize_scalar(lambda a:evaluate(a)[0],bounds=(grid[max(0,i-1)],grid[min(400,i+1)]),method='bounded')
        a=min([-2.,2.,float(opt.x)],key=lambda a:evaluate(a)[0])
    loss,amp,pred=evaluate(a)
    return dict(training_sum=loss,scale=10**a,amplitude=amp,boundary=abs(abs(a)-2)<1e-4,predictions=pred.tolist(),
                residual_sigma=((pred-obs)/errors).tolist(),all_sum=float(np.sum(((pred-obs)/errors)**2)))

out={'scope':'Post-comparison diagnostic, not blind prediction or a physical cluster transfer',
     'nfw_analytic_relative_error':analytic_error,'datasets':{}}
for label,obs,errors in [('published',y,err),('blank_subtracted_sensitivity',yc,ec)]:
    result={}
    for kind in functions:
        full=fit(kind,range(6),obs,errors);inner=fit(kind,range(3),obs,errors)
        loo=[fit(kind,[j for j in range(6) if j!=i],obs,errors) for i in range(6)]
        result[kind]=dict(full=full,inner=inner,leave_one_out=loo,
                         loo_sum=float(sum(m['residual_sigma'][i]**2 for i,m in enumerate(loo))),
                         inner_fit_outer_sum=float(sum(v*v for v in inner['residual_sigma'][3:])))
    out['datasets'][label]=result
# Verify prior fits reproduced rather than silently changing the models.
old=json.loads((HERE/'cluster-model-comparison-results.json').read_text())['models']
for m in old[:4]:
    new=out['datasets']['published'][m['model']]['inner']
    assert abs(new['training_sum']-m['training_diagonal_sum'])<1e-5
    assert np.max(abs(np.array(new['predictions'])-m['predicted_shear']))<1e-7
(HERE/'cluster-comparison-detail-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
lines=['# Coma comparison: why the NFW result looked worse','',
       'The earlier outer-only score is not a full-profile goodness-of-fit score. The authors fitted all six bins; we initially fitted only three. This diagnostic retains every measurement, including the negative outer shear, and uses the same two-parameter shapes. Point mass is a diagnostic control, not a physical cluster candidate.','']
for label,result in out['datasets'].items():
    lines += ['## '+label,'','| Shape | Fit all six | Predict outer three from inner three | Leave-one-bin-out sum |','|---|---:|---:|---:|']
    for kind,m in result.items():lines.append(f"| {kind} | {m['full']['all_sum']:.4f} | {m['inner_fit_outer_sum']:.4f} | {m['loo_sum']:.4f} |")
lines += ['','## Interpretation','',
          'NFW is competitive: its full-profile score is 3.8549 versus 3.6839/3.7245 for the companion shapes. Its leave-one-out score 6.8209 is essentially tied with transparent capture 6.8573 and better than strong interception 10.7919. These differences do not establish superiority. The point-mass control deteriorates to 32.1738 on leave-one-out prediction.','',
          'For the original inner-only NFW fit, outer bins 4 and 6 lie 1.28 and 1.74 plotted standard errors below its predictions; bin 5 differs by only 0.20. NFW fitted to all six prefers scale 0.431 first-radius units, versus 1.006 from the inner three. This is why extrapolating the inner fit looked worse: the three inner bins do not fix the scale tightly. It is not a failed NFW projection.','',
          'All sums use diagonal plotted errors. Lower is closer under these assumptions; leave-one-out predicts each omitted bin using five others and sums those six errors. It is an interpolation/stability diagnostic, not an independent next-cluster prediction. No probability or model-selection significance is claimed. The blank-subtracted rows assume an additive control offset and independent errors, and are sensitivity results rather than a mandatory correction.','',
          'The published NFW full-profile chi-square is 3.87 for four degrees of freedom. Our reconstructed-bin, bin-center fit is only an approximate reproduction; exact source weighting/bin averaging is unavailable. We independently checked the numerical NFW projection against the known analytic projected profile, with maximum relative error '+f'{analytic_error:.3g}'+'. The prior inner fits are reproduced to the checks recorded in the script.','',
          'A point mass places all lens mass inside the measured radii and has shear proportional to 1/R^2. It can be an exterior approximation, but these six points do not demonstrate that Coma is a point mass or that its extended gas and galaxy motions can be explained that way. Keep it out of the physical candidate ranking.','',
          'The companion alternatives and compact-baryon MOND remain freely normalized shapes. Neither is a baryon-calibrated cluster prediction; the one-third retention factor is still absorbed by the fitted amplitude. See the previous protocol.','',
          'Sources: [Kubo et al., section 3.2](https://arxiv.org/pdf/0709.0506); [known analytic NFW lensing, Wright and Brainerd](https://arxiv.org/abs/astro-ph/9908213). All six research goals remain open.']
(HERE/'cluster-comparison-detail-report.md').write_text('\n'.join(lines)+'\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
colors={'transparent':'#167d8d','strong_interception':'#c47d17','nfw':'#3657bd','mond_point_baryons':'#974d9d'}
labels={'transparent':'Companion: transparent','strong_interception':'Companion: interception','nfw':'NFW','mond_point_baryons':'MOND: compact baryons'}
rad=np.linspace(.85,8.2,300)
axes[0].errorbar(R,y*1000,yerr=err*1000,fmt='ko',capsize=3,label='Coma measurements')
for kind,color in colors.items():
    m=out['datasets']['published'][kind];p=m['full']
    axes[0].plot(rad,p['amplitude']*functions[kind](rad/p['scale'])*1000,color=color,label=labels[kind])
    axes[1].plot(range(1,7),[v['residual_sigma'][i] for i,v in enumerate(m['leave_one_out'])],'-o',color=color,label=labels[kind])
axes[0].axhline(0,color='gray',lw=.7);axes[0].set(xlabel='Radius / first measured radius',ylabel='Tangential shear × 1000',title='Fit all six measurements')
axes[0].legend(fontsize=8)
axes[1].axhline(0,color='gray',lw=.7);axes[1].axhspan(-1,1,color='gray',alpha=.1)
axes[1].set(xlabel='Omitted bin (fit the other five)',ylabel='Prediction minus observation / plotted error',title='Predict each omitted measurement')
fig.suptitle('Coma: exploratory shape comparison; free amplitudes and scales',fontsize=12)
fig.savefig(HERE/'cluster-comparison-detail.png',dpi=160)
plt.close(fig)
for label,result in out['datasets'].items():
    print(label)
    for kind,m in result.items():print(kind,round(m['full']['all_sum'],4),round(m['loo_sum'],4),m['full']['boundary'])
