"""Exploratory Coma shape comparison; not a frozen galaxy-to-cluster test."""
import json
import hashlib
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
source = HERE.parent/'cluster-observation-readiness/kubo-figure-data.json'
rows = json.loads(source.read_text())['rows']
R = np.array([r['published_radius_h_inverse_Mpc'] for r in rows])
R /= R[0]
y = np.array([r['shear_t'] for r in rows])
err = np.array([r['plotted_sigma_t'] for r in rows])

def shape(b, kind, n=160):
    # Equal-potential weak lensing: gamma proportional to alpha/b-alpha'.
    # With z=b*tan(theta), shape=2 integral cos(theta)*(g-r*g') dtheta.
    nodes, weights = leggauss(n)
    theta = (nodes+1)*np.pi/4
    weights = weights*np.pi/4
    r = np.asarray(b)[:,None]/np.cos(theta)
    if kind == 'nfw':
        m = np.log1p(r)-r/(1+r)
        # g=m/r^2, r*g'=1/(1+r)^2-2*m/r^2.
        kernel = 3*m/r**2-1/(1+r)**2
    elif kind == 'mond_point_baryons':
        v = 1/r**2
        root = np.sqrt(v*v/4+v)
        g = v/2+root
        kernel = g+v+(v*v+2*v)/(2*root)
    elif kind == 'newton_point':
        kernel = 3/r**2
    else:
        raise ValueError(kind)
    return 2*np.sum(np.cos(theta)*kernel*weights,axis=1)

def fit(kind):
    def evaluate(logscale):
        f = shape(R/10**logscale,kind)
        amp = max(0.,float(np.sum(f[:3]*y[:3]/err[:3]**2)/np.sum((f[:3]/err[:3])**2)))
        pred = amp*f
        return float(np.sum(((pred[:3]-y[:3])/err[:3])**2)),amp,pred
    grid = np.linspace(-2,2,161)
    i = int(np.argmin([evaluate(a)[0] for a in grid]))
    opt = minimize_scalar(lambda a:evaluate(a)[0],bounds=(grid[max(0,i-1)],grid[min(160,i+1)]),method='bounded')
    a = min([-2.,2.,float(opt.x)],key=lambda a:evaluate(a)[0])
    loss,amp,pred = evaluate(a)
    assert np.max(abs(shape(R/10**a,kind,320)/shape(R/10**a,kind)-1)) < 1e-6
    return dict(model=kind,scale_in_first_radius_units=10**a,positive_amplitude=amp,
                scale_boundary=abs(abs(a)-2)<1e-5,training_diagonal_sum=loss,
                outer_diagonal_sum=float(np.sum(((pred[3:]-y[3:])/err[3:])**2)),predicted_shear=pred.tolist())

assert np.allclose(shape(R,'newton_point'),4/R**2,rtol=1e-12)
old = json.loads((HERE/'coma-results.json').read_text())
assert old['observed'] == rows
models = old['models'][:2]+[fit('nfw'),fit('mond_point_baryons'),old['models'][2]]
out = dict(scope='One cluster, six reconstructed shear bins; exploratory shape-only comparison',
           input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),observed=rows,
           training_bins=[1,2,3],prediction_bins=[4,5,6],models=models,
           limitations=['No frozen one-third normalization transfer','MOND uses compact point baryons and stipulated equal-potential lensing; free amplitude and transition scale',
                        'No gas profile or source geometry; shear rather than reduced shear approximation',
                        'No covariance or bin averaging; outer bins previously explored; not an untouched test',
                        'NFW represents a fitted total shape, not baryons plus a calibrated cosmological halo'])
# Post-fit boundary diagnostic: compact-baryon deep-MOND shear tends to 1/R.
f = 1/R
amp = float(np.sum(f[:3]*y[:3]/err[:3]**2)/np.sum((f[:3]/err[:3])**2))
pred = amp*f
out['post_fit_deep_mond_limit'] = dict(positive_amplitude=amp,
    training_diagonal_sum=float(np.sum(((pred[:3]-y[:3])/err[:3])**2)),
    outer_diagonal_sum=float(np.sum(((pred[3:]-y[3:])/err[3:])**2)),predicted_shear=pred.tolist())
(HERE/'cluster-model-comparison-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
lines=['# Coma: companion, MOND and NFW shape comparison','',
       'This is one exploratory cluster comparison, not a completed multi-cluster test or a transfer of the frozen galaxy law. Fit the first three bins and predict the last three. Smaller diagonal residual sums mean closer agreement; these are not calibrated probabilities.','',
       '| Shape | Inner residual sum | Outer residual sum | Scale at bound |',
       '|---|---:|---:|---|']
for m in models:
    lines.append(f"| {m['model']} | {m['training_diagonal_sum']:.4f} | {m['outer_diagonal_sum']:.4f} | {m.get('scale_boundary',False)} |")
lines += ['', '## Predicted and observed shear', '', '| Bin | Observed | Error | Transparent | Interception | NFW | MOND compact baryons | Point mass |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
for i in range(6):
    lines.append('| '+str(i+1)+' | '+' | '.join(f'{v:.7f}' for v in [y[i],err[i]]+[m['predicted_shear'][i] for m in models])+' |')
lines += ['', 'The noisy outer bins permit several shapes. No result here establishes a physical winner. The MOND point-baryon approximation is especially restrictive for an extended, gas-rich cluster.', '',
          'The compact MOND fit hits the lower scale bound. Its post-fit deep-MOND 1/R limit gives inner/outer sums '+f"{out['post_fit_deep_mond_limit']['training_diagonal_sum']:.4f}/{out['post_fit_deep_mond_limit']['outer_diagonal_sum']:.4f}"+'. This is a boundary diagnostic, not an independently chosen new candidate.', '',
          'See [protocol](cluster-model-comparison-protocol.md) for equations, provenance, missing inputs and the reason the one-third exponent cannot be tested by a freely normalized shear profile. All six research objectives remain open.']
(HERE/'cluster-model-comparison-report.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(models,indent=2))
