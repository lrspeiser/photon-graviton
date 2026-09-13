"""Summarize the bounded depth experiment; no parameter fitting here."""
from pathlib import Path
import ast,json
import numpy as np
from scipy.integrate import cumulative_trapezoid

P=Path(__file__).resolve().parent
c=json.loads((P/'depth-capture-results.json').read_text())
f=json.loads((P/'depth-capture-refined.json').read_text())
saved=json.loads((P.parent/'isotropic-galaxy-transfer/model-comparison-predictions.json').read_text())
obs={r['galaxy']:r for r in saved if r['model']=='baryons'}
meta={}
for line in (P.parents[2]/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt').read_text().splitlines():
 t=line.split()
 if len(t)==19:
  try:meta[t[0]]=float(t[11])
  except ValueError:pass

# Independent analytical spherical-potential check, using the solver's function.
node=next(n for n in ast.parse((P/'depth-capture.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='grav')
x=np.geomspace(1e-6,3000,640)
ns={'x':x,'cumulative_trapezoid':cumulative_trapezoid}
exec(compile(ast.Module(body=[node],type_ignores=[]),'grav-check','exec'),ns)
_,B=ns['grav'](1/(x*(1+x)**3))
exact=1/(2*(1+x));use=(x>.001)&(x<100)
potential_error=float(np.max(abs(B[use]/exact[use]-1)))
assert potential_error<.003,potential_error
base={r['galaxy']:r for r in c['baseline']}
gate=c['branches'][3]
assert gate['family']=='gate' and gate['beta']==.25
groups=[]
for lo,hi in [(0,80),(80,160),(160,float('inf'))]:
 rows=[r for r in gate['rows'] if lo<=max(obs[r['galaxy']]['predicted_kms'])<hi]
 delta=[]
 for r in rows:
  y=np.array(obs[r['galaxy']]['observed_kms'])
  delta.append(float(np.mean(np.log10(np.array(r['predicted_kms'])/y)**2)-np.mean(np.log10(np.array(base[r['galaxy']]['predicted_kms'])/y)**2)))
 groups.append(dict(baryonic_speed_bin=f'{lo}-{hi}',n=len(rows),mean_mass_ratio=float(np.mean([r['mass_ratio'] for r in rows])),mean_log_loss_change=float(np.mean(delta))))

def radial(rows):
 out={}
 for name,lo,hi in [('inner',0,1),('middle',1,3),('outer',3,float('inf'))]:
  errors=[]
  for r in rows:
   o=obs[r['galaxy']];z=np.array(o['R_kpc'])/meta[r['galaxy']]
   errors.extend((np.array(r['predicted_kms'])-o['observed_kms'])[(z>=lo)&(z<hi)])
  out[name]=dict(n=len(errors),mean_error_kms=float(np.mean(errors)))
 return out

fg=f['branches'][0];cg={r['galaxy']:r for r in gate['rows']}
refinement=max(float(np.max(abs(np.array(r['predicted_kms'])-cg[r['galaxy']]['predicted_kms']))) for r in fg['rows'])
seed=max(r['alternate_seed_max_speed_difference'] for r in fg['rows'])
assert all(r['converged'] and r['alternate_seed_converged'] for r in fg['rows'])
summary=dict(selection=c['selected'],coarse_groups_exploratory=groups,radial_bias_fine={'reference':radial(f['baseline']),'gate_0.25':radial(fg['rows'])},verification=dict(spherical_potential_relative_error=potential_error,gate_max_refinement_kms=refinement,alternate_seed_max_difference_kms=seed,ray_relative_error=f['ray_integral_relative_error'],baseline_max_archived_speed_difference_kms=f['baseline_max_archived_speed_difference']))
(P/'depth-capture-summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8',newline='\n')
lines=['# Depth-dependent capture: results','',
'The original exact-third reference remains selected. Across 149 galaxies and 3150 radii, all three families choose zero depth dependence under the declared training logarithmic loss. This rejects these bounded replacements under that criterion, not every possible depth law. All data were previously exposed; the old splits are transfer checks, not new blind evidence.','',
'## Hypothesis and provenance','',
'Let B=-Phi be positive well depth, zero at infinity, and S=B/(B+(150 km/s)^2). Test opacity multipliers 1+beta*S (boost), 1-beta+beta*S (shallow suppression), and 1+beta*(2*S-1) (both). Beta is shared and scanned over 0, 0.25, 0.5, 1. The saturating function, Newtonian potential and attenuation integrals are known mathematics. Their use as a companion-capture law is a proposed phenomenological modification, not an established interaction or proven novelty. The original one-third retention exponent and source normalization stay fixed.','',
'The same multiplier changes both absorption at a location and interception along incoming rays. The deposited potential feeds back into capture until a stationary fixed point converges. Ordinary-matter midplane depth is integrated from the saved rotation baseline and extended spherically, with a finite-mass outer tail. This is not a full disk/bulge potential, physical formation history, energy-supply measurement or stability proof. Stronger capture may increase total stored inventory; it is not mass-preserving redistribution.','',
'## Same-resolution comparison','',
'RMSE is the square root of the mean of per-galaxy mean squared speed errors; log RMS weights proportional errors. Selection uses training log RMS. Every trial is reported.','',
'| Rule | beta | Train RMSE | Validation RMSE | Test RMSE | Train log RMS |','|---|---:|---:|---:|---:|---:|']
for label,beta,s in [('Reference',0,c['baseline_scores'])]+[(b['family'],b['beta'],b['scores']) for b in c['branches']]:
 lines.append(f"| {label} | {beta:g} | {s['train']['RMSE_kms']:.3f} | {s['validation']['RMSE_kms']:.3f} | {s['test']['RMSE_kms']:.3f} | {s['train']['log_RMS']:.6f} |")
lines+=['','Boosting deeper wells increases errors. Weak shallow suppression offers a small km/s improvement but worsens training and validation proportional errors. Picking the lowest test RMSE would change the selection rule after observing the outcome.','',
'An exploratory grouping by maximum ordinary-matter predicted speed explains the tradeoff: the weak gate retains about 83%, 91% and 93% of reference deposited mass in the <80, 80-160 and >=160 km/s groups. The weakest group loses the most inventory and its average proportional error grows. These group boundaries are post hoc diagnostics, not fitted or independently validated physics.','',
'## Refinement and radial behavior','',
'Because every family chose zero, the weak gate was refined only as a diagnostic of its metric-dependent gain; it does not replace the selected reference.','',
'| Fine calculation | Train RMSE | Validation RMSE | Test RMSE |','|---|---:|---:|---:|']
for label,s in [('Reference',f['baseline_scores']),('Gate beta=0.25',fg['scores'])]:
 lines.append(f"| {label} | {s['train']['RMSE_kms']:.3f} | {s['validation']['RMSE_kms']:.3f} | {s['test']['RMSE_kms']:.3f} |")
lines+=['','Point-weighted mean speed errors for R/Rd <1, 1-3, >=3:']
for label,s in summary['radial_bias_fine'].items():
 lines.append(f"- {label}: "+', '.join(f"{k} {v['mean_error_kms']:+.3f} km/s" for k,v in s.items())+'.')
lines+=['',f'All coarse branches converge. Fine gate solutions also converge from initial density factors 1 and 0.1; maximum speed difference is {seed:.6g} km/s. Maximum coarse/fine gate speed change is {refinement:.6g} km/s. Fine reference differs from archived reference by at most {f["baseline_max_archived_speed_difference"]:.6g} km/s. Analytic ray relative error is {f["ray_integral_relative_error"]:.3g}; an independent analytic spherical-potential profile agrees within {potential_error:.3g} relative error. Numerical convergence is not dynamical stability.','',
'## Interpretation and next step','',
'Well depth can enter the math, but deepening the inner well also changes the gravity felt farther out. Increasing capture is not guaranteed to correct opposite inner and outer residuals. This test does not establish a shared improvement or solve lensing. Keep the unchanged reference. A useful next bounded hypothesis would redistribute a fixed captured inventory toward intermediate radii while limiting central accumulation, with its capture, support and energy rules specified before fitting. Do not add separate correction factors for each galaxy. Direction-dependent depth feedback and uncertainty in the unmeasured outer potential remain untested.','',
'Reproduce: run depth-capture.py, depth-capture.py --refine, then depth-capture-report.py. Source hashes and all predicted curves are retained in the result JSON files.']
(P/'depth-capture-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summary,indent=2))
