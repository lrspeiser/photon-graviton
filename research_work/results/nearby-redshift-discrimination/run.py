"""Conditional discrimination sensitivity; no re-fit or new galaxy outcomes."""
from pathlib import Path
import json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
source=HERE.parent/'nearby-redshift-crosscheck/results.json';protocol=HERE.parent/'nearby-redshift-crosscheck/protocol.json'
r=json.loads(source.read_text());p=json.loads(protocol.read_text());rows=[]
c=p['c_kms'];alpha=p['alpha_per_mpc']
for x in r['rows']:
 d=x['distance_mpc'];sd=x['published_distance_error_mpc']
 delta=c*(x['exponential_predicted_z']-x['linear_predicted_z'])
 derivative=c*alpha*np.exp(alpha*d)
 sigma_distance=derivative*sd
 sigma_spectral=x['solar_to_cmb_factor']*x['statistical_center_error_kms']
 # Finite symmetric perturbation checks linear propagation, without changing adopted distances.
 finite=c*(np.expm1(alpha*(d+sd))-np.expm1(alpha*(d-sd)))/2
 assert abs(finite/sigma_distance-1)<1e-6
 curvature=c*(np.expm1(alpha*d)-alpha*d)
 rows.append({'galaxy':x['target_name'],'distance_mpc':d,'fixed_rule_difference_kms':delta,
    'exponential_curvature_above_its_own_tangent_kms':curvature,'distance_error_projected_kms':sigma_distance,
    'spectral_center_error_projected_kms':sigma_spectral,
    'rule_difference_over_distance_and_spectral_error':abs(delta)/np.hypot(sigma_distance,sigma_spectral)})
delta=np.array([x['fixed_rule_difference_kms'] for x in rows]);sd=np.array([x['distance_error_projected_kms'] for x in rows]);sv=np.array([x['spectral_center_error_projected_kms'] for x in rows])
scenarios=[]
for motion in [0.,100.,300.]:
 for rho in [0.,.5,.9]:
  covariance=np.diag((1-rho)*sd**2+sv**2+motion**2)+rho*np.outer(sd,sd)
  separation=float(delta@np.linalg.solve(covariance,delta))
  assert separation>=0 and np.linalg.eigvalsh(covariance).min()>0
  scenarios.append({'illustrative_independent_motion_sigma_kms':motion,'illustrative_distance_error_correlation':rho,
    'squared_mahalanobis_rule_separation':separation,'sqrt_separation':float(np.sqrt(separation)),
    'expected_log_likelihood_ratio_if_one_fixed_rule_true':separation/2})
# Independent simulation of the equal-covariance simple-hypothesis identity.
cov=np.diag(sd**2+sv**2);inv_delta=np.linalg.solve(cov,delta)
rng=np.random.default_rng(5301);noise=rng.normal(size=(100000,8))*np.sqrt(np.diag(cov))
llr=.5*(delta@inv_delta)+noise@inv_delta
expected=.5*(delta@inv_delta);se=float(llr.std(ddof=1)/np.sqrt(len(llr)))
assert abs(llr.mean()-expected)<5*se
out={'scope':'Post-evaluation test-design sensitivity using already-exposed eight-galaxy inputs; no new observations, parameter fitting, motion estimates or significance claim.',
 'formula_provenance':'Known first-order error propagation and Gaussian simple-hypothesis likelihood algebra, conditional on assumed covariance.',
 'rows':rows,'scenarios':scenarios,'gaussian_identity_check':{'draws':len(llr),'expected_llr':float(expected),'simulated_mean_llr':float(llr.mean()),'monte_carlo_standard_error':se},
 'source_hashes':{str(x.relative_to(HERE.parents[2])):hashlib.sha256(x.read_bytes()).hexdigest() for x in [source,protocol]}}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'individual_ratio_range':[min(x['rule_difference_over_distance_and_spectral_error'] for x in rows),max(x['rule_difference_over_distance_and_spectral_error'] for x in rows)],'scenarios':scenarios},indent=2))
