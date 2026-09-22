#!/usr/bin/env python3
"""Independent arithmetic readback of JR-3 results; does not refit models."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_triangular

def main():
 ap=argparse.ArgumentParser(__doc__);ap.add_argument('--jr2',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 target=a.output/'readback-audit.json'
 if target.exists():raise FileExistsError('Preserve earlier audit')
 def read(p):return json.loads(p.read_text())
 original=read(a.jr2/'original-JR1/evidence/R10_spheroid_stellar_population.json')['metrics']
 obs={r['name']:r for r in original['sparc_rows']};lobs={r['name']:r for r in original['lens_rows']}
 checks=[]
 def check(name,value,tolerance):
  checks.append(dict(name=name,value=float(value),tolerance=tolerance,passed=bool(np.isfinite(value) and value<=tolerance)))
 result=read(a.output/'disk-results.json')
 for row in result:
  d=obs[row['name']];y=np.asarray(d['observed_kms']);err=np.asarray(d['error_kms'])
  for kind,v in row['cases'].items():
   res=np.asarray(v['predicted_kms'])-y
   check(row['name']+':'+kind+':frms',abs(np.sqrt(np.mean((res/y)**2))-v['fractional_RMS']),1e-12)
   chi=np.sum((res/err)**2);check(row['name']+':'+kind+':chi',abs(chi-v['chi2'])/max(chi,1),1e-12)
   check(row['name']+':'+kind+':transmission',abs((1-v['net_dimming_fraction'])*v['intrinsic_light_factor']-1),1e-12)
 for row in read(a.output/'lens-results.json'):
  d=lobs[row['name']];y=np.asarray(d['observed_vrms_kms']);L=np.linalg.cholesky(d['covariance'])
  for kind,v in row['cases'].items():
   diff=np.asarray(v['predicted_vrms'])-y;w=solve_triangular(L,diff,lower=True)
   check(row['name']+':'+kind+':frms',abs(np.sqrt(np.mean((diff/y)**2))-v['stellar_fractional_RMS']),1e-12)
   check(row['name']+':'+kind+':chi',abs(w@w-v['stellar_chi2'])/max(w@w,1),1e-12)
   if kind.startswith('ring_'):check(row['name']+':'+kind+':ring',abs(v['theta_predicted']/v['theta_observed']-1),1e-8)
 for row in read(a.output/'reversible-shared.json')['rows']:
  d=obs[row['name']];y=np.asarray(d['observed_kms']);r=np.asarray(row['predicted_kms'])-y
  check(row['name']+':shared_reverse:frms',abs(np.sqrt(np.mean((r/y)**2))-row['fractional_RMS']),1e-12)
  T=1-row['net_dimming_fraction'];X=row['incoming_companion_flux_per_primary_photon_flux'];Z=row['outgoing_companion_flux_per_primary_photon_flux']
  check(row['name']+':shared_reverse:energy',abs(T+Z-1-X),1e-12)
 manifest=read(a.output/'input-manifest.json')['sha256']
 for name,sha in manifest.items():check('input hash '+name,0 if hashlib.sha256((a.jr2/name).read_bytes()).hexdigest()==sha else 1,0)
 o=dict(check_count=len(checks),failures=sum(not c['passed'] for c in checks),max_recorded_residual=max(c['value'] for c in checks),checks=checks)
 target.write_text(json.dumps(o,indent=2)+'\n');print({k:v for k,v in o.items() if k!='checks'})
 if o['failures']:raise SystemExit(1)
if __name__=='__main__':main()
