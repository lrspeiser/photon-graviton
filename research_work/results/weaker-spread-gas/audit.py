#!/usr/bin/env python3
"""Independent readback and analytic checks of JR-4; no refitting."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--jr1',type=Path,required=True);args=ap.parse_args();p=args.results
 checks=[]
 def ck(name,condition,error=None):checks.append(dict(name=name,passed=bool(condition),error=None if error is None else float(error)))
 allm=json.loads((p/'all-model-predictions.json').read_text());summary=json.loads((p/'summary.json').read_text());manifest=json.loads((p/'manifest.json').read_text())
 for name,m in allm.items():
  for row in m['galaxies']:
   y=np.array(row['observed_kms']);v=np.array(row['predicted_kms']);e=np.array(row['errors_kms']);f=(v-y)/y
   values={'fractional_RMS':np.sqrt(np.mean(f*f)),'RMSE_kms':np.sqrt(np.mean((v-y)**2)),'chi2':np.sum(((v-y)/e)**2),'mean_fractional':f.mean()}
   for key,value in values.items():
    er=abs(value-row[key])/max(1,abs(value));ck(name+'/'+row['name']+'/'+key,er<1e-12,er)
  rows=m['galaxies']
  for group in ('all','train','validation','test'):
   rr=rows if group=='all' else [r for r in rows if r['split']==group];s=m['stats'][group]
   ck(name+'/'+group+'/below20',sum(r['fractional_RMS']<.2 for r in rr)==s['below20'])
   ck(name+'/'+group+'/count',len(rr)==s['galaxies'])
  for row in m['lenses']:
   y=np.array(row['observed_vrms']);v=np.array(row['predicted_vrms']);cov=np.array(row['covariance']);e=v-y
   chi=e@np.linalg.solve(cov,e);ck(name+'/'+row['name']+'/lens_chi2',abs(chi-row['stellar_chi2'])<1e-8,chi-row['stellar_chi2'])
   rms=np.sqrt(np.mean((e/y)**2));ck(name+'/'+row['name']+'/lens_rms',abs(rms-row['stellar_fractional_RMS'])<1e-12,rms-row['stellar_fractional_RMS'])
 for f,h in manifest['inputs'].items():ck('unchanged/'+f,hashlib.sha256((args.jr1/'repository'/f).read_bytes()).hexdigest()==h)
 # The density identity is algebraic, independent of the implemented force sampler.
 # M(r)=A*r/G*[(r/rc)^q/(1+(r/rc)^q)]/sqrt(1+(r/rt)^2).
 # d ln M/d ln r = 1/(1+(r/rt)^2) + q/(1+(r/rc)^q)>0.
 for q in (1.6760087204080696,1.8906354271311694):
  r=np.geomspace(1e-6,1e6,501);slope=1/(1+(r/20)**2)+q/(1+(r/2)**q);ck('positive_mass_slope/'+str(q),np.all(slope>0))
  for lam in (1,1.5,3,10):
   mass=lambda r,A,rc,rt:A*r*((r/rc)**q/(1+(r/rc)**q))/np.sqrt(1+(r/rt)**2)
   left=mass(r,1/lam,2*lam,20*lam);right=mass(r/lam,1,2,20)
   err=np.max(abs(left/right-1));ck('dilation_identity/'+str(q)+'/'+str(lam),err<1e-12,err)
 gas=json.loads((p/'gas-physics-controls.json').read_text())
 for r in gas['rows']:ck('collapse/'+str(r['s']),abs(r['relative_error'])<1e-8,r['relative_error'])
 selected=allm[summary['selected']];base=allm['baseline'];inv=json.loads((p/'individual-diagnostics.json').read_text())
 corr={}
 for mode,key in [('total_strength','strength'),('companion_strength','strength'),('spreading','spread')]:
  x=np.array([r['models'][mode][key] for r in inv['rows']]);gas=np.array([r['fgas'] for r in inv['rows']]);sigma=np.array([r['Sigma'] for r in inv['rows']])
  corr[mode]=dict(Spearman_gas_fraction=float(spearmanr(x,gas).statistic),Spearman_surface_density=float(spearmanr(x,sigma).statistic))
 improved=sum(b['fractional_RMS']>=.2 and n['fractional_RMS']<.2 for b,n in zip(base['galaxies'],selected['galaxies']))
 harmed=sum(b['fractional_RMS']<.2 and n['fractional_RMS']>=.2 for b,n in zip(base['galaxies'],selected['galaxies']))
 equiv=max(abs(np.array(a['predicted_kms'])-np.array(b['predicted_kms'])).max() for a,b in zip(allm['spread_only']['galaxies'],allm['weaker_spread_companion']['galaxies']))
 result=dict(checks=len(checks),failures=sum(not c['passed'] for c in checks),rows=checks,correlations=corr,selected_threshold_improved=improved,selected_threshold_harmed=harmed,pure_vs_combined_spread_max_velocity_difference_kms=float(equiv),scope='Independent metric arithmetic, matrix covariance calculation, analytic mass/dilation checks, original-byte hashes and cloud collapse identity; not a physical validation')
 (p/'readback-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
 if result['failures']:raise SystemExit(1)
if __name__=='__main__':main()
