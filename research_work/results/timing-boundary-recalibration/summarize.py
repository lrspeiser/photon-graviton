"""Audit checkpoints; compare completed development calibration with frozen base."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.stats import beta,chi2
H=Path(__file__).resolve().parent;R=H.parents[2]
source=H/'results.json'
if not source.exists():source=R/'research_work/generated/timing-boundary-recalibration/checkpoint.json'
r=json.loads(source.read_text(encoding='utf-8'))
for name,value in r['hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==value,name
base=json.loads((H.parent/'timing-coverage-calibration/results.json').read_text(encoding='utf-8'))
old={v['label']:v for v in base['cases']}
assert len(old)==160
labels=[v['label'] for v in r['cases']];assert len(labels)==len(set(labels))
for v in r['cases']:
 o=old[v['label']];assert all(v[k]==o[k] for k in ['family','snr','truth_b','array_sha256'])
 p=R/'research_work/generated/timing-coverage-calibration'/(v['label']+'.npz')
 assert hashlib.sha256(p.read_bytes()).hexdigest()==v['array_sha256']
 q=v['fit']['best']['parameters'];a,b,s=q
 assert np.all(np.isfinite(q)) and np.log(5)-1e-10<=a<=np.log(100)+1e-10 and -2<=b<=3 and 0<=s<=.6
 assert v['accepted_95']==bool(v['valid'] and v['lr'] is not None and v['lr']<=chi2.ppf(.95,1))
 if v['valid']:assert v['raw_lr']>=-2e-6 and v['fit']['best']['maximum_outside']<=1.00001e-6
print(json.dumps(dict(audited_cases=len(labels),total=160,source_and_array_hashes_verified=True,complete=len(labels)==160)))
if len(labels)!=160:
 if '--audit' in sys.argv:sys.exit(0)
 raise SystemExit('Incomplete: no final frequency or comparison report written')
assert set(labels)==set(old)
p=json.loads((H/'protocol.json').read_text(encoding='utf-8'));gate=p['screens']
def ci(k,n):return [float(beta.ppf(.025,k,n-k+1)) if k else 0.,float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]
cells=[]
for key in sorted(set((v['family'],v['snr'],v['truth_b']) for v in r['cases'])):
 vs=[v for v in r['cases'] if (v['family'],v['snr'],v['truth_b'])==key];assert len(vs)==20
 errors=np.array([v['fit']['best']['parameters'][1]-v['truth_b'] for v in vs]);k=sum(v['accepted_95'] for v in vs);bad=sum(not v['valid'] for v in vs)
 bmean=float(errors.mean());screen=abs(bmean)<=gate['maximum_absolute_mean_bias'] and k/20>=gate['minimum_accepted_95_fraction'] and bad<=gate['maximum_invalid_per_cell']
 cells.append(dict(family=key[0],snr=key[1],truth_b=key[2],n=20,mean_bias=bmean,bias_monte_carlo_se=float(errors.std(ddof=1)/np.sqrt(20)),accepted_95=k,binomial_95=ci(k,20),invalid=bad,screen=bool(screen),original_accepted_95=sum(old[v['label']]['coverage']['0.95'] for v in vs),original_invalid=sum(not old[v['label']]['valid_fit'] for v in vs),near_zero_scatter=sum(v['fit']['best']['parameters'][2]<1e-6 for v in vs)))
comparison=[dict(label=v['label'],old_b=old[v['label']]['fit']['b'],new_b=v['fit']['best']['parameters'][1],old_valid=old[v['label']]['valid_fit'],new_valid=v['valid'],old_accepted=old[v['label']]['coverage']['0.95'],new_accepted=v['accepted_95']) for v in r['cases']]
out=dict(scope='Completed exposed-simulation development comparison; not fresh-data validation',result_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cells=cells,paired_case_comparison=comparison)
(H/'summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
lines=['# Continuous-scatter timing recalibration','','All 160 exposed artificial samples were retained. This is a development comparison, not a test on untouched observations and not evidence of photon conversion. All declared source hashes and event-array hashes were verified.','','| Shape | SNR | Truth b | Mean bias | Monte Carlo SE | Old accepted /20 | Revised accepted /20 | Revised exact 95% interval | Old/new invalid | Screen |','|---|---:|---:|---:|---:|---:|---:|---|---|---|']
for c in cells:
 lo,hi=c['binomial_95']
 lines.append(f"| {c['family']} | {c['snr']} | {c['truth_b']} | {c['mean_bias']:.5f} | {c['bias_monte_carlo_se']:.5f} | {c['original_accepted_95']} | {c['accepted_95']} | [{lo:.3f}, {hi:.3f}] | {c['original_invalid']}/{c['invalid']} | {c['screen']} |")
lines+=['','Acceptance means both valid fit and nominal likelihood-ratio interval membership at the injected exponent. Invalid fits count as unsuccessful acceptance, not necessarily a valid interval falsely excluding truth. The nominal chi-square threshold is a known asymptotic statistical approximation; this experiment checks its behavior rather than assuming it is exact. The normal width population and continuous interpolation are statistical tools, not new physical formulas.','','The screening limits are unchanged: absolute mean bias at most 0.1, accepted fraction at least 0.85, and no invalid fits in each cell. Twenty trials per cell leave wide uncertainty even if screens pass. Seeds are paired across cells; do not pool them as independent trials. Original results are preserved, not selectively replaced.','','Zero scatter is now permitted. Finite-support normalization, source priors and bounds still matter. Optimization start differences, width-grid sensitivity, independent simulation seeds, intrinsic source evolution, actual filters and calibrated brightness remain required. These results cannot establish a propagation law, absolute energy transfer or gravitational deposition. All six research goals remain open.','']
(H/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
