"""Check exact-third lens transfer and summarize paired galaxy errors, without refits."""
import json,hashlib
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent
def read(name):return json.loads((P/name).read_text())
new=read('third-retention-optics-results.json');old=read('retention-optics-results.json')
assert new['capture_input_sha256']==hashlib.sha256((P/'third-radiation-retention-results.json').read_bytes()).hexdigest()
assert read('third-radiation-retention-results.json')['models']['attenuated']['q']==1/3
assert new['photometric_input_sha256']==old['photometric_input_sha256']
assert new['optical_input_sha256']==old['optical_input_sha256']
previous={(r['Name'],r['model']):r for r in old['rows']}
for r in new['rows']:
    prior=previous[r['Name'],r['model']]
    assert r['observed_stellar_vrms']==prior['observed_stellar_vrms'] and r['geometry']==prior['geometry']
    assert 0<r['retention_mapping']['eta']<1
    assert r['optimizer_successes']>0
data=read('model-comparison-predictions.json')
pairs={}
for r in data:
    if r['split']=='train' or r['model'] not in ['companion_third','MOND_simple_fitted','baryons']:continue
    item=pairs.setdefault(r['galaxy'],dict(galaxy=r['galaxy'],split=r['split']))
    obs=np.array(r['observed_kms']);pred=np.array(r['predicted_kms'])
    item[r['model']]=dict(rms_kms=float(np.sqrt(np.mean((pred-obs)**2))),log_rms=float(np.sqrt(np.mean(np.log10(pred/obs)**2))))
summary=[]
for split in ['validation','test']:
    group=[r for r in pairs.values() if r['split']==split]
    for metric in ['rms_kms','log_rms']:
        summary.append(dict(split=split,metric=metric,n=len(group),
            companion_better_than_mond=sum(r['companion_third'][metric]<r['MOND_simple_fitted'][metric] for r in group),
            companion_better_than_baryons=sum(r['companion_third'][metric]<r['baryons'][metric] for r in group)))
out=dict(lens_summary=new['summary'],paired_galaxy_summary=summary,paired_galaxies=list(pairs.values()),
         scope='Descriptive, exposed samples; no new fit or significance claim',
         input_sha256={n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in ['third-retention-optics-results.json','retention-optics-results.json','model-comparison-predictions.json']})
(P/'cross-test-audit-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
lines=['# Exact one-third lens transfer and paired galaxy assessment','',
       'The exact-one-third law has now been run on the six existing lens systems, with all galaxy capture constants and optical parameters frozen. Stellar mass and orbital anisotropy are fitted only to inner stellar measurements, as before. This exposed-sample test uses conditional luminosity and morphology proxies.','',
       '| Population proxy | Outer motion residual-square sum | Lens-angle fractional RMS | Orbit boundary count |','|---|---:|---:|---:|']
for r in new['summary']:lines.append(f"| {r['model']} | {r['outer_conditional_residual_square_sum']:.5f} | {100*r['lens_fractional_rms']:.5f}% | {r['orbit_boundary_count']} |")
lines += ['', 'Exact one-third does not repair the existing lens mismatch. Its 13.923/14.094% errors remain close to the fitted-exponent values 13.938/14.110% and above the earlier intercepted-source value 13.112%. Better outer-motion errors than the older source model do not establish joint agreement.','',
          '## Paired galaxy results','', '| Split | Error measure | Galaxies | Companion closer than MOND | Companion closer than baryons |','|---|---|---:|---:|---:|']
for r in summary:lines.append(f"| {r['split']} | {r['metric']} | {r['n']} | {r['companion_better_than_mond']} | {r['companion_better_than_baryons']} |")
lines += ['', 'These counts describe how often the candidate helps, complementing the previously reported aggregate RMS. They are not discovery probabilities, and a win count does not weight the size of a miss. Do not select only winning galaxies or retune on comparison targets. Per-galaxy errors are retained in cross-test-audit-results.json.','',
          'Verification: exact q=1/3; matching capture-file hash; unchanged optical and photometric hashes; identical observed stellar velocities and geometry; probabilities inside (0,1); at least one successful optimizer start for every system/population.','',
          'See third-retention-lensing-protocol.md for the pre-run specification. No cluster-transfer or mechanism goal is completed.']
(P/'cross-test-audit-report.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(lens_summary=new['summary'],paired_summary=summary),indent=2))
