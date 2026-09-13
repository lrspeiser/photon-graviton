"""Verify frozen transfers and publish every redistribution outcome."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import quad
P=Path(__file__).resolve().parent
def read(n):return json.loads((P/n).read_text())
models=read('redistribution-results.json')['models'];models.update({k:v for k,v in read('redistribution-compact-results.json')['models'].items() if k!='baseline'})
baseline=read('third-retention-optics-results.json');prior={(r['Name'],r['model']):r for r in baseline['rows']}
lines=['# Conservative redistribution of the exact-one-third model','',
       'All variants preserve the original one-third retention, luminosity proxy and capture constants. Only shared redistribution parameters were fitted on 89 training galaxies. Predictions use the same 29 validation/31 test galaxies and six lens systems; all sets are previously exposed. This completes the bounded family and one boundary follow-up in redistribution-protocol.md.','',
       '| Variant | Validation RMS km/s | Test RMS km/s | Validation log RMS | Test log RMS | Chabrier lens RMS | Salpeter lens RMS |','|---|---:|---:|---:|---:|---:|---:|']
results={}
for name,m in models.items():
    lens=baseline if name=='baseline' else read('redistribution-'+name+'-optics-results.json')
    assert lens['capture_input_sha256']==baseline['capture_input_sha256']
    assert lens['optical_input_sha256']==baseline['optical_input_sha256']
    assert lens['photometric_input_sha256']==baseline['photometric_input_sha256']
    if name!='baseline':
        file='redistribution-compact-results.json' if name=='compact_partial' else 'redistribution-results.json'
        assert lens['redistribution_input_sha256']==hashlib.sha256((P/file).read_bytes()).hexdigest()
    for r in lens['rows']:
        old=prior[r['Name'],r['model']]
        assert r['geometry']==old['geometry'] and r['observed_stellar_vrms']==old['observed_stellar_vrms']
        assert r['retention_mapping']==old['retention_mapping'] and r['optimizer_successes']>0
    a,b=m['scores']['validation'],m['scores']['test'];l=lens['summary']
    lines.append(f"| {name} | {a['RMSE_kms']:.4f} | {b['RMSE_kms']:.4f} | {a['log_RMS']:.6f} | {b['log_RMS']:.6f} | {100*l[0]['lens_fractional_rms']:.5f}% | {100*l[1]['lens_fractional_rms']:.5f}% |")
    results[name]=dict(galaxy=m,lens_summary=l)
# Independent finite-volume conservation identity for the same radial scaling/mixing.
# Positive smooth profile, integrate density directly rather than its interpolated CDF.
checks=[]
for s,f in [(np.exp(models['shared']['parameters'][0]),1),
            (np.exp(models['partial']['parameters'][1]),models['partial']['parameters'][0]),
            (np.exp(models['compact_partial']['parameters'][1]),models['compact_partial']['parameters'][0])]:
    rho=lambda r:(1+r*r)**-2
    mass=lambda R:quad(lambda r:4*np.pi*r*r*rho(r),0,R,epsabs=1e-10)[0]
    for R in [.1,1,10]:
        direct=quad(lambda r:4*np.pi*r*r*((1-f)*rho(r)+f*rho(r/s)/s**3),0,R,epsabs=1e-10)[0]
        expected=(1-f)*mass(R)+f*mass(R/s);checks.append(abs(direct/expected-1))
assert max(checks)<1e-8
lines += ['', '## Outcome and iteration','',
          'The shared scale remains almost unchanged (s=0.99984), and the retention-conditioned scale also remains close to one. Neither gives a useful improvement. Partial redistribution initially selects f=0.01572 and s=0.25 at its allowed boundary: about 1.6% moves inward, not outward. It improves galaxy errors but worsens lens angles. An outward local optimum exists, but was not the best multi-start fit; it is not the selected result.','',
          'A separately declared extension to s>=1/16 finds an interior optimum: f=0.00369849 and s=0.10217473 (approximately). Thus about 0.37% of the existing deposit inventory is concentrated to about one-tenth its original radius. All six optimizer starts converge to essentially the same optimum, without a bound hit. This is an exploratory post-boundary extension and adds two shared parameters; it is not a newly derived law or proof of a central particle population.','',
          'The compact follow-up lowers validation/test speed RMS from 32.495/23.591 to 31.666/22.591 km/s. It remains worse than the matched MOND benchmark 26.876/16.398 km/s. Its lens errors remain about 14%; changes at the hundredth-of-a-percentage-point scale do not resolve the lens discrepancy or establish statistical improvement. The inner-star nuisance fits can compensate for redistribution, so lens response need not track a fixed-stellar-mass intuition.','',
          'Keep the compact variant as a modest rotation improvement candidate, not a replacement established across all tests. Do not continue adding radial parameters to these exposed objects. The bounded redistribution exercise has found no substantial joint rotation/lensing solution. Coma free-profile fits cannot distinguish a global scale from their already fitted scale, and do not validate this extension.','',
          '## Formulas and provenance','',
          'eta=X^(1/3)/(1+X^(1/3)) is unchanged. For each baseline deposited density rho0, rho_new(r)=(1-f)*rho0(r)+f*rho0(r/s)/s^3. Consequently M_new(<r)=(1-f)*M0(<r)+f*M0(<r/s). These are known conservative scaling/mixture identities; interpreting them as companion migration is our proposed hypothesis. They preserve total deposit inventory, not automatically total field-plus-kinetic energy. Work, release channels and support still require a physical model. No ages, external supply amplitudes or per-galaxy capture constants were changed.','',
          '## Verification','',
          'Input hashes, q=1/3, unchanged photometric/optical inputs, observed stellar data, geometry and retention mappings were checked. Every lens fit has a successful optimizer start. Galaxy radial/angular refinement changes aggregate RMS by less than 0.001 km/s; the unchanged baseline reproduces the prior result within 0.05 km/s. Direct finite-volume quadrature independently checks the conservative mixture identity to relative error '+f'{max(checks):.3g}'+'. These numerical checks do not establish a physical mechanism or significance of tiny improvements.','',
          'Commands: redistribution.py; redistribution.py --compact-followup; lensing.py --redistribution=shared (and partial, retention_conditioned, compact_partial); redistribution-report.py. Every candidate output is preserved separately. All six research goals remain open.']
(P/'redistribution-report.md').write_text('\n'.join(lines)+'\n')
(P/'redistribution-comparison.json').write_text(json.dumps(dict(models=results,conservation_identity_max_relative_error=max(checks)),indent=2)+'\n')
print('\n'.join(lines[:12]));print('compact scale',np.exp(models['compact_partial']['parameters'][1]))
