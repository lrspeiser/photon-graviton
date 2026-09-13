"""Conditional aperture ordering, not a lensing fit or cluster simulation."""
from pathlib import Path
import json
import hashlib

HERE = Path(__file__).resolve().parent
# Clowe et al. 2006 Table 2; masses in 1e12 solar masses under its calibration.
# B denotes the BCG-centered aperture, P the plasma-centered aperture.
inputs = [dict(name='main', gas_B=5.5, gas_P=6.6, star_B=.54, star_P=.23),
          dict(name='subcluster', gas_B=2.7, gas_P=5.8, star_B=.58, star_P=.12)]
rows = []
for d in inputs:
    dg, ds = d['gas_P']-d['gas_B'], d['star_B']-d['star_P']
    assert dg > 0 and ds > 0
    for stellar_mass_scale in [.25, 1., 1.5]:
        threshold = dg/(ds*stellar_mass_scale)
        # Effective mass is B_star * stellar_mass + B_gas * gas_mass.
        # Verify the crossing algebra directly, for several gas loadings.
        for bg in [1., 2., 10.]:
            def contrast(bs):
                return (bs*stellar_mass_scale*d['star_B']+bg*d['gas_B']
                        -bs*stellar_mass_scale*d['star_P']-bg*d['gas_P'])
            assert abs(contrast(bg*threshold)) < 1e-12
            assert contrast(bg*threshold*.99) < 0 < contrast(bg*threshold*1.01)
        rows.append(dict(name=d['name'], stellar_mass_scale=stellar_mass_scale,
                         required_Bstar_over_Bgas_strictly_greater_than=threshold,
                         additional_stellar_rest_energy_over_original_rest_energy_if_Bgas1=threshold-1,
                         universal_loading_preserves_plasma_larger=True))
out = dict(scope='Conditional additive receiver mass contrast; no lensing likelihood',
           source='https://arxiv.org/abs/astro-ph/0608407', source_table=2,
           reference_mass_inputs=inputs, sensitivity_rows=rows,
           uncertainty='No confidence intervals: mass covariance and calibration posterior unavailable',
           script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'results.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps(rows, indent=2))
