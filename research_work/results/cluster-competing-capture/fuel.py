"""Age-independent hydrogen-burning accounting from measured atomic masses."""
from pathlib import Path
import json
import numpy as np
OUT=Path(__file__).resolve().parent
mH=1.00782503223;mHe=4.00260325413
eps=(4*mH-mHe)/(4*mH)
assert 0<eps<.01
rows=[]
for fraction in [1.,.5,.1,.01]:
    burned=1/(eps*fraction)
    helium=burned*(1-eps)
    released=eps*burned
    assert abs(burned-helium-released)<1e-10*burned
    assert abs(released*fraction-1)<1e-14
    rows.append(dict(fraction_of_nuclear_release_retained_in_target=fraction,
        minimum_hydrogen_burned_per_unit_deposit_mass_equivalent=burned,
        helium_product_per_unit_deposit_mass_equivalent=helium,
        energy_not_in_target_as_mass_equivalent=released-1))
history=json.loads((OUT/'history-results.json').read_text())
costs=[dict(old_age_alpha_c_T=v['old_burst_age_alpha_c_T'],
    reference_minimum_burned_mass_in_energy_unit_over_c2=v['baseline_total_emitted']/eps,
    alternative_minimum_burned_mass_in_energy_unit_over_c2=v['alternative_total_emitted']/eps) for v in history['cases']]
result=dict(scope='Conditional hydrogen-to-helium photon funding; optimistic bound assuming all nuclear release can become photons. Not a total stellar/AGN energy limit or measured cosmic fuel inventory.',
    atomic_masses_u=dict(hydrogen1=mH,helium4=mHe),
    sources=['https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=H','https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=He'],
    energy_release_fraction=eps,rows=rows,previous_history_fuel_costs=costs)
(OUT/'fuel-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))
