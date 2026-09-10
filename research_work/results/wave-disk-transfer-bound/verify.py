from pathlib import Path
import json
import numpy as np
from astropy.constants import G
import astropy.units as u
H=Path(__file__).resolve().parent
rows=json.loads((H/'galaxy-bounds.json').read_text())
summary=json.loads((H/'results.json').read_text())
grav=G.to(u.kpc*(u.km/u.s)**2/u.Msun).value
assert abs(grav/4.30091727003628e-6-1)<1e-12
assert len(rows)==298 and len({r['galaxy'] for r in rows})==149
max_relation=0.
for r in rows:
    minimum=max(r['last_observed_speed_km_s']**2-r['last_baryonic_speed_km_s']**2,0)*r['last_radius_kpc']/grav/r['stellar_mass_upper_Msun']
    max_relation=max(max_relation,abs(minimum-r['minimum_f_from_last_radius'])/max(1,minimum))
    assert (r['last_maximum_speed_with_source_km_s']<r['last_observed_speed_km_s'])==r['fails_last_radius_necessary_bound']
assert max_relation<1e-12
verification={'astropy_G_kpc_kms2_per_Msun':grav,'max_speed_to_required_mass_identity_relative_error':max_relation,
              'bound_failures_agree_with_direct_maximum_speed_comparison':True,
              'all_149_galaxies_retained':True,'scope':'Unit and algebra checks, not full observational uncertainty'}
(H/'verification.json').write_text(json.dumps(verification,indent=2)+'\n',newline='\n')
print(json.dumps(verification,indent=2))
