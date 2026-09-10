from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from astropy.cosmology import FlatLambdaCDM
H=Path(__file__).resolve().parent
rows=json.loads((H/'normalization-sensitivity.json').read_text())
old=json.loads((H.parent/'lensing-data-readiness/lens-observations-and-image-models.json').read_text());old={r['Name']:r for r in old}
assert all(old[r['Name']]['role']=='training' for r in rows)
assert len(set(r['Name'] for r in rows))==32 and len(rows)==256
cosmo=FlatLambdaCDM(H0=70,Om0=.3,Tcmb0=0)
err=[]
for n in set(r['Name'] for r in rows):
    z=old[n]['zFG']
    independent=(1+z)*299792.458/70*quad(lambda zz:1/np.sqrt(.3*(1+zz)**3+.7),0,z,epsabs=1e-12)[0]
    err.append(abs(independent/cosmo.luminosity_distance(z).value-1))
assert max(err)<1e-11
for r in rows:
    assert abs(r['conditional_log10_stellar_mass']-r['published_log10_stellar_mass']-np.log10(r['normalization_factor']))<1e-12
single={r['Name']:r for r in rows}
sizes=np.array([r['new_I_Re_arcsec']/r['old_Re_arcsec'] for r in single.values() if r['new_I_Re_arcsec'] is not None])
out={'matched_mass_systems':32,'training_only':True,'normalization_scenarios':len(rows),'reference_distance_max_relative_difference':max(err),
    'new_to_old_I_size_ratio_median':float(np.median(sizes)),'new_to_old_I_size_ratio_range':[float(sizes.min()),float(sizes.max())],
    'note':'The size changes are not applied to the existing dynamical pilot. A consistent updated profile calculation is required before a final mass comparison.',
    'full_stellar_population_refit':False,'holdout_scores_opened':False}
(H/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
