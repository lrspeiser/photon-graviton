"""Angular benchmark only: no assumed distance, capture fit, or merger simulation."""
from pathlib import Path
import hashlib
import json
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = Path(__file__).resolve().parent
p = json.loads((HERE / 'inputs.json').read_text(encoding='utf-8'))
def coord(prefix):
    return SkyCoord(p[prefix+'_ra_hms'], p[prefix+'_dec_dms'], unit=(u.hourangle, u.deg))
bcg, gas = coord('bcg'), coord('plasma_aperture')
east, north = bcg.spherical_offsets_to(gas)
g = np.array([east.arcsec, north.arcsec])
l = np.array([p['lensing_peak_east_of_bcg_arcsec'], p['lensing_peak_north_of_bcg_arcsec']])
lens = bcg.spherical_offsets_by(l[0]*u.arcsec, l[1]*u.arcsec)
# Independent unit-vector spherical separation, stable at these angular distances.
def vector_sep(a, b):
    av, bv = a.cartesian.xyz.value, b.cartesian.xyz.value
    return np.arctan2(np.linalg.norm(np.cross(av, bv)), np.dot(av, bv))*180/np.pi*3600
separations = {}
for name, a, b in [('bcg_to_plasma', bcg, gas), ('bcg_to_lensing', bcg, lens), ('plasma_to_lensing', gas, lens)]:
    value = a.separation(b).arcsec
    assert abs(value-vector_sep(a, b)) < 1e-7
    separations[name+'_arcsec'] = float(value)
# Pure geometric projection, NOT a captured gas fraction or a parameter fit.
projection = float(np.dot(l, g)/np.dot(g, g))
residual = l-projection*g
out = dict(scope=p['scope'], gas_offset_east_north_arcsec=g.tolist(),
           lens_offset_east_north_arcsec=l.tolist(), separations=separations,
           lens_projection_on_bcg_gas_segment=projection,
           transverse_distance_from_segment_arcsec=float(np.linalg.norm(residual)),
           coordinate_crosscheck_tolerance_arcsec=1e-7,
           caveat='Peak positions are not reservoir centroids. No model likelihood or collision outcome is computed.',
           hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [HERE/'inputs.json', HERE/'run.py']})
(HERE/'results.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps(out, indent=2))
