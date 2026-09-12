"""Independent AB count-rate integration against archived calibration log."""
from pathlib import Path
import hashlib,json,re
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'research_work/generated/des-photometric-calibration/selected'
for filename in ['acquisition.json','dependencies.json']:
    for v in json.loads((HERE/filename).read_text(encoding='utf-8'))['selected']:
        assert hashlib.sha256((ROOT/v['cache_path']).read_bytes()).hexdigest()==v['sha256']
cal=(BASE/'kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.input').read_text()
log=(BASE/'kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.log').read_text()
rows=[]
for band in 'griz':
    wavelength,transmission=np.loadtxt(BASE/f'filters/DES/DES-SN3YR_DECam/DECam_{band}.dat',unpack=True)
    assert np.all(np.diff(wavelength)>0) and np.isfinite(transmission).all() and np.all(transmission>=0)
    h=6.62607015e-27;fnu=10**(-.4*48.6)
    counts=float(fnu/h*np.trapezoid(transmission/wavelength,wavelength))
    published=float(re.search(r'DES-'+band+r'\s+\(\s*AB\)\s+([0-9.e+-]+)',log).group(1))
    difference=abs(counts/published-1)
    expression=re.search(r'FILTER:\s+DES-'+band+r'\s+\S+\s+([^\s]+)',cal).group(1)
    offset=sum(float(v) for v in re.findall(r'[+-]?\d+(?:\.\d+)?',expression))
    pivot=float(np.sqrt(np.trapezoid(transmission*wavelength,wavelength)/np.trapezoid(transmission/wavelength,wavelength)))
    rows.append({'band':band,'points':len(wavelength),'pivot_angstrom':pivot,'calibration_offset_expression':expression,'calibration_offset_mag':offset,'calculated_AB0_photons_per_s_cm2':counts,'archived_log_AB0_photons_per_s_cm2':published,'relative_difference':difference,'agreement_within_0_1_percent':bool(difference<.001)})
result={'scope':'Independent count-weighted filter integration for a known flat-fnu AB reference. No SN flux fit, applied correction or cosmological/source model.', 'rows':rows,'all_agree_within_0_1_percent':all(v['agreement_within_0_1_percent'] for v in rows),'offset_sign_applied_to_observed_data':False}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result))
