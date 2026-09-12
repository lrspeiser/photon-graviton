"""Check conventions against archived FITS reference magnitudes and derivatives."""
from pathlib import Path
import hashlib,json
import numpy as np
from astropy.io import fits
from calibration import ab_magnitude_to_native_fluxcal,native_fluxcal_to_ab_fluxcal,calibration_covariance
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
manifest=json.loads((HERE/'source-manifest.json').read_text(encoding='utf-8'))
cache=ROOT/'research_work/generated/des-photometric-calibration/snana-source'
for f in manifest['files']:assert hashlib.sha256((cache/f['path']).read_bytes()).hexdigest()==f['sha256']
a=json.loads((HERE.parent/'des-photometric-calibration/acquisition.json').read_text(encoding='utf-8'))
r=next(v for v in a['selected'] if v['archive_path'].endswith('calib_DES-SN5YR_DES.fits.gz'))
path=ROOT/r['cache_path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==r['sha256']
with fits.open(path) as hd:table=hd['ZPoff'].data.copy()
rows=[]
for row in table:
    delta=float(row['ZPoff(Primary)']);reference=float(row['Primary Mag'])
    native=ab_magnitude_to_native_fluxcal(0.,delta)
    recovered_mag=27.5-2.5*np.log10(native)
    assert abs(recovered_mag-reference)<1e-12
    assert float(row['ZPoff(SNpot)'])==0.
    rows.append({'band':str(row['Filter Name']),'offset_mag':delta,'model_native_flux_factor':float(10**(-.4*delta)),'data_AB_equivalent_flux_factor':float(10**(.4*delta)),'reference_magnitude_error':float(abs(recovered_mag-reference))})
f=np.array([10.,20.,30.]);bands=np.array([0,0,1]);cov=np.array([[.01**2,.2*.01*.02],[.2*.01*.02,.02**2]])
analytic=calibration_covariance(f,bands,cov);eps=1e-6;jac=np.zeros((3,2))
for k in range(2):
    shift=np.zeros(2);shift[k]=eps
    jac[:,k]=(native_fluxcal_to_ab_fluxcal(f,shift[bands])-native_fluxcal_to_ab_fluxcal(f,-shift[bands]))/(2*eps)
finite=jac@cov@jac.T;assert np.allclose(analytic,finite,rtol=1e-8,atol=1e-12)
assert np.isclose(analytic[0,1],2*analytic[0,0])
assert np.isfinite(native_fluxcal_to_ab_fluxcal([-1,0,1],.01)).all()
out={'scope':'Checked transforms for the historical calibration system. No observed flux edited and no measured offset covariance supplied.','rows':rows,'finite_difference_covariance_max_error':float(np.max(abs(analytic-finite))),'shared_band_covariance_verified':True,'hypothetical_covariance_used_only_for_test':True,'release_compatibility_established':False}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out))
