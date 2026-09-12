"""Read published zero-point covariance without pickle, preserve full matrix."""
from pathlib import Path
import hashlib,json,sys,urllib.request
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CACHE=ROOT/'research_work/generated/des-calibration-covariance'
CACHE.mkdir(parents=True,exist_ok=True)
m=json.loads((HERE/'source-manifest.json').read_text(encoding='utf-8'))
for row in m['records']:
    path=CACHE/row['name']
    if not path.exists():path.write_bytes(urllib.request.urlopen(row['url'],timeout=60).read())
    data=path.read_bytes();assert len(data)==row['bytes']
    assert hashlib.sha256(data).hexdigest()==row['sha256']
    assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['git_blob_sha1']
with np.load(CACHE/'FRAGILISTIC_COVARIANCE.npz',allow_pickle=False) as a:cov=a['cov'];labels=np.char.strip(a['labels'])
assert cov.shape==(len(labels),len(labels)) and len(set(labels))==len(labels)
assert np.isfinite(cov).all() and np.allclose(cov,cov.T,atol=1e-15,rtol=0)
eigen=np.linalg.eigvalsh(cov);assert eigen.min()>0
selected=['DES5YR '+b for b in 'griz'];index=[int(np.flatnonzero(labels==name)[0]) for name in selected]
sub=cov[np.ix_(index,index)];sigma=np.sqrt(np.diag(sub));corr=sub/np.outer(sigma,sigma)
sys.path.insert(0,str(HERE.parent/'des-calibration-convention'))
from calibration import calibration_covariance
# Illustrative predictions, not observed supernova flux.
flux=np.array([10.,20.,30.,40.,50.]);bands=np.array([0,0,1,2,3]);propagated=calibration_covariance(flux,bands,sub)
assert np.isclose(propagated[0,1],2*propagated[0,0])
assert np.linalg.eigvalsh(propagated).min()>-1e-12
out={'scope':'Published DES5YR zero-point covariance subblock and first-order propagation. No observed SN fit or complete calibration budget.', 'full_matrix_shape':list(cov.shape),'full_min_eigenvalue':float(eigen.min()),'maximum_symmetry_error':float(abs(cov-cov.T).max()),'selected_labels':selected,'selected_zero_based_indices':index,'covariance_mag_squared':sub.tolist(),'correlation':corr.tolist(),'sigma_mag':sigma.tolist(),'first_order_fractional_flux_sigma':(.4*np.log(10)*sigma).tolist(),'illustrative_flux_covariance':propagated.tolist(),'full_matrix_retained_in_cache':True,'cross_survey_covariances_must_be_retained_for_joint_sample':True,'all_calibration_systematics_included':False}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
np.savez_compressed(HERE/'des5yr-zero-point-covariance.npz',cov=sub,labels=np.array(selected))
print(json.dumps({k:out[k] for k in ['full_matrix_shape','sigma_mag','first_order_fractional_flux_sigma','full_min_eigenvalue']}))
