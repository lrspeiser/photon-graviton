"""Historical DES synthetic-system transforms; observed data are not modified."""
import numpy as np

def ab_magnitude_to_native_fluxcal(magnitude_ab,offset_mag):
    return 10.**(.4*(27.5-np.asarray(magnitude_ab)-np.asarray(offset_mag)))

def native_fluxcal_to_ab_fluxcal(flux,offset_mag):
    # Linear conversion also supports negative noisy flux measurements.
    return np.asarray(flux)*10.**(.4*np.asarray(offset_mag))

def calibration_covariance(model_ab_fluxcal,band_index,offset_covariance):
    """First-order covariance in AB-equivalent FLUXCAL from shared offset errors.

    Evaluate the Jacobian at predicted fluxes, not at noisy measurements.
    """
    f=np.asarray(model_ab_fluxcal,dtype=float);index=np.asarray(band_index)
    covariance=np.asarray(offset_covariance,dtype=float)
    if f.ndim!=1 or index.shape!=f.shape or index.dtype.kind not in 'iu':raise ValueError('Invalid flux/band dimensions')
    if covariance.ndim!=2 or covariance.shape[0]!=covariance.shape[1]:raise ValueError('Invalid covariance dimensions')
    if np.any(index<0) or np.any(index>=len(covariance)):raise ValueError('Invalid band index')
    if not np.isfinite(f).all() or not np.isfinite(covariance).all():raise ValueError('Nonfinite input')
    if not np.allclose(covariance,covariance.T) or np.linalg.eigvalsh(covariance).min() < -1e-12:raise ValueError('Covariance must be symmetric positive semidefinite')
    jac=np.zeros((len(f),len(covariance)));jac[np.arange(len(f)),index]=.4*np.log(10)*f
    return jac@covariance@jac.T
