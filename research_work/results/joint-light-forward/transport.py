"""Conditional static-universe transport to historical DES native FLUXCAL."""
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
H_ERG_S=6.62607015e-27
C_ANGSTROM_S=2.99792458e18
MPC_CM=3.0856775814913673e24
AB0_FNU=10.**(-.4*48.6)
DEFAULT_ALPHA_PER_MPC=0.0002488993286382367

def load_filters():
    manifest=json.loads((ROOT/'research_work/results/des-photometric-calibration/dependencies.json').read_text(encoding='utf-8'))
    calibration=json.loads((ROOT/'research_work/results/des-calibration-convention/results.json').read_text(encoding='utf-8'))
    offsets={r['band'][-1]:r['offset_mag'] for r in calibration['rows']}
    result={}
    for r in manifest['selected']:
        if '/filters/' not in '/'+r['archive_path']:continue
        path=ROOT/r['cache_path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==r['sha256']
        band=path.stem[-1];lam,transmission=np.loadtxt(path,unpack=True)
        reference=AB0_FNU/H_ERG_S*np.trapezoid(transmission/lam,lam)
        result[band]=(lam,transmission,float(reference),offsets[band])
    if set(result)!=set('griz'):raise ValueError('Need all four calibrated filters')
    return result

def predict(observer_days,bands,distance_mpc,emitted_luminosity_lambda,*,alpha_per_mpc=DEFAULT_ALPHA_PER_MPC,stretch_exponent=1.,filters=None):
    """SED callable accepts emitted days and angstrom, returns erg/s/angstrom.

    Input time is observer elapsed time relative to the arrival of source time zero.
    Distance is independently supplied, not inferred from observed redshift/flux.
    No lensing, dust, Doppler corrections or receiver dynamics are implemented.
    """
    t=np.asarray(observer_days,dtype=float);band=np.asarray(bands,dtype=str)
    if t.ndim!=1 or band.shape!=t.shape or not np.isfinite(t).all():raise ValueError('Invalid observation arrays')
    if not np.isfinite([distance_mpc,alpha_per_mpc,stretch_exponent]).all() or distance_mpc<=0 or alpha_per_mpc<0:raise ValueError('Invalid distance or transfer parameters')
    filters=load_filters() if filters is None else filters
    if not set(band).issubset(filters):raise ValueError('Unsupported band label')
    optical_depth=alpha_per_mpc*distance_mpc
    S=float(np.exp(optical_depth));stretch=float(S**stretch_exponent)
    if not np.isfinite([S,stretch]).all() or stretch<=0:raise ValueError('Nonfinite transport')
    area=4*np.pi*(distance_mpc*MPC_CM)**2
    counts=np.empty_like(t);native=np.empty_like(t);ab_equivalent=np.empty_like(t)
    for name in np.unique(band):
        keep=band==name;lam,response,reference,offset=filters[name]
        source=np.asarray(emitted_luminosity_lambda(t[keep,None]/stretch,lam[None,:]/S),dtype=float)
        source=np.broadcast_to(source,(int(keep.sum()),len(lam)))
        if not np.isfinite(source).all() or np.any(source<0):raise ValueError('SED must be finite and nonnegative')
        flux_lambda=source/(area*S**(stretch_exponent+2))
        counts[keep]=np.trapezoid(flux_lambda*response[None,:]*lam[None,:]/(H_ERG_S*C_ANGSTROM_S),lam,axis=1)
        ab_equivalent[keep]=10**11*counts[keep]/reference
        native[keep]=ab_equivalent[keep]*10**(-.4*offset)
    return {'predicted_redshift':float(np.expm1(optical_depth)),'spectral_stretch':S,'event_stretch':stretch,'photon_rate_per_cm2_s':counts,'ab_equivalent_fluxcal':ab_equivalent,'native_fluxcal':native}
