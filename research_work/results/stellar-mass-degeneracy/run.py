from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq,minimize_scalar
H=Path(__file__).resolve().parent;R=H.parent
protocol=json.loads((H/'protocol.json').read_text())
geometryfile=R/'lensing-data-readiness/conditional-geometry.json'
geometry={r['Name']:r for r in json.loads(geometryfile.read_text()) if r['role'] in ['training','validation']}
trainfile=R/'wave-lens-training/parameter-grid/cell-4/predictions.json'
valfile=R/'wave-lens-validation/predictions.json'
massfile=R/'lens-photometric-audit/normalization-sensitivity-updated-profile.json'
valmassfile=R/'wave-lens-validation/mass-inputs.json'
phot=[r for r in json.loads(massfile.read_text()) if r['imf']=='Salpeter' and r['propagation_branch']=='energy_loss_and_event_stretch' and r['model']=='baryons']
photos={r['Name']:r for r in phot}
photos.update({r['Name']:r for r in json.loads(valmassfile.read_text())})
G=4.30091727003628e-6;C=299792.458;RAD=np.pi/(180*3600)

def projected(b):
    # Known projected Hernquist mass, with direct quadrature near the removable singularity.
    if abs(b-1)<1e-4:
        return b*quad(lambda t:(b/np.cos(t))**2/(1+b/np.cos(t))**2*np.cos(t)/b,0,np.pi/2,epsabs=1e-11)[0]
    X=np.arccosh(1/b)/np.sqrt(1-b*b) if b<1 else np.arccos(1/b)/np.sqrt(b*b-1)
    return b*b*(1-X)/(b*b-1)

def prediction(row,lam):
    n=row['Name'];re=photos[n]['new_I_Re_arcsec'];geo=geometry[n]
    a=geo['conditional_Dl_Mpc']*1000*re*RAD/1.8153
    factor=4*G*row['stellar_mass_Msun']/(a*C*C)*geo['conditional_Dls_over_Ds']/RAD
    def bend(angle):
        b=angle/re*1.8153
        return factor*projected(b)/b
    theta=np.exp(brentq(lambda log:lam*bend(np.exp(log))/np.exp(log)-1,-24,4))
    return dict(Name=n,role=row['role'],model='shared_stellar_mass',lambda_mass=lam,
        original_stellar_mass_Msun=row['stellar_mass_Msun'],stellar_mass_Msun=lam*row['stellar_mass_Msun'],
        sigma_observed_km_s=row['sigma_observed_km_s'],sigma_pred_km_s=row['sigma_pred_km_s']*np.sqrt(lam),
        theta_SIE_arcsec=row['theta_SIE_arcsec'],theta_pred_arcsec=float(theta),
        individually_required_lambda_from_sigma=(row['sigma_observed_km_s']/row['sigma_pred_km_s'])**2,
        individually_required_lambda_from_lens=row['theta_SIE_arcsec']/bend(row['theta_SIE_arcsec']))

def score(rows):
    return float(.5*np.mean([np.log(r['sigma_pred_km_s']/r['sigma_observed_km_s'])**2+
                            np.log(r['theta_pred_arcsec']/r['theta_SIE_arcsec'])**2 for r in rows]))

def metrics(rows):
    ds=np.array([r['sigma_pred_km_s']-r['sigma_observed_km_s'] for r in rows])
    da=np.array([r['theta_pred_arcsec']-r['theta_SIE_arcsec'] for r in rows])
    return {'score':score(rows),'sigma_rmse_km_s':float(np.sqrt(np.mean(ds**2))),
            'theta_rmse_arcsec':float(np.sqrt(np.mean(da**2))),
            'median_sigma_ratio':float(np.median([r['sigma_pred_km_s']/r['sigma_observed_km_s'] for r in rows])),
            'median_theta_ratio':float(np.median([r['theta_pred_arcsec']/r['theta_SIE_arcsec'] for r in rows]))}

training=json.loads(trainfile.read_text());base=[r for r in training if r['model']=='baryons']
assert len(base)==32 and all(r['role']=='training' for r in base)
evaluations=[]
def objective(loglam):
    lam=float(np.exp(loglam));value=score([prediction(r,lam) for r in base])
    evaluations.append({'lambda':lam,'score':value});return value
fit=minimize_scalar(objective,bounds=np.log(protocol['lambda_range']),method='bounded',options={'xatol':1e-8})
assert fit.success
lam=float(np.exp(fit.x))
frozen={'lambda_mass':lam,'log10_mass_shift':float(np.log10(lam)),
        'fit_role':'training','objective_value':float(fit.fun),
        'hashes':{str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [trainfile,massfile,geometryfile,H/'protocol.json']}}
(H/'frozen-calibration.json').write_text(json.dumps(frozen,indent=2)+'\n',newline='\n')
# Only now load validation predictions; no parameter update follows.
validation=json.loads(valfile.read_text())
assert len(validation)==14 and all(r['role']=='validation' for r in validation)
allpred=[];summary={'calibration':frozen,'scores':{},'objective_evaluations':evaluations}
for role,observed in [('training',training),('validation',validation)]:
    bb=[r for r in observed if r['model']=='baryons'];ww=[r for r in observed if r['model']=='stationary_wave']
    predicted=[prediction(r,lam) for r in bb];allpred.extend(predicted)
    summary['scores'][role]={'ordinary_fixed_mass':metrics(bb),'ordinary_shared_mass':metrics(predicted),'frozen_wave':metrics(ww)}
# Check projected analytic expression independently across radii used by root solver.
errors=[]
for b in np.geomspace(.001,100,100):
    reference=quad(lambda t:(b/np.cos(t))**2/(1+b/np.cos(t))**2*np.cos(t),0,np.pi/2,epsabs=1e-12)[0]
    errors.append(abs(projected(b)/reference-1))
repeat=max(abs(prediction(r,1.)['theta_pred_arcsec']/r['theta_pred_arcsec']-1) for r in base)
assert max(errors)<1e-7 and repeat<1e-7
# A dense coarse scan verifies the bounded solver did not select a worse local minimum.
scan=[{'lambda':float(t),'score':score([prediction(r,float(t)) for r in base])} for t in np.geomspace(.5,3.,121)]
assert fit.fun<=min(r['score'] for r in scan)+1e-10
summary['checks']={'projected_mass_max_relative_error':max(errors),'original_angle_repeat_max_relative_error':repeat,
                   'dense_scan_points':len(scan),'dense_scan_best_score':min(r['score'] for r in scan)}
summary['input_sha256']={str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [trainfile,valfile,massfile,valmassfile,geometryfile,H/'protocol.json']}
for file,obj in [('results.json',summary),('predictions.json',allpred),('training-scan.json',scan)]:
    (H/file).write_text(json.dumps(obj,indent=2)+'\n',newline='\n')
print(json.dumps({'lambda':lam,'scores':summary['scores'],'checks':summary['checks']},indent=2))
