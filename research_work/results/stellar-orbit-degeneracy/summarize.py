from pathlib import Path
import json,hashlib
import numpy as np
H=Path(__file__).resolve().parent;R=H.parent
files=[R/'wave-lens-training/parameter-grid/cell-4/orbital-predictions.json',
       R/'wave-lens-validation/orbital-predictions.json',R/'stellar-mass-degeneracy/predictions.json',
       R/'stellar-mass-degeneracy/frozen-calibration.json']
train,val,benchmark,calibration=[json.loads(f.read_text()) for f in files]
bench={(r['Name'],r['role']):r for r in benchmark};lam=calibration['lambda_mass']
allrows=[];scores={}
for role,rows in [('training',train),('validation',val)]:
    assert len(rows)==(32 if role=='training' else 7)*6 and all(r['role']==role for r in rows)
    for beta in [-.3,0.,.3]:
        for model in ['shared_stellar_mass','stationary_wave']:
            origin='baryons' if model=='shared_stellar_mass' else model
            sub=[]
            for r in rows:
                if r['model']!=origin or r['beta']!=beta:continue
                out=dict(r);out['model']=model
                if model=='shared_stellar_mass':
                    out['sigma_pred_km_s']*=np.sqrt(lam)
                    out['theta_pred_arcsec']=bench[r['Name'],role]['theta_pred_arcsec']
                    out['stellar_mass_Msun']*=lam
                sub.append(out)
            allrows.extend(sub)
            ds=np.array([r['sigma_pred_km_s']-r['sigma_observed_km_s'] for r in sub])
            da=np.array([r['theta_pred_arcsec']-r['theta_SIE_arcsec'] for r in sub])
            scores[f'{role}/{model}/beta{beta:+.1f}']={
                'n':len(sub),'sigma_rmse_km_s':float(np.sqrt(np.mean(ds**2))),
                'theta_rmse_arcsec':float(np.sqrt(np.mean(da**2))),
                'score':float(.5*np.mean([np.log(r['sigma_pred_km_s']/r['sigma_observed_km_s'])**2+
                                         np.log(r['theta_pred_arcsec']/r['theta_SIE_arcsec'])**2 for r in sub]))}
for role in ['training','validation']:
    for model in ['shared_stellar_mass','stationary_wave']:
        angles={beta:{r['Name']:r['theta_pred_arcsec'] for r in allrows if r['role']==role and r['model']==model and r['beta']==beta} for beta in [-.3,0.,.3]}
        assert angles[-.3]==angles[0.]==angles[.3]
result={'classification':'Shared orbital sensitivity, no new fitted beta or holdout selection',
        'fixed_stellar_mass_multiplier':lam,'scores':scores,
        'input_sha256':{str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}
max_pressure=max(r['independent_pressure_variance_relative_error'] for r in train+val)
assert max_pressure<1e-4
zero_error=0.
for role,filename in [('training',R/'wave-lens-training/parameter-grid/cell-4/predictions.json'),('validation',R/'wave-lens-validation/predictions.json')]:
    reference={(r['Name'],r['model']):r for r in json.loads(filename.read_text())}
    for row in (train if role=='training' else val):
        if row['beta']==0:
            zero_error=max(zero_error,abs(row['sigma_pred_km_s']/reference[row['Name'],row['model']]['sigma_pred_km_s']-1))
assert zero_error<1e-10
result['checks']={'max_independent_pressure_variance_relative_error':max_pressure,
                  'max_original_isotropic_dispersion_relative_change':zero_error,
                  'lensing_identical_across_beta':True}
(H/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
(H/'predictions.json').write_text(json.dumps(allrows,indent=2)+'\n',newline='\n')
print(json.dumps(scores,indent=2))
