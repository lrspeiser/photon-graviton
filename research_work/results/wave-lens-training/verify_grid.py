from pathlib import Path
import json
import numpy as np
H=Path(__file__).resolve().parent;D=H/'parameter-grid'
summary=json.loads((D/'summary.json').read_text());assert summary['grid_complete']
original=json.loads((H/'predictions.json').read_text())
original={(r['Name'],r['model']):r for r in original}
score_baseline=.5*np.mean([np.log(r['sigma_pred_km_s']/r['sigma_observed_km_s'])**2+
                         np.log(r['theta_pred_arcsec']/r['theta_SIE_arcsec'])**2
                         for r in original.values() if r['model']=='baryons'])
max_repeat=max_projection=max_norm=max_virial=0.
for cell in summary['cells']:
    path=D/f"cell-{cell['cell']}"
    rows=json.loads((path/'predictions.json').read_text())
    result=json.loads((path/'results.json').read_text())
    assert {(r['Name'],r['model']) for r in rows}==set(original)
    assert all(r['role']=='training' for r in rows)
    for row in rows:
        ref=original[row['Name'],row['model']]
        assert row['stellar_mass_Msun']==ref['stellar_mass_Msun']
        assert row['sigma_observed_km_s']==ref['sigma_observed_km_s'] and row['theta_SIE_arcsec']==ref['theta_SIE_arcsec']
        max_projection=max(max_projection,row['projected_mass_lensing_identity_relative_error'])
        if cell['cell']==5 or row['model']=='baryons':
            for k in ['sigma_pred_km_s','theta_pred_arcsec']:
                max_repeat=max(max_repeat,abs(row[k]/ref[k]-1))
    max_norm=max(max_norm,max(r['normalization_relative_error'] for r in result['equilibrium_checks']))
    max_virial=max(max_virial,max(r['virial_relative_residual'] for r in result['equilibrium_checks']))
assert max_repeat<1e-10 and max_projection<1e-6 and max_norm<1e-6 and max_virial<2e-5
result={'all_nine_cells_complete':True,'same_32_systems_and_fixed_stellar_masses':True,
    'baseline_log_error_score':float(score_baseline),
    'max_baseline_or_original_point_repeat_relative_error':max_repeat,
    'max_projected_lensing_relative_error':max_projection,'max_source_normalization_error':max_norm,
    'max_source_virial_residual':max_virial,
    'best_joint_cell':summary['best_successful_cell'],
    'best_dispersion_rms_cell':min(summary['cells'],key=lambda r:r['scores']['stationary_wave']['sigma_rmse_km_s'])['cell'],
    'best_lensing_rms_cell':min(summary['cells'],key=lambda r:r['scores']['stationary_wave']['theta_rmse_arcsec'])['cell'],
    'scope':'Complete local training grid; not global optimum, uncertainty likelihood or holdout validation'}
(D/'verification.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
