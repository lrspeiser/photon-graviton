"""Transfer brightness-required event exponent to spectral aging observations."""
import csv,hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
opacity=json.loads((HERE/'opacity-results.json').read_text())
timing=json.loads((HERE.parent/'conversion-first/channel-budget.json').read_text())
source=HERE.parent/'electromagnetic-audit/spectral-aging-predictions.csv'
with source.open(newline='',encoding='utf-8') as f:rows=list(csv.DictReader(f))
assert len(rows)==35 and hashlib.sha256(source.read_bytes()).hexdigest()==timing['input_sha256']
b_required=1+opacity['epsilon'];b_des=timing['calibration']['published_b']
models=[];predictions=[]
for name,b in [('baseline',1.),('DES_central',b_des),('brightness_time_only',b_required)]:
    scores={'all':0.,'low_z':0.,'high_z':0.}
    for row in rows:
        z=float(row['z']);obs=float(row['observed_aging_rate']);sigma=float(row['sigma'])
        pred=(1+z)**(-b);e=(obs-pred)/sigma
        scores['all']+=e*e;scores['low_z' if z<.04 else 'high_z']+=e*e
        predictions.append(dict(model=name,b=b,object=row['object'],z=z,observed_aging_rate=obs,sigma=sigma,predicted_aging_rate=pred,standardized_residual=e))
    models.append(dict(model=name,b=b,spectral_aging_diagonal_scores=scores,stretch_at_z1=2**b))
assert abs(models[0]['spectral_aging_diagonal_scores']['all']-timing['comparisons'][0]['diagonal_chi_square']['all'])<1e-10
geometry=[dict(z=z,required_beam_distance_over_Euclidean=(1+z)**((b_required-b_des)/2),required_beam_area_ratio=(1+z)**(b_required-b_des)) for z in [.1,.5,1.,1.2]]
out=dict(brightness_required_b=b_required,brightness_formal_sigma=opacity['epsilon_formal_sigma'],DES_b=b_des,DES_statistical_scale=timing['calibration']['statistical_scale'],DES_systematic_estimate=timing['calibration']['systematic_estimate'],exponent_difference=b_required-b_des,models=models,required_geometry_if_timing_fixed=geometry,brightness_transfer_baseline_score=opacity['transfer_baseline_chi2'],brightness_transfer_time_only_score_on_same_covariance=opacity['transfer_revised_chi2_on_baseline_covariance'],predictions=predictions)
(HERE/'time-only-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='predictions'},indent=2))
