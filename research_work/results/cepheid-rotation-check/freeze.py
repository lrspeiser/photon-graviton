"""Freeze an external-catalog check before reading its numerical table."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
paths=[ROOT/'research_work/results/conservative-field-completion/run.py',
       ROOT/'research_work/results/conservative-field-completion/results-refined.json',
       ROOT/'research_work/results/joint-galaxy-audit/results.json',
       ROOT/'research_work/results/milky-way-mass-response/results.json',
       ROOT/'research_work/data-cache/bar-field/bar-L64.npz',
       ROOT/'research_work/data-cache/bar-field/nuclei-L16.npz']
mass_result=json.loads(paths[3].read_text())
rotation=next(row for row in mass_result['cases'] if row['calibration_objective']=='rotation')
protocol=dict(status='Frozen before numerical table acquisition; verify commit order in git',
    target_title='The rotation curve of the Milky Way measured by classical Cepheids from Gaia DR3',
    target_url='https://academic.oup.com/mnras/article/546/2/stag011/8416425',
    input_target='Published final circular-velocity Table 1, including its stated asymmetric-drift correction and uncertainties',
    prior_exposure='Search-result title, sample size 903, stated range 6<R<18 kpc, and table availability seen; numerical table not yet read in this check. Same Galaxy and Gaia astrometry as earlier work; individual-star overlap not established.',
    eligibility='All finite positive-radius, positive-speed rows in the published final Table 1 within 6<=R<=18 kpc; no residual-based exclusions. Missing uncertainties retain a row for unweighted metrics but prevent uncertainty scoring.',
    baseline='Same axisymmetrized bar, nuclear components, two stellar and two gas disks, central mass as archived conservative completion; no halo',
    prediction_resolution=dict(refined=True,response_outer_radius_kpc=200),
    models=[dict(name='ordinary',lambda_mass=1.,role='fixed ordinary-matter reference'),
            dict(name='conservative_completion',lambda_mass=1.,role='primary frozen empirical-field prediction'),
            dict(name='rotation_mass_adjusted_completion',lambda_mass=rotation['lambda_mass'],role='secondary normalization previously fitted to Eilers rotation rows; no Cepheid fitting')],
    metrics=['Unweighted circular-speed RMS in km/s','Mean signed prediction-minus-observation residual','Median predicted/observed ratio','Per-row residual with published error when available; any error-standardized summary is diagnostic only because full covariance/systematics are not established'],
    prohibitions=['No parameter refit on this table','No choosing the winning mass or geometry after opening the table','No calling shared Gaia systematics or model-derived circular speeds fully independent raw-force measurements','No treating rotation agreement as proof of photon origin or vertical/lensing agreement'],
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths})
path=HERE/'protocol.json'
if path.exists():raise RuntimeError('Do not overwrite a frozen protocol')
path.write_text(json.dumps(protocol,indent=2)+'\n',newline='\n')
print(json.dumps(protocol,indent=2))
