"""Freeze all prediction inputs before evaluating reserved validation motions."""
from pathlib import Path
import ast
import hashlib
import json

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
target=HERE/'protocol.json'
assert not target.exists(),'Refuse to replace a frozen protocol.'
assert not (HERE/'results.json').exists()
source=HERE.parent/'baryon-component-response/run.py'
adapter=HERE/'components.py'
def classes(path):return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(path.read_text()).body if isinstance(n,ast.ClassDef)}
assert classes(source)==classes(adapter)
mass=json.loads((HERE.parent/'baryon-component-response/results.json').read_text())
paths=[HERE/'evaluate.py',adapter,HERE.parent/'cepheid-common-frame/run.py',HERE.parent/'cepheid-common-frame/protocol.json',
    source,HERE.parent/'baryon-component-response/results.json',HERE.parent/'conservative-field-completion/run.py',
    HERE.parent/'joint-galaxy-audit/results.json',ROOT/'research_work/data-cache/bar-field/bar-L64.npz',ROOT/'research_work/data-cache/bar-field/nuclei-L16.npz',
    ROOT/'research_work/data-cache/cepheid-stars/gaia-dr3-dcep-parent-with-flags.parquet',ROOT/'research_work/data-cache/cepheid-stars/cepheid-common-frame-split.parquet']
protocol=dict(stage='Freeze before first inferred validation-star outcomes',
    selection='Use existing validation role and measurement/mode flag; apply identical distance calibration, common frame, R=6..18 kpc, abs(phi)<=30 degrees, abs(z)<=0.5 kpc, abs(vz)<=100 km/s.',
    bins='Twelve fixed one-kpc bins. Score only n>=5 and valid corrected moments, per frozen training routine. Record every unscored bin. No residual exclusions or merging.',
    estimator='Same frozen simplified Jeans moments and independent 7-percent distance-error scenario as training.',
    primary='Original frozen conservative completion; all ordinary component multipliers unity.',
    secondary='Previously fitted balanced ordinary-component scales, no q deformation.',
    balanced_component_scales=mass['balanced_candidate_scales'],
    controls=['Original ordinary matter','Balanced ordinary matter'],
    metrics='Unweighted bin RMS, bias, median predicted/proxy ratio, count below proxy; no significance threshold or complete physical-validation claim.',
    numerical_check='Refined vs coarse per-bin prediction difference below 0.1 km/s.',
    final_test='Do not compute or inspect final-test inferred positions, velocities or model residuals.',
    component_adapter_classes_identical=True,
    frozen_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    limitations=['Same Gaia calibration and Galactic population as training; published aggregate curves already exposed.',
        'Conditional validation of a frozen approximate pipeline, not a full selection/orbit likelihood or a new independent catalog.',
        'Published calibration and adopted radial profiles are shared assumptions.',
        'Mass adjustment is exploratory and lacks independent allowed-mass priors.',
        'No photon production, event timing, transport, storage, lensing or source-energy validation.',
        'Validation becomes exposed once evaluated; later retuning is not fresh validation.'])
target.write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
print('Protocol frozen; validation outcomes not processed.')
