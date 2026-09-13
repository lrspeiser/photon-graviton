"""Audit absolute-capacity correction against the original fit and prior outputs."""
from pathlib import Path
import hashlib,json,subprocess
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
fitfile=P.parent/'isotropic-galaxy-transfer/third-radiation-retention-results.json'
fit=json.loads(fitfile.read_text())['models']['attenuated']
assert abs(fit['source_C_before_retention_Msun_kpc3']/(2*fit['C_Msun_kpc3'])-1)<1e-12
old_revision='9abdd90'
def old(name):
    path=str((P/name).relative_to(ROOT)).replace('\\','/')
    return json.loads(subprocess.check_output(['git','show',old_revision+':'+path],cwd=ROOT))
cap=json.loads((P/'capture-capacity-results.json').read_text());oldcap=old('capture-capacity-results.json')
recycle=json.loads((P/'conservative-recycling-results.json').read_text());oldrecycle=old('conservative-recycling-results.json')
for new,prior in zip(cap['rows'],oldcap['rows']):
    assert new['galaxy']==prior['galaxy']
    for key in ['reference_capacity_mass_Msun','unattenuated_capacity_mass_Msun']:
        assert abs(new[key]/prior[key]-2)<1e-12
    for key in ['area_per_reference_capacity_kpc2_per_Msun','area_per_unattenuated_capacity_kpc2_per_Msun']:
        assert abs(new[key]/prior[key]-.5)<1e-12
    assert new['capture_area_kpc2']==prior['capture_area_kpc2']
    assert new['attenuated_to_unattenuated_capacity']==prior['attenuated_to_unattenuated_capacity']
for new,prior in zip(recycle['rows'],oldrecycle['rows']):
    assert new['galaxy']==prior['galaxy']
    assert abs(new['total_deposit_mass_Msun']/prior['total_deposit_mass_Msun']-2)<1e-12
    assert new['predicted_kms']==prior['predicted_kms'] and new['extra_mass_ratio']==prior['extra_mass_ratio']
assert recycle['scores']==oldrecycle['scores']
assert len(cap['rows'])==len(recycle['rows'])==149
out=dict(scope='Correction of fit-parameter versus density-amplitude convention; no physics change or new fit',
    prior_revision=old_revision,fit_C0_Msun_kpc3=fit['C_Msun_kpc3'],source_A_Msun_kpc3=fit['source_C_before_retention_Msun_kpc3'],
    corrected_capacity_factor=2,corrected_area_per_capacity_factor=.5,galaxies_checked=149,
    rotation_predictions_and_scores_unchanged=True,
    artifact_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [fitfile,P/'capture-capacity-results.json',P/'conservative-recycling-results.json']})
(P/'capacity-normalization-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))
