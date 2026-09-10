from pathlib import Path
import json,hashlib
H=Path(__file__).resolve().parent
load=lambda p:json.loads(p.read_text())
coarse=load(H/'predictions-updated-profile-validation.json')
fine=load(H/'predictions-refined-updated-profile-validation.json')
result=load(H/'results-refined-updated-profile-validation.json')
training=load(H/'results-refined-updated-profile.json')
catalog=load(H.parent/'lensing-data-readiness/lens-observations-and-image-models.json')
eligible={r['Name'] for r in catalog if r['role']=='validation' and r['Mph']=='E' and r['spectroscopic_dispersion_available']}
assert {r['Name'] for r in fine}==eligible and len(eligible)==7
assert len(fine)==56 and all(r['role']=='validation' for r in fine)
key=lambda r:(r['Name'],r['model'],r['seeing_fwhm_arcsec'],r['cutoff_over_Re'])
assert [key(r) for r in fine]==[key(r) for r in coarse]
assert result['parameters']==training['parameters']
assert result['protocol_sha256']==hashlib.sha256((H/'validation-protocol.json').read_bytes()).hexdigest()
delta=max(abs(r['theta_pred_arcsec']-s['theta_pred_arcsec']) for r,s in zip(fine,coarse));assert delta<.001
primary={m:[r for r in fine if r['model']==m and r['seeing_fwhm_arcsec']==1.5 and r['cutoff_over_Re']==20] for m in ['baryons','empirical_companion']}
out={'eligible_systems':7,'all_eligible_systems_retained':True,'shared_parameters_unchanged':True,
     'max_refinement_difference_arcsec':delta,'validation_scores_opened':True,'test_scores_opened':False,
     'primary_overpredicted_systems':{m:sum(r['residual_arcsec']>0 for r in s) for m,s in primary.items()},
     'primary_companion_absolute_residual_smaller_count':sum(abs(a['residual_arcsec'])<abs(b['residual_arcsec']) for a,b in zip(primary['empirical_companion'],primary['baryons']))}
(H/'validation-verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
