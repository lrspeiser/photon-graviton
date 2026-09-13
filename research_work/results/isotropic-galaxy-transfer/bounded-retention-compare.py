"""Compare frozen metrics and verify identical observed inputs."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
prefixes=['','radiative-flux-','bounded-radiation-retention-','bounded-depth-retention-']
rows=[]
baseline=json.loads((HERE/'predictions.json').read_text())
original={r['galaxy']:r for r in baseline if r['model']=='attenuated'}
for prefix in prefixes:
    model=json.loads((HERE/(prefix+'results.json')).read_text())['models']['attenuated']
    predictions=json.loads((HERE/(prefix+'predictions.json')).read_text())
    for row in predictions:
        if row['model']!='attenuated':continue
        before=original[row['galaxy']]
        assert before['observed_kms']==row['observed_kms'] and before['R_kpc']==row['R_kpc'] and before['split']==row['split']
    rows.append(dict(candidate=prefix.rstrip('-') or 'original',q=model.get('q'),C_Msun_kpc3=model['C_Msun_kpc3'],k0_per_kpc=model['k0_per_kpc'],scale_to_disk=model['scale_to_disk'],scores=model['finer_scores'],boundary=model['boundary'],optimizer_success=model['optimizer_success'],retention_range=[model.get('minimum_retention'),model.get('maximum_retention')]))
(HERE/'bounded-retention-comparison.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
for r in rows:
    print(r['candidate'],'q',r['q'],'retention',r['retention_range'])
    for s,v in r['scores'].items(): print(s,v['RMSE_kms'],v['log_RMS'])
