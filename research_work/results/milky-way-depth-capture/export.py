from pathlib import Path
import json,csv
HERE=Path(__file__).resolve().parent
p=json.loads((HERE/'results.json').read_text())
p['vertical_admissibility']='Excluded from accepted observational validation: published Kz inferences use a potential family including a dark halo. Retained only as a historical alternate-model diagnostic; not fitted.'
(HERE/'results.json').write_text(json.dumps(p,indent=2)+'\n',encoding='utf8',newline='\n')
bins=json.loads((HERE.parent/'cepheid-common-frame/training-bins.json').read_text())
rows=[]
for i,b in enumerate(bins):
    row=dict(bin=i,R_mean_kpc=b['R_mean_kpc'],stars=b['n'],inferred_Jeans_proxy_km_s=p['radial_observed'][i],ordinary_prediction_km_s=p['ordinary_radial_prediction'][i])
    for r in p['runs']:
        if r['refined']:row[f"p{r['p']}_outer{r['outer_kpc']:g}_km_s"]=r['radial_predictions'][i]
    rows.append(row)
with (HERE/'radial-comparison.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
