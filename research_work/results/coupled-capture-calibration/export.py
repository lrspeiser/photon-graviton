from pathlib import Path
import json,csv
HERE=Path(__file__).resolve().parent
p=json.loads((HERE/'results.json').read_text())
bins=json.loads((HERE.parent/'cepheid-common-frame/training-bins.json').read_text())
rows=[]
for i,b in enumerate(bins):
    row=dict(bin=i,R_mean_kpc=b['R_mean_kpc'],stars=b['n'],inferred_Jeans_proxy_km_s=p['observed_training_proxy'][i],ordinary_prediction_km_s=p['ordinary_prediction'][i])
    for fit in p['fits']:
        row[f"p{fit['p']}_coupled_prediction_km_s"]=fit['refined']['predictions'][i]
        row[f"p{fit['p']}_residual_km_s"]=fit['refined']['residuals'][i]
    rows.append(row)
with (HERE/'radial-comparison.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
