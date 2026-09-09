"""Fixed-parameter external diagnostic; neither fitting nor blind validation."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
C=299792.458
ALPHA=.0002488993286382367
LINEAR_K=75.17651383428588
# Pesce et al., arXiv:2001.09213v2, Table 1, visually checked against PDF page 4.
# name, distance, lower/upper 16-84 percentile errors, optical CMB c*z, error
data=[('UGC 3789',51.5,4.,4.5,3319.9,.8),('NGC 6264',132.1,17.,21.,10192.6,.8),
      ('NGC 6323',109.4,23.,34.,7801.5,1.5),('NGC 5765b',112.2,5.1,5.4,8525.7,.7),
      ('CGCG 074-064',87.6,7.2,7.9,7172.2,1.9),('NGC 4258',7.58,.11,.11,679.3,.4)]
features=[dict(name=n,distance_mpc=d,distance_error_minus=lo,distance_error_plus=hi) for n,d,lo,hi,_,_ in data]
(HERE/'maser-features.json').write_text(json.dumps(features,indent=2)+'\n',newline='\n')
predictions=[dict(name=f['name'],exponential_z=float(np.expm1(ALPHA*f['distance_mpc'])),linear_z=LINEAR_K*f['distance_mpc']/C) for f in features]
(HERE/'maser-frozen-predictions.json').write_text(json.dumps(predictions,indent=2)+'\n',newline='\n')
seal=dict(prediction_sha256=hashlib.sha256((HERE/'maser-frozen-predictions.json').read_bytes()).hexdigest(),protocol_sha256=hashlib.sha256((HERE/'external-maser-protocol.md').read_bytes()).hexdigest(),source_url='https://arxiv.org/pdf/2001.09213v2',source_pdf_sha256='c9f2d2b626fee60d77928856d217a2cbe546f06c69cd9e7ac804f4ab652e731a',status='Protocol/rates fixed before table extraction, but labels viewed during extraction before prediction-file seal; not blinded. Full alias/group overlap unresolved.')
(HERE/'maser-seal.json').write_text(json.dumps(seal,indent=2)+'\n',newline='\n')
rows=[]
for f,p,(_,d,lo,hi,cz,e) in zip(features,predictions,data):
    rows.append(dict(**f,observed_cmb_z=cz/C,cz_error_km_s=e,**{k:v for k,v in p.items() if k!='name'},exponential_residual_km_s=C*p['exponential_z']-cz,linear_residual_km_s=C*p['linear_z']-cz,range_status='within' if 10.2046943194<=d<=93.1965852566 else 'extrapolation'))
def metrics(selected,key):
    r=np.array([v[key+'_residual_km_s'] for v in selected])
    return dict(n=len(r),rms=float(np.sqrt(np.mean(r*r))),mae=float(np.mean(abs(r))),median_absolute=float(np.median(abs(r))),bias=float(r.mean()))
result={part:{name:metrics([r for r in rows if part=='all' or r['range_status']==part],name) for name in ['exponential','linear']} for part in ['all','within','extrapolation']}
with (HERE/'maser-comparison.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
(HERE/'maser-results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
