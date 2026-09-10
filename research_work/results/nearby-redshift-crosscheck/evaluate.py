"""Evaluate only the eight identities fixed in protocol.json; never refit."""
from pathlib import Path
import json,hashlib,csv
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CACHE=ROOT/'research_work/data-cache/nearby-redshift-crosscheck'
P=json.loads((HERE/'protocol.json').read_text())
for path,digest in P['input_hashes'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
ids={t['agc'] for t in P['targets']};names={t['target_name'] for t in P['targets']}
hi={int(line[:6]):line for line in (CACHE/'haynes-table2.dat').read_text().splitlines() if line.strip() and int(line[:6]) in ids}
body=(CACHE/'elves-table.txt').read_text().rsplit('-'*80,1)[1]
elves={line[:13].strip():line for line in body.splitlines() if line[:13].strip() in names}
frames={r['target_name']:float(r['solar_to_cmb_spectral_factor']) for r in csv.DictReader((HERE.parent/'redshift-priority/elves-frame-factors.csv').open())}
rows=[]
for t in P['targets']:
 a=hi[t['agc']];b=elves[t['target_name']]
 ra=15*(float(a[32:34])+float(a[34:36])/60+float(a[36:40])/3600)
 dec=(1 if a[40:41]=='+' else -1)*(float(a[41:43])+float(a[43:45])/60+float(a[45:47])/3600)
 direction=SkyCoord(float(b[14:22])*u.deg,float(b[23:30])*u.deg)
 sep=direction.separation(SkyCoord(ra*u.deg,dec*u.deg)).arcsec
 assert sep<30 and a[112:113]=='1' and b[31:32]==t['distance_method']
 d=float(b[33:40]);derr=float(b[41:47]);vh=float(a[48:53]);ve=float(a[58:61])/2
 observed=(1+vh/P['c_kms'])*frames[t['target_name']]-1
 exp=np.expm1(P['alpha_per_mpc']*d);linear=P['linear_alpha_per_mpc']*d
 rows.append({**t,'distance_mpc':d,'published_distance_error_mpc':derr,'heliocentric_optical_velocity_kms':vh,'statistical_center_error_kms':ve,
   'source_position_separation_arcsec':float(sep),'solar_to_cmb_factor':frames[t['target_name']], 'observed_cmb_z':observed,
   'exponential_predicted_z':exp,'linear_predicted_z':linear,'exponential_residual_kms':P['c_kms']*(exp-observed),'linear_residual_kms':P['c_kms']*(linear-observed),
   'below_original_calibration_distance_range':d<10.204694319421995})
assert len(rows)==8
def score(group,key):
 residual=np.array([r[key] for r in group]);return {'n':len(group),'rms_kms':float(np.sqrt(np.mean(residual**2))),'mae_kms':float(np.mean(abs(residual))),'bias_kms':float(residual.mean())}
scores={model:score(rows,model+'_residual_kms') for model in ['exponential','linear']}
out={'scope':P['scope'],'rows':rows,'scores':scores,'by_distance_method':{m:{model:score([r for r in rows if r['distance_method']==m],model+'_residual_kms') for model in ['exponential','linear']} for m in sorted({r['distance_method'] for r in rows})},
 'protocol_sha256':hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),'parameters_refitted':False,'outcomes_now_exposed':True}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
with (HERE/'predictions.csv').open('w',newline='',encoding='utf-8') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
print(json.dumps(scores,indent=2))
