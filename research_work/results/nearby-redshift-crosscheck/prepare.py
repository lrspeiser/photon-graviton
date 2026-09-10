"""Freeze identities and analysis before extracting the selected outcomes."""
from pathlib import Path
import json,hashlib,csv,gzip
import numpy as np
import requests
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PRIOR=HERE.parent/'redshift-priority'
CACHE=ROOT/'research_work/data-cache/nearby-redshift-crosscheck';CACHE.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if (HERE/'protocol.json').exists():raise RuntimeError('Do not overwrite frozen protocol')
registry=json.loads((PRIOR/'elves-quarantine-screen.json').read_text())
pending={r['target_name'] for r in registry['records'] if r['decision']=='pending_freshness_and_host_audit'}
radio=json.loads((PRIOR/'elves-radio-quality.json').read_text())
features={x['target_name']:x for x in json.loads((PRIOR/'elves-field-feature-audit.json').read_text())['records']}
targets=[]
for r in radio['records']:
 if r['target_name'] not in pending:continue
 hits=[x for x in r['matches'] if x['explicit_supported_alias'] and x['notes']=='1']
 if len(hits)==1:targets.append({'target_name':r['target_name'],'agc':hits[0]['agc'],'distance_method':features[r['target_name']]['distance_method']})
assert len(targets)==8 and len(pending)==24
p=CACHE/'haynes-table2.dat'
if not p.exists():
 response=requests.get('https://cdsarc.cds.unistra.fr/ftp/J/ApJ/861/49/table2.dat.gz',timeout=45);response.raise_for_status()
 p.write_bytes(gzip.decompress(response.content))
# Only ID, optical counterpart position and quality are read before freezing.
catalog={int(line[:6]):line for line in p.read_text().splitlines() if line.strip()}
for target in targets:
 line=catalog[target['agc']]
 assert line[112:113]=='1'
 target['haynes_class']=1
 target['haynes_optical_position_present']=bool(line[32:47].strip())
 assert target['haynes_optical_position_present']
source=Path(r'C:/Users/henry/Documents/Codex/2026-09-09/cr/work/catalog-audit/elves-source-0.txt')
assert sha(source)=='890b1d178feb505cb2d9105926d13b77bc639d882e18776a3b3561b0650d8903'
(CACHE/'elves-table.txt').write_bytes(source.read_bytes())
# Record accidental display of first 65 lines during header inspection, without reprinting values.
names=set(features)
displayed=[line[:13].strip() for line in source.read_text().splitlines()[:65] if line[:13].strip() in names]
selected_names={t['target_name'] for t in targets}
exposure={'header_inspection_first_65_lines_displayed_before_freeze':displayed,'selected_targets_in_that_display':sorted(selected_names&set(displayed)),
 'scope':'Historical distances and cross-task object exposure may remain. This is not certified fresh validation.'}
assert not exposure['selected_targets_in_that_display']
old=list(csv.DictReader((ROOT/'redshift_paper/all_164_groups.csv').open(encoding='utf-8-sig')))
train=[x for x in old if x['split']=='train'];d=np.array([float(x['catalog_distance_mpc']) for x in train]);z=np.array([float(x['observed_cmb_z']) for x in train])
linear=float(d@z/(d@d));alpha=json.loads((ROOT/'papers/cumulative-time-companions/analysis/metrics.json').read_text())['alpha_per_mpc']
inputs=[p,CACHE/'elves-table.txt',PRIOR/'elves-quarantine-screen.json',PRIOR/'elves-radio-quality.json',PRIOR/'elves-frame-factors.csv',ROOT/'redshift_paper/all_164_groups.csv',HERE/'evaluate.py']
protocol={'targets':targets,'alpha_per_mpc':alpha,'linear_alpha_per_mpc':linear,'c_kms':299792.458,
 'scope':'Fixed-formula conditional cross-catalog diagnostic, not certified independent/fresh or a complete motion/selection likelihood.',
 'distance':'ELVES publisher adopted TRGB/SBF distances, unchanged. No ALFALFA flow distances imported.',
 'velocity':'Corrected Haynes ALFALFA catalog heliocentric optical midpoint velocity, not presumed identical to the ELVES adopted velocity. Center statistical uncertainty is half e_W50.',
 'frame':'Use frozen exact solar-to-CMB spectral factors. z_helio=Vhel/c; 1+z_CMB=(1+z_helio)*F.',
 'models':'Frozen exponential and historical training-fitted linear control; no intercept or per-object velocity fitted.',
 'outputs':'All eight rows, residual c*(predicted_z-observed_z), RMS/MAE/bias, distance-method summaries. No residual exclusions or significance/acceptance threshold. Below original distance calibration range is flagged as extrapolation.',
 'identity_check':'Haynes optical position must be within 30 arcsec of frozen ELVES direction; stop on mismatch, do not choose by residual.',
 'exposure':exposure,'input_hashes':{str(x.relative_to(ROOT)):sha(x) for x in inputs}}
(HERE/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'targets':targets,'exposure':exposure},indent=2))
