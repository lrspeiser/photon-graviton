"""Carry the Cepheid amplitude into unrefitted vertical force predictions."""
from pathlib import Path
import json,hashlib
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[2]
source=H.parent/'component-loading-forces/results.json';scale=H.parent/'component-energy-scale/results.json'
f=json.loads(source.read_text(encoding='utf-8'));s=json.loads(scale.read_text(encoding='utf-8'));k=s['results'][1]['best_rms_scale']
for d in [f,s]:
 for name,value in d['hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==value
rows=[]
for v in f['rows']:
 loc=v['location'];ordinary=np.array(v['ordinary_acceleration']);target=np.array(v['target_extra_acceleration']);extra=np.array(v['fine_extra_acceleration']);low=np.array(v['coarse_extra_acceleration'])
 if loc['z_kpc']==0:continue
 reference=ordinary+target;predicted=ordinary+k*extra;coarse=ordinary+k*low
 assert reference[2]<0 and predicted[2]<0
 radial=np.array([np.cos(loc['phi_rad']),np.sin(loc['phi_rad']),0])
 rows.append(dict(R_kpc=loc['R_kpc'],z_kpc=loc['z_kpc'],phi_rad=loc['phi_rad'],ordinary_vertical_pull=-ordinary[2],target_total_vertical_pull=-reference[2],predicted_total_vertical_pull=-predicted[2],vertical_ratio=float(predicted[2]/reference[2]),radial_ratio=float(predicted@radial/(reference@radial)),vertical_refinement_relative=float(abs(coarse[2]-predicted[2])/abs(reference[2])) ))
summary=[]
for name,subset in [('inner',[v for v in rows if v['R_kpc']<=8 and v['z_kpc']<=1]),('all_offplane',rows),('solar_R8',[v for v in rows if v['R_kpc']==8])]:
 ratios=np.array([v['vertical_ratio'] for v in subset]);summary.append(dict(subset=name,n=len(subset),minimum_vertical_ratio=float(ratios.min()),median_vertical_ratio=float(np.median(ratios)),maximum_vertical_ratio=float(ratios.max()),maximum_refinement_relative=max(v['vertical_refinement_relative'] for v in subset)))
out=dict(scope='Unrefitted vertical force transfer of an exposed Cepheid amplitude; empirical target comparison, not measured vertical gravity',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),scale_sha256=hashlib.sha256(scale.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),amplitude=k,units='(km/s)^2/kpc for accelerations',summaries=summary,rows=rows)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(summary))
