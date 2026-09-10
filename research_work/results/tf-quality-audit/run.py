from pathlib import Path
import csv,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
raw=ROOT/'research_work/data-cache/tf-quality/features.tsv'
lines=raw.read_text().splitlines(); start=next(i for i,s in enumerate(lines) if s.startswith('PGC\t'))
quality=list(csv.DictReader([lines[start]]+[s for s in lines[start+3:] if s.strip() and not s.startswith('#')],delimiter='\t'))
quality=[{k:v.strip() for k,v in r.items()} for r in quality]
assert len(quality)==10737 and len({r['PGC'] for r in quality})==10737
lookup={int(r['PGC']):r for r in quality}
cal=ROOT/'research_work/results/indicator-tf-calibration'
source=list(csv.DictReader((cal/'cross-validation.csv').open()))
source=[r for r in source if r['photometry_p']=='1']
slope=json.loads((cal/'results.json').read_text())['scenarios'][1]['slope']
rows=[]
for r in source:
    q=lookup[int(r['pgc'])]
    vals={k:float(q[k]) if q[k] else np.nan for k in ['e_Ai','logWmxi','e_logWmxi']}
    scale=np.sqrt(float(r['sigma_mu'])**2+0.05**2+vals['e_Ai']**2+slope**2*vals['e_logWmxi']**2)
    rows.append(dict(pgc=r['pgc'],name=r['name'],QSflag=int(q['QSflag']),QWflag=int(q['QWflag']),
        **vals,reconstructed_logW=float(r['x'])+2.5,reconstructed_sigma_logW=float(r['sigma_x']),
        residual_mag=float(r['residual_mag']),distance_ratio=float(r['distance_ratio']),measurement_only_scale_mag=scale))
assert len(rows)==73
with (HERE/'joined.csv').open('w',encoding='utf-8',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
def stats(items):
    if not items: return dict(n=0)
    return dict(n=len(items),rms_mag=float(np.sqrt(np.mean([r['residual_mag']**2 for r in items]))),
        median_abs_fractional_distance_error=float(np.median([abs(r['distance_ratio']-1) for r in items])),
        median_measurement_only_scale_mag=float(np.nanmedian([r['measurement_only_scale_mag'] for r in items])),
        finite_measurement_scales=sum(bool(np.isfinite(r['measurement_only_scale_mag'])) for r in items),
        beyond_three_measurement_only_scales=sum(bool(abs(r['residual_mag'])>3*r['measurement_only_scale_mag']) for r in items))
result=dict(source_rows=len(quality),joined=len(rows),by_QSflag={str(k):stats([r for r in rows if r['QSflag']==k]) for k in range(6)},
    quality_at_least_4=stats([r for r in rows if r['QSflag']>=4]),all=stats(rows),
    max_abs_logW_reconstruction_difference=max(abs(r['logWmxi']-r['reconstructed_logW']) for r in rows),
    median_abs_sigma_logW_difference=float(np.median([abs(r['e_logWmxi']-r['reconstructed_sigma_logW']) for r in rows])),
    max_abs_sigma_logW_difference=max(abs(r['e_logWmxi']-r['reconstructed_sigma_logW']) for r in rows),
    ngc4424=next(r for r in rows if int(r['pgc'])==40809),
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw,cal/'cross-validation.csv',cal/'results.json']})
(HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
