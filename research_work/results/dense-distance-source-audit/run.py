from pathlib import Path
import csv,json,hashlib
from collections import Counter
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
raw=ROOT/'research_work/data-cache/dense-distance-audit/tf-observer-features.tsv'
metadata=ROOT/'research_work/data-cache/dense-distance-audit/tf-ReadMe.txt'
cal=ROOT/'research_work/results/extended-environment-geometry/selected-tracers.csv'
lines=raw.read_text().splitlines();start=next(i for i,s in enumerate(lines) if s.startswith('PGC\t'))
rows=list(csv.DictReader([lines[start]]+[s for s in lines[start+3:] if s.strip() and not s.startswith('#')],delimiter='\t'))
assert len(rows)==10737
allowed=['PGC','Name','Inc','e_Inc','Wmx','e_Wmx','gmag','rmag','imag','W1mag','W2mag','Rei','ReW1','icmag','Ai']
assert list(rows[0])==allowed
clean=[{k:v.strip() for k,v in r.items()} for r in rows]
assert len({int(r['PGC']) for r in clean})==len(clean)
def write(name,data):
    with (HERE/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
write('observer-features.csv',clean)
calibrators={int(r['pgc']):r for r in csv.DictReader(cal.open())}
matched=[]
for r in clean:
    if int(r['PGC']) in calibrators:
        c=calibrators[int(r['PGC'])]
        matched.append(dict(pgc=r['PGC'],name=r['Name'],method=c['method'],
            distance_mpc=c['distance_mpc'],sigma_mag=c['sigma_mag'],
            imag=r['imag'],W1mag=r['W1mag'],inclination_deg=r['Inc'],
            adjusted_linewidth_kms=r['Wmx'],adjusted_linewidth_error_kms=r['e_Wmx'],
            role='possible_calibration; shared_history_and_corrections_pending; not_holdout'))
write('possible-calibrators.csv',matched)
def present(r,k):return bool(r[k]) and float(r[k])!=-100
result=dict(source_rows=len(rows),unique_pgc=len(clean),
    photometry_present={k:sum(present(r,k) for r in clean) for k in ['gmag','rmag','imag','W1mag','W2mag']},
    exact_overlap_with_previously_selected_indicator_tracers=len(matched),
    overlap_methods=dict(Counter(r['method'] for r in matched)),
    overlap_with_i_photometry=sum(present(r,'imag') for r in matched),
    overlap_with_W1_photometry=sum(present(r,'W1mag') for r in matched),
    source_query='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJ/902/145/table1&-out=PGC,Name,Inc,e_Inc,Wmx,e_Wmx,gmag,rmag,imag,W1mag,W2mag,Rei,ReW1,icmag,Ai&-out.max=unlimited',
    interpretation='Measurement-stage inputs only. Wmx is already adjusted, icmag includes corrections. No velocities or distance-outcome columns queried; features are not certified fresh observations.',
    raw_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw,metadata,cal]},
    normalized_feature_sha256=hashlib.sha256((HERE/'observer-features.csv').read_bytes()).hexdigest())
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
