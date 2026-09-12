"""Compare fixed-law angular quadratures without retuning acceptance gates."""
from pathlib import Path
import json,csv
HERE=Path(__file__).resolve().parent
rows=[];checks=[]
old=json.loads((HERE.parent/'bar-source-persistence/audit-fine-L40.json').read_text())
datasets=[(6,old['records'],lambda r:r['extra_refined'] or r['refined'])]
for n in (8,12):
    d=json.loads((HERE/f'grid{n}-L40.json').read_text())
    assert len(d['records'])==d['expected']==n*n
    assert d['source_hash']==old['source_hash']
    datasets.append((n,d['records'],lambda r:r['refined'] or r['primary']))
    checks.extend(dict(n=n,probe=i+1,passes=r['passes'],comparison=r['comparison'],
                       invariant_error=(r['refined'] or r['primary'])['invariant_error']) for i,r in enumerate(d['records']))
for n,records,get in datasets:
    for R in (1.,3.):
        group=[r for r in records if r['R']==R];norm=sum(r['raw_weight'] for r in group)
        for j,T in enumerate((.05,.1,.25)):
            cohort=history=inside=late=0.
            for r in group:
                a=get(r);w=r['raw_weight']/norm;t=a['entry_time']
                if t is not None and t<=T:
                    cohort+=w;history+=w*(T-t)/T
                    late+=w*bool(a['entry_after_turn'])
                inside+=w*a['inside_durations'][j]/T
            assert 0<=inside<=history<=cohort<=1+1e-12
            rows.append(dict(n=n,R=R,T=T,cohort_ever_entered=cohort,continuous_ever_entered=history,
                             continuous_inside=inside,late_entry=late,source_norm=norm,
                             maximum_normalized_direction_weight=max(r['raw_weight']/norm for r in group)))
comparisons=[]
for na,nb in ((6,8),(8,12)):
    for a in [r for r in rows if r['n']==na]:
        b=next(r for r in rows if r['n']==nb and r['R']==a['R'] and r['T']==a['T'])
        for metric,limit in [('cohort_ever_entered',.05),('continuous_ever_entered',.05),('continuous_inside',.005)]:
            delta=abs(a[metric]-b[metric]);relative=delta/max(abs(b[metric]),1e-30)
            comparisons.append(dict(grid_a=na,grid_b=nb,R=a['R'],T=a['T'],metric=metric,
                                    absolute_change=delta,relative_to_finer=relative,limit=limit,passes=delta<limit))
out=dict(rows=rows,checks=checks,comparisons=comparisons,all_orbit_checks_pass=all(r['passes'] for r in checks))
(HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
with (HERE/'summary.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
print(json.dumps(dict(orbit_checks=out['all_orbit_checks_pass'],failed=[r for r in comparisons if not r['passes']],
                     final_epoch=[r for r in rows if r['T']==.25]),indent=2))
