"""Export corrected orbit/source diagnostics; retain failed convergence gates."""
from pathlib import Path
import json,csv

HERE=Path(__file__).resolve().parent
rows=[]; orbit_checks=[]; missed=[]
for grid,L,N in [('coarse',64,16),('fine',64,36),('fine',40,36)]:
    data=json.loads((HERE/f'audit-{grid}-L{L}.json').read_text())
    assert len(data['records'])==N, 'Audit unfinished'
    for i,r in enumerate(data['records']):
        orbit_checks.append(dict(grid=grid,L=L,probe=i+1,**r['final_check']))
        for resolution in ('coarse','refined','extra_refined'):
            a=r[resolution]
            if a and len(a['crossings'])!=a['naive_crossing_count']:
                missed.append(dict(grid=grid,L=L,probe=i+1,resolution=resolution,
                                   bracketed=len(a['crossings']),naive=a['naive_crossing_count']))
    for R in (1.,3.):
        group=[r for r in data['records'] if r['R']==R]
        norm=sum(r['raw_weight'] for r in group)
        for j,T in enumerate((.05,.1,.25)):
            cohort=history=residence=late=0.
            for r in group:
                a=r['extra_refined'] or r['refined'];w=r['raw_weight']/norm;t=a['entry_time']
                if t is not None and t<=T:
                    cohort+=w;history+=w*(T-t)/T
                    if a['entry_after_turn']:late+=w
                residence+=w*a['inside_durations'][j]/T
            rows.append(dict(grid=grid,L=L,R_kpc=R,T=T,T_Myr=T*977.79222168,
                             cohort_ever_entered=cohort,continuous_ever_entered=history,
                             continuous_inside=residence,cohort_entry_after_turn=late))
with (HERE/'summary.csv').open('w',newline='',encoding='utf8') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
gates=[]
for name,grid,L,thresholds in [('angular','coarse',64,(.05,.05,.005)),('field','fine',40,(.02,.02,.005))]:
    for a in [r for r in rows if r['grid']==grid and r['L']==L]:
        b=next(r for r in rows if r['grid']=='fine' and r['L']==64 and r['R_kpc']==a['R_kpc'] and r['T']==a['T'])
        for metric,limit in zip(('cohort_ever_entered','continuous_ever_entered','continuous_inside'),thresholds):
            delta=abs(a[metric]-b[metric])
            gates.append(dict(kind=name,R_kpc=a['R_kpc'],T=a['T'],metric=metric,
                              difference=delta,limit=limit,passes=delta<limit))
out=dict(rows=rows,orbit_checks=orbit_checks,statistical_gates=gates,
         naive_crossing_disagreements=missed,
         all_orbit_gates_pass=all(r['passes'] for r in orbit_checks),
         all_statistical_gates_pass=all(r['passes'] for r in gates))
(HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(dict(orbit_pass=out['all_orbit_gates_pass'],failed_gates=[r for r in gates if not r['passes']],
                     naive_disagreements=missed,fine64=[r for r in rows if r['grid']=='fine' and r['L']==64]),indent=2))
