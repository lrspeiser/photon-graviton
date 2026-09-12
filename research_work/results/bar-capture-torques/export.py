from pathlib import Path
import json,csv
import numpy as np

HERE=Path(__file__).resolve().parent
data={mode:json.loads((HERE/f'{mode}.json').read_text()) for mode in ('spherical','static','rotating')}
rows=[];field_checks=[]
for mode,d in data.items():
    assert all(r['check']['passes'] for r in d['records'])
    for r in d['records']:
        if r['order']!=64:continue
        x=r['refined']
        rows.append(dict(mode=mode,launch_radius=r['launch_radius'],phi=r['phi'],z_over_r=r['z_over_r'],event=x['event'],
                         event_radius=x['radius'],event_time=x['time'],Lz=x['L_final'][2],L_norm=float(np.linalg.norm(x['L_final'])),
                         delta_E=x['delta_E'],rotation_work=x['rotation_work'],invariant_error=x['invariant_error_over_220_squared']))
        if mode!='spherical':
            other=next(o for o in d['records'] if o['order']==40 and all(o[k]==r[k] for k in ('launch_radius','phi','z_over_r')))
            error=abs(x['radius']-other['refined']['radius'])
            field_checks.append(dict(mode=mode,launch_radius=r['launch_radius'],phi=r['phi'],z_over_r=r['z_over_r'],
                                     radius_difference=error,event_agreement=x['event']==other['refined']['event'],passes=error<.01 and x['event']==other['refined']['event']))
assert all(c['passes'] for c in field_checks)
(HERE/'field-checks.json').write_text(json.dumps(field_checks,indent=2)+'\n',encoding='utf8',newline='\n')
with (HERE/'comparison.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
lines=['## Final checkpoint','',
'| Field | Probes entering0.1kpc before first turn | First-turn radius range (kpc) |',
'|---|---:|---:|']
for mode in data:
    selected=[r for r in rows if r['mode']==mode];turns=[r['event_radius'] for r in selected if r['event']=='first_turn']
    interval=f'{min(turns):.6f}–{max(turns):.6f}' if turns else 'Not measured: stopped on entry'
    lines.append(f"| {mode} | {sum(r['event']=='core_entry' for r in selected)}/8 | {interval} |")
max_radius=max(c['radius_difference'] for c in field_checks)
max_j=max(r['invariant_error'] for r in rows)
lines+=['',f'All integration and harmonic-order event checks pass. The largest L40/L64 first-turn radius change is {max_radius:.6g}kpc; the largest retained invariant drift divided by220^2 is {max_j:.4g}. These gates concern the numerical bar-component experiment, not its astrophysical normalization.', '',
'Representative probe: launch radius3kpc, azimuth pi/6, z/r=.2:','',
'| Field | Event | Radius (kpc) | Final Lz (kpc km/s) | Energy change ((km/s)^2) | Omega delta Lz ((km/s)^2) |',
'|---|---|---:|---:|---:|---:|']
for mode in data:
    r=next(r for r in rows if r['mode']==mode and r['launch_radius']==3 and r['phi']<.6 and r['z_over_r']==.2)
    lines.append(f"| {mode} | {r['event']} | {r['event_radius']:.6f} | {r['Lz']:.6f} | {r['delta_E']:.6f} | {r['rotation_work']:.6f} |")
lines+=['','The eight directions are deliberately chosen probes, not an isotropic sample or a measured capture fraction. The result establishes first-pass deflection and angular-momentum exchange in this conditional bar field. Later passages, symmetry-axis trajectories and a continuously supplied population remain untested.', '',
'[comparison.csv](comparison.csv) retains all24 preferred-order results; the JSON files preserve both field orders, tolerances, events and source-cache hashes.','']
p=HERE/'report.md';s=p.read_text(encoding='utf8').split('## Final checkpoint')[0]
p.write_text(s+'\n'+'\n'.join(lines),encoding='utf8',newline='\n')
print('\n'.join(lines))
