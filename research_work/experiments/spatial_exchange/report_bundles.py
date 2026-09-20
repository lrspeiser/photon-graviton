"""Render completed evolving-bundle cases; never label model diagnostics observations."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
directory=ROOT/'bundle-v1'
manifest=json.loads((directory/'manifest.json').read_text())
rows=[]
for cfg in manifest['configurations']:
    path=directory/(cfg['name']+'.json')
    if path.exists():rows.append(json.loads(path.read_text())['summary'])
if not rows:raise RuntimeError('No completed bundle case to report')
fig,axes=plt.subplots(1,2,figsize=(11,4.3))
for row in rows:
    name=row['config']['name'];raw=np.load(directory/(name+'.npz'));t=raw['time'];state=raw['probe_states']
    axes[0].plot(t,1000*(state[:,0,1]-state[0,0,1]),label=name)
    radius=np.linalg.norm(state[:,9:,:3],axis=2)
    # Freely moving bodies start at r=2 with tangential canonical p=.02,m=1.
    free=np.sqrt(4+(.02/np.sqrt(1+.02**2)*t)**2)
    axes[1].plot(t,1000*(radius.mean(axis=1)-free),label=name)
axes[0].set(title='Central ray: transverse displacement',xlabel='Model time',ylabel='1000 × (y − initial y)')
axes[1].set(title='Massive bodies: mean radial response',xlabel='Model time',ylabel='1000 × (radius − free-flight radius)')
for ax in axes:ax.grid(alpha=.25);ax.legend(fontsize=8)
fig.suptitle(f'Evolving candidate field — {len(rows)}/7 cases complete; no observational data',fontsize=12)
fig.tight_layout();fig.savefig(ROOT/'bundles.png',dpi=160);plt.close(fig)
lines=['# SE-B evolving bundle checkpoint','',f'{len(rows)} of seven declared cases have complete archives.', '',
       '| Case | Central bend (model radians) | Arrival offset (model time) | Area gain | Bundle derivative error | Numerical gates |',
       '|---|---:|---:|---:|---:|---|']
for r in rows:
    if r['all_crossed']:
        lines.append(f"| {r['config']['name']} | {r['central_bend']:.9g} | {r['arrival_offset']:.9g} | {r['transport']['area_gain']:.9g} | {r['transport']['bundle_error']:.3g} | {r['passed']} |")
    else:lines.append(f"| {r['config']['name']} | incomplete crossing | — | — | — | {r['passed']} |")
lines += ['', '![Candidate trajectories](bundles.png)', '',
          'The detector is x=1.5; rays launch at x=-1.5 near y=1 with parallel +x',
          'momentum. The detector map includes the evolving background, not a frozen',
          'lens snapshot. Arrival offsets are coordinate times relative to flat flight.',
          'Area gain refers to the parallel-ray transport map, not yet an astronomical',
          'magnification. No galaxy or cluster measurements appear in this chart.', '',
          'The massive-body chart subtracts free flight from the same initial position',
          'and velocity, derived from the zero-field candidate Hamiltonian. These short',
          'test-body paths are not stable orbits, circular speeds or a rotation curve.', '',
          'Independent audit: audit_bundles.py reconstructs detector crossings and maps',
          'from every-step traces and background energy from raw endpoint states.',
          'Probe cone extrema are recorded by the runner; the archive does not contain',
          'the full field at every step for independent reconstruction of those extrema.', '',
          'Comparative emitter effects must wait for matching completed cases. Time,',
          'space and probe-radius comparisons are still required before interpreting',
          'a point-ray result. Positive source enhancement alone is not lens enhancement.',
          'Hamiltonian optics and numerical bundle derivatives are established methods;',
          'the shared response law is a candidate assumption. All twelve goals remain active.']
if (directory/'summary.json').exists():
    summary=json.loads((directory/'summary.json').read_text())
    lines += ['', '## Completed campaign verdict', '',
              f"All seven executions finished. Overall declared campaign pass: {summary['passed']}.",
              'The time comparison passes, but the spatial comparison fails. Thus the',
              'small differences between emitter choices are not resolved accurately.', '',
              '| Comparison | Bend relative difference | Arrival difference | Matrix difference | Passed |',
              '|---|---:|---:|---:|---|']
    for c in summary['comparisons']:
        lines.append(f"| {c['name']} | {c['bend_relative_difference']:.9g} | {c['arrival_difference']:.9g} | {c['matrix_difference']:.9g} | {c['passed']} |")
    base=next(r for r in rows if r['config']['name']=='Y');radius=next(r for r in rows if r['config']['name']=='probe-radius')
    difference=abs(radius['central_bend']-base['central_bend'])/max(abs(base['central_bend']),1e-8)
    lines += ['', f'Changing the probe kernel radius changes the central bend by {100*difference:.6g}%.',
              'This radius dependence is a modeling/regularization sensitivity, not an',
              'observed photon-size effect. A point-ray interpretation is not established.']
(ROOT/'bundle-report.md').write_text('\n'.join(lines)+'\n')
print('Reported',len(rows),'completed cases')
