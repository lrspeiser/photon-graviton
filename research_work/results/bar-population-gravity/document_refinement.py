from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'source-refinement.json').read_text())
pfailed=sum(not r['potential_pass'] for r in d['comparisons'])
afailed=sum(not r['acceleration_pass'] for r in d['comparisons'])
lines=['# Direction-sampling test of the gravitational response','',
'This holds the source law and ordinary bar fixed and compares6x12 with8x16 angular sampling. The source is still a negligible test population, separately normalized on1 and3kpc launch spheres. No energy supply, self-gravity, observed star fit or new physical interaction is introduced.','',
f"All {len(d['orbit_checks'])} combined orbit/age status checks pass: **{d['all_orbit_checks_pass']}**. Of30 position/time/source comparisons, {pfailed} fail the2 percent potential gate and {afailed} fail the5 percent vector-force gate. These development thresholds are not measurement error bars. A passing pair is useful evidence but not proof of convergence.",
'','## Final-epoch comparison','',
'At approximately244Myr, values below are per G times total injected rest mass. Potential units are1/kpc; acceleration units are1/kpc². Changes use the refined magnitude with the predeclared0.01 floor.','',
'| Source R | Position | Refined potential / GM | Refined acceleration / GM | Potential change (%) | Force change (%) |',
'|---:|---|---:|---|---:|---:|']
for c in d['comparisons']:
    if c['T']==.25:
        r=c['refined'];pos=', '.join(f'{x:g}' for x in c['point']);acc=', '.join(f'{x:.5g}' for x in r['acceleration'])
        lines.append(f"| {c['R']:g} | ({pos}) | {r['potential']:.6g} | ({acc}) | {100*c['potential_scaled_change']:.3f} | {100*c['acceleration_scaled_change']:.3f} |")
lines+=['','## Interpretation and next step','',
'Potential averages contributions over the moving population. Local force is more sensitive to where that material sits relative to a test point, so a stable potential or enclosed central fraction cannot substitute for a force check. Resolve failed locations before interpreting directional force differences as a bulge prediction. Keep actual source anisotropy distinct from sampling error.','',
'The equations are the known Newtonian Green functions applied conditionally to the project cold-rest-mass source. Their derivation, units, age quadrature and analytic checks remain in [the original report](report.md). The refined force calculations do not include a fresh L40/L64 comparison, full ordinary Galactic components or collective evolution.','',
'All raw coarse and refined outputs are retained. The next required step is further force-source refinement or a demonstrably convergent integration method, followed by source normalization and evolving capture/self-gravity. A metric/stress law is needed for lensing. Joint redshift/timing predictions and all nine overall goals remain incomplete.']
(HERE/'source-refinement.md').write_text('\n'.join(lines)+'\n',encoding='utf8',newline='\n')
