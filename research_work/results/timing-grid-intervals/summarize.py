from pathlib import Path
import json
H=Path(__file__).resolve().parent;r=json.loads((H/'results.json').read_text(encoding='utf-8'))
assert len(r['cases'])==4
lines=['# Timing interval sensitivity to event-grid resolution','','For both exposed artificial samples, we profiled mean duration and intrinsic scatter at fixed timing exponent, using the same bounds and finite-support constraint as the revised estimator. Three optimizer starts were retained per evaluation. We located lower and upper crossings of the nominal chi-square(1) 95% likelihood-ratio threshold on both 321- and 641-node grids.','','The threshold and profile likelihood are known statistical methods, not new physics. These endpoints are nominal local profile crossings, not empirically calibrated confidence limits. The search does not establish absence of distant disconnected likelihood regions or prove global nuisance optimization.','','| Sample | Grid nodes | Lower crossing | Upper crossing |','|---|---:|---:|---:|']
for c in r['cases']:
 lo,hi=c['nominal_95_endpoints'];lines.append(f"| {c['label']} | {c['nodes']} | {lo:.8f} | {hi:.8f} |")
changes=[]
for label in sorted(set(c['label'] for c in r['cases'])):
 pair={c['grid']:c for c in r['cases'] if c['label']==label}
 changes.extend(abs(a-b) for a,b in zip(pair['coarse']['nominal_95_endpoints'],pair['fine']['nominal_95_endpoints']))
lines+=['',f'Maximum absolute crossing change: {max(changes):.8g}.','','All input arrays and source hashes are recorded. Every optimizer trial is retained, including unsuccessful or inferior fits. These two numerical checks do not establish population-wide coverage, nor do they validate a physical interpretation of redshift. The full 160-case development recalibration remains separate and unchanged. Wider grid checks and fresh-seed calibration remain necessary. All six goals remain open.','']
(H/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n');print({'maximum_endpoint_change':max(changes)})
