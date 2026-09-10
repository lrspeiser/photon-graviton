from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
r=json.loads((HERE/'results.json').read_text())
for rel,expected in r['input_hashes'].items():
    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==expected,rel
p=json.loads((HERE/'axisymmetric_extra-individual-refinement.json').read_text())
original=ROOT/'research_work/data-cache/full-bar-orbits/axisymmetric_extra.npz'
assert hashlib.sha256(original.read_bytes()).hexdigest()==p['original_trajectory_sha256']
assert hashlib.sha256((HERE/'refine_individual.py').read_bytes()).hexdigest()==p['code_sha256']
assert p['failing_paths']==[0] and p['all_selected_path_checks_pass']
assert r['models']['ordinary']['integration_refinement_pass']
assert r['models']['full']['integration_refinement_pass']
assert r['frame_check_pass'] and r['actual_trajectory_field_gate_pass']
deps=list(HERE.glob('*.py'))+[ROOT/'research_work/results/bar-field-foundation/field.py',
    ROOT/'research_work/results/rotating-bar-orbits/fast_multipole.py']
out=dict(current_input_hashes_verified=len(r['input_hashes']),original_axis_failure_preserved=not r['models']['axisymmetric_extra']['integration_refinement_pass'],
         repaired_probe=0,all_declared_checks_pass_after_repair=True,
         dependency_hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(deps)},
         observational_fit_completed=False)
(HERE/'integrity.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))
