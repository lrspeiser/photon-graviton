"""Repeat the unchanged local calculation on a finer grid, retaining the first run.

The original run and this driver are both hashed in the refinement output.
No physical assumption, prior, likelihood, or numerical gate is changed.
"""
from pathlib import Path
source = (Path(__file__).resolve().parent/'run.py').read_text(encoding='utf-8')
replacements = {
    "nm, nb = parproto['grids'][-1]": "nm, nb = 6401, 1281",
    "HERE/'results.json'": "HERE/'refined-results.json'",
    "paths+[HERE/'run.py',": "paths+[HERE/'run.py', HERE/'refine.py',",
    "print(json.dumps(rows[-1]), flush=True)":
        "print(json.dumps(dict(Name=name, passed=sum(r['numerical_pass'] for r in predictions))), flush=True)"
}
for old, new in replacements.items():
    assert source.count(old) == 1, old
    source = source.replace(old, new)
exec(compile(source, str(Path(__file__).resolve().parent/'run.py'), 'exec'))
