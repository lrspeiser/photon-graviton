"""Reuse the recorded inner-bin fitting algorithm with one explicit mass-profile change."""
from pathlib import Path
import json
import sys

folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder))
template = folder.parent/'slacs-outer-bin-check/run.py'
original = template.read_text(encoding='utf-8')
for eta in json.loads((folder/'protocol.json').read_text(encoding='utf-8'))['eta_values']:
    tag = str(eta).replace('-', 'minus').replace('.', 'p')
    source = original
    replacements = {
        'from model import ComponentModel':
            f'from gradient_model import GradientModel\nfrom functools import partial\nComponentModel=partial(GradientModel, eta={eta!r})',
        "HERE/'results.json'": f"HERE/'eta-{tag}.json'",
        "paths+[HERE/'run.py',":
            "paths+[HERE/'run.py', HERE/'gradient_model.py', HERE.parent/'slacs-outer-bin-check/run.py',"
    }
    for old, new in replacements.items():
        assert source.count(old) == 1, old
        source = source.replace(old, new)
    exec(compile(source, str(template), 'exec'), {'__file__':str(folder/'run.py'), '__name__':'__main__'})
    print(json.dumps(dict(eta_completed=eta)), flush=True)
