"""Run one declared TF-1 stage into a new immutable evidence directory (after CWC-1's run.py)."""
import argparse
from datetime import datetime, timezone
import importlib
import importlib.metadata
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
from common import HERE, ROOT, CWC, PM, hashes, save

INPUTS = [CWC/'evidence/spatial2d-v1/primary-final-field.npz', CWC/'evidence/spatial2d-v1/results.json',
          CWC/'evidence/spatial1d-v1/results.json', PM/'cl1-results.json', PM/'cl2-results.json',
          CWC/'model.py', CWC/'spatial.py', CWC/'two_dimensional.py', PM/'cl2_sources.py', PM/'cl2_response.py', PM/'steady_field.py']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=['exact', 'scan', 'lenses'], required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise ValueError('Output already exists; preserve previous evidence')
    a.output.mkdir(parents=True)
    source = sorted(list(HERE.glob('*.py')) + list(HERE.glob('*.md')))
    meta = dict(stage=a.stage, started_utc=datetime.now(timezone.utc).isoformat(),
                commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                source_sha256=hashes(source), consumed_inputs_sha256=hashes([x for x in INPUTS if x.exists()]),
                python=platform.python_version(), packages={k: importlib.metadata.version(k) for k in ['numpy', 'scipy']})
    save(a.output/'manifest.json', meta)
    start = time.monotonic()
    try:
        module = importlib.import_module(dict(exact='exact', scan='lattice', lenses='lenses')[a.stage])
        result = module.run(a.output)
        save(a.output/'results.json', result)
        print({k: v for k, v in result.items() if k in ['stage', 'numerical_pass', 'gates', 'outcome']}, flush=True)
    except Exception:
        meta['exception'] = traceback.format_exc()
        raise
    finally:
        meta['elapsed_seconds'] = time.monotonic() - start
        meta['source_changed_during_run'] = meta['source_sha256'] != hashes(source)
        save(a.output/'manifest.json', meta)
        if meta['source_changed_during_run']:
            raise RuntimeError('Source changed during the run')


if __name__ == '__main__':
    main()
