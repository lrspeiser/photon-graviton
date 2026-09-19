"""Run one predeclared CWC-1 stage into a new immutable directory."""
import argparse
from datetime import datetime,timezone
import importlib
import importlib.metadata
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
from common import HERE,ROOT,hashes,save

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--stage',choices=['small','one_dimensional','two_dimensional','readiness'],required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise ValueError('Output already exists; preserve previous evidence')
    a.output.mkdir(parents=True)
    source=sorted(list(HERE.glob('*.py'))+list(HERE.glob('*.md')))
    metadata=dict(stage=a.stage,started_utc=datetime.now(timezone.utc).isoformat(),
                  commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  source_sha256=hashes(source),python=platform.python_version(),
                  packages={k:importlib.metadata.version(k) for k in ['numpy','scipy']})
    save(a.output/'manifest.json',metadata)
    start=time.monotonic()
    try:
        result=importlib.import_module(a.stage).run(a.output) if a.stage!='small' else importlib.import_module(a.stage).run()
        save(a.output/'results.json',result)
        print({k:v for k,v in result.items() if k in ['stage','numerical_pass','mechanism_pass','observational_ready','gates']},flush=True)
    except Exception:
        metadata['exception']=traceback.format_exc()
        raise
    finally:
        metadata['elapsed_seconds']=time.monotonic()-start
        metadata['source_changed_during_run']=metadata['source_sha256']!=hashes(source)
        metadata['excluded_modules_loaded']=sorted(set(sys.modules)&{'cl1','cl2','cl2_sources','astropy.cosmology'})
        save(a.output/'manifest.json',metadata)
        if metadata['source_changed_during_run'] or metadata['excluded_modules_loaded']:
            raise RuntimeError('Run provenance/source-scope guard failed')

if __name__=='__main__':main()
