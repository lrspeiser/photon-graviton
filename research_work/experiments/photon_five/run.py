"""Run an immutable PF5 stage. Scientific gate failures are retained, not exceptions."""
import argparse
from datetime import datetime, timezone
import importlib
import importlib.metadata
import sys
from pathlib import Path
import platform
import subprocess
import time
import traceback
from common import HERE, ROOT, save, provenance

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['radiation','transport','fluctuations','fluctuations_colored','shock','strain','audit'],required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise ValueError('Evidence directory already exists; choose a new name')
    args.output.mkdir(parents=True)
    start=time.monotonic()
    metadata=dict(stage=args.stage,started_utc=datetime.now(timezone.utc).isoformat(),
                  commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  python=platform.python_version(),
                  packages={p:importlib.metadata.version(p) for p in ['numpy','scipy']})
    try:
        result=importlib.import_module(args.stage).run()
        save(args.output/'results.json',result)
        print({k:v for k,v in result.items() if k in ['experiment','gates','numerical_pass','mechanism_pass','observational_status']},flush=True)
    except Exception:
        metadata['exception']=traceback.format_exc()
        raise
    finally:
        metadata['elapsed_seconds']=time.monotonic()-start
        blocked={'cl1','cl2','cl2_sources','astropy.cosmology'}
        metadata['excluded_modules_loaded']=sorted(blocked.intersection(sys.modules))
        metadata['repository_modules']={name:str(Path(mod.__file__).resolve().relative_to(ROOT))
                                       for name,mod in list(sys.modules.items())
                                       if getattr(mod,'__file__',None) and Path(mod.__file__).resolve().is_relative_to(ROOT)}
        if metadata['excluded_modules_loaded']: raise RuntimeError('Excluded physics module was loaded')
        metadata['sha256']=provenance()
        save(args.output/'manifest.json',metadata)

if __name__=='__main__': main()
