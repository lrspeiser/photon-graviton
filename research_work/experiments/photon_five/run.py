"""Run an immutable PF5 stage. Scientific gate failures are retained, not exceptions."""
import argparse
from datetime import datetime, timezone
import importlib
from pathlib import Path
import platform
import subprocess
import time
import traceback
from common import HERE, ROOT, save, provenance

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['radiation','transport','fluctuations','shock','strain'],required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise ValueError('Evidence directory already exists; choose a new name')
    args.output.mkdir(parents=True)
    start=time.monotonic()
    metadata=dict(stage=args.stage,started_utc=datetime.now(timezone.utc).isoformat(),
                  commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  python=platform.python_version())
    try:
        result=importlib.import_module(args.stage).run()
        save(args.output/'results.json',result)
        print({k:v for k,v in result.items() if k in ['experiment','gates','numerical_pass','mechanism_pass','observational_status']},flush=True)
    except Exception:
        metadata['exception']=traceback.format_exc()
        raise
    finally:
        metadata['elapsed_seconds']=time.monotonic()-start
        metadata['sha256']=provenance()
        save(args.output/'manifest.json',metadata)

if __name__=='__main__': main()
