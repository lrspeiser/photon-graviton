"""Immutable runner for the declared source-symmetry and observational-readiness audits."""
import argparse,importlib,importlib.metadata,platform,subprocess,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path
from common import HERE,ROOT,hashes,save

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['mirror_audit','readiness'],required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise ValueError('Preserve previous evidence; use new directory')
    a.output.mkdir(parents=True)
    source=sorted(list(HERE.glob('*.py'))+list(HERE.glob('*.md')))
    meta=dict(stage=a.stage,started_utc=datetime.now(timezone.utc).isoformat(),
              commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              source_sha256=hashes(source),python=platform.python_version(),
              packages={k:importlib.metadata.version(k) for k in ['numpy','scipy']})
    save(a.output/'manifest.json',meta);start=time.monotonic()
    try:
        result=importlib.import_module(a.stage).run(a.output);save(a.output/'results.json',result)
        print({k:v for k,v in result.items() if k in ['stage','gates','numerical_pass','mechanism_pass']},flush=True)
    except Exception:
        meta['exception']=traceback.format_exc();raise
    finally:
        meta['elapsed_seconds']=time.monotonic()-start
        meta['source_changed_during_run']=meta['source_sha256']!=hashes(source)
        meta['excluded_modules_loaded']=sorted(set(sys.modules)&{'cl1','cl2','cl2_sources','astropy.cosmology'})
        save(a.output/'manifest.json',meta)
        if meta['source_changed_during_run'] or meta['excluded_modules_loaded']:raise RuntimeError('Provenance guard failed')
if __name__=='__main__':main()
