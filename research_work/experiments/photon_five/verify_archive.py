"""Read-only archive verification; optional JSON output for publication."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FINAL=['e1-v2','e2-v2','e3-v1','e3-colored-v1','e4-v1','e5-v1','audit-v2']

def git(*args):
    return subprocess.run(['git',*args],cwd=ROOT,capture_output=True,check=False)

def blob_matches(ref,path,expected):
    result=git('show',ref+':'+path)
    return result.returncode==0 and hashlib.sha256(result.stdout).hexdigest()==expected

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    ast_count=0
    for file in HERE.glob('*.py'):
        ast.parse(file.read_text(encoding='utf-8'));ast_count+=1
    resolved=[];unresolved=[];later=[]
    for file in sorted((HERE/'evidence').glob('*/manifest.json')):
        manifest=json.loads(file.read_text(encoding='utf-8'))
        assert not manifest.get('excluded_modules_loaded',[])
        for path,expected in manifest['sha256'].items():
            path=path.replace('\\','/')
            commit=manifest['commit']
            if blob_matches(commit,path,expected):
                resolved.append((file.parent.name,path,commit));continue
            commits=git('log','--format=%H','--all','--',path).stdout.decode().splitlines()
            found=next((c for c in commits if blob_matches(c,path,expected)),None)
            if found:
                resolved.append((file.parent.name,path,found))
                later.append(dict(run=file.parent.name,path=path,resolved_commit=found,
                                  explanation='Ancillary file snapshot differs from starting commit; hash resolves in history.'))
            else:unresolved.append(dict(run=file.parent.name,path=path,sha256=expected))
    statuses=[]
    for name in FINAL:
        d=json.loads((HERE/'evidence'/name/'results.json').read_text())
        assert all(type(v) is bool for v in d['gates'].values())
        assert d['numerical_pass']==all(d['gates'].values())
        assert d['numerical_pass']
        if 'mechanism_pass' in d:assert d['mechanism_pass'] is False
        statuses.append(dict(run=name,numerical_pass=d['numerical_pass'],gates=len(d['gates'])))
    old=json.loads((HERE/'evidence/e1-v1/results.json').read_text())
    assert old['numerical_pass'] is False
    failed=json.loads((HERE/'evidence/e2-v1/manifest.json').read_text())
    assert 'exception' in failed
    assert not unresolved,unresolved
    result=dict(verified=True,python_modules_parsed=ast_count,manifest_entries_resolved=len(resolved),
                ancillary_history_resolutions=later,unresolved=unresolved,statuses=statuses,
                original_failure_records_retained=True,
                scope='Archive/AST/gate consistency and exact source/input hash recovery; scientific evidence is not rerun.')
    text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text,encoding='utf-8',newline='\n')
    print(text)

if __name__=='__main__':main()
