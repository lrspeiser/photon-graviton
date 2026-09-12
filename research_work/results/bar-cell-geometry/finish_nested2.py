"""Wait for specified existing Windows preparation processes, then run dependent audits.
No preparation is restarted. Incomplete or failed preparations stop the pipeline.
"""
from pathlib import Path
import subprocess,sys,json,concurrent.futures
H=Path(__file__).resolve().parent
pids=[int(x) for x in sys.argv[1:]]
assert len(pids)==2 and all(x>0 for x in pids)
command='Wait-Process -Id '+','.join(str(x) for x in pids)+' -ErrorAction SilentlyContinue'
subprocess.run(['powershell','-NoProfile','-Command',command],check=False)
for R in (1,3):
 d=json.loads((H/f'prepared2-R{R}.json').read_text(encoding='utf8'))
 if len(d['records'])!=d['expected'] or not all(r['passes'] for r in d['records']):
  raise RuntimeError(f'Preparation R={R} incomplete or failed; dependent work stopped')
subprocess.run([sys.executable,str(H/'check_nested2_geometry.py')],check=True)
def calculate(job):
 layers,R=job
 logfile=H/f'nested2-force-t{layers}-R{R}.log'
 with logfile.open('w',encoding='utf8',newline='\n') as f:
  p=subprocess.run([sys.executable,str(H/'nested2_volumes.py'),str(layers),str(R)],stdout=f,stderr=subprocess.STDOUT)
 print('Volume run',layers,R,'exit',p.returncode,flush=True)
 return p.returncode
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 statuses=list(pool.map(calculate,[(256,1),(256,3),(512,1),(512,3)]))
if any(statuses):raise RuntimeError('At least one volume job failed; inspect retained logs')
subprocess.run([sys.executable,str(H/'export_nested2.py')],check=True)
print('Primary nested2 geometry and force comparisons complete; physical validation remains open.',flush=True)
