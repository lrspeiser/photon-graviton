from __future__ import annotations
import argparse,json,time
from pathlib import Path
from dataclasses import asdict
import numpy as np
from atlas import Grid,Mesh,Case
ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--max-nr',type=int,default=16);args=ap.parse_args()
p=args.results;out=p/'refinement'
if out.exists():raise FileExistsError('refinement already exists')
out.mkdir();primary=[json.loads(line) for line in (p/'primary.jsonl').read_text().splitlines()]
ids=json.loads((p/'extended-selection-before-refinement.json').read_text())['ids'];rows=[];t=time.monotonic()
for mesh in [Mesh(),Mesh(12,16,48),Mesh(16,20,64)]:
 if mesh.nr>args.max_nr:continue
 g=Grid(mesh)
 for i in ids:
  c=Case(**primary[i]['case']);d,f,A=g.solve(c);r=dict(id=i,mesh=asdict(mesh),case=asdict(c),diagnostics=d,readout=g.readout(f));rows.append(r)
  (out/f'case_{i}_{mesh.nr}.json').write_text(json.dumps(r,indent=2))
  if mesh.nr==args.max_nr:np.savez_compressed(out/f'fields_{i}.npz',P=f[0],C=f[1],gas=f[2],J=f[3],weight=g.w,R=g.R,Z=g.Z,phi=g.F)
  print('REFINE',i,mesh.nr,'Ec',d['Ec'],'seconds',round(time.monotonic()-t,2),flush=True)
(out/'all.json').write_text(json.dumps(rows,indent=2));print('DONE',time.monotonic()-t,flush=True)
