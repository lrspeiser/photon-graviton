"""Compare selected actual thin cells using 60-digit evaluation of the formula."""
from pathlib import Path
import json,sys
import numpy as np
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'bar-volume-gravity'))
from tetra import gravity
HERE=Path(__file__).resolve().parent;mp.mp.dps=60
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def sub(a,b):return [x-y for x,y in zip(a,b)]
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def norm(a):return mp.sqrt(dot(a,a))
def evaluate(vertices,mass,point):
    v=[[mp.mpf(x) for x in r] for r in vertices];point=[mp.mpf(x) for x in point]
    volume=abs(dot(sub(v[1],v[0]),cross(sub(v[2],v[0]),sub(v[3],v[0]))))/6
    rho=mp.mpf(mass)/volume;potential=mp.mpf(0);force=[mp.mpf(0)]*3
    for opposite in range(4):
        face=[v[i] for i in range(4) if i!=opposite]
        n=cross(sub(face[1],face[0]),sub(face[2],face[0]))
        if dot(n,sub(v[opposite],face[0]))>0:face[1],face[2]=face[2],face[1]
        n=cross(sub(face[1],face[0]),sub(face[2],face[0]));n=[x/norm(n) for x in n]
        r=[sub(x,point) for x in face];length=[norm(x) for x in r];d=dot(n,r[0])
        num=dot(r[0],cross(r[1],r[2]));den=mp.fprod(length)
        for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):den+=dot(r[i],r[j])*length[k]
        integral=-d*2*mp.atan2(num,den)
        for i,j in ((0,1),(1,2),(2,0)):
            edge=sub(face[j],face[i]);ell=norm(edge);outward=cross([x/ell for x in edge],n)
            h=dot(outward,r[i]);dd=length[i]+length[j]-ell
            integral+=h*mp.log1p(2*ell/dd)
        potential-=rho*d*integral/2
        force=[x-rho*y*integral for x,y in zip(force,n)]
    return float(potential),np.array([float(x) for x in force])

rows=[]
for path in sorted(HERE.glob('volume-adaptive-t*.json')):
    d=json.loads(path.read_text(encoding='utf8'))
    if len(d['records']) not in (1,3):continue
    for record in d['records']:
        cell=record.get('thinnest_cell')
        if cell is None:continue
        fp,fa=gravity([cell['vertices']],[cell['mass']],d['points'])
        for i,point in enumerate(d['points']):
            p,a=evaluate(cell['vertices'],cell['mass'],point)
            ep=abs(fp[i]-p);ea=float(np.linalg.norm(fa[i]-a))
            rows.append(dict(file=path.name,R=record['R'],T=record['T'],point=point,
                             potential_absolute_error=ep,force_absolute_error=ea,passes=bool(max(ep,ea)<1e-7)))
assert rows and all(r['passes'] for r in rows)
(HERE/'roundoff-adaptive.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n')
print('Checks',len(rows),'max potential',max(r['potential_absolute_error'] for r in rows),'max force',max(r['force_absolute_error'] for r in rows))
