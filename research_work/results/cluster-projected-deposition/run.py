"""Projected point-emitter deposit power; flags unresolved source cusps."""
from pathlib import Path
import json,ast,hashlib
import numpy as np
from scipy.special import roots_legendre
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
reader=ROOT/'research_work/results/cluster-catalog-pilot/run.py';tree=ast.parse(reader.read_text());ns={}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom)) or isinstance(n,ast.FunctionDef) and n.name=='read'],type_ignores=[]),str(reader),'exec'),ns)
file=ROOT/'companion_causal_test/data/themis.dat';sources,_=ns['read'](file)
center=next(s['xyz'] for s in sources if s['name']=='NGC4486');sight=center/np.linalg.norm(center)
east=np.array([-sight[1],sight[0],0]);east/=np.linalg.norm(east);north=np.cross(sight,east);basis=np.array([east,north,sight])
positions=np.array([basis@(s['xyz']-center) for s in sources]);lum=np.array([s['lum'] for s in sources]);a=.0002488993265191759;k=10.

def density(points):
    out=np.zeros(len(points))
    for src,L in zip(positions,lum):
        delta=points-src;length=np.linalg.norm(delta,axis=1)
        if np.any(length==0):raise ValueError('Point-emitter cusp; need finite source extent')
        if np.linalg.norm(src)<1: entry=np.zeros(len(points))
        else:
            proj=delta@src/length
            entry=-proj-np.sqrt(np.maximum(proj*proj-(src@src-1),0))
        path=np.maximum(length-entry,0);p0=np.exp(-a*entry);c0=-np.expm1(-a*entry)
        comp=c0*np.exp(-k*path)+p0*a*(np.exp(-a*path)-np.exp(-k*path))/(k-a)
        out+=k*L*comp/(4*np.pi*length*length)
    return out

def project(xy,n):
    z,w=roots_legendre(n);height=np.sqrt(np.maximum(1-np.sum(xy*xy,axis=1),0))
    points=np.column_stack([np.repeat(xy[:,0],n),np.repeat(xy[:,1],n),(height[:,None]*z).ravel()])
    return height*(density(points).reshape(len(xy),n)@w)

axis=np.linspace(-.975,.975,40);xx,yy=np.meshgrid(axis,axis);inside=(xx*xx+yy*yy)<1
xy=np.column_stack([xx[inside],yy[inside]])
coarse=project(xy,64);fine=project(xy,128)
assert np.all(fine>=0) and np.all(np.isfinite(fine))
internal=positions[np.linalg.norm(positions,axis=1)<1,:2]
separation=np.min(np.linalg.norm(xy[:,None,:]-internal[None,:,:],axis=2),axis=1)
safe=separation>.05
relative=np.abs(fine-coarse)/np.maximum(fine,1e-100)
center_values=[float(project(np.array([[0.,0.]]),n)[0]) for n in [32,64,128,256]]
rows=[dict(east_Mpc=float(x),north_Mpc=float(y),surface_power_Lsun_Mpc2=float(v),relative_refinement=float(e),near_projected_source=bool(d<=.05)) for (x,y),v,e,d in zip(xy,fine,relative,separation)]
(OUT/'map.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
summary=dict(input_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),scope='Point-sampled projected power, 1 Mpc fixed receiver and THEMIS sources; no finite source sizes, pixel averages, accumulation history or gravity/lens map.',samples=len(rows),off_source_samples=int(sum(safe)),
 off_source_max_relative_refinement=float(np.max(relative[safe])),all_max_relative_refinement=float(max(relative)),center_projection_by_nodes=dict(zip([32,64,128,256],center_values)),finite_peak_prediction=False)
(OUT/'results.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary,indent=2))
