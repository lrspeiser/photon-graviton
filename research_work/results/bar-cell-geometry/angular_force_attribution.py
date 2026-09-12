"""Mass-preserving affine source-age cells, with exact tetrahedral gravity."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.spatial import ConvexHull
from scipy.interpolate import CubicHermiteSpline
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'bar-volume-gravity'))
from tetra import gravity
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
POINTS=np.array([[.1,0,0]])
SYMMETRIES=np.array([[1,1,1],[-1,-1,1],[1,1,-1],[-1,-1,-1]])
level=int(sys.argv[1]);assert level in (1,2)
layers=512;launch_R=3.;n=f'nested{level}'
meshfile=HERE/('mesh.json' if level==1 else 'mesh2.json')
prepfile=HERE/('prepared-R3.json' if level==1 else 'prepared2-R3.json')
mesh=json.loads(meshfile.read_text(encoding='utf8'));d=json.loads(prepfile.read_text(encoding='utf8'))
assert len(d['records'])==d['expected']==len(mesh['directions']) and all(r['passes'] for r in d['records'])
assert d['mesh_hash']==hashlib.sha256(meshfile.read_bytes()).hexdigest()
out=[]
for R in (launch_R,):
    group=[r for r in d['records'] if r['R']==R]
    captures=[];splines=[]
    for r in group:
        p=ROOT/r['cache'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['cache_sha256']
        c=np.load(p);splines.append(CubicHermiteSpline(c['ages'],c['position'],c['velocity']))
        captures.append(r['capture'])
    directions=np.array(mesh['vertices']);lookup=np.array(mesh['lookup']);signs=np.array(mesh['signs'])
    captures=np.array(captures)[lookup]
    faces=np.sort(np.array(mesh['faces']),axis=1)
    va,vb,vc=np.moveaxis(directions[faces],1,0)
    omega=2*np.arctan2(abs(np.einsum('ij,ij->i',va,np.cross(vb,vc))),1+np.sum(va*vb+vb*vc+vc*va,axis=1))
    assert abs(sum(omega)-4*np.pi)<1e-10
    raw=omega*np.mean(captures[faces],axis=1)/(4*np.pi);weights=raw/sum(raw)
    # Explicit symmetrization removes the arbitrary diagonal choices of coplanar hull faces.
    targets,inverse=np.unique((POINTS[:,None,:]*SYMMETRIES).reshape(-1,3),axis=0,return_inverse=True)
    T=.1;ages=np.linspace(0,T,layers+1)
    x=np.stack([s(ages) for s in splines],axis=1)[:,lookup,:]*signs[None,:,:]
    potential=np.zeros((len(targets),len(faces)));force=np.zeros((len(targets),len(faces),3))
    for j in range(320,384):
        lo=x[j,faces];hi=x[j+1,faces]
        cells=np.concatenate([np.stack([lo[:,0],lo[:,1],lo[:,2],hi[:,2]],axis=1),
                              np.stack([lo[:,0],lo[:,1],hi[:,1],hi[:,2]],axis=1),
                              np.stack([lo[:,0],hi[:,0],hi[:,1],hi[:,2]],axis=1)])
        mass=np.tile(weights/(3*layers),3)
        p,a=gravity(cells,mass,targets,per_cell=True)
        potential+=p.reshape(len(targets),3,len(faces)).sum(axis=1)
        force+=a.reshape(len(targets),3,len(faces),3).sum(axis=1)
    p=potential[inverse].mean(axis=0)
    a=(force[inverse]*SYMMETRIES[:,None,:]).mean(axis=0)
    masses=weights/8
    if level==2:
        p=p.reshape(-1,4).sum(axis=1);a=a.reshape(-1,4,3).sum(axis=1);masses=masses.reshape(-1,4).sum(axis=1)
    ref=json.loads((HERE/f'age-attribution-n{level}.json').read_text(encoding='utf8'))['records'][0]['bins'][5]
    ep=abs(p.sum()-ref['potential']);ea=np.linalg.norm(a.sum(axis=0)-ref['acceleration'])
    assert max(ep,ea)<1e-10 and abs(masses.sum()-.125)<1e-12
    rows=[dict(parent_face=i,mass=float(masses[i]),potential=float(p[i]),acceleration=a[i].tolist()) for i in range(len(p))]
    result=dict(level=level,age_interval=[.0625,.075],R=3,T=.1,point=POINTS[0].tolist(),records=rows,
                potential_reconstruction_error=float(ep),force_reconstruction_error=float(ea))
    (HERE/f'angular-attribution-n{level}.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
    print(level,ep,ea,flush=True)
