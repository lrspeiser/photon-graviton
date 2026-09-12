"""Mass-preserving affine source-age cells, with exact tetrahedral gravity."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.spatial import ConvexHull
from scipy.interpolate import CubicHermiteSpline
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'bar-volume-gravity'))
from tetra import gravity
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
POINTS=np.array([[.1,0,0],[1,0,0],[1,0,1],[0,1,1],[3,0,0]])
SYMMETRIES=np.array([[1,1,1],[-1,-1,1],[1,1,-1],[-1,-1,-1]])
layers=int(sys.argv[1]);assert layers in (256,512,1024)
launch_R=float(sys.argv[2]);assert launch_R in (1.,3.)
n='nested1'
mesh=json.loads((HERE/'mesh.json').read_text(encoding='utf8'))
d=json.loads((HERE/f'prepared-R{launch_R:g}.json').read_text(encoding='utf8'))
assert len(d['records'])==d['expected']==len(mesh['directions'])
assert all(r['passes'] for r in d['records'])
assert d['mesh_hash']==hashlib.sha256((HERE/'mesh.json').read_bytes()).hexdigest()
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
    selected_times=(float(sys.argv[3]),) if len(sys.argv)>3 else (.05,.1,.25)
    assert all(T in (.05,.1,.25) for T in selected_times)
    for T in selected_times:
        ages=np.linspace(0,T,layers+1)
        x=np.stack([s(ages) for s in splines],axis=1)[:,lookup,:]*signs[None,:,:]
        potential=np.zeros(len(targets));force=np.zeros((len(targets),3));mass_sum=0.;minimum_ratio=np.inf;thinnest=None
        for j in range(layers):
            lo=x[j,faces];hi=x[j+1,faces]
            # Sorted shared vertices give conforming diagonals on adjacent prism faces.
            cells=np.concatenate([np.stack([lo[:,0],lo[:,1],lo[:,2],hi[:,2]],axis=1),
                                  np.stack([lo[:,0],lo[:,1],hi[:,1],hi[:,2]],axis=1),
                                  np.stack([lo[:,0],hi[:,0],hi[:,1],hi[:,2]],axis=1)])
            mass=np.tile(weights/(3*layers),3);mass_sum+=sum(mass)
            scale=np.max(np.linalg.norm(cells-cells[:,:1],axis=2),axis=1)
            volume=abs(np.einsum('ij,ij->i',cells[:,1]-cells[:,0],np.cross(cells[:,2]-cells[:,0],cells[:,3]-cells[:,0])))/6
            ratios=volume/scale**3;index=int(np.argmin(ratios))
            if ratios[index]<minimum_ratio:
                minimum_ratio=float(ratios[index]);thinnest=dict(vertices=cells[index].tolist(),mass=float(mass[index]))
            p,a=gravity(cells,mass,targets);potential+=p;force+=a
        assert abs(mass_sum-1)<1e-12
        p=potential[inverse].reshape(len(POINTS),4).mean(axis=1)
        a=(force[inverse].reshape(len(POINTS),4,3)*SYMMETRIES).mean(axis=1)
        result=dict(R=R,T=T,layers=layers,source_n=n,triangles=len(faces),mass=mass_sum,
                    angular_source_norm=float(sum(raw)),minimum_volume_scale_ratio=minimum_ratio,thinnest_cell=thinnest,
                    potential=p.tolist(),acceleration=a.tolist())
        out.append(result);print(n,layers,R,T,'mass',mass_sum,'potential',p.tolist(),flush=True)
        (HERE/f'volume-t{layers}-R{R:g}.json').write_text(json.dumps(dict(points=POINTS.tolist(),records=out),indent=2)+'\n',encoding='utf8',newline='\n')
