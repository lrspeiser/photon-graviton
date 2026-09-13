"""Directional moments of fixed catalog-fed deposit power, not lens peaks."""
from pathlib import Path
import json,ast,hashlib
import numpy as np
from scipy.special import roots_legendre
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
reader=ROOT/'research_work/results/cluster-catalog-pilot/run.py'
tree=ast.parse(reader.read_text());ns={}
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom)) or isinstance(n,ast.FunctionDef) and n.name=='read'],type_ignores=[]),str(reader),'exec'),ns)
a=.0002488993265191759;k=10.

def moments(d,n):
    x,w=roots_legendre(n); v,wv=roots_legendre(n)
    if d>=1:
        y=(x+1)/2;wa=w/2*y/(2*d*np.sqrt(d*d-1+y*y))
        entry=np.sqrt(d*d-1+y*y);mu=entry/d;entry=entry-y;L=2*y
        p0=np.exp(-a*entry);c0=-np.expm1(-a*entry)
    else:
        mu=x;wa=w/2;entry=np.zeros_like(mu)
        L=-d*mu+np.sqrt(1-d*d*(1-mu*mu));p0=np.ones_like(mu);c0=np.zeros_like(mu)
    s=L[:,None]*(v+1)/2
    comp=c0[:,None]*np.exp(-k*s)+p0[:,None]*a*(np.exp(-a*s)-np.exp(-k*s))/(k-a)
    weights=wa[:,None]*wv[None,:]*L[:,None]/2*k*comp
    travel=entry[:,None]+s
    z=d-travel*mu[:,None] if d>=1 else d+s*mu[:,None]
    transverse2=travel**2*(1-mu[:,None]**2)
    return np.array([np.sum(weights),np.sum(weights*z),np.sum(weights*z*z),np.sum(weights*transverse2)/2])

def calculate(stars,center,n):
    total=0.;first=np.zeros(3);second=np.zeros((3,3))
    for src in stars:
        xyz=src['xyz']-center;d=float(np.linalg.norm(xyz));u=xyz/d if d else np.array([0.,0.,1.])
        m=moments(d,n)*src['lum'];total+=m[0];first+=m[1]*u
        second+=m[3]*np.eye(3)+(m[2]-m[3])*np.outer(u,u)
    centroid=first/total;cov=second/total-np.outer(centroid,centroid)
    assert np.min(np.linalg.eigvalsh(cov))>=-1e-12
    return total,centroid,cov

central=moments(0.,128)
assert abs(central[1])<1e-14 and abs(central[2]-central[3])<1e-14
runs=[]
pilot=json.loads((ROOT/'research_work/results/cluster-catalog-pilot/results.json').read_text())
for model in ['themis','dl14']:
    file=ROOT/f'companion_causal_test/data/{model}.dat';stars,_=ns['read'](file)
    center=next(x['xyz'] for x in stars if x['name']=='NGC4486')
    coarse=calculate(stars,center,64);power,centroid,cov=calculate(stars,center,128)
    expected=next(r for r in pilot['runs'] if r['model']==model)
    ref=expected['internal_power_Lsun']+expected['external_power_Lsun']
    assert abs(power/ref-1)<1e-7
    sight=center/np.linalg.norm(center);east=np.array([-sight[1],sight[0],0]);east/=np.linalg.norm(east)
    north=np.cross(sight,east);basis=np.array([east,north])
    pc=basis@centroid;pv=basis@cov@basis.T;ev=np.linalg.eigvalsh(pv)
    delta=float(np.linalg.norm(centroid-coarse[1]));assert delta<1e-5
    runs.append(dict(model=model,catalog_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        deposited_power_Lsun=power,centroid_xyz_Mpc=centroid.tolist(),centroid_offset_Mpc=float(np.linalg.norm(centroid)),
        covariance_Mpc2=cov.tolist(),projected_centroid_east_north_Mpc=pc.tolist(),
        projected_offset_Mpc=float(np.linalg.norm(pc)),projected_rms_minor_major_Mpc=np.sqrt(ev).tolist(),
        projected_second_moment_axis_ratio=float(np.sqrt(ev[0]/ev[1])),
        centroid_refinement_Mpc=delta,covariance_refinement_max_Mpc2=float(np.max(np.abs(cov-coarse[2]))),
        total_power_relative_error=abs(power/ref-1)))
(OUT/'results.json').write_text(json.dumps(dict(scope='Spatial moments of deposition power in stipulated static 1-Mpc receiver; orthographic sky projection, not a lensing peak or accumulated reservoir.',runs=runs),indent=2)+'\n',encoding='utf-8')
print(json.dumps(runs,indent=2))
