"""Conditional snapshot supply from cached DustPedia; no cluster mass fit."""
from pathlib import Path
import json,hashlib,ast
import numpy as np
from scipy.special import roots_legendre
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
# Reuse only definitions/imports, without executing the external solver's driver.
source=ROOT/'research_work/results/cluster-external-photon-supply/run.py'
tree=ast.parse(source.read_text()); nodes=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom,ast.FunctionDef))]
ns={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),ns)
a=.0002488993265191759; radius=1.; k=10.

def read(file):
    good=[];bad=[]
    for line in file.read_text().splitlines():
        name=line[:23].strip()
        try:
            ra,dec,dist,lum=map(float,[line[24:33],line[34:43],line[49:61],line[100:109]])
            assert np.all(np.isfinite([ra,dec,dist,lum])) and dist>0 and lum>0
        except (ValueError,AssertionError):bad.append(name);continue
        ra,dec=np.deg2rad([ra,dec]);xyz=dist*np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)])
        good.append(dict(name=name,xyz=xyz,lum=lum))
    assert len({v['name'] for v in good})==len(good)
    return good,bad

def internal(d,n):
    mu,w=roots_legendre(n)
    length=-d*mu+np.sqrt(1-d*d+d*d*mu*mu)
    photon=np.exp(-a*length);comp=a*(photon-np.exp(-k*length))/(k-a)
    return float(np.dot(w,1-photon-comp)/2)

runs=[]
for model in ['themis','dl14']:
    file=ROOT/f'companion_causal_test/data/{model}.dat';stars,bad=read(file)
    center=next(v for v in stars if v['name']=='NGC4486')['xyz']
    values=[]
    for v in stars:
        d=float(np.linalg.norm(v['xyz']-center))/radius
        if d<1:
            f=internal(d,128); refined=internal(d,256);kind='internal'
        else:
            _,parts=ns['source'](d,a,k,128);_,fine=ns['source'](d,a,k,256)
            f=float(parts[2]);refined=float(fine[2]);kind='external'
        assert abs(f-refined)<1e-9
        values.append(dict(name=v['name'],separation_Mpc=d*radius,source_luminosity_Lsun=v['lum'],
            kind=kind,deposit_fraction=refined,deposited_power_Lsun=v['lum']*refined))
    values.sort(key=lambda x:x['deposited_power_Lsun'],reverse=True)
    runs.append(dict(model=model,input_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        usable_sources=len(stars),excluded_sources=bad,center_xyz_Mpc=center.tolist(),
        internal_count=sum(v['kind']=='internal' for v in values),
        internal_power_Lsun=sum(v['deposited_power_Lsun'] for v in values if v['kind']=='internal'),
        external_power_Lsun=sum(v['deposited_power_Lsun'] for v in values if v['kind']=='external'),sources=values))
result=dict(scope='Conditional established-illumination snapshot, NGC4486-centered R=1 Mpc sphere, alpha archived empirical, kappa R=10; catalog luminosities are SED inferences. No source history, completeness correction, universe age, cluster mass or lensing fit.',kernel_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),runs=runs)
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
for r in runs:print(json.dumps({k:v for k,v in r.items() if k not in ['sources','excluded_sources']}));print(json.dumps(r['sources'][:5]))
