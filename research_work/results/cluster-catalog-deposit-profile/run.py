"""Radial cumulative capture from the frozen catalog supply pilot."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.special import roots_legendre
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
INPUT=ROOT/'research_work/results/cluster-catalog-pilot/results.json'
pilot=json.loads(INPUT.read_text()); a=.0002488993265191759;k=10.
radii=np.linspace(0,1,21)

def cumulative(s,p0,c0):
    ep=np.exp(-a*s);ec=np.exp(-k*s)
    return c0*(-np.expm1(-k*s))+p0*(1-ep-a*(ep-ec)/(k-a))

def profile(source,n):
    x,w=roots_legendre(n);d=source['separation_Mpc']
    if source['kind']=='external':
        y=(x+1)/2;weight=w/2*y/(2*d*np.sqrt(d*d-1+y*y))
        entry=np.sqrt(d*d-1+y*y)-y
        p0=np.exp(-a*entry);c0=-np.expm1(-a*entry)
        closest=y;impact2=1-y*y;length=2*y
    else:
        weight=w/2;mu=x;closest=-d*mu;impact2=d*d*(1-mu*mu)
        length=closest+np.sqrt(1-impact2);p0=1.;c0=0.
    output=[]
    for radius in radii:
        half=np.sqrt(np.maximum(radius*radius-impact2,0))
        lo=np.clip(closest-half,0,length);hi=np.clip(closest+half,0,length)
        deposited=np.where(radius*radius>=impact2,cumulative(hi,p0,c0)-cumulative(lo,p0,c0),0)
        output.append(float(weight@deposited)*source['source_luminosity_Lsun'])
    return np.array(output)

runs=[]
for row in pilot['runs']:
    coarse={key:np.zeros(len(radii)) for key in ['internal','external']}
    fine={key:np.zeros(len(radii)) for key in coarse}
    for src in row['sources']:
        coarse[src['kind']]+=profile(src,256)
        fine[src['kind']]+=profile(src,512)
    total=fine['internal']+fine['external'];power=row['internal_power_Lsun']+row['external_power_Lsun']
    discrepancy=float(np.max(np.abs(total-coarse['internal']-coarse['external']))/power)
    assert discrepancy<.002
    assert abs(total[-1]/power-1)<1e-8
    assert np.all(np.diff(total)>=-1e-5) and total[0]==0
    runs.append(dict(model=row['model'],radii_Mpc=radii.tolist(),internal_enclosed_power_Lsun=fine['internal'].tolist(),
        external_enclosed_power_Lsun=fine['external'].tolist(),total_enclosed_power_Lsun=total.tolist(),
        enclosed_fraction=(total/power).tolist(),outer_half_fraction=1-total[10]/power,
        internal_outer_half_fraction=1-fine['internal'][10]/fine['internal'][-1],
        external_outer_half_fraction=1-fine['external'][10]/fine['external'][-1],
        quadrature_max_difference_over_total=discrepancy,total_recovery_relative_error=abs(total[-1]/power-1)))
(OUT/'results.json').write_text(json.dumps(dict(input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
    scope='Radial cumulative deposition POWER of same fixed catalog, alpha, kappa and radius; not a full spatial density, accumulated mass or lensing map.',runs=runs),indent=2)+'\n',encoding='utf-8')
for r in runs: print(json.dumps({key:value for key,value in r.items() if not isinstance(value,list)}))
