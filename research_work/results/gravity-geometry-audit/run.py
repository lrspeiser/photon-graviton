"""Frozen empirical relation, two conservative geometric completions.

No fitted parameters and no new stellar holdout evaluation. Published force
summaries remain model-dependent, previously exposed diagnostic inputs.
"""
from pathlib import Path
import sys,json,hashlib,argparse
import numpy as np
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'rotating-bar-orbits'))
from run import DiskGrid,Field
parser=argparse.ArgumentParser()
parser.add_argument('--refined',action='store_true')
parser.add_argument('--bar-order',type=int,default=40,choices=[40,64])
args=parser.parse_args()
suffix=('-refined' if args.refined else '')+('-L64' if args.bar_order==64 else '')

def save(name,x):
    f=Path(name)
    (HERE/(f.stem+suffix+f.suffix)).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
inputs=HERE.parent/'milky-way-capture/inputs.json'
fitfile=HERE.parent/'joint-galaxy-audit/results.json'
par=json.loads(fitfile.read_text())['sparc']['parameters']
KPC=3.085677581491367e19;G=4.30091727003628e-6
A,p,astar=par['A'],par['p'],par['a_star_m_s2']*KPC/1e6
def extra(g):return A*astar*(g/astar)**p
disk=DiskGrid(256,129,None)
field=Field(disk,order=args.bar_order)
na=64 if args.refined else 32
angles=np.arange(na)*2*np.pi/na
def averaged(R,z):
    R,z=np.broadcast_arrays(np.atleast_1d(R),np.atleast_1d(z))
    xyz=np.stack([R[:,None]*np.cos(angles),R[:,None]*np.sin(angles),np.broadcast_to(z[:,None],(len(R),len(angles)))],axis=-1)
    pot,a=field.evaluate(xyz.reshape(-1,3));a=a.reshape(-1,len(angles),3)
    ar=(a[:,:,0]*np.cos(angles)+a[:,:,1]*np.sin(angles)).mean(axis=1)
    return pot.reshape(-1,len(angles)).mean(axis=1),np.c_[ar,a[:,:,2].mean(axis=1)]

# Angle averaging of a scalar potential remains conservative. Calibrate the
# potential composition along its equatorial curve, without observed refitting.
R=np.geomspace(.2,32,960 if args.refined else 480)
phi,acc=averaged(R,0)
gb=-acc[:,0]
assert np.all(gb>0) and np.all(np.diff(phi)>0)
radial=PchipInterpolator(R,extra(gb),extrapolate=False)
ratio=PchipInterpolator(phi,extra(gb)/gb,extrapolate=False)
radial_potential=radial.antiderivative()
composition_potential=ratio.antiderivative()
def models(R,z):
    phi,a=averaged(R,z);R,z=np.broadcast_arrays(np.atleast_1d(R),np.atleast_1d(z));r=np.hypot(R,z)
    shell=-radial(r)[:,None]*np.c_[R/r,z/r]
    comp=ratio(phi)[:,None]*a
    direct=(extra(np.linalg.norm(a,axis=1))/np.linalg.norm(a,axis=1))[:,None]*a
    return {'baryons':a,'spherical':a+shell,'potential_composition':a+comp,'direct_multiplier':a+direct},phi

data=json.loads(inputs.read_text());rows=[];scores={}
for obs,src in [('rotation',data['eilers']['rows']),('vertical',data['bovy']['rows'])]:
    rr=np.array([x['R_kpc'] for x in src]);zz=np.zeros(len(rr)) if obs=='rotation' else np.full(len(rr),1.1)
    forces,_=models(rr,zz)
    for name,a in forces.items():
        pred=np.sqrt(-rr*a[:,0]) if obs=='rotation' else np.abs(a[:,1])/(2*np.pi*G*1e6)
        seen=np.array([x['vc_kms'] if obs=='rotation' else x['Kz_over_2piG_Msun_pc2'] for x in src])
        assert np.isfinite(pred).all()
        scores[f'{obs}/{name}']={'rows':len(rr),'rmse':float(np.sqrt(np.mean((pred-seen)**2))),'bias':float(np.mean(pred-seen))}
        for i in range(len(rr)):rows.append(dict(observable=obs,model=name,R_kpc=float(rr[i]),z_kpc=float(zz[i]),observed=float(seen[i]),predicted=float(pred[i])))

# Work around rectangular meridional loops: a static scalar well gives zero.
# This is a constrained diagnostic path, not an assumed freely falling orbit.
def work(bounds,n):
    r0,r1,z0,z1=bounds;q,w=leggauss(n);out={k:0. for k in ['baryons','spherical','potential_composition','direct_multiplier']}
    for start,end in [(np.array([r0,z0]),np.array([r1,z0])),(np.array([r1,z0]),np.array([r1,z1])),(np.array([r1,z1]),np.array([r0,z1])),(np.array([r0,z1]),np.array([r0,z0]))]:
        pts=(start+end)/2+q[:,None]*(end-start)/2
        fs,_=models(pts[:,0],pts[:,1])
        for name,a in fs.items():out[name]+=float(np.sum(w*(a@(end-start)/2)))
    return out
loops=[]
for bounds in [(1,3,.1,1.5),(3,5,.1,1.5),(5,9,.1,1.5),(9,15,.1,1.5)]:
    low=work(bounds,64 if args.refined else 32);high=work(bounds,128 if args.refined else 64)
    loops.append(dict(bounds_R0_R1_z0_z1_kpc=bounds,work_kms2=high,quadrature_change_kms2={k:abs(high[k]-low[k]) for k in high}))
    print('loop',bounds,high,flush=True)
save('predictions.json',rows)
sources=[inputs,fitfile,disk.path,ROOT/f'research_work/data-cache/bar-field/bar-L{args.bar_order}.npz',ROOT/'research_work/data-cache/bar-field/nuclei-L16.npz',HERE.parent/'bar-field-foundation/field.py',HERE.parent/'rotating-bar-orbits/run.py',HERE.parent/'rotating-bar-orbits/fast_multipole.py']
save('results.json',dict(parameters=par,bar_order=args.bar_order,angles=na,radial_nodes=len(R),field='Azimuth average of existing ordinary-matter bar, L16 nuclei, conservative disk grid, central mass; no halo',
    scores=scores,loops=loops,input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources},
    holdout_scores_opened=False,new_parameters_fitted=False,
    caveats=['Previously exposed, model-dependent summary inputs; not independent raw-force measurements.',
             'Completions are geometry assumptions, not a derived companion mechanism.',
             'Direct multiplier is diagnostic only and is not admitted as a conservative static well.',
             'Finite radial domain 0.2–32 kpc; no outer boundary, positive deposit density, source history or lensing closure established.']))
print(json.dumps(scores,indent=2))
