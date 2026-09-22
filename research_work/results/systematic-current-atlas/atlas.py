#!/usr/bin/env python3
"""JR-7 positive two-channel, fully spatial current atlas; not a gravity validation.

Cylindrical 3D finite volumes. Source and gas patterns have a fixed relative
orientation in the chosen frame; relative angular currents are prescribed.
Run: OPENBLAS_NUM_THREADS=1 python atlas.py --output fresh_directory
"""
from __future__ import annotations
from dataclasses import dataclass, asdict, replace
from functools import lru_cache
from itertools import product
import argparse, hashlib, json, math, time
from pathlib import Path
import numpy as np
from scipy.sparse import bmat, coo_matrix, diags, eye, kron
from scipy.sparse.linalg import spsolve, gmres, LinearOperator, splu, expm_multiply

@dataclass(frozen=True)
class Mesh:
    nr: int=8
    nz: int=10
    nphi: int=24
    rmax: float=6.
    zmax: float=3.

@dataclass(frozen=True)
class Case:
    family: str='mixed'
    phase: float=0.
    omega_c: float=0.
    structured: int=1
    height: float=.35
    loss_c: float=.05
    dz_ratio: float=1.
    source_tilt: float=0.
    source_offset: float=0.
    two_source: int=0
    omega_p: float=-.6

FAMILIES={'constant':(0.,0.),'enhance':(2.,0.),'return':(0.,4.),'mixed':(2.,4.)}

class Grid:
    def __init__(self,m):
        self.m=m; self.shape=(m.nr,m.nz,m.nphi); self.n=np.prod(self.shape)
        self.dr=m.rmax/m.nr;self.dz=2*m.zmax/m.nz;self.df=2*np.pi/m.nphi
        self.rf=np.linspace(0,m.rmax,m.nr+1);self.r=(self.rf[:-1]+self.rf[1:])/2
        self.z=(np.arange(m.nz)+.5)*self.dz-m.zmax
        self.phi=(np.arange(m.nphi)+.5)*self.df
        self.R,self.Z,self.F=np.meshgrid(self.r,self.z,self.phi,indexing='ij')
        self.X=self.R*np.cos(self.F);self.Y=self.R*np.sin(self.F)
        self.wr=np.diff(self.rf**2)/2
        self.w=np.broadcast_to(self.wr[:,None,None]*self.dz*self.df,self.shape).ravel().copy()
        up=self.rf[1:-1]/self.dr/self.wr[:-1];down=self.rf[1:-1]/self.dr/self.wr[1:]
        lr=diags([down,-np.r_[up,0.]-np.r_[0.,down],up],[-1,0,1],shape=(m.nr,m.nr))
        cz=np.full(m.nz,-2.);cz[[0,-1]]=-1
        lz=diags([np.ones(m.nz-1),cz,np.ones(m.nz-1)],[-1,0,1])/self.dz**2
        self.lr=kron(kron(lr,eye(m.nz)),eye(m.nphi),format='csc')
        self.lz=kron(kron(eye(m.nr),lz),eye(m.nphi),format='csc')
        self.idx=np.arange(self.n).reshape(self.shape)
        self._transport={}
        self.ref_gas_mass=float((np.exp(-self.R/2)/np.cosh(self.Z/.35)**2).ravel()@self.w)
    def transport(self,d,ratio,omega):
        key=(d,ratio,omega)
        if key in self._transport:return self._transport[key]
        eta=d/self.R**2;pe=omega*self.df/eta
        def B(x):
            out=np.ones_like(x);keep=abs(x)>1e-8
            out[keep]=x[keep]/np.expm1(x[keep]);out[~keep]=1-x[~keep]/2+x[~keep]**2/12
            return out
        prev=eta/self.df**2*B(-pe);nxt=eta/self.df**2*B(pe)
        rr=self.idx.ravel()
        angular=coo_matrix((np.r_[prev.ravel(),nxt.ravel(),-(prev+nxt).ravel()],
            (np.tile(rr,3),np.r_[np.roll(self.idx,1,axis=2).ravel(),np.roll(self.idx,-1,axis=2).ravel(),rr])),shape=(self.n,self.n)).tocsc()
        L=d*self.lr+d*ratio*self.lz+angular
        self._transport[key]=L
        return L
    def fields(self,c):
        # Positivity and equal gas inventory; no velocity data determine a shape.
        gas=np.exp(-self.R/2)/np.cosh(self.Z/c.height)**2
        if c.structured:
            f=self.R**2/(self.R**2+.8**2)
            gas*=1+.7*f*np.cos(2*(self.F-c.phase))+.2*f*np.cos(4*(self.F-c.phase))
        gas*=self.ref_gas_mass/(gas.ravel()@self.w)
        a,b=FAMILIES[c.family];h=gas/(gas+.15)
        kp=.5*np.exp(a*h);km=.5*np.exp(b*h*h)
        tilt=c.source_tilt
        xx=(self.X-c.source_offset)*np.cos(tilt)+self.Z*np.sin(tilt)
        zz=-(self.X-c.source_offset)*np.sin(tilt)+self.Z*np.cos(tilt)
        rr2=xx*xx+self.Y*self.Y
        J=np.exp(-rr2/(2*1.5**2)-zz*zz/(2*.35**2))
        # A smooth quadrupolar source avoids a coordinate singularity.
        J*=1+.6*(xx*xx-self.Y*self.Y)/(rr2+.8**2)
        if c.two_source:
            J=np.exp(-((self.X-1.3)**2+self.Y**2+self.Z**2)/(2*.55**2))
            J+=np.exp(-((self.X+1.3)**2+self.Y**2+self.Z**2)/(2*.55**2))
            # Two hot, broad gas maxima displaced perpendicular to source pair.
            gas=np.exp(-((self.X-.6)**2+(self.Y-.7)**2+self.Z**2)/(2*.95**2))
            gas+=np.exp(-((self.X+.6)**2+(self.Y+.7)**2+self.Z**2)/(2*.95**2))
            gas*=self.ref_gas_mass/(gas.ravel()@self.w)
            h=gas/(gas+.15);kp=.5*np.exp(a*h);km=.5*np.exp(b*h*h)
        J/=J.ravel()@self.w
        return kp.ravel(),km.ravel(),J.ravel(),gas.ravel()
    def solve(self,c,save_fields=False):
        kp,km,J,gas=self.fields(c)
        LP=self.transport(.5,1.,c.omega_p);LC=self.transport(.05,c.dz_ratio,c.omega_c)
        A=bmat([[-LP+diags(kp+.5),diags(-km)],[diags(-kp),-LC+diags(km+c.loss_c)]],format='csc')
        rhs=np.r_[J,np.zeros(self.n)]
        t0=time.monotonic();y=spsolve(A,rhs);P,C=np.split(y,2)
        if not np.all(np.isfinite(y)):raise ArithmeticError('Nonfinite state')
        ep=float(P@self.w);ec=float(C@self.w)
        forward=float((kp*P)@self.w);reverse=float((km*C)@self.w)
        av=lambda v:float(v@self.w/self.w.sum())
        corr_forward=forward-self.w.sum()*av(kp)*av(P)
        corr_reverse=reverse-self.w.sum()*av(km)*av(C)
        ang=np.exp(2j*self.F.ravel());r2=self.R.ravel()**2
        q=complex((C*r2*ang)@self.w)/max(float((C*r2)@self.w),1e-300)
        center=[float((C*v.ravel())@self.w/ec) for v in (self.X,self.Y,self.Z)]
        d=dict(Ep=ep,Ec=ec,total_energy=ep+ec,outgoing=.5*ep+c.loss_c*ec,
            source=float(J@self.w),gas_inventory=float(gas@self.w),
            min_state=float(y.min()),residual=float(np.linalg.norm(A@y-rhs)/np.linalg.norm(rhs)),
            energy_balance_error=abs(.5*ep+c.loss_c*ec-1.),
            companion_balance_error=abs(forward-reverse-c.loss_c*ec),
            forward=forward,reverse=reverse,cov_forward=corr_forward,cov_reverse=corr_reverse,
            normalized_quadrupole=abs(q),major_axis_deg=float(np.angle(q)*90/np.pi),
            mean_R2=float((C*r2)@self.w/ec),mean_Z2=float((C*self.Z.ravel()**2)@self.w/ec),
            centroid=center,solve_seconds=time.monotonic()-t0)
        return d,(P,C,gas,J),A
    def readout(self,fields,softening=.35):
        P,C,gas,J=fields
        xyz=np.array([self.X.ravel(),self.Y.ravel(),self.Z.ravel()]).T
        weights=np.array([C*self.w,(P+C)*self.w])
        radius=2.5;angles=np.arange(16)*np.pi/8;inclinations=(0.,np.pi/4,np.pi/2)
        force=[];lens=[]
        for a in angles:
            er=np.array([np.cos(a),np.sin(a),0.]);et=np.array([-np.sin(a),np.cos(a),0.])
            delta=xyz-radius*er
            accel=(weights/((delta*delta).sum(1)+softening**2)**1.5)@delta
            force.append(dict(angle=float(a),inward=(-accel@er).tolist(),tangent=(accel@et).tolist()))
        for tilt in inclinations:
            los=np.array([0.,np.sin(tilt),np.cos(tilt)])
            e1=np.array([1.,0.,0.]);e2=np.cross(los,e1)
            proj=np.column_stack([xyz@e1,xyz@e2]);rows=[]
            for a in angles:
                er=np.array([np.cos(a),np.sin(a)]);delta=proj-radius*er
                bend=(4*weights/((delta*delta).sum(1)+softening**2))@delta
                rows.append((-bend@er).tolist())
            z=np.array(rows)
            lens.append(dict(inclination_deg=float(tilt*180/np.pi),mean_inward=z.mean(0).tolist(),half_range=(np.ptp(z,axis=0)/2).tolist()))
        fr=np.array([a['inward'] for a in force]);ft=np.array([a['tangent'] for a in force])
        return dict(radius=radius,softening=softening,channels=['companion','both'],mean_inward=fr.mean(0).tolist(),
            tangent_max_over_mean=(abs(ft).max(0)/abs(fr.mean(0))).tolist(),force_angles=force,lens= lens)

def factorial():
    for vals in product(FAMILIES,(0.,np.pi/4,np.pi/2),(-.6,0.,.6),(0,1),(.35,.7),(.05,.3),(.25,4.)):
        yield Case(*vals)

def analyze(rows):
    # Matched phase contrasts, rather than selecting unrelated extremes.
    groups={}
    for r in rows:
        c=r['case'];key=tuple((k,v) for k,v in c.items() if k!='phase')
        groups.setdefault(key,[]).append(r)
    contrasts=[]
    for group in groups.values():
        group=sorted(group,key=lambda r:r['case']['phase'])
        if len(group)!=3:continue
        a,b=group[0],group[-1]
        contrasts.append(dict(base_case={k:v for k,v in a['case'].items() if k!='phase'},
            aligned=a['id'],quarter=group[1]['id'],perpendicular=b['id'],
            aligned_Ec=a['diagnostics']['Ec'],perpendicular_Ec=b['diagnostics']['Ec'],
            aligned_over_perpendicular_minus1=a['diagnostics']['Ec']/b['diagnostics']['Ec']-1))
    variable=[r for r in contrasts if r['base_case']['structured'] and r['base_case']['family']!='constant']
    families={}
    for fam in FAMILIES:
        v=[r['aligned_over_perpendicular_minus1'] for r in contrasts if r['base_case']['structured'] and r['base_case']['family']==fam]
        families[fam]=dict(min=min(v),max=max(v),median=float(np.median(v)),positive=sum(x>1e-5 for x in v),negative=sum(x< -1e-5 for x in v))
    # Freeze representative IDs before finer runs. Retain both signs and constant control.
    selected=[min(variable,key=lambda r:r['aligned_over_perpendicular_minus1']),max(variable,key=lambda r:r['aligned_over_perpendicular_minus1'])]
    picked=sorted(set(v for x in selected for v in [x['aligned'],x['perpendicular']]))
    return dict(phase_contrasts=contrasts,families=families,representative_ids=picked,
        max_balance=max(r['diagnostics']['energy_balance_error'] for r in rows),
        max_linear_residual=max(r['diagnostics']['residual'] for r in rows),
        min_state=min(r['diagnostics']['min_state'] for r in rows),
        axisymmetric_phase_control=max(abs(r['aligned_over_perpendicular_minus1']) for r in contrasts if not r['base_case']['structured']),
        no_gas_rate_phase_control=max(abs(r['aligned_over_perpendicular_minus1']) for r in contrasts if r['base_case']['family']=='constant'))

def main():
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--limit',type=int,default=0);ap.add_argument('--skip-refinement',action='store_true');args=ap.parse_args()
    out=args.output
    if out.exists():raise FileExistsError('Use a new output directory')
    out.mkdir(parents=True);t0=time.monotonic();mesh=Mesh();grid=Grid(mesh);cases=list(factorial())
    if args.limit:cases=cases[:args.limit]
    (out/'declaration.json').write_text(json.dumps(dict(mesh=asdict(mesh),cases=[asdict(c) for c in cases],code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
    rows=[]
    with (out/'primary.jsonl').open('w') as f:
        for i,c in enumerate(cases):
            d,_,_=grid.solve(c);r=dict(id=i,case=asdict(c),diagnostics=d);rows.append(r);f.write(json.dumps(r)+'\n');f.flush()
            if i%24==0:print('PRIMARY',i,'of',len(cases),'seconds',round(time.monotonic()-t0,2),flush=True)
    if args.limit:return
    analysis=analyze(rows);(out/'selection-before-refinement.json').write_text(json.dumps(analysis,indent=2));print('PHASE',json.dumps(analysis['families']),flush=True)
    refined=[]
    if not args.skip_refinement:
        for m in [mesh,Mesh(12,16,48),Mesh(16,20,72)]:
            g=Grid(m)
            for i in analysis['representative_ids']:
                d,f,A=g.solve(cases[i]);r=dict(id=i,mesh=asdict(m),case=asdict(cases[i]),diagnostics=d,readout=g.readout(f));refined.append(r)
                np.savez_compressed(out/f'field_{i}_{m.nr}.npz',P=f[0],C=f[1],gas=f[2],J=f[3],weight=g.w,R=g.R,Z=g.Z,phi=g.F)
                print('REFINED',i,m.nr,'Ec',d['Ec'],'seconds',round(time.monotonic()-t0,2),flush=True)
    controls=[]
    g=Grid(Mesh(10,12,32))
    for label,c in [('tilt0',Case()),('tilt45',Case(source_tilt=np.pi/4)),('tilt90',Case(source_tilt=np.pi/2)),
                    ('offset',Case(source_offset=1.)),('two_source_cluster_analogue',Case(two_source=1)),
                    ('mirror_reference',Case(phase=np.pi/4,omega_c=.6)),('mirror',Case(phase=-np.pi/4,omega_c=-.6,omega_p=.6))]:
        d,f,A=g.solve(c);controls.append(dict(label=label,case=asdict(c),diagnostics=d,readout=g.readout(f)))
        np.savez_compressed(out/(label+'.npz'),P=f[0],C=f[1],gas=f[2],J=f[3],weight=g.w,R=g.R,Z=g.Z,phi=g.F)
    summary=dict(experiment='JR-7',scope='576 coarse 3D positive spatial transport solutions; not a real-galaxy fit',
        primary_count=len(rows),analysis=analysis,refinements=refined,controls=controls,seconds=time.monotonic()-t0)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False));print('DONE',time.monotonic()-t0,flush=True)
if __name__=='__main__':main()
