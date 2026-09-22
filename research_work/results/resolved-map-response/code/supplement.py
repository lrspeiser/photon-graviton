"""Post-primary force/lensing adjoint and finite-domain control."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import splu
from scipy.optimize import brentq
from response import Mapped,Setting
from jr7_atlas import Grid,Mesh,FAMILIES
from maps import shape


def observable_vectors(M):
    g=M.g;xyz=np.c_[g.X.ravel(),g.Y.ravel(),g.Z.ravel()];r=.6*g.m.rmax;soft=.15;fr=0.;lr=0.
    inc=np.arccos(M.record['assumed_axis_ratio']);los=np.array([0.,np.sin(inc),np.cos(inc)]);e1=np.array([1.,0.,0.]);e2=np.cross(los,e1);projected=np.c_[xyz@e1,xyz@e2]
    for angle in np.arange(24)*2*np.pi/24:
        er=np.array([np.cos(angle),np.sin(angle),0.]);d=xyz-r*er
        fr+=-(d@er)/((d*d).sum(1)+soft**2)**1.5*g.w/24
        er2=er[:2];d=projected-r*er2;lr+=-4*(d@er2)/((d*d).sum(1)+soft**2)*g.w/24
    return np.r_[fr,fr],np.r_[lr,lr]


def adjoint_observables(record):
    M=Mapped(record);g=M.g;s=Setting();A,rhs,gas,h,kp,km=M.system(s);lu=splu(A);y=lu.solve(rhs);P,C=np.split(y,2)
    fv,lv=observable_vectors(M);F=fv@y;L=lv@y
    a,b=FAMILIES[s.family];q=(a*kp*P-2*b*h*km*C)*h*(1-h)
    gradients=[];verifications=[]
    for label,target,value in [('force',fv,F),('lensing',lv,L)]:
        adj=lu.solve(target,trans='T');sensitivity=(adj[g.n:]-adj[:g.n])*q;gradients.append(sensitivity)
        k=int(np.argmax(abs(sensitivity)));values=[];eps=1e-4
        for sign in [1,-1]:
            gg=gas.copy();gg[k]*=np.exp(sign*eps);AA,rr,*_=M.system(s,gg);values.append(float(target@splu(AA).solve(rr)))
        num=(values[0]-values[1])/(2*eps)
        verifications.append(dict(observable=label,analytic=float(sensitivity[k]),finite_difference=num,relative_error=float(abs(num-sensitivity[k])/max(abs(sensitivity[k]),1e-12))))
    mass=gas*g.w;dF=np.divide(gradients[0],mass,out=np.zeros(g.n),where=mass>0);dL=np.divide(gradients[1],mass,out=np.zeros(g.n),where=mass>0)
    # Project lensing gradient off total-mass and force directions in the
    # density-weighted inner product. This is a constrained feasibility probe.
    B=np.c_[np.ones(g.n),dF];coef=np.linalg.solve(B.T@(mass[:,None]*B),B.T@(mass*dL));u=dL-B@coef;u/=max(abs(u))
    normal=dF-float(mass@dF);normal/=max(abs(normal))
    tests=[]
    for epsilon in [.1,.3,.6]:
        def state(beta):
            gg=gas*np.exp(epsilon*u+beta*normal);gg/=gg@g.w
            AA,rr,*_=M.system(s,gg);yy=splu(AA).solve(rr);return gg,yy
        f=lambda beta:float(fv@state(beta)[1]-F)
        fminus,fplus=f(-.5),f(.5)
        if fminus*fplus>0:
            tests.append(dict(epsilon=epsilon,status='constant-force correction not bracketed'));continue
        beta=brentq(f,-.5,.5,xtol=1e-9);gg,yy=state(beta)
        tests.append(dict(epsilon=epsilon,status='solved',force_relative_change=float(fv@yy/F-1),lensing_relative_change=float(lv@yy/L-1),gas_inventory_error=float(gg@g.w-1),gas_fraction_moved=float(.5*g.w@abs(gg-gas)),min_gas_multiplier=float((gg/gas).min()),max_gas_multiplier=float((gg/gas).max()),force_correction_coefficient=beta))
    return dict(object=record['id'],force_both_channels=float(F),lensing_both_channels=float(L),adjoint_checks=verifications,linear_constraint_mass=float(mass@u),linear_constraint_force=float(gradients[0]@u),counterfactuals=tests,scope='Common-source observables in model units; constructed gas rearrangement, not actual galaxy or lens-image fitting')

class Expanded(Mapped):
    def __init__(self,record):
        self.record=record;self.limit=record['aperture_Re'];self.g=Grid(Mesh(20,12,24,2*self.limit,1.5))
    def inputs(self,s):
        g=self.g;mask=g.R<=self.limit
        J=shape(np.array(self.record['source_coeff']),g.R,g.F,self.limit)*np.exp(-.5*(g.Z/.2)**2)*mask;J/=J.ravel()@g.w
        H=shape(np.array(self.record['gas_coeff']),g.R,g.F-np.deg2rad(s.phase_deg),self.limit)
        gas=H**s.proxy_power*np.exp(-.5*(g.Z/s.gas_height)**2)*mask;gas/=gas.ravel()@g.w
        return J.ravel(),gas.ravel()

if __name__=='__main__':
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--maps',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    if a.output.exists():raise FileExistsError('Use a new output')
    records=json.loads(a.maps.read_text());out=dict(adjoint_observables=[],boundary_controls=[])
    for rec in records:
        row=adjoint_observables(rec);out['adjoint_observables'].append(row);print('TANGENT',rec['id'],row['counterfactuals'],flush=True)
        M=Expanded(rec)
        for phase in [0,90]:
            s=Setting(phase_deg=phase);d,_=M.solve(s);out['boundary_controls'].append(dict(object=rec['id'],phase=phase,mesh=M.g.shape,diagnostics=d))
        a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
