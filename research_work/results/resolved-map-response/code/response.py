"""Map-conditioned positive 3D exchange, phase controls and exact adjoints.

Dimensionless shape-only construction. Not calibrated galaxy gravity.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
from pathlib import Path
from itertools import product
import argparse,json,hashlib,time
import numpy as np
from scipy.sparse import bmat,diags
from scipy.sparse.linalg import splu
from jr7_atlas import Mesh,Grid,FAMILIES
from maps import shape

@dataclass(frozen=True)
class Setting:
    family:str='mixed'
    phase_deg:float=0.
    omega:float=0.
    gas_height:float=.20
    proxy_power:float=1.

class Mapped:
    def __init__(self,record,nr=10,nz=12,nphi=24):
        self.record=record;self.g=Grid(Mesh(nr,nz,nphi,record['aperture_Re'],1.5))
        self.Jxy=shape(np.array(record['source_coeff']),self.g.R,self.g.F,record['aperture_Re'])
    def inputs(self,s):
        g=self.g
        J=self.Jxy*np.exp(-.5*(g.Z/.20)**2);J/=J.ravel()@g.w
        hx=shape(np.array(self.record['gas_coeff']),g.R,g.F-np.deg2rad(s.phase_deg),self.record['aperture_Re'])
        gas=hx**s.proxy_power*np.exp(-.5*(g.Z/s.gas_height)**2);gas/=gas.ravel()@g.w
        return J.ravel(),gas.ravel()
    def system(self,s,gas_override=None):
        g=self.g;J,gas=self.inputs(s)
        if gas_override is not None:gas=np.asarray(gas_override)
        h=gas/(gas+.15);a,b=FAMILIES[s.family];kp=.5*np.exp(a*h);km=.5*np.exp(b*h*h)
        LP=g.transport(.5,1.,0.);LC=g.transport(.05,1.,s.omega)
        A=bmat([[-LP+diags(kp+.5),-diags(km)],[-diags(kp),-LC+diags(km+.05)]],format='csc')
        return A,np.r_[J,np.zeros(g.n)],gas,h,kp,km
    def solve(self,s,adjoint=False):
        g=self.g;A,rhs,gas,h,kp,km=self.system(s);lu=splu(A);y=lu.solve(rhs);P,C=np.split(y,2)
        ep=float(g.w@P);ec=float(g.w@C)
        fwd=float(g.w@(kp*P));rev=float(g.w@(km*C));diag=dict(Ep=ep,Ec=ec,energy=ep+ec,source=float(g.w@rhs[:g.n]),outgoing=.5*ep+.05*ec,energy_error=abs(.5*ep+.05*ec-1),min_state=float(y.min()),linear_residual=float(np.linalg.norm(A@y-rhs)/np.linalg.norm(rhs)),forward=fwd,reverse=rev)
        sens={};arrays=dict(P=P,C=C,gas=gas,J=rhs[:g.n],weight=g.w,R=g.R,Z=g.Z,phi=g.F)
        if adjoint:
            target=np.r_[np.zeros(g.n),g.w];psi=lu.solve(target,trans='T');pp,pc=np.split(psi,2)
            a,b=FAMILIES[s.family];q=(a*kp*P-2*b*h*km*C)*h*(1-h);S=(pc-pp)*q
            Sn=S-gas*g.w/(gas@g.w)*S.sum()
            arrays.update(adjoint_P=pp,adjoint_C=pc,sensitivity_loggas=S,sensitivity_fixed_inventory=Sn)
            sens=dict(positive_response_sum=float(S[S>0].sum()),negative_response_sum=float(S[S<0].sum()),global_loggas_derivative=float(S.sum()),normalized_redistribution_sum=float(Sn.sum()),positive_adjoint_difference_fraction=float(np.mean(pc>pp)),largest_positive_density=float(gas[np.argmax(S)]),largest_negative_density=float(gas[np.argmin(S)]))
            checks=[]
            for kind,der in [('local',S),('fixed_inventory',Sn)]:
                for label,k in [('positive',int(np.argmax(der))),('negative',int(np.argmin(der)))]:
                    eps=1e-4;vs=[]
                    for sign in [1,-1]:
                        gg=gas.copy();gg[k]*=np.exp(sign*eps)
                        if kind=='fixed_inventory':gg*=float(gas@g.w)/float(gg@g.w)
                        AA,rr,*_=self.system(s,gg);yy=splu(AA).solve(rr);vs.append(float(g.w@yy[g.n:]))
                    numerical=(vs[0]-vs[1])/(2*eps);scale=max(abs(der[k]),1e-12)
                    checks.append(dict(kind=kind,sign=label,index=k,analytic=float(der[k]),finite_difference=numerical,relative_error=float(abs(numerical-der[k])/scale)))
            sens['finite_difference_checks']=checks
            # Observable finite redistribution: move 5% proxy mass FROM the
            # low-susceptibility half TO high-susceptibility half. This is an
            # operator-directed counterfactual, never a changed observed map.
            susceptibility=np.divide(S,gas*g.w,out=np.zeros_like(S),where=gas>0)
            mass=gas*g.w;ids=np.argsort(susceptibility);cumulative=np.cumsum(mass[ids]);low=ids[cumulative<=.5];high=ids[cumulative>.5]
            fraction=.05;gg=gas.copy();gg[low]*=1-fraction/mass[low].sum();gg[high]*=1+fraction/mass[high].sum()
            AA,rr,*_=self.system(s,gg);yy=splu(AA).solve(rr);delta=gg-gas
            linear=float(np.sum(S*delta/np.maximum(gas,1e-300)));actual=float(g.w@yy[g.n:]-ec)
            sens['finite_mass_redistribution']=dict(mass_fraction=.05,inventory_error=float(gg@g.w-1),linear_delta_Ec=linear,actual_delta_Ec=actual,fractional_Ec_change=actual/ec,scope='Computed counterfactual from adjoint, NOT a measured change or independent prediction')
        diag['sensitivity']=sens
        return diag,arrays
    def readout(self,a):
        g=self.g;xyz=np.c_[g.X.ravel(),g.Y.ravel(),g.Z.ravel()];mass=np.array([a['C']*g.w,(a['P']+a['C'])*g.w]);soft=.15
        radii=[.6*g.m.rmax,.9*g.m.rmax];theta=np.arange(24)*2*np.pi/24;res=[]
        for r in radii:
            fr=[];ft=[];lens=[]
            for ang in theta:
                er=np.array([np.cos(ang),np.sin(ang),0.]);et=np.array([-np.sin(ang),np.cos(ang),0.]);d=xyz-r*er
                acceleration=(mass/((d*d).sum(1)+soft**2)**1.5)@d;fr.append(-acceleration@er);ft.append(acceleration@et)
            for incl in [0.,np.arccos(self.record['assumed_axis_ratio'])]:
                los=np.array([0.,np.sin(incl),np.cos(incl)]);e1=np.array([1.,0.,0.]);e2=np.cross(los,e1);projected=np.c_[xyz@e1,xyz@e2];bends=[]
                for ang in theta:
                    er=np.array([np.cos(ang),np.sin(ang)]);d=projected-r*er;bend=(4*mass/((d*d).sum(1)+soft**2))@d;bends.append(-bend@er)
                lens.append(dict(inclination_deg=float(incl*180/np.pi),mean_radial_deflection=np.mean(bends,0).tolist()))
            res.append(dict(radius_Re=r,mean_inward=np.mean(fr,0).tolist(),tangent_RMS=np.sqrt(np.mean(np.array(ft)**2,0)).tolist(),force_angle_inward=np.array(fr).tolist(),lens=lens))
        return dict(channels=['C','P+C'],softening_Re=soft,radii=res)

def primary(records,out):
    out.mkdir(parents=True,exist_ok=False);rows=[];t0=time.monotonic()
    for rec in records:
        M=Mapped(rec)
        settings=[Setting(f,p,o) for f,p,o in product(FAMILIES,[0,90,180,270],[0.,-.6,.6])]
        settings += [Setting(f,p,0.,h,power) for f,p,(h,power) in product(FAMILIES,[0,90,180,270],[(.45,1.),(.20,.5),(.45,.5)])]
        for k,s in enumerate(settings):
            adj=s==Setting();d,a=M.solve(s,adjoint=adj);row=dict(object=rec['id'],setting=asdict(s),mesh=asdict(M.g.m),diagnostics=d)
            if s.phase_deg in [0,90] and s.omega==0 and s.gas_height==.2 and s.proxy_power==1:
                row['readout']=M.readout(a)
            rows.append(row)
            if adj:np.savez_compressed(out/(rec['id']+'-fields.npz'),**a)
            if k%16==0:print('SOLVE',rec['id'],k,len(settings),'seconds',round(time.monotonic()-t0,1),flush=True)
        (out/'primary.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    print('PRIMARY_DONE',len(rows),time.monotonic()-t0,flush=True)
    return rows

def refine(records,out):
    rows=[]
    for nr,nz,nphi in [(14,16,48),(18,24,64)]:
        for rec in records:
            M=Mapped(rec,nr,nz,nphi)
            for phase in [0,90]:
                s=Setting(phase_deg=phase);d,a=M.solve(s);row=dict(object=rec['id'],setting=asdict(s),mesh=asdict(M.g.m),diagnostics=d,readout=M.readout(a));rows.append(row)
                print('REFINE',rec['id'],nr,phase,d['Ec'],flush=True)
                if nr==18 and phase==0:np.savez_compressed(out/(rec['id']+'-fine-fields.npz'),**a)
                (out/'refinement.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    return rows

if __name__=='__main__':
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--maps',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--refine-only',action='store_true');a=ap.parse_args();records=json.loads(a.maps.read_text())
    if a.refine_only:refine(records,a.output)
    else:primary(records,a.output)
