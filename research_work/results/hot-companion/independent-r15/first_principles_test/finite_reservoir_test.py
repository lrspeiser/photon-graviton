#!/usr/bin/env python3
"""Finite-fuel refinement: no fixed pump and no nonradiative sink.
A dark internal mode is the finite energy store of the matter; motion linearly
mixes it into radiating modes. Energy lost by the modes is counted as radiation.
This is a local envelope model, NOT spatial gravity or a test-body force model.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
from numba import njit
from scipy.linalg import expm

@njit(cache=True)
def rate(x,y0,y1,y2,a0,a1,a2,g0):
    return (-g0*x-a0*y0-a1*y1-a2*y2,
            -y0+a0*x,-y1+a1*x,-y2+a2*x,
            2*g0*x*x+2*(y0*y0+y1*y1+y2*y2),
            2*(y0*y0+y1*y1+y2*y2))

@njit(cache=True)
def sim(q,nu,seed,cells=64,T=80.,burn=10.,dt=.005,g0=1e-6,stop=-1.):
    np.random.seed(seed)
    rho=math.exp(-nu*dt);noise=q*math.sqrt(1-rho*rho)
    out=0.;hot=0.;min_final=1.;max_res=0.;final_hot=0.;store_mean=0.
    nstep=int(round(T/dt));n_burn=int(round(burn/dt))
    for j in range(cells):
        a0=q*np.random.normal();a1=q*np.random.normal();a2=q*np.random.normal()
        x=1.;y0=0.;y1=0.;y2=0.;esc=0.;o=0.;h=0.;u=0.
        for k in range(nstep):
            if stop>=0 and k*dt>=stop:a0=0.;a1=0.;a2=0.
            elif nu>0:
                a0=rho*a0+noise*np.random.normal();a1=rho*a1+noise*np.random.normal();a2=rho*a2+noise*np.random.normal()
            r1=rate(x,y0,y1,y2,a0,a1,a2,g0)
            r2=rate(x+.5*dt*r1[0],y0+.5*dt*r1[1],y1+.5*dt*r1[2],y2+.5*dt*r1[3],a0,a1,a2,g0)
            r3=rate(x+.5*dt*r2[0],y0+.5*dt*r2[1],y1+.5*dt*r2[2],y2+.5*dt*r2[3],a0,a1,a2,g0)
            r4=rate(x+dt*r3[0],y0+dt*r3[1],y1+dt*r3[2],y2+dt*r3[3],a0,a1,a2,g0)
            x+=dt*(r1[0]+2*r2[0]+2*r3[0]+r4[0])/6
            y0+=dt*(r1[1]+2*r2[1]+2*r3[1]+r4[1])/6
            y1+=dt*(r1[2]+2*r2[2]+2*r3[2]+r4[2])/6
            y2+=dt*(r1[3]+2*r2[3]+2*r3[3]+r4[3])/6
            de=dt*(r1[4]+2*r2[4]+2*r3[4]+r4[4])/6
            esc+=de
            if k>=n_burn:
                o+=de;h+=dt*(r1[5]+2*r2[5]+2*r3[5]+r4[5])/6;u+=dt*x*x
        e=x*x+y0*y0+y1*y1+y2*y2
        max_res=max(max_res,abs(e+esc-1.))
        min_final=min(min_final,e)
        out+=o/(T-burn);hot+=h/(T-burn);store_mean+=u/(T-burn)
        final_hot+=2*(y0*y0+y1*y1+y2*y2)
    return out/cells,hot/cells,max_res,min_final,final_hot/cells,store_mean/cells


def group(q,nu,**kw):
    a=np.array([sim(q,nu,s,**kw) for s in range(1,9)])
    return dict(q=q,nu=nu,total=float(a[:,0].mean()),hot=float(a[:,1].mean()),
                hot_se=float(a[:,1].std(ddof=1)/math.sqrt(8)),
                max_energy_residual=float(a[:,2].max()),lowest_remaining_fraction=float(a[:,3].min()),
                final_hot=float(a[:,4].mean()),mean_store=float(a[:,5].mean()),per_seed_hot=a[:,1].tolist())


def exact_control():
    delta=np.array([.01,.02,.03]);g0=1e-6
    A=np.zeros((4,4));A[0,0]=-g0;A[1:,1:]=-np.eye(3);A[0,1:]=-delta;A[1:,0]=delta
    y=np.array([1.,0.,0.,0.]);dt=.005;T=80.
    for k in range(round(T/dt)):
        k1=A@y;k2=A@(y+dt*k1/2);k3=A@(y+dt*k2/2);k4=A@(y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
    exact=expm(A*T)@np.array([1.,0.,0.,0.])
    return dict(max_state_error=float(np.max(np.abs(y-exact))))


def main():
    p=Path(__file__).parent/'results';p.mkdir(exist_ok=True)
    free=[group(q,0.) for q in [0.,.0005,.001,.002,.004,.008,.016,.032]]
    cold=free[0]['total']
    for i,r in enumerate(free):
        r['excess']=r['total']-cold
        if i>1:
            r['doubling_ratio']=r['excess']/free[i-1]['excess']
    collision=[group(.004,nu) for nu in [1.,10.,20.]]
    base=next(r for r in free if r['q']==.004)
    for r in collision:
        r['suppression']=r['hot']/base['hot'];r['weak_coupling_prediction']=1/(1+r['nu'])
    stopped=[group(.004,nu,T=100.,burn=70.,stop=40.) for nu in [0.,20.]]
    refinement=[dict(group(.004,20.,dt=dt),dt=dt) for dt in [.01,.0025]]
    result=dict(parameters={'gamma_0':1e-6,'gamma_b':1.,'initial_store_energy':1.,'pump':0.,'internal_dissipation':0.,
                             'seeds':8,'cells_per_seed':64,'T':80.,'burn':10.,'dt':.005},
                caveat='Local finite-fuel leakage, not test-body force, spatial gravity, or the astrophysical suite.',
                free=free,collisions=collision,stopped=stopped,refinement=refinement,exact_control=exact_control())
    (p/'finite_reservoir.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
