#!/usr/bin/env python3
"""Independent, dimensionless first-principles candidate test (not the astrophysical suite).

A finite-supply gain reservoir, an approximately non-radiating mode x,
and three radiating quadrature modes y_j:
 dx = (g*n-gamma_i-gamma_0)*x - sum(delta_j*y_j)
 dy_j = -gamma_b*y_j + delta_j*x
 dn = pump - gamma_R*n - 2*g*n*x*x.

Here x and y are real coordinates of complex envelopes d=x and b_j=-i*y_j.
Their common overall phase is arbitrary. Squared envelopes are normalized modal
energies. Pump input is debited from an upstream supply; all losses are accumulated.
The simulation uses only linear-in-detuning couplings, never a sigma**2 source law.
Detunings represent prescribed relative-motion disturbances, not evolved particle
trajectories. This is NOT a force, lensing, transport, or cosmology derivation.

Needs numpy, scipy and numba. Run: python test_dark_bright.py --out results
"""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PARAMS = dict(g=1., pump=4., gamma_R=1., gamma_i=1., gamma_0=1e-4, gamma_b=1.)

@njit(cache=True)
def rhs(x, y0, y1, y2, n, a0, a1, a2, pump, gi, g0, gb, gr, g):
    xx = x*x
    yy = y0*y0+y1*y1+y2*y2
    return ((g*n-gi-g0)*x-a0*y0-a1*y1-a2*y2,
            -gb*y0+a0*x, -gb*y1+a1*x, -gb*y2+a2*x,
            pump-gr*n-2*g*n*xx,
            2*g0*xx+2*gb*yy, 2*gi*xx+gr*n, 2*gb*yy)

@njit(cache=True)
def ensemble(q, nu, seed=1, cells=64, T=120., burn=40., dt=0.005,
             stop_at=-1., init_scale=1., pump=4., gi=1., g0=1e-4,
             gb=1., gr=1., g=1.):
    np.random.seed(seed)
    gd=gi+g0
    n0=gd/g
    x0=math.sqrt((pump-gr*n0)/(2*gd))
    rho=math.exp(-nu*dt)
    noise=q*math.sqrt(max(0.,1-rho*rho))
    steps=int(round(T/dt)); nb=int(round(burn/dt))
    sum_hot=0.; sum_out=0.; sum_dark=0.; max_resid=0.; final_out=0.; min_n=1e100
    # Independent cells; same seed reused between q values for paired speed comparisons.
    for j in range(cells):
        a0=q*np.random.normal(); a1=q*np.random.normal(); a2=q*np.random.normal()
        x=x0*init_scale; y0=0.;y1=0.;y2=0.;n=n0
        e_initial=x*x+n
        escaped=0.; internal=0.; cell_hot=0.;cell_out=0.;cell_dark=0.;count=0
        for k in range(steps):
            time=k*dt
            if stop_at>=0. and time>=stop_at:
                a0=0.;a1=0.;a2=0.
            elif nu>0.:
                a0=rho*a0+noise*np.random.normal()
                a1=rho*a1+noise*np.random.normal()
                a2=rho*a2+noise*np.random.normal()
            r1=rhs(x,y0,y1,y2,n,a0,a1,a2,pump,gi,g0,gb,gr,g)
            r2=rhs(x+.5*dt*r1[0], y0+.5*dt*r1[1], y1+.5*dt*r1[2], y2+.5*dt*r1[3],
                   n+.5*dt*r1[4],a0,a1,a2,pump,gi,g0,gb,gr,g)
            r3=rhs(x+.5*dt*r2[0], y0+.5*dt*r2[1], y1+.5*dt*r2[2], y2+.5*dt*r2[3],
                   n+.5*dt*r2[4],a0,a1,a2,pump,gi,g0,gb,gr,g)
            r4=rhs(x+dt*r3[0], y0+dt*r3[1], y1+dt*r3[2], y2+dt*r3[3],
                   n+dt*r3[4],a0,a1,a2,pump,gi,g0,gb,gr,g)
            x+=dt*(r1[0]+2*r2[0]+2*r3[0]+r4[0])/6
            y0+=dt*(r1[1]+2*r2[1]+2*r3[1]+r4[1])/6
            y1+=dt*(r1[2]+2*r2[2]+2*r3[2]+r4[2])/6
            y2+=dt*(r1[3]+2*r2[3]+2*r3[3]+r4[3])/6
            n+=dt*(r1[4]+2*r2[4]+2*r3[4]+r4[4])/6
            escaped+=dt*(r1[5]+2*r2[5]+2*r3[5]+r4[5])/6
            internal+=dt*(r1[6]+2*r2[6]+2*r3[6]+r4[6])/6
            if k>=nb:
                cell_hot+=dt*(r1[7]+2*r2[7]+2*r3[7]+r4[7])/6
                cell_out+=dt*(r1[5]+2*r2[5]+2*r3[5]+r4[5])/6
                cell_dark+=dt*x*x
                count+=1
            if n<min_n: min_n=n
        e_final=x*x+y0*y0+y1*y1+y2*y2+n
        resid=abs(e_final+escaped+internal-e_initial-pump*T)/(e_initial+pump*T)
        max_resid=max(max_resid,resid)
        sum_hot+=cell_hot/(count*dt)
        sum_out+=cell_out/(count*dt)
        sum_dark+=cell_dark/(count*dt)
        final_out+=2*g0*x*x+2*gb*(y0*y0+y1*y1+y2*y2)
    return (sum_hot/cells, sum_out/cells, sum_dark/cells, max_resid, final_out/cells, min_n)


def run_group(q:float, nu:float, seeds=range(1,9), **kwargs):
    rows=[ensemble(q,nu,seed=int(s),**kwargs) for s in seeds]
    a=np.array(rows)
    return dict(q=float(q), nu=float(nu), seeds=len(rows),
                cells_per_seed=kwargs.get('cells',64),
                hot_mean=float(a[:,0].mean()), hot_se=float(a[:,0].std(ddof=1)/math.sqrt(len(rows))),
                total_mean=float(a[:,1].mean()), dark_energy=float(a[:,2].mean()),
                max_energy_residual=float(a[:,3].max()), final_output=float(a[:,4].mean()),
                min_reservoir=float(a[:,5].min()), per_seed_hot=a[:,0].tolist())


def generator_motor():
    """Independent active one-port oscillator; input-output coupling balances power exactly.
    a'=(g*n-gi-ge)a+sqrt(2ge)s; n'=P-gr*n-2gn a^2.
    s_out=s-sqrt(2ge)a. No target phase/arcsin or feed/absorb switch is used.
    """
    pump,gi,ge,gr,g=4.,.1,1.,1.,1.
    gt=gi+ge; kap=math.sqrt(2*ge)
    a0=math.sqrt((pump/gt-gr/g)/2); n0=gt/g
    result=[]
    for s in [0.,.1,1.,2.,3.,4.,10.]:
        def fn(t,y):
            a,n,loss,out_net=y
            return [(g*n-gt)*a+kap*s,pump-gr*n-2*g*n*a*a,
                    gr*n+2*gi*a*a,(s-kap*a)**2-s*s]
        sol=solve_ivp(fn,(0,180),[a0,n0,0,0],rtol=1e-10,atol=1e-12)
        assert sol.success
        a,n,loss,out_net=sol.y[:,-1]
        balance=a*a+n+loss+out_net-a0*a0-n0-pump*180
        def steady(aa):
            return (g*pump/(gr+2*g*aa*aa)-gt)*aa+kap*s
        ar=a0 if s==0 else brentq(steady,a0+1e-10, max(100.,20*s))
        result.append(dict(input_amplitude=s,input_power=s*s,amplitude=float(a),reservoir=float(n),
                           net_to_wave=float((s-kap*a)**2-s*s),
                           energy_residual=float(abs(balance)/(pump*180+a0*a0+n0)),
                           steady_amplitude_error=float(abs(a-ar))))
    return dict(parameters=dict(pump=pump,gi=gi,ge=ge,gr=gr,g=g),
                crossover_input_amplitude=math.sqrt((pump-gr*gi/g)*ge/(4*gi)),runs=result)


def symmetry_controls():
    rng=np.random.default_rng(712)
    x1=rng.normal(size=(100000,3));x2=rng.normal(size=(100000,3))
    n=(x1-x2);n/=np.linalg.norm(n,axis=1)[:,None]
    omega=np.array([.7,-.3,1.1]);V=np.array([3.,-7.,2.])
    v1=np.cross(omega,x1)+V;v2=np.cross(omega,x2)+V
    radial=np.sum(n*(v1-v2),axis=1)
    v1r=rng.normal(size=n.shape);v2r=rng.normal(size=n.shape)
    random=np.sum(n*(v1r-v2r),axis=1)
    boosted=np.sum(n*((v1r+V)-(v2r+V)),axis=1)
    return dict(pairs=len(n),rigid_rotation_max=float(np.max(np.abs(radial))),
                uniform_boost_max_change=float(np.max(np.abs(random-boosted))),
                independent_random_radial_variance=float(np.var(random)), expected_variance=2.,
                caveat='A radial relative-velocity coupling is postulated; differential rotation/shear is not eliminated.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=Path('results'))
    ap.add_argument('--part',choices=['static','collisions','controls','no-sink','all'],default='all')
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    result=dict(model='gain-reservoir dark/bright modal test, not a force calculation',parameters=PARAMS)
    if args.part in ('static','all'):
        result['ballistic']=[run_group(q,0.) for q in [0.,.005,.01,.02,.04,.08,.16,.32,.64]]
        for i,r in enumerate(result['ballistic']):
            if i>1:
                prev=result['ballistic'][i-1]
                pair=np.array(r['per_seed_hot'])/np.array(prev['per_seed_hot'])
                r['doubling_ratio']=float(pair.mean());r['doubling_ratio_se']=float(pair.std(ddof=1)/math.sqrt(len(pair)))
        result['generator_motor']=generator_motor()
        result['symmetry']=symmetry_controls()
    if args.part in ('collisions','all'):
        result['collisions']=[run_group(q,nu) for q in [.02,.04,.08] for nu in [1.,10.,20.]]
    if args.part in ('controls','all'):
        result['stop']=[dict(run_group(.08,nu,stop_at=60.,T=160.,burn=120.),initial_nu=nu) for nu in [0.,20.]]
        result['initial_states']=[dict(run_group(.04,0.,init_scale=s),initial_amplitude_scale=s) for s in [.1,1.,3.]]
        result['timestep']=[dict(run_group(.04,nu,dt=dt,seeds=range(1,5)),dt=dt) for nu in [0.,20.] for dt in [.01,.005,.0025]]
    if args.part in ('no-sink','all'):
        result['no_internal_sink']=[run_group(q,0.,seeds=range(1,5),cells=32,T=1000.,burn=600.,dt=.01,pump=.04,gi=0.,g0=.1,gr=0.) for q in [0.,.02,.04,.08,.16]]
    file=args.out/(args.part+'.json');file.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
