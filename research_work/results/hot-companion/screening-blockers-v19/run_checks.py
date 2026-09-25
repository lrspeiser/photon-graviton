#!/usr/bin/env python3
"""Rate-based screening: independent blockers and reset-gate controls.

Exploratory reduced stochastic kinetics, NOT a gravity/astronomical fit.
All times in the stochastic checks are s=gamma*t; x=lambda/gamma.
The event simulator samples creation and survival, never an exp(-x) gate.
A full time-homogeneous master equation is an independent verification.
"""
from __future__ import annotations
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import csv
import json
import math
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from scipy.stats import poisson

PC_M = 3.085677581491367e16
YEAR_S = 365.25 * 86400.0


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def constants() -> dict:
    u, gd, Lpc = 169000.0, 2.03e-10, 0.15
    L = Lpc * PC_M
    tau = L/u
    gamma = 1.0/tau
    vb = gd/gamma
    eta = 1.0/vb
    G, a = 6.67430e-11, 6.30e-11
    ell = a*u/2
    packet_energy_per_cross_section = eta**(-1)*ell/(4*math.pi*G)
    return dict(u_m_s=u, g_d_m_s2=gd, L_pc=Lpc, L_m=L, tau_s=tau,
                tau_years=tau/YEAR_S, gamma_per_s=gamma, eta_s_m=eta,
                v_b_m_s=vb, v_b_over_u=vb/u, gd_L_m2_s2=gd*L,
                single_speed_L_pc=u*u/gd/PC_M,
                single_speed_length_mismatch=u/vb,
                response_corner_per_s=gamma,
                a_m_s2=a,ell_W_kg=ell,
                capture_packet_energy_over_area_J_m2=packet_energy_per_cross_section,
                L019_v_b_m_s=gd*0.19*PC_M/u)


def generator(x: float, cap: int) -> np.ndarray:
    """Column-probability generator of n->n+1 at x, n->n-1 at n.
    Creation at the artificial upper cap is suppressed. Truncation checked.
    """
    Q = np.zeros((cap+1, cap+1))
    for n in range(cap+1):
        if n < cap:
            Q[n+1, n] = x
            Q[n, n] -= x
        if n:
            Q[n-1, n] = n
            Q[n, n] -= n
    return Q


def master_checks() -> tuple[list[dict], dict]:
    rows = []
    for cap in (80, 120):
        for x in (0., .1, .3, 1., 3., 5., 10.):
            # Ambient equilibrium plus exactly one extra launch blocker.
            p = np.zeros(cap+1)
            p[1:] = poisson.pmf(np.arange(cap), x)
            p /= p.sum()
            Q = generator(x, cap)
            for s in (.1, .5, 1., 3., 5.):
                pn = expm(Q*s) @ p
                target = math.exp(-x)*(-math.expm1(-s))
                rows.append(dict(cap=cap,x=x,time_tau=s,master_open=float(pn[0]),
                                 predicted_open=target,error=float(pn[0]-target),
                                 normalization_error=float(pn.sum()-1.),
                                 upper_tail=float(pn[-1])))
    summary = dict(max_abs_open_error=max(abs(r['error']) for r in rows),
                   max_abs_normalization_error=max(abs(r['normalization_error']) for r in rows),
                   max_upper_probability=max(r['upper_tail'] for r in rows))
    assert summary['max_abs_open_error'] < 2e-12, summary
    assert summary['max_abs_normalization_error'] < 2e-12, summary
    return rows, summary


def event_step(n: np.ndarray, x: float, ds: float, rng: np.random.Generator):
    """Exact interval transition using independent Poisson births and lifetimes.

    Actual births: Poisson(x*ds). Each birth time is uniform in the interval;
    averaging its exponential survival gives (1-exp(-ds))/ds. This aggregation
    is exact for interval endpoints, not an Euler approximation. Every created
    blocker is counted, including those that disappear before the endpoint.
    """
    if ds <= 0:
        raise ValueError('ds must be positive')
    keep = math.exp(-ds)
    old_survivors = rng.binomial(n, keep)
    births = rng.poisson(x*ds, size=n.size)
    birth_survival = -math.expm1(-ds)/ds
    new_survivors = rng.binomial(births, birth_survival)
    deaths = (n-old_survivors)+(births-new_survivors)
    new = old_survivors+new_survivors
    assert np.array_equal(new-n, births-deaths)
    return new, int(births.sum()), int(deaths.sum())


def stochastic_checks(n_per_seed: int, seeds: list[int]) -> tuple[list[dict], list[dict], dict]:
    seed_rows, aggregate = [], []
    clocks = [.1, .5, 1., 3., 5.]
    max_ledger_error = 0
    # Start the ambient process empty and evolve for 25 lifetimes, rather than
    # supplying a Poisson equilibrium distribution by hand.
    burn = 25.
    for x in (0., .1, .3, 1., 3., 5., 10.):
        for seed in seeds:
            rng=np.random.default_rng(seed + int(1000*x)*101)
            n=np.zeros(n_per_seed, dtype=np.int64)
            n, births, deaths=event_step(n,x,burn,rng)
            max_ledger_error=max(max_ledger_error,abs(int(n.sum())-(births-deaths)))
            seed_rows.append(dict(x=x,seed=seed,time_tau=0.,stage='ambient',
                                  open_count=int((n==0).sum()),n=n.size,
                                  mean_blockers=float(n.mean()),births=births,deaths=deaths,
                                  predicted_open=math.exp(-x*(1-math.exp(-burn)))))
            n=n+1  # emission creates one additional identical blocker
            initial_count=int(n.sum()); btot=dtot=0
            prev=0.
            for s in clocks:
                n, b, d=event_step(n,x,s-prev,rng)
                btot+=b;dtot+=d
                ledger=int(n.sum())-initial_count-(btot-dtot)
                max_ledger_error=max(max_ledger_error,abs(ledger))
                ambient_mean=x*(1-math.exp(-(burn+s)))
                prediction=math.exp(-ambient_mean)*(-math.expm1(-s))
                seed_rows.append(dict(x=x,seed=seed,time_tau=s,stage='launched',
                                      open_count=int((n==0).sum()),n=n.size,
                                      mean_blockers=float(n.mean()),births=btot,deaths=dtot,
                                      predicted_open=prediction))
                prev=s
        for stage,s in [('ambient',0.)]+[('launched',s) for s in clocks]:
            rr=[r for r in seed_rows if r['x']==x and r['stage']==stage and r['time_tau']==s]
            total=sum(r['n'] for r in rr); count=sum(r['open_count'] for r in rr)
            p=count/total; pred=rr[0]['predicted_open']
            se=math.sqrt(pred*(1-pred)/total)
            aggregate.append(dict(x=x,stage=stage,time_tau=s,n=total,open_count=count,
                                  simulated_open=p,predicted_open=pred,
                                  standard_error=se,z=(p-pred)/se if se else 0.,
                                  mean_blockers=float(np.mean([r['mean_blockers'] for r in rr]))))
    summary=dict(realizations_per_condition=n_per_seed*len(seeds),
                 n_seeds=len(seeds),burn_lifetimes=burn,
                 max_abs_z=max(abs(r['z']) for r in aggregate),
                 max_count_balance_error=max_ledger_error,
                 note='Exact count balance is NOT a complete physical energy/momentum audit.')
    return seed_rows,aggregate,summary


def step_checks():
    rows=[]
    for x0,x1 in ((10.,0.),(0.,10.),(3.,.3),(.3,3.)):
        p=poisson.pmf(np.arange(121),x0)
        Q=generator(x1,120)
        for s in (.1,.5,1.,2.,3.,5.,10.):
            m=x1+(x0-x1)*math.exp(-s)
            blockers=math.exp(-m)
            reset=math.exp(-x1)+(math.exp(-x0)-math.exp(-x1))*math.exp(-s)
            direct=float((expm(Q*s)@p)[0])
            rows.append(dict(x_initial=x0,x_final=x1,time_tau=s,
                             blocker_open=blockers,reset_open=reset,
                             master_open=direct,master_error=direct-blockers))
    assert max(abs(r['master_error']) for r in rows)<2e-12
    return rows


def preparation_controls():
    rows=[]
    for x in (0.,.3,1.,3.,10.):
        Q=generator(x,120)
        for s in (.05,.1,.3,1.,3.,5.):
            y=-math.expm1(-s)
            target=math.exp(-x)*y
            alone=y*math.exp(-x*y)
            p=np.zeros(121);p[1]=1.
            direct=float((expm(Q*s)@p)[0])
            two=math.exp(-x)*y*y
            rows.append(dict(x=x,time_tau=s,ambient_plus_one=target,
                             one_only=alone,one_only_master=direct,
                             master_error=direct-alone,ambient_plus_two=two))
    assert max(abs(r['master_error']) for r in rows)<2e-12
    return rows


def spatial_checks():
    """Dimensionless prescribed radial field: x(r)=(0.3/r)^2, distances in L.
    Not an astrophysical force or Cassini calculation. No trajectories evolved.
    Both models start closed; blocker ambient population starts at local equilibrium.
    """
    rows=[]
    for r0 in (.01,.1,.3):
        def fun(r,y):
            m,launch,reset=y
            x=(.3/r)**2
            return [x-m,-launch,math.exp(-x)-reset]
        for tol in (1e-9,1e-11):
            sol=solve_ivp(fun,(r0,30.),[(.3/r0)**2,1.,0.],method='DOP853',
                          rtol=tol,atol=tol*.01,dense_output=True,max_step=.2)
            assert sol.success
            for r in (1.,3.,5.,10.,20.):
                m,launch,reset=sol.sol(r)
                gate=(1-launch)*math.exp(-m)
                endpoint=math.exp(-(.3/r)**2)*(-math.expm1(-(r-r0)))
                rows.append(dict(r_start_L=r0,rtol=tol,r_L=r,
                                 blocker_open=gate,reset_open=reset,
                                 endpoint_product=endpoint,blocker_mean=m))
    return rows


def covariance_checks():
    rows=[]
    for x in (.1,1.,3.,10.):
        Q=generator(x,120)
        p0=math.exp(-x);p=np.zeros(121);p[0]=1.
        for s in (0.,.1,.5,1.,3.):
            joint=p0*float((expm(Q*s)@p)[0])
            cov=joint-p0*p0
            ana=p0*p0*math.expm1(x*math.exp(-s))
            variance=p0*(1-p0)
            rows.append(dict(x=x,time_tau=s,blocker_cov=cov,analytical_cov=ana,
                             covariance_error=cov-ana,
                             blocker_normalized=cov/variance,reset_normalized=math.exp(-s)))
    assert max(abs(r['covariance_error']) for r in rows)<2e-12
    return rows


def small_signal_checks():
    """Drive both mean equations with x=x0+eps*cos(w t), fit amplitude/phase.
    Same linear susceptibility, differences only higher order in amplitude.
    """
    rows=[]
    x0=1.;eps=1e-3
    for w in (.1,1.,10.):
        period=2*math.pi/w
        tend=max(40.,10*period)
        def f(t,y):
            x=x0+eps*math.cos(w*t)
            return [x-y[0],math.exp(-x)-y[1]]
        sol=solve_ivp(f,(0.,tend),[x0,math.exp(-x0)],method='DOP853',
                      rtol=1e-11,atol=1e-13,dense_output=True,max_step=period/40.)
        tt=np.linspace(tend-5*period,tend,5000,endpoint=False)
        m,rr=sol.sol(tt); rb=np.exp(-m)
        mat=np.column_stack([np.ones(len(tt)),np.cos(w*tt),np.sin(w*tt)])
        cb=np.linalg.lstsq(mat,rb,rcond=None)[0]
        cr=np.linalg.lstsq(mat,rr,rcond=None)[0]
        pred=-math.exp(-x0)/(1+1j*w)
        # y=Re(H eps exp(iwt)) => coefficient_cos=eps ReH, coefficient_sin=-eps ImH
        hb=complex(cb[1],-cb[2])/eps;hr=complex(cr[1],-cr[2])/eps
        rows.append(dict(omega_over_gamma=w,eps=eps,
                         blocker_re=hb.real,blocker_im=hb.imag,
                         reset_re=hr.real,reset_im=hr.imag,
                         linear_re=pred.real,linear_im=pred.imag,
                         blocker_relative_error=abs(hb-pred)/abs(pred),
                         reset_relative_error=abs(hr-pred)/abs(pred)))
    return rows



def reset_gate_checks(n_per_seed: int, seeds: list[int]):
    """A literal two-state gate refreshed by two Gaussian energy quadratures.
    At Poisson clock events (rate 1), draw E/E*=0.5*(X^2+Y^2); open iff E>x.
    It never uses exp(-x) to set a gate. The exponential survival of the
    refresh clock is sampled exactly over finite intervals.
    """
    rows=[]
    for x in (0.,.1,.3,1.,3.,5.,10.):
        totals={s:0 for s in (.1,.5,1.,3.,5.)}
        for seed in seeds:
            rng=np.random.default_rng(seed+404+int(x*1000)*17)
            state=np.zeros(n_per_seed,dtype=bool)
            prev=0.
            for ss in totals:
                refreshed=rng.random(n_per_seed)<-math.expm1(-(ss-prev))
                q=rng.standard_normal((int(refreshed.sum()),2))
                state[refreshed]=.5*np.sum(q*q,axis=1)>x
                totals[ss]+=int(state.sum());prev=ss
        n=n_per_seed*len(seeds)
        for ss,count in totals.items():
            pred=math.exp(-x)*(-math.expm1(-ss));se=math.sqrt(pred*(1-pred)/n)
            rows.append(dict(x=x,time_tau=ss,n=n,open_count=count,simulated_open=count/n,
                             predicted_open=pred,standard_error=se,
                             z=(count/n-pred)/se if se else 0.))
    return rows



def capacity_checks():
    """Analytic independent finite-site control, not fitted to any data."""
    rows=[]
    for x in (.1,1.,3.,10.):
        for K in (10,100,1000):
            rr=(1+x/K)**(-K)
            rows.append(dict(x=x,available_sites=K,finite_site_open=rr,
                             poisson_open=math.exp(-x),
                             relative_difference=rr/math.exp(-x)-1.))
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'results')
    ap.add_argument('--n-per-seed',type=int,default=200000)
    ap.add_argument('--seeds',type=int,nargs='+',default=[11,23,37,53,71])
    a=ap.parse_args()
    if a.n_per_seed<1 or len(set(a.seeds))!=len(a.seeds):
        ap.error('Require positive ensemble size and distinct seeds')
    a.output_dir.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();summary=dict(status='candidate stochastic mechanism, not validated gravity',
         python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
         config=dict(n_per_seed=a.n_per_seed,seeds=a.seeds),constants=constants())
    rows,s=master_checks();write_csv(a.output_dir/'master_equation.csv',rows);summary['master']=s
    print('Master equation:',s,flush=True)
    seeds,rows,s=stochastic_checks(a.n_per_seed,a.seeds)
    write_csv(a.output_dir/'stochastic_by_seed.csv',seeds);write_csv(a.output_dir/'stochastic_aggregate.csv',rows)
    summary['stochastic']=s;print('Stochastic:',s,flush=True)
    resetrows=reset_gate_checks(a.n_per_seed,a.seeds)
    write_csv(a.output_dir/'reset_gate_stochastic.csv',resetrows)
    summary['reset_gate']=dict(max_abs_z=max(abs(r['z']) for r in resetrows),realizations_per_condition=a.n_per_seed*len(a.seeds))
    print('Reset gate:',summary['reset_gate'],flush=True)
    for name,f in [('field_steps',step_checks),('initial_preparation',preparation_controls),
                   ('spatial_history',spatial_checks),('gate_covariance',covariance_checks),
                   ('small_signal',small_signal_checks),('finite_capacity_control',capacity_checks)]:
        rows=f();write_csv(a.output_dir/f'{name}.csv',rows)
        print(name,len(rows),'rows',flush=True)
    summary['seconds']=time.monotonic()-start
    (a.output_dir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    main()
