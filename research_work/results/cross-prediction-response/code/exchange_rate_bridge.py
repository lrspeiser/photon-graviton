#!/usr/bin/env python3
"""Conditional two-state conversion/return bridge; NOT a galaxy fit.

u, Delta and Gamma have dimensions inverse time. The worked examples use a
single arbitrary reference time, not a fitted astronomical clock.

A two-state Hamiltonian has H/hbar = [[Delta/2,u],[u,-Delta/2]].
With specified coherence decay Gamma, define z=p_R-p_X and c=x+i*y:
 z_dot=-4*u*y, x_dot=-Gamma*x+Delta*y,
 y_dot=u*z-Delta*x-Gamma*y.
Adiabatic elimination gives p_X_dot=k*(p_R-p_X),
 k=2*u**2*Gamma/(Gamma**2+Delta**2).

The damping is a prescribed reduced open-system assumption, not a derived
material bath or a closed total-energy ledger. Exact below means exact for
THIS reduced linear system, not exact microscopic photon-companion physics.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.linalg import eigvals, expm


def rate(u: float, detuning: float, gamma: float) -> float:
    if not np.isfinite([u, detuning, gamma]).all() or u < 0 or gamma <= 0:
        raise ValueError('Finite u>=0 and gamma>0 are required')
    return float(2*u*u*gamma/(gamma*gamma+detuning*detuning))


def calculate() -> dict:
    u, delta = .02, 1.
    rows=[]
    for gamma in (.1, 1., 10.):
        matrix=np.array([[0.,0.,-4*u],[0.,-gamma,delta],[u,-delta,-gamma]])
        poles=eigvals(matrix)
        slow=poles[np.argmax(poles.real)]
        k=rate(u,delta,gamma)
        duration=1/(2*k)
        times=np.linspace(0,duration,251)
        probabilities=np.array([(1-(expm(matrix*t)@np.array([1.,0.,0.]))[0])/2 for t in times])
        approximate=.5*(1-np.exp(-2*k*times))
        rows.append(dict(u=u,detuning=delta,coherence_decay=gamma,
            effective_forward_and_return_rate=k,
            exact_reduced_slow_population_difference_decay=float(-slow.real),
            approximate_slow_decay=2*k,
            relative_slow_decay_error=float(abs(-slow.real/(2*k)-1)),
            max_population_error_over_one_slow_time=float(np.max(abs(probabilities-approximate))),
            minimum_population=float(probabilities.min()),maximum_population=float(probabilities.max())))
    grid=np.geomspace(.01,100.,2001)
    values=2*u*u*grid/(grid*grid+delta*delta)
    analytic_gamma_peak=abs(delta)
    derivative=lambda g:2*u*u*(delta*delta-g*g)/(g*g+delta*delta)**2
    errors=[]
    for g in (.1,.2,2.,10.):
        eps=g*1e-5
        numeric=(rate(u,delta,g+eps)-rate(u,delta,g-eps))/(2*eps)
        errors.append(abs(numeric/derivative(g)-1))
    return dict(status='Executed conditional mechanism example; no observed gravity prediction',
       conventions='All numerical rates are inverse arbitrary reference time',
       formula='k = 2*u^2*Gamma/(Gamma^2+Delta^2)',
       rows=rows,analytic_peak_gamma=analytic_gamma_peak,
       sampled_peak_gamma=float(grid[np.argmax(values)]),
       maximum_derivative_relative_error=float(max(errors)),
       assumptions=['Two effective channels that must eventually include photon, companion and material recoil states',
         'Prescribed coherence decay, not a derived energy-conserving bath',
         'Weak coupling and separation of coherence/population timescales for the rate approximation',
         'State populations are not automatically fractions of radiation versus gravitational energy'],
       open_physical_links=['Actual transition matrix element and gas response spectrum',
         'Absolute rates, outgoing spectral/angular distribution and momentum exchange',
         'Travel, binding, capacity, self-gravity and finite source fuel',
         'Shared force/light coupling and sufficient gravitational amplitude'],
       checks=dict(nonnegative_populations=all(r['minimum_population']>=-1e-12 for r in rows),
         probabilities_below_one=all(r['maximum_population']<=1+1e-12 for r in rows),
         reduced_slow_rate_accuracy_below_0p2percent=all(r['relative_slow_decay_error']<.002 for r in rows),
         derivative_relative_accuracy_below_1e_7=max(errors)<1e-7))


def main():
    ap=argparse.ArgumentParser(__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    if args.output.exists():raise FileExistsError('Preserve old results: choose a new filename')
    result=calculate()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))
    if not all(result['checks'].values()):raise SystemExit('Reduced-model check failed')

if __name__=='__main__':main()
