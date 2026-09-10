"""Autonomous homogeneous receiver with universal matter coupling.

Dimensionless synthetic controls, not astronomical observations. Existing
photon-only receiver: H = Pi**2/(2*K) + E_r/n. Here E_m/n is included under
the already proposed universal-clock prescription. No fitted parameters.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent


def analytic_time(n, K, U):
    return np.sqrt(K/(2*U))*(np.sqrt(n*(n-1)) + np.arcsinh(np.sqrt(n-1)))


def background(matter, inertia, tol):
    radiation = 1.0
    total = radiation + matter
    scale = np.sqrt(inertia/total)
    # n, canonical Pi, accumulated proper time tau; c=1.
    def rhs(t, y):
        n, pi, tau = y
        return [pi/inertia, total/n**2, 1/n]
    sol = solve_ivp(rhs, (0, 5*scale), [1., 0., 0.],
                    rtol=tol, atol=tol*.01, max_step=.025*scale,
                    dense_output=True)
    assert sol.success
    t = np.linspace(0, 5*scale, 501)
    n, pi, tau = sol.sol(t)
    energy = pi*pi/(2*inertia) + total/n
    energy_error = float(np.max(np.abs(energy/total-1)))
    inverse_time_error = float(np.max(np.abs(analytic_time(n, inertia, total)-t))/scale)
    rows = []
    for target in [1.1, 1.5, 2., 3.]:
        time = brentq(lambda t: sol.sol(t)[0]-target, 0, 5*scale,
                      xtol=1e-13*scale)
        n1, pi1, _ = sol.sol(time)
        photon_loss = radiation*(1-1/n1)
        matter_loss = matter*(1-1/n1)
        field_gain = pi1*pi1/(2*inertia)
        rows.append(dict(matter_to_radiation=matter, inertia=inertia,
                         n=target, reference_time=time,
                         analytic_reference_time=float(analytic_time(target,inertia,total)),
                         photon_reference_energy_loss=photon_loss,
                         matter_reference_energy_loss=matter_loss,
                         field_energy_gain=field_gain,
                         fraction_gain_supplied_by_radiation=photon_loss/field_gain,
                         normalized_budget_error=(field_gain-photon_loss-matter_loss)/total))
    signals=[]
    for emission in [0., .2*scale, .6*scale]:
        for distance in [.2*scale, .7*scale, 1.4*scale]:
            tau_e = sol.sol(emission)[2]
            def arrival(te):
                target_tau=sol.sol(te)[2]+distance
                return brentq(lambda t: sol.sol(t)[2]-target_tau,
                              te, 5*scale, xtol=1e-13*scale)
            arrival1 = arrival(emission)
            # Finite emission separation defined by the emitter's own clock.
            delta_tau = .001*scale
            emission2 = brentq(lambda t: sol.sol(t)[2]-tau_e-delta_tau,
                              emission,5*scale,xtol=1e-13*scale)
            arrival2 = arrival(emission2)
            n_e = sol.sol(emission)[0]
            n_o = sol.sol(arrival1)[0]
            # Independent neighboring-ray derivative in reference time.
            dt = 1e-5*scale
            derivative = (arrival(emission+dt)-arrival(emission))/dt
            s_reference = n_o/n_e
            # Local photon E=H/q=P and identical local atomic standards.
            s_local_frequency = s_reference*(n_e/n_o)
            s_local_event = (sol.sol(arrival2)[2]-sol.sol(arrival1)[2])/delta_tau
            signals.append(dict(matter_to_radiation=matter,inertia=inertia,
                                emission_reference_time=emission,distance=distance,
                                reference_redshift=s_reference-1,
                                local_observed_redshift=s_local_frequency-1,
                                local_event_stretch=s_local_event,
                                reference_arrival_derivative=derivative,
                                reference_frequency_stretch=s_reference,
                                clock_travel_time=float(sol.sol(arrival1)[2]-tau_e)))
    return rows, signals, dict(matter_to_radiation=matter,inertia=inertia,
                              normalized_energy_error=energy_error,
                              analytic_time_error_in_scaled_units=inverse_time_error)


def main():
    rows=[]; signals=[]; checks=[]; refinements=[]
    for matter in [0., .1, 1., 10., 1000.]:
        for inertia in [.5, 2.]:
            a,b,c=background(matter,inertia,2e-10)
            aa,bb,cc=background(matter,inertia,2e-12)
            rows.extend(aa); signals.extend(bb); checks.append(cc)
            refinements.append(max(abs(x['reference_redshift']-y['reference_redshift'])
                                   for x,y in zip(b,bb)))
    summary=dict(
        status='Synthetic universal homogeneous closure fails to produce observed redshift or event stretch',
        backgrounds=10,energy_samples=len(rows),signal_pairs=len(signals),
        max_normalized_energy_error=max(c['normalized_energy_error'] for c in checks),
        max_analytic_time_error=max(c['analytic_time_error_in_scaled_units'] for c in checks),
        max_local_frequency_redshift=max(abs(s['local_observed_redshift']) for s in signals),
        max_local_event_stretch_error=max(abs(s['local_event_stretch']-1) for s in signals),
        reference_redshift_range=[min(s['reference_redshift'] for s in signals),
                                 max(s['reference_redshift'] for s in signals)],
        max_reference_derivative_relative_error=max(abs(s['reference_arrival_derivative']/s['reference_frequency_stretch']-1) for s in signals),
        max_tolerance_refinement_redshift_difference=max(refinements),
        astronomical_data_fitted=False,heldout_scores_opened=False,
        total_astrophysical_source_budget_tested=False)
    assert summary['max_normalized_energy_error']<1e-9
    assert summary['max_analytic_time_error']<1e-8
    assert summary['max_local_event_stretch_error']<1e-7
    assert summary['max_reference_derivative_relative_error']<1e-5
    assert max(abs(r['fraction_gain_supplied_by_radiation']-1/(1+r['matter_to_radiation'])) for r in rows)<1e-8
    assert summary['max_tolerance_refinement_redshift_difference']<1e-8
    # Partial clock response is a separate stipulated diagnostic, not a change
    # to the autonomous runs above: q_clock=n**(-beta), fixed spatial rulers.
    # Can reducing universality recover a shift while retaining local c?
    alpha = .0002488993286382367
    distance_Mpc = 30.660139
    target_stretch = float(np.exp(alpha*distance_Mpc))
    partial=[]
    for beta in [0., .25, .5, .75, .9, .99, 1.]:
        if beta == 1:
            partial.append(dict(beta=beta,target_redshift=target_stretch-1,
                                finite_solution=False,n_ratio=None,
                                local_speed_ratio=None))
        else:
            ratio = target_stretch**(1/(1-beta))
            c_ratio = ratio**(beta-1)
            assert abs(c_ratio-1/target_stretch)<1e-14
            partial.append(dict(beta=beta,target_redshift=target_stretch-1,
                                finite_solution=True,n_ratio=ratio,
                                local_speed_ratio=c_ratio))
    result=dict(summary=summary,energy_samples=rows,signal_pairs=signals,checks=checks,
                partial_clock_response=dict(
                    status='Conditional algebraic diagnostic; no microscopic atomic model',
                    distance_Mpc=distance_Mpc,alpha_per_Mpc=alpha,
                    target_origin='Earlier fitted conversion law at approximately 100 million light-years; not a new observation',
                    fixed_spatial_rulers=True,rows=partial))
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
