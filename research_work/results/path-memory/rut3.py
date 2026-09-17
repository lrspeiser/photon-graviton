"""RUT-1 stage 3: formation from an empty field, and the two-stage candidate (protocol-rut1.md).

    python rut3.py [--output-dir DIR] [--canonical]

Bodies and field are evolved together from zero, with no prescribed trajectory, no velocity resets, no
target speed and no writing rate adjusted to reach a chosen orbit. Two constitutive models are compared on
matched initial conditions -- same bodies, same total writing rate, same retention, same width, same
positions and velocities, ONLY the maturation process changing:

    one-stage:   dC/dt = S - C/tau_keep                      (the preserved baseline)
    two-stage:   tau_form dE/dt = S - E,  dC/dt = E - C/tau_keep

with Phi = -C and a_mem = grad C in both. The one-stage model is exactly tau_form = 0, and both have the
same steady state C = tau_keep S, so the comparison does not weaken gravity.

Support, stability, settling, heating and formation cost are reported separately, three termination
labels are distinguished, and unbinding is asserted only from a total energy that includes the memory
potential.

Regenerates rut3-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))
import formation as FM  # noqa: E402
import rut1 as U  # noqa: E402

GM = 1.
R0, PERIOD = 1., 2*np.pi
OMEGA = 1.


def _label_rate(w, tau, fraction):
    """The writing rate the mature-ring calculation predicts for `fraction` of Newtonian: a LABEL only."""
    sup, _ = U.self_force(R0, w, tau, OMEGA)
    return float(fraction*(GM/R0**2)/abs(sup))


# ---------------------------------------------------------------- gates
def gate_close_passage(pericentres=(.11, .063, .03), steps=(.01, .002, .0005), eta=.01):
    """A fixed step mishandles close passages while conserving angular momentum; the adaptive step does
    not. The no-memory circular control could never have caught this, which is why it is replaced."""
    def kepler(rp, h=None, eta_=None, n_orb=10.):
        a, e = 1., 1 - rp
        r = a*(1 + e)
        x, v = np.array([r, 0.]), np.array([0., np.sqrt(GM*(2/r - 1/a))])
        acc = lambda p: -GM*p/np.linalg.norm(p)**3
        E0, L0 = .5*v@v - GM/r, x[0]*v[1] - x[1]*v[0]
        A, t, T = acc(x), 0., n_orb*2*np.pi
        while t < T:
            rr = np.linalg.norm(x)
            step = h if h is not None else min(.01, eta_*np.sqrt(rr**3/GM))
            step = min(step, T - t)
            vh = v + .5*step*A
            x = x + step*vh
            A = acc(x)
            v = vh + .5*step*A
            t += step
        rr = np.linalg.norm(x)
        return dict(energy_error=float(abs((.5*v@v - GM/rr)/E0 - 1)),
                    angular_momentum_error=float(abs((x[0]*v[1] - x[1]*v[0])/L0 - 1)),
                    spuriously_unbound=bool(.5*v@v - GM/rr > 0))
    fixed = [dict(pericentre=rp, h=h, **kepler(rp, h=h)) for rp in pericentres for h in steps]
    adaptive = [dict(pericentre=rp, eta=eta, **kepler(rp, eta_=eta)) for rp in pericentres]
    worst_fixed = max(r['energy_error'] for r in fixed if r['h'] == max(steps))
    worst_adaptive = max(r['energy_error'] for r in adaptive)
    return dict(fixed_step=fixed, adaptive_step=adaptive,
                worst_energy_error_at_coarsest_fixed=worst_fixed,
                worst_energy_error_adaptive=worst_adaptive,
                angular_momentum_is_no_guide=float(max(r['angular_momentum_error'] for r in fixed)),
                improvement_over_coarsest_fixed=float(worst_fixed/max(worst_adaptive, 1e-300)),
                tolerance=1e-3, passed=bool(worst_adaptive < 1e-3 and worst_fixed/worst_adaptive > 1e3),
                statement='at pericentre 0.063, the radius RUT-1 plunges reached, fixed h = 0.01 carries a '
                          '4.4% energy error while its angular momentum stays good to 1e-14; at 0.03 it '
                          'throws a bound orbit out and calls it unbound. The adaptive step removes that')


def gate_two_state_ode(tau_keep=10*PERIOD, tau_form=3*PERIOD, S=.7, h=.05, n=400):
    """The exact-per-step two-state update against an independently integrated ODE at constant source."""
    # drive the REAL update with a single writer parked on the grid centre, so the gate tests the
    # shipped code path rather than a re-implementation of it
    w = .1
    fld = FM.MemoryField(.4, .05, w, tau_keep, tau_form)
    q = S/1.                                        # the Gaussian peaks at q where the writer sits
    for _ in range(n):
        fld.advance([[0., 0.]], [q], h)
    got = float(fld.C[0, fld.n//2, fld.n//2])
    sol = solve_ivp(lambda t, y: [(S - y[0])/tau_form, y[0] - y[1]/tau_keep], (0, n*h), [0., 0.],
                    rtol=1e-12, atol=1e-14)
    ref = float(sol.y[1, -1])
    return dict(stepped=got, ode=ref, absolute_difference=float(abs(got - ref)),
                relative_difference=float(abs(got/ref - 1)),
                steady_state_C=float(tau_keep*S), tolerance=1e-10,
                passed=bool(abs(got/ref - 1) < 1e-10),
                note='the steady state tau_keep*S is the same as the one-stage model, so the two-stage '
                     'response does not weaken the equilibrium field')


def gate_transfer_function(tau_keep=10*PERIOD, tau_form=3*PERIOD):
    """H_new(w) = H_old(w)/(1 + i w tau_form): zero frequency untouched, fast variation attenuated."""
    w = np.array([0., .1, 1., 10.])/PERIOD
    old = tau_keep/(1 + 1j*w*tau_keep)
    new = old/(1 + 1j*w*tau_form)
    return dict(omega_per_period=(w*PERIOD).tolist(),
                magnitude_ratio=np.abs(new/old).tolist(),
                zero_frequency_unchanged=float(abs(new[0]/old[0] - 1)),
                passed=bool(abs(new[0]/old[0] - 1) < 1e-14),
                caveat='a statement about one harmonic of a prescribed writing pattern, NOT a formula for '
                       'the drag of a freely evolving collective, where transient torques of either sign '
                       'stay in the accounting')


# ---------------------------------------------------------------- the matched comparison
def _run(n, wr, tau_form, orbits, fraction=.1, phase_jitter=0., speed_jitter=0., seed=0,
         h_max=.01, spacing=None, memory=True, prime=False, write=True, freeze=False, eta=.01):
    w, tau = wr*R0, 10*PERIOD
    q = _label_rate(w, tau, fraction)
    extra = U.ring_inward_acceleration(R0, R0, w, q*tau) if prime else 0.
    pos, vel, rates = FM.ring_start(n, R0, q, phase_jitter, speed_jitter, seed, extra_inward=extra)
    out = FM.run(pos, vel, rates, w, tau, t_end=orbits*PERIOD, h_max=h_max, tau_form=tau_form,
                 spacing=spacing, memory=memory, prime_ring_R=R0 if prime else None,
                 write=write, freeze_field=freeze, eta=eta)
    s = FM.summarize(out, R0, PERIOD)
    s.update(writers=n, width_ratio=wr, tau_form_periods=tau_form/PERIOD, seed=seed,
             phase_jitter=phase_jitter, speed_jitter=speed_jitter, model='two_stage' if tau_form else
             'one_stage', nominal_support_label=fraction)
    return s


def matched_comparison(orbits=20., jitter=.02, seeds=(1, 2), tau_form_periods=3.):
    """Same bodies, same total rate, same retention, same width, same ICs. Only maturation changes."""
    tf = tau_form_periods*PERIOD
    rows = []
    for n, wr, use in ((16, .1, seeds[:1]), (32, .2, seeds)):
        for seed in use:
            common = dict(n=n, wr=wr, orbits=orbits, phase_jitter=jitter, speed_jitter=jitter, seed=seed)
            one = _run(tau_form=0., **common)
            two = _run(tau_form=tf, **common)
            rows.append(dict(
                writers=n, width_ratio=wr, seed=seed, jitter=jitter,
                one_stage={k: one[k] for k in ('status', 'orbits', 'support_fraction_mean_late',
                                               'angular_momentum_final_over_initial',
                                               'radius_change_fraction', 'radial_velocity_rms_late',
                                               'per_body_L_ratio_spread', 'any_body_positive_energy')},
                two_stage={k: two[k] for k in ('status', 'orbits', 'support_fraction_mean_late',
                                               'angular_momentum_final_over_initial',
                                               'radius_change_fraction', 'radial_velocity_rms_late',
                                               'per_body_L_ratio_spread', 'any_body_positive_energy')}))
    n_one = sum(r['one_stage']['status'] == FM.COMPLETED for r in rows)
    n_two = sum(r['two_stage']['status'] == FM.COMPLETED for r in rows)
    return dict(rows=rows, orbits=orbits, jitter=jitter, tau_form_periods=tau_form_periods,
                completed_one_stage=n_one, completed_two_stage=n_two, cases=len(rows),
                two_stage_completes_more=bool(n_two > n_one),
                note='identical initial conditions under both equations; the writing rate, retention and '
                     'width are held fixed, so any difference is the maturation process')


def split_primed_controls(orbits=12.):
    """Frozen field / decaying without writing / decaying with writing, for one and for sixteen writers.

    Correction 10: the earlier single primed control could not separate field shape from writing, because
    it let the old field decay while the body added new uneven structure.
    """
    out = []
    for n in (1, 16):
        base = dict(n=n, wr=.1, tau_form=0., orbits=orbits, prime=True)
        out.append(dict(writers=n, control='frozen_field', **{k: v for k, v in
                                                              _run(freeze=True, **base).items()
                                                              if k in ('status', 'orbits',
                                                                       'angular_momentum_final_over_initial',
                                                                       'radius_change_fraction')}))
        out.append(dict(writers=n, control='decaying_no_writing', **{k: v for k, v in
                                                                     _run(write=False, **base).items()
                                                                     if k in ('status', 'orbits',
                                                                              'angular_momentum_final_over_initial',
                                                                              'radius_change_fraction')}))
        out.append(dict(writers=n, control='decaying_with_writing', **{k: v for k, v in
                                                                       _run(**base).items()
                                                                       if k in ('status', 'orbits',
                                                                                'angular_momentum_final_over_initial',
                                                                                'radius_change_fraction')}))
    return dict(rows=out,
                note='comparing the three isolates whether the trouble is the field already present, its '
                     'decay, or the new uneven structure a moving body adds')


def convergence_on_a_disturbed_run(orbits=10., jitter=.02, seed=1, tau_form_periods=3.):
    """Refined on the DISTURBED case that carries the claim, not on the symmetric control."""
    tf = tau_form_periods*PERIOD
    base = dict(n=16, wr=.1, tau_form=tf, orbits=orbits, phase_jitter=jitter, speed_jitter=jitter,
                seed=seed)
    ref = _run(**base)
    fine_t = _run(h_max=.005, eta=.0025, **base)
    fine_x = _run(spacing=.1/10, **base)
    key = ('radius_change_fraction', 'angular_momentum_final_over_initial', 'support_fraction_mean_late')
    shifts = {k: dict(half_timestep=float(abs(ref[k] - fine_t[k])),
                      half_spacing=float(abs(ref[k] - fine_x[k]))) for k in key}
    worst = max(max(v.values()) for v in shifts.values())
    return dict(baseline={k: ref[k] for k in key}, shifts=shifts, worst_shift=worst,
                tolerance=5e-3, passed=bool(worst < 5e-3),
                note='the earlier convergence check used the unperturbed ring for ten periods, which never '
                     'approaches the centre and so could not test the close passages the breakup claim '
                     'depends on')


def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    close = gate_close_passage()
    ode = gate_two_state_ode()
    transfer = gate_transfer_function()
    matched = matched_comparison()
    primed = split_primed_controls()
    conv = convergence_on_a_disturbed_run()
    no_mem = _run(n=16, wr=.1, tau_form=0., orbits=20., phase_jitter=.02, speed_jitter=.02, seed=1,
                  memory=False)
    result = dict(
        experiment='RUT-1 stage 3: formation from an empty field, and the two-stage candidate',
        protocol='protocol-rut1.md',
        models=dict(
            one_stage='dC/dt = S - C/tau_keep (the preserved baseline)',
            two_stage='tau_form dE/dt = S - E,  dC/dt = E - C/tau_keep',
            shared='Phi = -C, a_mem = grad C; the one-stage model is exactly tau_form = 0 and both have '
                   'steady state C = tau_keep S, so the comparison does not weaken the equilibrium field',
            owed='no field energy or momentum budget, and a formation delay is not spatial causality'),
        gates=dict(close_passage=close, two_state_ode=ode, transfer_function=transfer,
                   convergence_on_a_disturbed_run=conv),
        matched_comparison=matched,
        split_primed_controls=primed,
        no_memory_control={k: no_mem[k] for k in ('status', 'angular_momentum_final_over_initial',
                                                  'radius_change_fraction')},
        termination_labels=dict(
            completed=FM.COMPLETED, left_domain=FM.LEFT_DOMAIN, unresolved=FM.UNRESOLVED_CENTRE,
            note='never interchanged. Unbinding is asserted only from a total energy including the memory '
                 'potential, which is why the scalar field is carried and not only its gradient'),
        passed=bool(close['passed'] and ode['passed'] and transfer['passed'] and conv['passed']
                    and no_mem['status'] == FM.COMPLETED
                    and abs(no_mem['angular_momentum_final_over_initial'] - 1) < 1e-9),
        what_this_is_not='one planar ring around a fixed centre, an instantaneous Gaussian kernel with no '
                         'spatial causality, and no energy or momentum accounting for the field. A run '
                         'that has not crossed a boundary but is heating is not a success',
        input_sha256={'protocol-rut1.md': hashlib.sha256((HERE/'protocol-rut1.md').read_bytes()).hexdigest()},
        checks_short_run=dict(
            close_passage_adaptive=close['worst_energy_error_adaptive'],
            two_state_ode=ode['absolute_difference'],
            convergence_worst=conv['worst_shift'],
            completed_one_stage=matched['completed_one_stage'],
            completed_two_stage=matched['completed_two_stage']),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:1600])
    status = evidence_io.finish(args, 'path-memory-rut3', text, HERE/'rut3-results.json',
                                ignore={'/runtime_seconds'})
    return status if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
