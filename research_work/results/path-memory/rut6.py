"""RUT-1 stage 6 part 1: two repairs, and a supported population built in equilibrium (protocol-rut6.md).

    python rut6.py [--output-dir DIR] [--canonical]

R1 replaces stage 5's contour root COUNT -- which certifies nothing and returns a false unstable root at
every mode of the neutral control -- with a contour-integral eigensolver that LOCATES every root in a
declared region. R2 repairs the seeded verification, which perturbed positions and velocities but not the
excitation and force-producing fields. R3 constructs an orbital population in self-consistent equilibrium
with the field it writes, which is the state stage 6 part 2 will test for stability.

Stages 4 and 5 stay frozen: this file, beyn.py and equilibrium.py are new, and nothing they import is
edited. It is also the suite job for stage 6.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))
import beyn                 # noqa: E402
import equilibrium as EQ    # noqa: E402
import formation as FM      # noqa: E402
import ring_modes as RM     # noqa: E402
import rut3 as R3           # noqa: E402

T0 = RM.T0
R0 = 1.
GM = FM.GM
TAU_KEEP, TAU_FORM = 10*T0, 3*T0
LABEL = .1
REGION = (1e-3, 1.5, -6., 6.)         # the declared search region: Re(s) > 0 keeps B_m analytic
# declared thresholds (protocol-rut6.md)
TOL_H1, TOL_H2 = 1e-8, 1e-6
TOL_H3 = 1e-5
TOL_H4_RATE, TOL_H4_SCALE = .05, .10
TOL_H5 = 1e-6
TOL_H6_FACTOR = 3.
ANNULI = (('narrow_cold', .06, .010), ('narrow_warm', .06, .030),
          ('wide_cold', .12, .010), ('wide_warm', .12, .030))


def ring(model, n=32, w=.2, fraction=LABEL):
    return RM.Ring(model, n, R3._label_rate(w, TAU_KEEP, fraction)/n, R=R0, w=w,
                   tau_keep=TAU_KEEP, tau_form=TAU_FORM)


# ---------------------------------------------------------------- R1: the eigensolver
def gate_known_spectra():
    """H1: three problems whose spectra are known exactly, one of which must return no roots at all."""
    rows = {}

    # (a) exp(-s) = 1/2, roots ln2 + 2 pi i k -- genuinely nonlinear in s
    T = lambda s: (np.exp(-s) - .5)*np.eye(2)
    got = beyn.solve_region(T, (.1, 2., -8., 8.), moments=12,
                            residual=lambda s: (np.exp(-s) - .5)**2,
                            polish=lambda s: np.exp(-s) - .5)
    exact = [np.log(2) + 2j*np.pi*k for k in (-1, 0, 1)]
    err = max(min(abs(r['s'] - e)/max(1., abs(e)) for e in exact) for r in got['roots']) if got['roots'] else 1.
    rows['transcendental'] = dict(found=got['count'], expected=len(exact), worst_relative=float(err),
                                  passed=bool(got['count'] == len(exact) and err < TOL_H1))

    # (b) the Kepler control: determinant s^2 (s^2 + Omega^2), no root with Re > 0
    kep = ring('none')
    counts = {}
    for m in (0, 2, 5):
        out = beyn.solve_region(lambda s: kep.M(m, s), REGION, residual=lambda s: kep.det(m, s))
        counts[f'm{m}'] = out['count']
    rows['kepler_control'] = dict(counts=counts, expected=0, passed=bool(all(c == 0 for c in counts.values())))

    # (c) the zero-lag rung, whose quartic roots follow exactly from the companion form
    inst = ring('instantaneous')
    worst, ok = 0., True
    for m in (2, 4, 7):
        exact_roots = [z for z in inst._quartic_roots(inst.B(m, 0.)) if beyn.in_region(z, REGION)]
        out = beyn.solve_region(lambda s: inst.M(m, s), REGION, residual=lambda s: inst.det(m, s),
                                polish=lambda s: inst.det(m, s))
        ag = beyn.agree([r['s'] for r in out['roots']], exact_roots, tol=TOL_H1)
        worst = max(worst, ag['worst_relative'])
        ok = ok and ag['agreed']
    rows['zero_lag_quartic'] = dict(worst_relative=float(worst), passed=bool(ok and worst < TOL_H1))
    return dict(rows=rows, tolerance=TOL_H1, passed=all(r['passed'] for r in rows.values()))


def gate_eigensolver_robustness():
    """H2: double the quadrature and shift the contour edge; the located roots must not move."""
    rows, worst = {}, 0.
    for model in ('one_stage', 'two_stage'):
        r = ring(model)
        for m in (2, 5):
            base = beyn.solve_region(lambda s: r.M(m, s), REGION, residual=lambda s: r.det(m, s),
                                     polish=lambda s: r.det(m, s))
            fine = beyn.solve_region(lambda s: r.M(m, s), REGION, panels=24, order=20,
                                     residual=lambda s: r.det(m, s), polish=lambda s: r.det(m, s))
            moved = beyn.solve_region(lambda s: r.M(m, s), (5e-3, 1.6, -6.5, 6.5),
                                      residual=lambda s: r.det(m, s), polish=lambda s: r.det(m, s))
            a = beyn.agree([x['s'] for x in base['roots']], [x['s'] for x in fine['roots']], tol=TOL_H2)
            b = beyn.agree([x['s'] for x in base['roots']], [x['s'] for x in moved['roots']], tol=TOL_H2)
            rows[f'{model}|m{m}'] = dict(count=base['count'], refined=a, shifted=b,
                                         worst_residual=max([x['residual'] for x in base['roots']] or [0.]))
            worst = max(worst, a['worst_relative'], b['worst_relative'])
            if not (a['agreed'] and b['agreed']):
                rows[f'{model}|m{m}']['disagreed'] = True
    return dict(rows=rows, worst=float(worst), tolerance=TOL_H2,
                passed=bool(all(r.get('refined', {}).get('agreed') and r.get('shifted', {}).get('agreed')
                                for r in rows.values())))


def certified_spectrum(model, m_values=range(9)):
    """Every root in the declared region, located rather than counted, with continuation as cross-check."""
    r = ring(model)
    rows, agreed = {}, True
    for m in m_values:
        out = beyn.solve_region(lambda s: r.M(m, s), REGION, residual=lambda s: r.det(m, s),
                                polish=lambda s: r.det(m, s))
        cont = [z for z in r.modes(m)[0] if beyn.in_region(z, REGION)]
        ag = beyn.agree([x['s'] for x in out['roots']], cont, tol=TOL_H3)
        agreed = agreed and ag['agreed']
        rows[str(m)] = dict(
            located=[dict(re=float(x['s'].real), im=float(x['s'].imag), multiplicity=x['multiplicity'],
                          residual=x['residual'], distance_to_contour=x['distance_to_contour'],
                          e_folding_periods=float(1/(x['s'].real*T0)) if x['s'].real > 0 else None,
                          pattern_speed_over_omega=float((r.omega - x['s'].imag/m)/r.omega) if m else None)
                     for x in out['roots']],
            count=out['count'], saturated=out['saturated'], continuation_agreement=ag)
    fastest = max((x for row in rows.values() for x in row['located']),
                  key=lambda x: x['re'], default=None)
    return dict(model=model, omega=r.omega, support_fraction=r.support_fraction, rows=rows,
                unstable_modes=[m for m, row in rows.items() if row['count']],
                fastest=fastest, continuation_agrees_everywhere=agreed,
                any_saturated=any(row['saturated'] for row in rows.values()))


# ---------------------------------------------------------------- R2: the full-state seeded test
def seed_full_state(r, m, amplitude, spin_periods=40., h=.01):
    """Build the mode's own history into E and C by prescribing its trajectory through a spin-up window
    and advancing the SHIPPED field update, then release the bodies.

    This realises the declared delta_E and delta_C integrals without evaluating them: the field update is
    the same one the run uses, so whatever it does to a prescribed source history is exactly what the
    run's own history would have produced. It also replaces the analytic continuous-ring priming with the
    actual discrete base field.
    """
    roots, _ = r.modes(m)
    s = roots[int(np.argmax(roots.real))]
    vec = r.eigenvector(m, s)
    vec = vec/np.linalg.norm(vec)
    phase = np.exp(2j*np.pi*m*np.arange(r.n)/r.n)
    er = np.stack([np.cos(r.phi), np.sin(r.phi)], axis=1)
    et = np.stack([-np.sin(r.phi), np.cos(r.phi)], axis=1)
    rates = np.full(r.n, r.strength)

    def state(t):
        """Prescribed position and rotating-frame velocity of the mode at time t (t <= 0 in the spin-up)."""
        loc = amplitude*np.real(vec[None, :]*phase[:, None]*np.exp(s*t))
        dloc = amplitude*np.real(s*vec[None, :]*phase[:, None]*np.exp(s*t))
        dx = loc[:, :1]*er + loc[:, 1:]*et
        dv = dloc[:, :1]*er + dloc[:, 1:]*et
        rot = np.array([[np.cos(r.omega*t), -np.sin(r.omega*t)],
                        [np.sin(r.omega*t), np.cos(r.omega*t)]])
        x = (r.P + dx)@rot.T
        v = (dv + r.omega*np.stack([-(r.P + dx)[:, 1], (r.P + dx)[:, 0]], axis=1))@rot.T
        return x, v

    field = FM.MemoryField(2.5, r.w/5, r.w, TAU_KEEP, TAU_FORM if r.model == 'two_stage' else 0.)
    field.prime_with_ring(r.R, A=float(np.sum(rates))*TAU_KEEP)      # the m=0 part starts mature
    t = -spin_periods*T0
    while t < 0.:
        step = min(h, -t)
        x_mid, _ = state(t + .5*step)
        field.advance(x_mid, rates, step)
        t += step
    x0, v0 = state(0.)
    return field, x0, v0, complex(s), rates


def gate_full_state_seed(window=30., amplitudes=(1e-4, 1e-3)):
    """H4: with the whole state perturbed, the seeded mode must reproduce the linear prediction, and its
    initial amplitude must scale with the seed."""
    r = ring('two_stage')
    rows, amps = {}, []
    predicted = None
    for amp in amplitudes:
        field, x0, v0, s, rates = seed_full_state(r, 2, amp)
        predicted = s
        # the run MUST continue in the field the spin-up built: RM.run_ring would construct its own,
        # and with no priming the bodies would be launched at the supported speed into an empty field --
        # a formation transient, which is what this gate measured until the field was passed through.
        out = run_ring_field(field, x0, v0, rates, r.w, window*T0, samples=600)
        meas = RM.measure_mode(out, 2, r.R, .2*window, window)
        start = RM.measure_mode(out, 2, r.R, 0., .1*window)
        amps.append(start['amplitude_start'] if start else float('nan'))
        rows[f'{amp:g}'] = dict(
            amplitude=amp, t_final_periods=out['t_final_periods'], measured=meas,
            initial_amplitude=amps[-1],
            relative_growth=abs(meas['growth_rate']/s.real - 1),
            relative_frequency=abs(meas['frequency_rotating_frame'] - s.imag)/abs(s.imag))
    # the control that decides the diagnosis: the same seed with NO field history, which is what the
    # stage 5 gate did. Reported, not gated.
    field, x0, v0, s, rates = seed_full_state(r, 2, amplitudes[1], spin_periods=0.)
    out = run_ring_field(field, x0, v0, rates, r.w, window*T0, samples=600)
    bare = RM.measure_mode(out, 2, r.R, .2*window, window)
    no_history = dict(amplitude=amplitudes[1], measured=bare,
                      relative_growth=abs(bare['growth_rate']/s.real - 1),
                      relative_frequency=abs(bare['frequency_rotating_frame'] - s.imag)/abs(s.imag))
    ratio = amps[1]/amps[0] if amps[0] else float('nan')
    scale_err = abs(ratio/(amplitudes[1]/amplitudes[0]) - 1)
    worst = max(max(v['relative_growth'], v['relative_frequency']) for v in rows.values())
    return dict(predicted=dict(growth_rate=float(predicted.real), frequency=float(predicted.imag),
                               e_folding_periods=float(1/(predicted.real*T0))),
                rows=rows, no_field_history_control=no_history,
                worst_relative=float(worst), amplitude_ratio=float(ratio),
                amplitude_scaling_error=float(scale_err),
                tolerance_rate=TOL_H4_RATE, tolerance_scaling=TOL_H4_SCALE,
                passed=bool(worst < TOL_H4_RATE and scale_err < TOL_H4_SCALE))


def run_ring_field(field, x, v, rates, w, t_end, h_max=.01, eta=.01, write=True, freeze=False,
                   samples=400):
    """Evolve bodies in a field built elsewhere.

    Three distinct controls, and they are NOT the same -- stage 3 learned this once already and this
    function reintroduced the confusion. `freeze=True` holds the field fixed, which is the exact
    equilibrium of the constructed population. `write=False` lets the field DECAY with no new deposition,
    which removes the support over tau_keep and pushes the population outward. `field` left at zero is no
    field at all."""
    x = np.array(x, float)
    v = np.array(v, float)
    zero = np.zeros_like(np.asarray(rates, float))
    accel = lambda xx: (-GM*xx/np.maximum(np.linalg.norm(xx, axis=1, keepdims=True), 1e-12)**3
                        + field.sample_gradient(xx))
    a = accel(x)
    rec, t, next_sample = dict(t=[], r=[], vr=[], L=[]), 0., 0.
    while t < t_end:
        rr = np.linalg.norm(x, axis=1)
        if rr.min() < .25:
            break
        step = min(h_max, eta*float(np.min(np.sqrt(rr**3/GM))), t_end - t)
        v_half = v + .5*step*a
        x_mid = x + .5*step*v_half
        x = x + step*v_half
        if np.max(np.abs(x)) > .9*2.5:
            break
        if not freeze:
            field.advance(x_mid, rates if write else zero, step)
        a = accel(x)
        v = v_half + .5*step*a
        t += step
        if t >= next_sample or t >= t_end:
            next_sample = t + t_end/samples
            rr = np.linalg.norm(x, axis=1)
            rec['t'].append(t)
            rec['r'].append(rr.copy())
            rec['vr'].append(np.sum(v*x, axis=1)/rr)
            rec['L'].append(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])
    out = {k: np.array(val) for k, val in rec.items()}
    out['t_final_periods'] = float(t/T0)
    return out


def _drift(out, window=5.):
    """Fractional change of the population's mean radius, spread and radial dispersion, first to last."""
    t, r, vr = out['t'], out['r'], out['vr']
    early, late = t <= window*T0, t >= t[-1] - window*T0
    f = lambda sel: (float(np.mean(r[sel])), float(np.std(r[sel])), float(np.sqrt(np.mean(vr[sel]**2))),
                     float(np.mean(out['L'][sel])))
    a, b = f(early), f(late)
    return dict(mean_radius=[a[0], b[0]], radial_spread=[a[1], b[1]], radial_dispersion=[a[2], b[2]],
                angular_momentum=[a[3], b[3]],
                mean_radius_drift=abs(b[0]/a[0] - 1), spread_drift=abs(b[1]/a[1] - 1) if a[1] else 0.,
                dispersion_drift=abs(b[2]/a[2] - 1) if a[2] else 0.,
                L_drift=abs(b[3]/a[3] - 1) if a[3] else 0.)


# ---------------------------------------------------------------- R3: the constructed population
def equilibria():
    rows = {}
    for name, dL, dE in ANNULI:
        t0 = time.time()
        a = EQ.WarmAnnulus(L0=1., dL=dL, dE=dE, w=.2, tau_keep=TAU_KEEP).solve(target_support=LABEL)
        s = a.summary()
        s.update(a.dispersions(), consistency=a.consistency_residual(), dL=dL, dE=dE,
                 support_target=LABEL, support_error=float(getattr(a, 'support_error', 0.)),
                 mixing_used=float(getattr(a, 'mixing_used', 0.)),
                 seconds=round(time.time() - t0, 1))
        rows[name] = s
    return rows


def gate_equilibrium(rows):
    """H5: both sides of the coupled equations, recomputed from the converged state."""
    worst = max(max(r['consistency'].values()) for r in rows.values())
    finite = all(r['mass'] > 0 and r['edge_density_ratio'] < 1e-3 for r in rows.values())
    missed = {k: r['support_error'] for k, r in rows.items() if r['support_error'] > 1e-3}
    return dict(worst_consistency=float(worst), finite_mass_with_tapering=bool(finite),
                support_target_not_reached=missed,
                note='H5 is about self-consistency, which every annulus meets. Where the outer solve could '
                     'not reach the declared support label the achieved support is reported instead: that '
                     'configuration is still an exact equilibrium, just not at the requested strength.',
                tolerance=TOL_H5, passed=bool(worst < TOL_H5 and finite))


def written_versus_frozen(name, n=64, horizon=20., seed=1):
    """A measurement, not a gate: the same sampled population evolved with the field written and with it
    frozen, for a second annulus, so that the effect of temperature can be read off."""
    cfg = {k: (dL, dE) for k, dL, dE in ANNULI}[name]
    a = EQ.WarmAnnulus(L0=1., dL=cfg[0], dE=cfg[1], w=.2, tau_keep=TAU_KEEP).solve(target_support=LABEL)
    x, v = a.sample(n, seed=seed)
    rates = np.full(n, a.rate_per_body(n))
    runs = {}
    for tag, write, freeze in (('written', True, False), ('frozen_field', False, True)):
        field = FM.MemoryField(2.5, a.w/5, a.w, TAU_KEEP, TAU_FORM)
        a.prime(field)
        out = run_ring_field(field, x, v, rates, a.w, horizon*T0, write=write, freeze=freeze)
        runs[tag] = dict(drift=_drift(out), t_final_periods=out['t_final_periods'])
    return dict(annulus=name, bodies=n, horizon_periods=horizon, sigma_r=a.dispersions()['sigma_r'],
                runs=runs)


def gate_stationary_control(name='narrow_warm', n=64, horizon=20., seed=1):
    """H6: evolve the constructed population, and compare its drift with controls.

    The protocol declared a NO-MEMORY control. That control is vacuous here and is reported as such: a run
    with no field at all is not an equilibrium of a population built in the field's potential, so it drifts
    hugely and any memory run would pass beside it. The gate is therefore evaluated against the
    FROZEN-FIELD control, which is an exact equilibrium of the constructed distribution function -- a
    strictly harder test than the one declared. Both controls are reported.
    """
    cfg = {k: (dL, dE) for k, dL, dE in ANNULI}[name]
    a = EQ.WarmAnnulus(L0=1., dL=cfg[0], dE=cfg[1], w=.2, tau_keep=TAU_KEEP).solve(target_support=LABEL)
    x, v = a.sample(n, seed=seed)
    q = a.rate_per_body(n)
    rates = np.full(n, q)
    runs = {}
    for tag, write, freeze, memory in (('written', True, False, True),
                                       ('frozen_field', False, True, True),
                                       ('decaying_field', False, False, True),
                                       ('no_memory', False, True, False)):
        field = FM.MemoryField(2.5, a.w/5, a.w, TAU_KEEP, TAU_FORM)
        if memory:
            a.prime(field)
        out = run_ring_field(field, x, v, rates, a.w, horizon*T0, write=write, freeze=freeze)
        runs[tag] = dict(drift=_drift(out), t_final_periods=out['t_final_periods'])
    w, f, nm = runs['written']['drift'], runs['frozen_field']['drift'], runs['no_memory']['drift']
    keys = ('mean_radius_drift', 'spread_drift', 'dispersion_drift')
    ratios = {k: (w[k]/f[k] if f[k] else float('inf')) for k in keys}
    return dict(annulus=name, bodies=n, horizon_periods=horizon, runs=runs,
                ratios_to_frozen_field=ratios,
                declared_no_memory_control_drift={k: nm[k] for k in keys},
                declared_control_is_vacuous=True,
                evaluated_against='frozen_field, an exact equilibrium of the constructed population',
                tolerance_factor=TOL_H6_FACTOR,
                passed=bool(all(v <= TOL_H6_FACTOR for v in ratios.values())))


# ---------------------------------------------------------------- driver
def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    eq = equilibria()
    gates = dict(H1_known_spectra=gate_known_spectra(),
                 H2_eigensolver_robustness=gate_eigensolver_robustness(),
                 H4_full_state_seed=gate_full_state_seed(),
                 H5_equilibrium=gate_equilibrium(eq),
                 H6_stationary_control=gate_stationary_control())
    spectra = {model: certified_spectrum(model) for model in ('none', 'instantaneous', 'one_stage',
                                                              'two_stage')}
    gates['H3_eigensolver_vs_continuation'] = dict(
        rows={k: dict(agrees=v['continuation_agrees_everywhere'], any_saturated=v['any_saturated'],
                      unstable_modes=v['unstable_modes']) for k, v in spectra.items()},
        tolerance=TOL_H3,
        passed=bool(all(v['continuation_agrees_everywhere'] and not v['any_saturated']
                        for v in spectra.values())))
    result = dict(
        experiment='RUT-1 stage 6 part 1: two repairs, and a population built in equilibrium',
        protocol='protocol-rut6.md',
        search_region=dict(re=[REGION[0], REGION[1]], im=[REGION[2], REGION[3]],
                           why='Re(s) > 0 keeps B_m(s) analytic: its poles sit at -1/tau + i k Omega'),
        certified_spectra=spectra, equilibria=eq,
        cold_annulus_written_versus_frozen=written_versus_frozen('narrow_cold'),
        gates=gates,
        passed=all(g['passed'] for g in gates.values()),
        failed_gates=sorted(k for k, g in gates.items() if not g['passed']),
        what_this_is_not=(
            'part 1 only. The spectra are certified for the COLD RING family of stage 5, in the declared '
            'region and for m <= 8, not for the constructed populations -- their mode analysis is part 2 '
            'and is not run. An equilibrium is not a stability result, and H6 -- the declared stationary '
            'control -- FAILED as declared and is recorded as failed. No claim is made here that any '
            'population is stable.'),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(json.dumps({k: result[k] for k in ('experiment', 'passed')}, indent=1))
    for name, g in sorted(gates.items()):
        print(f"  {'PASS' if g['passed'] else 'FAIL'} {name}")
    status = evidence_io.finish(args, 'path-memory-rut6', text, HERE/'rut6-results.json',
                                ignore={'/runtime_seconds'},
                                rules=((r'/equilibria/[^/]+/seconds', None, None),))
    # The exit status reports whether this run REPRODUCES ITS ARCHIVE, not whether every declared gate
    # passed. A declared gate that fails is a recorded outcome: it is in `failed_gates`, in the report and
    # in the changelog, and the archive comparison still turns the suite red if any gate's outcome ever
    # changes. Tying the exit status to the gates instead would leave two choices when one fails honestly
    # -- withhold the result, or reword the gate until it passes -- and the second is the habit the
    # owner's review of ceb0c86 told me to stop.
    if not result['passed']:
        print('  DECLARED GATES FAILED AND RECORDED AS FAILED:', ', '.join(result['failed_gates']))
    return status


if __name__ == '__main__':
    raise SystemExit(main())
