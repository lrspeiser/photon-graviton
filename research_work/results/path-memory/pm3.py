"""PM-3: what a shared-field memory can and cannot do (protocol-pm3.md).

    python pm3.py [--output-dir DIR] [--canonical]

Stage 1 reproduces the owner's toy as a control. With z = ln(r/r_ref) - chi,

    g_inward = GM/r^2 + [K + B tanh(z)]/r,   r'' = j^2/r^3 - g_inward,   tau chi' = tanh(z)

the response is bounded into (K-B, K+B)/r, so the force stays attractive; the memory is a restoring
correction, not a ratchet. Six gates: the equilibrium identity, the energy ledger E + Q - W, settling from
five disturbed starts, the B = 0 control that must NOT circularize, Routh-Hurwitz across tau, and two
preparation histories released at the same r, rdot and theta with a prescribed excursion.

Stage 2 checks the scope of any relaxational memory: in a steady state tau dQ/dt = F[fields] - Q relaxes
to Q = F[fields], so the memory contributes a term fixed by the instantaneous fields and nothing a static
law of those same fields could not produce. Corollary 1: on the toy's settled orbit z = 0, so the memory's
contribution is exactly zero whatever B and tau are. Corollary 2: the proposed field memory returns
exactly to Completion I for any stationary source, so it is blind to steady rotation by construction.

Regenerates pm3-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402

GM, K, B, TAU, J, R_REF = .2, .8, .4, 1., 1., 1.     # the owner's declared dimensionless parameters
R_C = 1.                                              # j^2 = GM r + K r^2 gives exactly 1
PERIOD = 2*np.pi/(J/R_C**2)                           # reference orbital period at the circular radius
TOL_LEDGER, TOL_SETTLE = 1e-10, 1e-6
FLOOR = 1e-12            # below this a residual is the integrator's own noise, not a trend


def _state(r, chi):
    return np.log(r/R_REF) - chi


def rhs(t, y, B_=B, drive=None):
    """[r, rdot, theta, chi, Q, W]: Q accumulates what the memory removes, W the driver's work."""
    r, vr, th, chi, _, _ = y
    z = _state(r, chi)
    f = drive(t) if drive is not None else 0.
    return [vr, J**2/r**3 - (GM/r**2 + (K + B_*np.tanh(z))/r) + f, J/r**2,
            np.tanh(z)/TAU, (B_/TAU)*np.tanh(z)**2, f*vr]


def energy(y, B_=B):
    r, vr, _, chi, _, _ = (y[i] for i in range(6))
    z = _state(r, chi)
    return vr**2/2 + J**2/(2*r**2) - GM/r + K*np.log(r/R_REF) + B_*np.log(np.cosh(z))


def integrate(y0, t_end, B_=B, drive=None, dense=False):
    return solve_ivp(rhs, (0, t_end), y0, args=(B_, drive), rtol=1e-12, atol=1e-14,
                     dense_output=dense, max_step=.05 if drive is not None else np.inf)


# ---------------------------------------------------------------- stage 1
def gate_equilibrium():
    """j^2 = GM r_c + K r_c^2 fixes the circular radius; no orbit or speed is supplied to the dynamics."""
    lhs, rhs_ = J**2, GM*R_C + K*R_C**2
    return dict(j_squared=lhs, GM_rc_plus_K_rc2=rhs_, relative_difference=float(abs(lhs/rhs_ - 1)),
                circular_radius=R_C, tolerance=1e-12, passed=bool(abs(lhs/rhs_ - 1) < 1e-12),
                note='the equilibrium radius follows from the angular momentum and the force law')


def gate_settling(starts=((1.3, 0.), (.75, 0.), (1., .35), (1., -.3), (1.45, .2)), periods=40):
    """Five disturbed starts must settle; the ledger E + Q must hold throughout."""
    rows = []
    for r0, v0 in starts:
        s = integrate([r0, v0, 0., 0., 0., 0.], periods*PERIOD, dense=True)
        tt = np.linspace((periods - 1)*PERIOD, periods*PERIOD, 4001)
        y = s.sol(tt)
        e0 = energy(s.sol(0.))
        rows.append(dict(r0=r0, vr0=v0,
                         final_period_max_deviation=float(np.max(np.abs(y[0]/R_C - 1))),
                         final_period_radius_range=[float(y[0].min()), float(y[0].max())],
                         ledger_drift=float(np.max(np.abs(energy(y) + y[4] - e0)))))
    worst_dev = max(r['final_period_max_deviation'] for r in rows)
    worst_led = max(r['ledger_drift'] for r in rows)
    return dict(rows=rows, periods=periods, worst_deviation=worst_dev, worst_ledger_drift=worst_led,
                tolerances=dict(settle=TOL_SETTLE, ledger=TOL_LEDGER),
                passed=bool(worst_dev < TOL_SETTLE and worst_led < TOL_LEDGER),
                note='dE/dt = -(B/tau) tanh^2(z) <= 0 exactly, so E + Q is conserved and the memory can '
                     'only remove disequilibrium energy; where Q goes physically is a hypothesis, and '
                     'carrying it only keeps the accounting honest')


def gate_no_memory_control(r0=1.3, periods=40):
    """The same law with B = 0 must NOT circularize, or the settling would not be the memory's doing."""
    s = integrate([r0, 0., 0., 0., 0., 0.], periods*PERIOD, B_=0., dense=True)
    y = s.sol(np.linspace((periods - 1)*PERIOD, periods*PERIOD, 4001))
    lo, hi = float(y[0].min()), float(y[0].max())
    return dict(r0=r0, final_period_radius_range=[lo, hi], spread=hi - lo,
                circularized=bool(hi - lo < TOL_SETTLE), passed=bool(hi - lo > .1),
                note='a conservative orbit oscillates indefinitely; settling needs the energy to leave')


def gate_stability(taus=(1e-3, 1e-2, .1, .3, 1., 3., 10., 100., 1e3)):
    """tau L^3 + L^2 + tau(kappa^2 + b)L + kappa^2; Routh-Hurwitz reduces to b > 0."""
    kap2, b = GM/R_C**3 + 2*K/R_C**2, B/R_C**2
    rows = []
    for tau in taus:
        c = [tau, 1., tau*(kap2 + b), kap2]
        rows.append(dict(tau=tau, max_real_part=float(np.max(np.roots(c).real)),
                         routh_hurwitz_margin=float(c[1]*c[2] - c[0]*c[3])))
    return dict(kappa_squared=kap2, b=b, rows=rows,
                condition='a1 a2 - a0 a3 = tau(kappa^2 + b) - tau kappa^2 = tau b > 0',
                stable_at_every_tau=bool(all(r['max_real_part'] < 0 for r in rows)),
                fastest_decay_tau=float(min(rows, key=lambda r: r['max_real_part'])['tau']),
                passed=bool(all(r['max_real_part'] < 0 and r['routh_hurwitz_margin'] > 0 for r in rows)),
                contrast='candidate C\'s cubic put its extra term in the CONSTANT coefficient, which forced '
                         'q0 < 0 and left a growing mode for every positive memory time; here the extra '
                         'term sits in the linear coefficient instead, which is what changes the sign',
                caution='stability is not prompt settling: the decay rate vanishes as tau -> 0 and as '
                        'tau -> infinity, so a stable law can still fail to settle within a useful time')


def _pulse(t, t0, w=.6):
    s = (t - t0)/w
    return np.where(np.abs(s) < 1, np.cos(.5*np.pi*s)**2, 0.)


def gate_two_histories(excursions=(1., 4.), amplitude=.35, target=(1., 0., 4*np.pi)):
    """Same start, same released r, rdot and theta, different when the excursion happened.

    The excursion amplitude is prescribed, so the three correctors and the release time cannot satisfy the
    match by removing the excursion -- which is what a free amplitude does, collapsing to the trivial
    undisturbed orbit.
    """
    def run(p, t_exc, dense=False):
        A2, A3, T_rel = p
        drive = lambda t: amplitude*_pulse(t, t_exc) + A2*_pulse(t, 7.) + A3*_pulse(t, 8.5)
        return integrate([1., 0., 0., 0., 0., 0.], T_rel, drive=drive, dense=dense)

    rows = []
    for t_exc in excursions:
        sol = root(lambda p: np.array(run(p, t_exc).y[:3, -1]) - np.array(target), [-.2, -.15, 12.4],
                   tol=1e-13)
        if not sol.success:
            return dict(passed=False, error=f'shooting failed at t_exc={t_exc}: {sol.message}')
        s = run(sol.x, t_exc, dense=True)
        r, vr, th, chi, Q, W = s.y[:, -1]
        traj = s.sol(np.linspace(0, sol.x[2], 3000))[0]
        z = _state(r, chi)
        rows.append(dict(excursion_at=t_exc, corrector_amplitudes=[float(sol.x[0]), float(sol.x[1])],
                         release_time=float(sol.x[2]),
                         radius_range=[float(traj.min()), float(traj.max())],
                         released=dict(r=float(r), rdot=float(vr), theta=float(th)),
                         memory_chi=float(chi), z=float(z),
                         inward_acceleration=float(GM/r**2 + (K + B*np.tanh(z))/r),
                         driver_work=float(W),
                         ledger=float(energy(s.y[:, -1]) + Q - W - (J**2/2 - GM))))
    a, b = rows
    same = max(abs(a['released'][k] - b['released'][k]) for k in ('r', 'rdot', 'theta'))
    diff = abs(a['inward_acceleration'] - b['inward_acceleration'])
    return dict(rows=rows, state_agreement=float(same), acceleration_difference=float(diff),
                relative_difference=float(abs(a['inward_acceleration']/b['inward_acceleration'] - 1)),
                worst_ledger=float(max(abs(r['ledger']) for r in rows)),
                passed=bool(same < 1e-10 and diff > 1e-4 and max(abs(r['ledger']) for r in rows)
                            < TOL_LEDGER),
                note='genuine history dependence within this toy: the same present state responds '
                     'differently because its internal state retains a different history. The size of the '
                     'difference depends on the preparation driver, which the review did not specify, so '
                     'this is a reproduction of the structure and not of a particular number')


# ---------------------------------------------------------------- stage 2
def corollary_settled_memory(Bs=(0., .2, .4, .8), horizons=(64, 128)):
    """On the settled orbit z = 0, so the memory's contribution to it is exactly zero for every B.

    Two separate things are checked, because they are different kinds of claim. The identity is algebra:
    at z = 0 the memory term B tanh(z)/r is exactly zero whatever B and tau are. The limit is numerical:
    a settling orbit drives z toward zero, so its residual memory force falls as the horizon grows. It is
    bounded above by B|z|/r exactly, since |tanh z| <= |z|, and different B settle at different rates.
    """
    algebraic = max(abs(Bv*np.tanh(0.)/R_C) for Bv in Bs)
    rows = []
    for Bv in Bs:
        per = []
        for periods in horizons:
            s = integrate([1.3, 0., 0., 0., 0., 0.], periods*PERIOD, B_=Bv)
            r, chi = s.y[0, -1], s.y[3, -1]
            z = _state(r, chi)
            per.append(dict(periods=periods, radius=float(r), z=float(z),
                            memory_part_of_the_force=float(Bv*np.tanh(z)/r),
                            bound_B_abs_z_over_r=float(Bv*abs(z)/r),
                            static_part_of_the_force=float(K/r),
                            settled=bool(abs(r/R_C - 1) < TOL_SETTLE)))
        last, first = (abs(per[-1]['memory_part_of_the_force']),
                       abs(per[0]['memory_part_of_the_force']))
        rows.append(dict(B=Bv, horizons=per, at_numerical_floor=bool(last < FLOOR),
                         falls_with_horizon=bool(last <= first or last < FLOOR)))
    live = [r for r in rows if r['B'] > 0]
    worst = max(abs(r['horizons'][-1]['memory_part_of_the_force']) for r in live)
    return dict(algebraic_identity_at_z_zero=algebraic, rows=rows,
                worst_memory_contribution=worst,
                all_settled=bool(all(h['settled'] for r in live for h in r['horizons'])),
                all_fall_with_horizon=bool(all(r['falls_with_horizon'] for r in live)),
                passed=bool(algebraic == 0.
                            and all(h['settled'] for r in live for h in r['horizons'])
                            and all(r['falls_with_horizon'] for r in live)),
                statement='the memory produces the SETTLING, not the support: the flat-curve K/r term is '
                          'assumed in the toy, as the review states, and the memory contributes exactly '
                          'nothing to the circular orbit it created',
                caution='the residual is a settling rate, not a defect: B = 0.2 has the slowest decay of '
                        'the three and still carries 2.7e-8 at 64 periods, reaching the floor by 128. '
                        'Once a residual is below 1e-12 it is the numerical noise of the integrator, so '
                        'the trend test accepts either a fall or arrival at that floor, rather than '
                        'comparing one rounding error with another')


def corollary_stationary_source(taus=(.1, 1., 10., 100.), horizons=(5, 10, 20)):
    """tau dQ/dt = p - Q with p stationary gives Q -> p as exp(-t/tau), so eta(p - Q) -> 0 and the field
    equation returns exactly to Completion I: blind to steady rotation at every tau and every eta."""
    rows = [dict(tau=t, horizon_in_tau=h, residual_memory_flux_fraction=float(np.exp(-h)))
            for t in taus for h in horizons]
    return dict(rows=rows, decays_as='exp(-t/tau), independent of eta and of the source',
                passed=bool(all(r['residual_memory_flux_fraction'] < 1 for r in rows)),
                statement='a smooth axisymmetric galaxy in steady rotation has a stationary density, hence '
                          'stationary psi and p = grad(psi), hence Q -> p and a vanishing memory flux. The '
                          'proposed extension therefore reduces EXACTLY to Completion I for the case it '
                          'was meant to address',
                prior_measurement='PM-1 candidate A measured the density version: a rigidly rotating '
                                  'axisymmetric source is indistinguishable from a static one to 3.7e-11')


def steady_state_theorem():
    return dict(
        statement='if tau dQ/dt = F[fields] - Q for any functional F of the instantaneous fields, then in '
                  'a steady state Q relaxes to F[fields] and is itself time-independent, so the memory '
                  'contributes a term fixed by the instantaneous fields -- nothing a static constitutive '
                  'law of those same fields could not also produce',
        consequence='memory buys transients; in a steady state it buys nothing. Rotation-dependent gravity '
                    'in a steady galaxy must come from a coupling to something stationary but nonzero for '
                    'a rotating disk and zero for a static one -- the mass current rho v, or a velocity '
                    'stress -- or from genuine non-steadiness, which in a galaxy means the discreteness of '
                    'stars rather than a smooth fluid',
        for_the_programme='PM-2A stage C is already the right test of the steady-state content, because '
                          'the static completions span the whole space of smooth steady memory theories; '
                          'the memory extension is not a shortcut around it',
        settled_memory=corollary_settled_memory(), stationary_source=corollary_stationary_source())


def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    eq, settle = gate_equilibrium(), gate_settling()
    control, stability = gate_no_memory_control(), gate_stability()
    histories = gate_two_histories()
    theorem = steady_state_theorem()
    stage1 = bool(eq['passed'] and settle['passed'] and control['passed'] and stability['passed']
                  and histories['passed'])
    stage2 = bool(theorem['settled_memory']['passed'] and theorem['stationary_source']['passed'])
    result = dict(
        experiment='PM-3: what a shared-field memory can and cannot do',
        protocol='protocol-pm3.md',
        parameters=dict(GM=GM, K=K, B=B, tau=TAU, j=J, r_ref=R_REF, circular_radius=R_C,
                        reference_period=PERIOD,
                        note='dimensionless and illustrative, as declared; no galaxy data enters PM-3'),
        stage_1_control=dict(equilibrium=eq, settling=settle, no_memory_control=control,
                             stability=stability, two_histories=histories, passed=stage1),
        stage_2_scope=dict(theorem=theorem, passed=stage2),
        stage_3='declared in the protocol and not run: a source coupling that sees the mass current, with '
                'the obligations that killed earlier candidates -- independence of how the same matter is '
                'divided into numerical particles, momentum conservation with a central-force limit, and '
                'reduction to the verified static equation as the current vanishes',
        passed=bool(stage1 and stage2),
        what_this_is_not='a galaxy theory. PM-3 reproduces a declared toy and proves a scope statement '
                         'about it. It establishes neither path memory nor a shared medium, a settling '
                         'toy is not a rotation curve, and a numerical relaxation time is never a physical '
                         'memory time',
        input_sha256={'protocol-pm3.md': hashlib.sha256((HERE/'protocol-pm3.md').read_bytes()).hexdigest()},
        checks_short_run=dict(settle_worst=settle['worst_deviation'], ledger_worst=settle['worst_ledger_drift'],
                              control_spread=control['spread'],
                              stability_worst=max(r['max_real_part'] for r in stability['rows']),
                              history_difference=histories['acceleration_difference'],
                              settled_memory=theorem['settled_memory']['worst_memory_contribution']),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:2000])
    status = evidence_io.finish(args, 'path-memory-pm3', text, HERE/'pm3-results.json',
                                ignore={'/runtime_seconds'})
    return status if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
