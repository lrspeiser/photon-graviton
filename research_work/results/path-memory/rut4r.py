"""RUT-1 stage 4R: the reciprocal, energy-accounted benchmark (protocol-rut1.md, stage 4R).

    python rut4r.py [--output-dir DIR] [--canonical]

Three checks on the reciprocal form of the two-stage response, C = W * h with W * W = K:

A. It IS the same response. Along prescribed circular trajectories, the spectral h-field and the grid
   two-stage field that every formation run uses must give the same C and grad C, at the bodies and away
   from them -- two independent spatial representations of one equation.
B. Its energy closes. H(t) - H(0) + int (gamma/alpha) int h_t^2 dt must vanish for freely moving bodies, and
   converge as the timestep shrinks, with H = sum m[v^2/2 - GM/r - C(X)] + tau_form/(2 alpha) int h_t^2
   + 1/(2 alpha tau_keep) int h^2. This is a derived dissipation rate with a stated destination, not a
   bookkeeping variable fitted afterwards.
C. The grid formation runs are corroborated. From the same initial conditions and the same fixed step, the
   free trajectories under the spectral field and under the grid field must agree.

Limits, stated in the protocol: writing proportional to mass and a symmetric Gaussian kernel; no identified
physical reservoir for the dissipated power; no spatial causality; and no momentum accounting for that
reservoir and the fixed centre.

Regenerates rut4r-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))
import formation as FM  # noqa: E402
import longrun as LR  # noqa: E402
import reciprocal as RC  # noqa: E402
import rut3 as R3  # noqa: E402

T0 = 2*np.pi
TAU_KEEP, TAU_FORM = 10*T0, 3*T0
N, WR, SEED, JITTER = 32, .2, 1, .02
BOX, MODES = 5., 64


def _setup():
    w = WR
    q = R3._label_rate(w, TAU_KEEP, .1)
    pos, vel, rates = FM.ring_start(N, 1., q, JITTER, JITTER, SEED)
    return w, q, pos, vel, rates


def check_same_response(orbits=3., h=.005):
    """A: prescribed circular trajectories, both fields driven identically, compared at bodies and probes."""
    w, q, pos, _, rates = _setup()
    th0 = np.arctan2(pos[:, 1], pos[:, 0])
    grid = FM.MemoryField(2.5, w/5, w, TAU_KEEP, TAU_FORM)
    spec = RC.ReciprocalField(BOX, MODES, w, TAU_KEEP, TAU_FORM, alpha=q/N)
    ones = np.ones(N)
    steps = int(round(orbits*T0/h))
    for k in range(steps):
        tm = (k + .5)*h
        mid = np.stack([np.cos(th0 + tm), np.sin(th0 + tm)], axis=1)
        grid.advance(mid, rates, h)
        spec.advance(mid, ones, h)
    t = steps*h
    X = np.stack([np.cos(th0 + t), np.sin(th0 + t)], axis=1)
    phi = np.linspace(0, 2*np.pi, 181)[:-1]
    probes = np.concatenate([np.stack([rr*np.cos(phi), rr*np.sin(phi)], axis=1) for rr in (.9, 1.1)])
    rows = {}
    for label, pts in (('bodies', X), ('probe_rings_0.9_and_1.1', probes)):
        Cg, gg = grid.sample_value(pts), grid.sample_gradient(pts)
        Cs, gs = spec.C_and_grad(pts)
        rows[label] = dict(C_relative=float(np.max(np.abs(Cg - Cs))/np.max(np.abs(Cs))),
                           grad_relative=float(np.max(np.linalg.norm(gg - gs, axis=1))/np.max(np.linalg.norm(gs, axis=1))))
    worst = max(max(r.values()) for r in rows.values())
    return dict(orbits=orbits, step=h, box=BOX, modes=MODES, rows=rows, worst_relative=worst, tolerance=1e-3,
                passed=bool(worst < 1e-3),
                note='the grid model interpolates a cubic spline of carried value and gradient arrays; the spectral '
                     'model evaluates its truncated Fourier series exactly, so their agreement bounds the grid '
                     'representation error of the field every formation run uses')


def check_energy_balance(orbits=3., steps=(.02, .01, .005)):
    """B: H(t) - H(0) + dissipated energy, for free motion, converging as the step halves."""
    w, q, pos, vel, _ = _setup()
    rows = []
    for h in steps:
        f = RC.ReciprocalField(BOX, MODES, w, TAU_KEEP, TAU_FORM, alpha=q/N)
        o = RC.run_free(pos, vel, np.ones(N), f, orbits*T0, h, samples=60)
        rows.append(dict(step=h, worst_balance_relative_to_kinetic=o['worst_balance_relative'],
                         dissipated=o['total_dissipated'], kinetic_scale=o['kinetic_scale']))
    ratios = [rows[i]['worst_balance_relative_to_kinetic']/rows[i + 1]['worst_balance_relative_to_kinetic']
              for i in range(len(rows) - 1)]
    f = RC.ReciprocalField(BOX, MODES, w, TAU_KEEP, TAU_FORM, alpha=q/N)
    return dict(orbits=orbits, rows=rows, halving_ratios=ratios,
                decay_rates_per_T0=dict(slow=float(-f.lp.real*T0), fast=float(-f.lm.real*T0),
                                        expected_slow=float(T0/TAU_KEEP), expected_fast=float(T0/TAU_FORM)),
                tolerance=dict(worst_at_finest=1e-5, minimum_halving_ratio=3.),
                passed=bool(rows[-1]['worst_balance_relative_to_kinetic'] < 1e-5 and min(ratios) > 3.),
                note='the dissipated power has a derived rate and a named destination -- the field damping -- '
                     'but the physical reservoir that receives it is not identified')


def check_trajectories(orbits=3., h=.005):
    """C: free motion from the same state under the spectral field and under the grid field."""
    w, q, pos, vel, rates = _setup()
    g = LR.run_instrumented(pos, vel, rates, w, TAU_KEEP, orbits*T0, h, tau_form=TAU_FORM, eta=1.,
                            band_diagnostics=False, checkpoints=(1., 2., 3.))
    f = RC.ReciprocalField(BOX, MODES, w, TAU_KEEP, TAU_FORM, alpha=q/N)
    s = RC.run_free(pos, vel, np.ones(N), f, orbits*T0, h, samples=30, checkpoints=(1., 2., 3.))
    rows = []
    for key in ('1', '2', '3'):
        xg, xs = np.array(g['checkpoints'][key]['positions']), np.array(s['checkpoints'][key]['positions'])
        vg, vs = np.array(g['checkpoints'][key]['velocities']), np.array(s['checkpoints'][key]['velocities'])
        rows.append(dict(t_T0=float(key), max_position_difference=float(np.max(np.linalg.norm(xg - xs, axis=1))),
                         max_velocity_difference=float(np.max(np.linalg.norm(vg - vs, axis=1)))))
    worst = max(r['max_position_difference'] for r in rows)
    return dict(orbits=orbits, step=h, rows=rows, worst_position_difference=worst, tolerance=1e-3,
                passed=bool(worst < 1e-3),
                note='an independent spatial representation of the same equation reproduces the trajectories '
                     'the grid formation runs compute')


def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    A, B, C = check_same_response(), check_energy_balance(), check_trajectories()
    result = dict(
        experiment='RUT-1 stage 4R: the reciprocal, energy-accounted benchmark', protocol='protocol-rut1.md',
        model=dict(field='tau_form h_tt + gamma h_t + h/tau_keep = alpha W * rho',
                   force='a = -GM x/r^3 + grad(W * h)',
                   energy='H = sum m[v^2/2 - GM/r - C(X)] + tau_form/(2 alpha) int h_t^2 + 1/(2 alpha tau_keep) int h^2',
                   balance='dH/dt = -(gamma/alpha) int h_t^2', gamma=float(1 + TAU_FORM/TAU_KEEP)),
        same_response=A, energy_balance=B, trajectories=C,
        passed=bool(A['passed'] and B['passed'] and C['passed']),
        limits='writing proportional to mass and a symmetric Gaussian kernel; no identified physical reservoir; '
               'no spatial causality; no momentum accounting for the reservoir and the fixed centre',
        input_sha256={'protocol-rut1.md': hashlib.sha256((HERE/'protocol-rut1.md').read_bytes()).hexdigest()},
        checks_short_run=dict(same_response=A['worst_relative'],
                              balance_finest=B['rows'][-1]['worst_balance_relative_to_kinetic'],
                              trajectory=C['worst_position_difference']),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:1500])
    status = evidence_io.finish(args, 'path-memory-rut4r', text, HERE/'rut4r-results.json',
                                ignore={'/runtime_seconds'})
    return status if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
