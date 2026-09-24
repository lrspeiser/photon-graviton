"""Round 15: which kinds of test body are pulled by a passing wave, and how hard -- from their own equations.

One test body at r = 6 wavelengths in the steady wave of a monopole of strength A at the origin (companion wavelength 1,
wave speed 1, k = omega = 2 pi; G = e^{ikR}/(4 pi R)). The body carries an oscillation a (its source strength) and feels
only the local wave E:  pull = -(1/2) Re(conj(a) grad E) . rhat  (toward the source: positive),  power it feeds into
the wave  fed = (omega/2) Im(conj(a) E).  lead = sin(arg E - arg a): +1 is a quarter cycle ahead (feeding), -1 behind.
The bodies (all in the emitters' own units, radiation rate 1, coupling to the wave i E):
  passive     an ordinary oscillator driven by the wave: da/dt = -a/2 + i E.
  amplifier   an inverted one below its own threshold (stimulated emission; the drive's sign flips): da/dt = -a/2 - i E.
  pll         rounds 10-14: fixed amplitude, held a quarter cycle ahead (the rule, not derived here).
  laser       the independent calculation's gain-reservoir oscillator (independent-r15/, part 1), written with a complex
              amplitude and the physical coupling: da/dt = (g n - g_i - 1/2 + i Delta) a + i E,
              dn/dt = P - g_R n - 2 g n |a|^2 (P = 4, g = 1, g_R = 1, g_i = 0.1). Self-sustained; its phase is free.
  bloch       a self-sustained, inverted emitter: a pumped ensemble whose coherence s is kept up by its own collective
              emission ('superradiant laser'), mean-field Bloch equations in the wave's frame, the wave entering only
              as a drive whose torque carries the inversion w:
                ds/dt = (-g_perp + (Gc/2) w + i Delta) s - (i/2) b E w
                dw/dt = W (1 - w) - g_par (1 + w) - 2 Gc |s|^2 + i b (E conj(s) - conj(E) s)
              (Gc = 1, W = 0.4, g_par = 0.01, g_perp = (W + g_par)/2, b = 1). Its phase is free too.
  bloch_ordinary  the same ensemble with the wave's torque given the ordinary (absorbing) sign: the control.
Each self-sustained body starts at 8 random phases and runs 4000 time units; means over the second half.
strong: the bloch body alone at zero detuning in ever stronger waves (A up to 10^4, 2000 time units, dt = 0.002): what
happens when the wave drives it harder than its own pump can keep up with.

    python code/receivers_v15.py --output run-reservoir-force-v15/receivers_v15.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

K = 2 * np.pi; OMEGA = K; R_TB = 6.0
G = np.exp(1j * K * R_TB) / (4 * np.pi * R_TB); G1 = (1j * K - 1 / R_TB) * G          # dG/dR at the test body


def measures(a, A):
    E = A * G; gE = A * G1
    return -0.5 * np.real(np.conj(a) * gE), 0.5 * OMEGA * np.imag(np.conj(a) * E), np.sin(np.angle(E) - np.angle(a))


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + 0.5 * dt * k1); k3 = f(y + 0.5 * dt * k2); k4 = f(y + dt * k3)
    return y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def self_sustained(kind, A, deltas, seeds=8, T=4000.0, dt=0.01, detail=False):
    E = A * G
    rng = np.random.default_rng(1)
    D = np.repeat(np.asarray(deltas, float), seeds); m = D.size
    ph = np.exp(1j * rng.uniform(0, 2 * np.pi, m))
    if kind == 'laser':
        P, g, gR, gi = 4.0, 1.0, 1.0, 0.1
        y = np.concatenate([np.sqrt((P / (gi + 0.5) - gR / g) / 2) * ph, np.full(m, (gi + 0.5) / g) + 0j])
        def f(y):
            a, n = y[:m], y[m:].real
            return np.concatenate([(g * n - gi - 0.5 + 1j * D) * a + 1j * E, (P - gR * n - 2 * g * n * np.abs(a) ** 2) + 0j])
    else:
        Gc, W, gpar = 1.0, 0.4, 0.01; gperp = (W + gpar) / 2; b = 1.0
        w0 = 2 * gperp / Gc; sgn = 1.0 if kind == 'bloch' else -1.0
        y = np.concatenate([np.sqrt((W * (1 - w0) - gpar * (1 + w0)) / (2 * Gc)) * ph, np.full(m, w0) + 0j])
        def f(y):
            s, w = y[:m], y[m:].real
            Om = b * E
            ds = (-gperp + 0.5 * Gc * w + 1j * D) * s - sgn * 0.5j * Om * (w if sgn > 0 else w0)
            dw = W * (1 - w) - gpar * (1 + w) - 2 * Gc * np.abs(s) ** 2 + sgn * np.real(1j * (Om * np.conj(s) - np.conj(Om) * s))
            return np.concatenate([ds, dw + 0j])
    steps = int(T / dt); acc = np.zeros((5, m)); cnt = 0
    for st in range(steps):
        y = rk4(f, y, dt)
        if st > steps // 2 and st % 10 == 0:
            a = y[:m]; pl, fd, ld = measures(a, A)
            acc += np.array([pl, fd, ld, np.abs(a), y[m:].real]); cnt += 1
    acc /= cnt
    out = []
    for i, d in enumerate(deltas):
        sl = slice(i * seeds, (i + 1) * seeds)
        out.append(dict(Delta=float(d), pull=float(acc[0, sl].mean()), fed=float(acc[1, sl].mean()), lead=float(acc[2, sl].mean()),
                        amplitude=float(acc[3, sl].mean()), store=float(acc[4, sl].mean()), pull_spread_over_starts=float(acc[0, sl].std())))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    res = dict(experiment='round 15: which test bodies a passing wave pulls, from their own equations', r=R_TB, bodies={})
    amps = [1.0, 10.0, 100.0]
    for name, a_of in (('passive', lambda A: 2j * A * G), ('amplifier', lambda A: -2j * A * G),
                       ('pll', lambda A: np.exp(1j * (np.angle(A * G) - np.pi / 2)))):
        rows = []
        for A in amps:
            pl, fd, ld = measures(a_of(A), A)
            rows.append(dict(A=A, E=float(abs(A * G)), pull=float(pl), fed=float(fd), lead=float(ld), pull_over_fed=float(pl / fd)))
        res['bodies'][name] = rows
        print(name, [(round(r['E'], 5), f"{r['pull']:.4e}", round(r['lead'], 3)) for r in rows], flush=True)
    for name in ('laser', 'bloch', 'bloch_ordinary'):
        rows = []
        for A in amps:
            E = abs(A * G)
            amp0 = 1.68 if name == 'laser' else 0.333
            wscale = 0.41 if name != 'laser' else 1.0
            lock = (E / amp0) * (0.5 * wscale if name != 'laser' else 1.0)     # the locking range's rough size
            deltas = [0.0, 0.5 * lock, -0.5 * lock, 2 * lock]
            for r_ in self_sustained(name, A, deltas):
                rows.append(dict(A=A, E=float(E), locking_range=float(lock), **r_))
        res['bodies'][name] = rows
        print(name, [(round(r['E'], 5), round(r['Delta'], 5), f"{r['pull']:.4e}", round(r['lead'], 3), round(r['amplitude'], 3)) for r in rows], flush=True)
    rows = []
    for A in (10.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0):
        r_ = self_sustained('bloch', A, [0.0], seeds=4, T=2000.0, dt=0.002)[0]
        rows.append(dict(A=A, E=float(abs(A * G)), **r_))
        print('strong', round(abs(A * G), 4), f"pull {r_['pull']:.4e} lead {r_['lead']:.3f} amplitude {r_['amplitude']:.4f} inversion {r_['store']:.4f}", flush=True)
    res['strong'] = rows
    args.output.write_text(json.dumps(res, indent=1) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
