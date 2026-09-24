"""Round 16, step 3b diagnostic: does a source keep a steady beat, given time? code/one_matter_v16.py's matter and main-run
constants (seed 1), run for 24,000 time units, tracking the unwrapped phase of every piece's quiet oscillation and of
the sources' wave at each receiver. Every 1,000 time units: the sources' mean rhythm and its spread, the receivers'
rhythm, the rhythm of the sources' wave at the receivers and its jitter, and the receivers' lead on that wave.

    python code/one_matter_rhythm_v16.py 0.0 24000 --output run-one-matter-v16/rhythm_cold.json     (a cold source)
    python code/one_matter_rhythm_v16.py 0.41 24000 --output run-one-matter-v16/rhythm_warm.json    (k = 2)
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import sys, json, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import one_matter_v16 as om
from shared_wave_v16 import ball_min_sep


def run(q, seed=1, T=8000.0, dt=0.02, gi=4.0, g0=0.01, piece=dict(G=4.0, W0=1.6, g_par=0.01), Rb=3.0, rp=6.0, Ns=48, Nr=8,
        window=250.0):
    rng = np.random.default_rng(seed)
    xs = ball_min_sep(Ns, Rb, 0.15, rng); xr = om.sphere_points(Nr, rp); x = np.vstack([xs, xr]); N = Ns + Nr
    p = dict(om.PIECE, **piece); mat = om.Matter(x, g0, 1.0, p, np.zeros(N), gi=gi)
    w0 = 2 * mat.lam / p['G']; amp0 = np.sqrt((p['W0'] * (1 - w0) - p['g_par'] * (1 + w0)) / (2 * p['G']))
    z = np.zeros((N, 4), complex); z[:, 0] = amp0 * np.exp(1j * rng.uniform(0, 2 * np.pi, N)); z = z.ravel()
    w = np.full(N, w0); n0 = 1e5; n = np.full(N, n0)
    e_dir = rng.normal(size=(Ns, 3)); dl = np.zeros((N, 3)); dl[:Ns] = q * e_dir
    Kb = om.mixing_block(dl)
    src4 = np.repeat(np.arange(N) < Ns, 4); Mrs = mat.M[4 * Ns::4, :][:, src4]
    nstep = int(T / dt); samp = 10; per = int(window / dt)
    ph = np.zeros(N); phE = np.zeros(Nr)
    s_prev = z.reshape(N, 4)[:, 0].copy(); E_prev = Mrs @ z[src4]
    mark_ph = ph.copy(); mark_phE = phE.copy(); jit = []; leads = []
    out = []
    for st in range(1, nstep + 1):
        f = lambda z_, w_, n_: mat.rhs(z_, w_, n_, Kb, n0)
        k1 = f(z, w, n); k2 = f(z + .5 * dt * k1[0], w + .5 * dt * k1[1], n + .5 * dt * k1[2])
        k3 = f(z + .5 * dt * k2[0], w + .5 * dt * k2[1], n + .5 * dt * k2[2]); k4 = f(z + dt * k3[0], w + dt * k3[1], n + dt * k3[2])
        z = z + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6; w = w + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
        n = n + dt * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
        if st % samp == 0:
            s = z.reshape(N, 4)[:, 0]; E = Mrs @ z[src4]
            ph += np.angle(s * np.conj(s_prev)); phE += np.angle(E * np.conj(E_prev))
            jit.append(np.angle(E * np.conj(E_prev)) / (samp * dt))
            leads.append(np.sin(np.angle(E) - np.angle(s[Ns:])))
            s_prev = s.copy(); E_prev = E.copy()
        if st % per == 0:
            r = (ph - mark_ph) / window; rE = (phE - mark_phE) / window
            J = np.array(jit)
            out.append(dict(t=st * dt, src_rhythm=float(r[:Ns].mean()), src_spread=float(r[:Ns].std()),
                            rec_rhythm=float(r[Ns:].mean()), rec_spread=float(r[Ns:].std()),
                            wave_rhythm=float(rE.mean()), wave_rhythm_spread=float(rE.std()),
                            wave_jitter=float(np.sqrt(np.mean((J - J.mean(0)) ** 2))), lead=float(np.mean(leads)),
                            sync=float(abs(np.mean(np.exp(1j * np.angle(z.reshape(N, 4)[:Ns, 0])))))))
            r_ = out[-1]
            print('  t {t:6.0f}  sources {src_rhythm:+.2e} (spread {src_spread:.1e})  receivers {rec_rhythm:+.2e}  wave at receivers {wave_rhythm:+.2e} (spread {wave_rhythm_spread:.1e})  lead {lead:+.2f}  sync {sync:.2f}'.format(**r_), flush=True)
            mark_ph = ph.copy(); mark_phE = phE.copy(); jit = []; leads = []
    return out


if __name__ == '__main__':
    q = float(sys.argv[1]); T = float(sys.argv[2]) if len(sys.argv) > 2 else 8000.0; t0 = time.time()
    out_path = sys.argv[sys.argv.index('--output') + 1] if '--output' in sys.argv else f'rhythm_q{q}.json'
    rows = run(q, T=T, window=1000.0)
    print(f'q {q}: {time.time() - t0:.0f} s')
    for r in rows:
        print('  t {t:6.0f}  sources {src_rhythm:+.2e} (spread {src_spread:.1e})  receivers {rec_rhythm:+.2e} (spread {rec_spread:.1e})  '
              'wave at receivers {wave_rhythm:+.2e} (spread {wave_rhythm_spread:.1e}, jitter {wave_jitter:.1e})  lead {lead:+.2f}  sync {sync:.2f}'.format(**r))
    json.dump(dict(q=q, T=T, rows=rows), open(out_path, 'w'), indent=1)
