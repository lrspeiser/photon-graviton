"""Round 16, step 3b: one kind of matter, every piece both sender and receiver, the shared wave fed back exactly.

An independent audit of round 15 asked for this as the next construction: "Keep the quiet and bright internal modes.
Give every piece the same finite-energy reservoir and the same local coupling to the wave. Derive its outgoing
radiation, its energy loss, its response to incoming radiation and its mechanical force from that coupling. Include
mutual feedback from the beginning. Then repeat the cold, freely moving and collisionally interrupted cases ... Do not
add a separate strong-field switch." Round 15 used two descriptions of matter: sources with a finite quiet store, and
receivers with a continuously applied pump. Here there is one.

Every piece (sources and receivers alike, same equations, same constants) holds:
  s      its quiet oscillation: the coherence of an inverted ensemble (it holds energy it is ready to give), Bloch
         normalization 4|s|^2 + w^2 <= 1; it radiates into the companion wave as a monopole (rate gamma_0)
  w      the ensemble's inversion; its energy is (1 + w)/2
  n      a finite store that refills the inversion (rate W0 n/n0), debited exactly
  B_m    three ordinary radiators (m = x, y, z), passive oscillators that radiate as dipoles (rate gamma = 1) and lose
         energy inside at gi; f = gamma/(gamma + gi) is the share of their energy that goes into the companion
and motion mixes s into the B's with delta = chi w (the independent calculation's postulate, still an input).
The piece keeps its quiet oscillation going with an internal collective gain G (an internal channel whose output
leaves the companion; energy tracked), as a superradiant laser does. A receiver that feeds a passing wave is pulled by
the fed power over the wave's speed, less 2f of it: its own radiators, stirred by the same wave, recoil against it
(README §26.4). Run sets: 'main' (f = 0.2), 'lowf' (f = 0.05), 'disp' (a cold source put out of tune by hand).

The wave coupling is ONE complex symmetric matrix M (code/shared_wave_v16.coupling: the exact Green's function, near
and far field, validated there: z^dagger Im(M) z is the power through a distant sphere, and the forces balance the
momentum carried away). With z the list of every piece's (s, B), the equations are
    ds_j/dt = -(i/2) w_j [M z + K_j z]_s + (G/2) w_j s_j - lambda s_j - i Delta_j s_j
    dB_j/dt = +(i/2)     [M z + K_j z]_B                  - (gi + i Delta_j) B_j
    dw_j/dt = W_j (1 - w_j) - g_par (1 + w_j) - 2 G |s_j|^2 - 2 Im(conj(s_j) [M z + K_j z]_s)
    dn_j/dt = -W_j (1 - w_j)/2,     W_j = W0 n_j/n0,     lambda = (W0 + g_par)/2
K_j is the motion's mixing, the 4x4 block [[0, -2 delta^T], [-2 delta, 0]]. The factor -w_j on the quiet channel is
what "inverted" means: an ordinary (w = -1) ensemble would respond like a passive oscillator. Energy: each channel's
loss equals the power it feeds the wave, Im(conj(z_a) (M z)_a), so the sum is z^dagger Im(M) z exactly -- the
interference between pieces is inside every piece's own budget. Force on piece j: (1/omega) Re sum_l z_j^dagger
dM_jl z_l (the same M, differentiated). Nothing converts power into force; no piece's timing is imposed.

Set-up: Ns source pieces in a ball (radius Rb, at least 0.15 wavelength apart) and Nr receivers on a sphere of radius
rp, at rest (cold). The sources are at rest (cold), move freely (delta fixed, Gaussian with rms q per component), or
collide (delta an Ornstein-Uhlenbeck process, rate nu). A run first lets every piece settle, then measures:
  the sources' companion output (from Im M restricted to the sources' own channels), the receivers' pull toward the
  cloud (split into the quiet channel's part and the radiators' part), their lead sin(arg E - arg s), the power they
  feed, the energy budget (stores + ensembles + radiators + internal output + radiated, against the start) and the
  momentum budget (all forces against the far-field momentum flux).

    python code/one_matter_v16.py --output run-one-matter-v16/one_matter_v16.json [--set main|lowf|disp] [--processes 4]
    python code/one_matter_summary_v16.py run-one-matter-v16/one_matter_v16.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from shared_wave_v16 import coupling, far_field, forces, ball_min_sep, OMEGA   # noqa: E402

PIECE = dict(G=1.0, W0=0.4, g_par=0.01)        # the internal constants (round 15's inverted receiver)


def sphere_points(n, r):
    i = np.arange(n) + 0.5; phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return r * np.vstack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]).T


def mixing_block(dl):
    """K as a (N, 4, 4) stack: the motion's mixing of s into the B's, entering like a drive -2 delta.B on s."""
    N = len(dl); Kb = np.zeros((N, 4, 4))
    Kb[:, 0, 1:] = -2 * dl; Kb[:, 1:, 0] = -2 * dl
    return Kb


class Matter:
    def __init__(self, x, g0, g, piece, Delta, gi=0.0):
        self.x = x; self.N = len(x); self.gi = gi               # gi: the radiators' internal loss (to the internal channel)
        self.M, self.dM = coupling(x, g0, g)
        self.Gam = self.M.imag
        self.p = piece; self.lam = (piece['W0'] + piece['g_par']) / 2
        self.Delta = Delta

    def rhs(self, z, w, n, Kb, n0):
        N, p = self.N, self.p
        zz = z.reshape(N, 4)
        Mf = self.M @ z
        P = float(np.imag(np.vdot(z, Mf)))                                   # = z^dagger Im(M) z, to infinity
        Mz = Mf.reshape(N, 4) + np.einsum('jab,jb->ja', Kb, zz)
        s = zz[:, 0]
        dz = np.empty_like(zz)
        dz[:, 0] = -0.5j * w * Mz[:, 0] + (0.5 * p['G'] * w - self.lam - 1j * self.Delta) * s
        dz[:, 1:] = 0.5j * Mz[:, 1:] - (self.gi + 1j * self.Delta[:, None]) * zz[:, 1:]
        W = p['W0'] * n / n0
        dw = W * (1 - w) - p['g_par'] * (1 + w) - 2 * p['G'] * np.abs(s) ** 2 - 2 * np.imag(np.conj(s) * Mz[:, 0])
        dn = -0.5 * W * (1 - w)
        sink = 0.5 * p['g_par'] * (1 + w) + p['G'] * np.abs(s) ** 2 + 2 * self.gi * np.sum(np.abs(zz[:, 1:]) ** 2, 1)
        return dz.ravel(), dw, dn, float(np.sum(sink)), P


def energy(z, w, n):
    zz = z.reshape(-1, 4)
    return float(np.sum(0.5 * (1 + w)) + np.sum(np.abs(zz[:, 1:]) ** 2) + np.sum(n))


def simulate(cfg):
    rng = np.random.default_rng(cfg['seed'])
    Ns, Nr = cfg['Ns'], cfg['Nr']
    g0, g = cfg['gamma0'], 1.0
    xs = ball_min_sep(Ns, cfg['Rb'], 0.15, rng)
    xr = sphere_points(Nr, cfg['rp'])
    x = np.vstack([xs, xr]); N = Ns + Nr
    Delta = np.zeros(N)
    if cfg.get('Delta0', 0.0) > 0:
        Delta[:Ns] = cfg['Delta0'] * rng.normal(size=Ns)
    p = dict(PIECE, **cfg.get('piece', {}))
    mat = Matter(x, g0, g, p, Delta, gi=cfg.get('gi', 0.0))
    w0 = 2 * mat.lam / p['G']
    amp0 = np.sqrt((p['W0'] * (1 - w0) - p['g_par'] * (1 + w0)) / (2 * p['G']))
    z = np.zeros((N, 4), complex); z[:, 0] = amp0 * np.exp(1j * rng.uniform(0, 2 * np.pi, N)); z = z.ravel()
    w = np.full(N, w0); n0 = cfg.get('n0', 1e5); n = np.full(N, n0)
    e_dir = rng.normal(size=(Ns, 3))
    dl = np.zeros((N, 3))
    kind, q, nu = cfg['kind'], cfg.get('q', 0.0), cfg.get('nu', 0.0)
    if kind in ('free', 'collisional'):
        dl[:Ns] = q * g * e_dir
    dt = cfg['dt']; nstep = int(round(cfg['T'] / dt)); nb = int(round(cfg['burn'] / dt))
    rho = np.exp(-nu * dt) if nu > 0 else 1.0; kick = q * g * np.sqrt(max(0.0, 1 - rho * rho))
    rhat = xr / cfg['rp']
    E_start = energy(z, w, n); sink_tot = 0.0; rad_tot = 0.0
    # indices of the sources' own channels, for their output
    src = np.zeros(4 * N, bool); src[:4 * Ns] = True
    Gss = mat.Gam[np.ix_(src, src)]
    acc = dict(src_out=0.0, src_out_cold=0.0, src_bright=0.0, pull=0.0, pull_s=0.0, pull_B=0.0, lead=0.0, fed=0.0,
               amp_r=0.0, amp_s=0.0, w_r=0.0, w_s=0.0, E_r=0.0, sync=0.0, lead_src=0.0, fed_src=0.0, E_src=0.0, pull_s_src=0.0)
    cnt = 0; mom = []; pulls = np.zeros(Nr); s_prev = None; t_prev = None; freq = np.zeros(N)
    src4 = np.repeat(np.arange(N) < Ns, 4)
    Mrs = mat.M[4 * Ns::4, :][:, src4]                       # receivers' quiet channels <- the sources' channels
    dMrs = mat.dM[:, 4 * Ns::4, :][:, :, src4]                # and its gradient: the pull from the sources' wave alone
    Kb = mixing_block(dl)
    for st in range(nstep):
        # RK4 on (z, w, n) and the two tallies, the mixing frozen over the step
        def f(z_, w_, n_):
            return mat.rhs(z_, w_, n_, Kb, n0)
        k1 = f(z, w, n)
        k2 = f(z + 0.5 * dt * k1[0], w + 0.5 * dt * k1[1], n + 0.5 * dt * k1[2])
        k3 = f(z + 0.5 * dt * k2[0], w + 0.5 * dt * k2[1], n + 0.5 * dt * k2[2])
        k4 = f(z + dt * k3[0], w + dt * k3[1], n + dt * k3[2])
        z = z + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
        w = w + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
        n = n + dt * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
        sink_tot += dt * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]) / 6          # the tallies ride along with RK4
        rad_tot += dt * (k1[4] + 2 * k2[4] + 2 * k3[4] + k4[4]) / 6
        if kind == 'collisional':
            dl[:Ns] = rho * dl[:Ns] + kick * rng.normal(size=(Ns, 3)); Kb = mixing_block(dl)
        if st >= nb and st % cfg.get('every', 5) == 0:
            zz = z.reshape(N, 4)
            # the force on each piece, split by which of its own channels feels it (quiet oscillation or radiators)
            v = np.einsum('cab,b->ca', mat.dM, z)
            parts = (np.real(np.conj(z)[None, :] * v) / OMEGA).reshape(3, N, 4)
            F = parts.sum(2).T; Fs = parts[:, :, 0].T
            pl = -(F[Ns:] * rhat).sum(1); pls = -(Fs[Ns:] * rhat).sum(1)
            acc['pull'] += pl.mean(); acc['pull_s'] += pls.mean(); acc['pull_B'] += (pl - pls).mean(); pulls += pl
            zsrc = z * src
            acc['src_out'] += float(np.real(np.conj(zsrc) @ (mat.Gam @ zsrc)))
            zb = zsrc.reshape(N, 4).copy(); zb[:, 0] = 0; zb = zb.ravel()
            acc['src_bright'] += float(np.real(np.conj(zb) @ (mat.Gam @ zb)))
            acc['src_out_cold'] += float(np.sum(2 * g0 * np.abs(zz[:Ns, 0]) ** 2))
            Mz = (mat.M @ z).reshape(N, 4)
            # the wave at each receiver from everything else (its own self-part removed)
            Eext = Mz[Ns:, 0] - 2j * g0 * zz[Ns:, 0]
            acc['lead'] += float(np.mean(np.sin(np.angle(Eext) - np.angle(zz[Ns:, 0]))))
            acc['fed'] += float(np.mean(np.imag(np.conj(zz[Ns:, 0]) * Eext)))
            acc['E_r'] += float(np.mean(np.abs(Eext)))
            acc['amp_r'] += float(np.mean(np.abs(zz[Ns:, 0]))); acc['amp_s'] += float(np.mean(np.abs(zz[:Ns, 0])))
            acc['w_r'] += float(np.mean(w[Ns:])); acc['w_s'] += float(np.mean(w[:Ns]))
            acc['sync'] += float(np.abs(np.mean(np.exp(1j * np.angle(zz[:Ns, 0])))))
            Esrc = Mrs @ z[src4]                                    # the wave at each receiver from the sources alone
            acc['lead_src'] += float(np.mean(np.sin(np.angle(Esrc) - np.angle(zz[Ns:, 0]))))
            acc['fed_src'] += float(np.mean(np.imag(np.conj(zz[Ns:, 0]) * Esrc)))
            acc['E_src'] += float(np.mean(np.abs(Esrc)))
            gsrc = np.einsum('crb,b->rc', dMrs, z[src4])
            acc['pull_s_src'] += float(np.mean(-(np.real(np.conj(zz[Ns:, 0])[:, None] * gsrc) / OMEGA * rhat).sum(1)))
            tnow = st * dt
            if s_prev is not None:
                freq += np.angle(zz[:, 0] * np.conj(s_prev))           # unwrapped phase advance (small per sample)
            s_prev = zz[:, 0].copy(); t_prev = tnow
            if cnt % 50 == 0:
                Pf, pf = far_field(x, z, g0, g, ndir=2000)
                mom.append(dict(force_sum=F.sum(0).tolist(), momentum_out=pf.tolist(), power_far=Pf,
                                power_M=float(np.real(np.conj(z) @ (mat.Gam @ z)))))
            cnt += 1
    res = {k: v / cnt for k, v in acc.items()}
    span = (cnt - 1) * cfg.get('every', 5) * dt
    res['rhythm_sources'] = float(np.mean(freq[:Ns]) / span); res['rhythm_receivers'] = float(np.mean(freq[Ns:]) / span)
    res['rhythm_spread_sources'] = float(np.std(freq[:Ns]) / span); res['rhythm_spread_receivers'] = float(np.std(freq[Ns:]) / span)
    E_end = energy(z, w, n)
    res.update(cfg=cfg, pull_se=float((pulls / cnt).std(ddof=1) / np.sqrt(Nr)),
               energy=dict(start=E_start, end=E_end, internal_out=sink_tot, radiated=rad_tot,
                           residual=float((E_start - E_end - sink_tot - rad_tot) / (E_start - E_end))),
               store_used=float(1 - n.mean() / n0),
               momentum=[dict(m, balance=float(np.linalg.norm(np.array(m['force_sum']) + np.array(m['momentum_out']))
                                               / max(np.linalg.norm(m['momentum_out']), 1e-30))) for m in mom[:4]],
               bloch_max=float(np.max(4 * np.abs(z.reshape(N, 4)[:, 0]) ** 2 + w ** 2)))
    return res


def k_nominal(cfg):
    """The heat weight of an isolated piece: 3 q^2 gamma / (gamma_0 (gamma + gi)^2)."""
    gt = 1.0 + cfg.get('gi', 0.0)
    return 3 * cfg.get('q', 0.0) ** 2 / (cfg['gamma0'] * gt * gt)


def configs(quick=False, which='main'):
    if which == 'disp':
        # a cold source whose pieces are put out of tune by hand (rms Delta0 in their own rhythms), by the amounts that
        # warming produced in the main run: is being out of tune, by itself, what loses the receivers?
        base = dict(Ns=48, Nr=8, Rb=3.0, rp=6.0, gamma0=0.01, gi=4.0, piece=dict(G=4.0, W0=1.6, g_par=0.01), dt=0.02,
                    T=8000.0, burn=4000.0, every=10, set='disp')
        return [dict(base, tag='cold_detuned', kind='rest', Delta0=d0, seed=sd) for sd in (1, 2) for d0 in (0.001, 0.0025, 0.004)]
    if which == 'lowf':
        # the same matter with radiators that send a twentieth of their energy into the companion (f = 0.05), and a
        # stronger internal throughput (G, W0 ten times larger) so the heat's extra drain stays small beside it
        base = dict(Ns=48, Nr=8, Rb=3.0, rp=6.0, gamma0=0.01, gi=19.0, piece=dict(G=40.0, W0=16.0, g_par=0.01), dt=0.01,
                    T=8000.0, burn=4000.0, every=20, set='lowf')
        runs = []
        for sd in (1, 2):
            runs.append(dict(base, tag='cold', kind='rest', seed=sd))
            for q in (1.63, 3.27):                        # k = 2, 8 for an isolated piece
                runs.append(dict(base, tag='free', kind='free', q=q, seed=sd))
            runs.append(dict(base, tag='collisional', kind='collisional', q=3.27, nu=200.0, seed=sd))   # 20/(20+200) = 0.09
        return runs
    # dilute source (48 pieces in radius 3), 8 receivers of the same matter at radius 6; radiators that send a fifth of
    # their energy into the companion (f = 1/(1 + gi) = 0.2): see the README for why f < 1/2 is needed
    base = dict(Ns=48, Nr=8, Rb=3.0, rp=6.0, gamma0=0.01, gi=4.0, piece=dict(G=4.0, W0=1.6, g_par=0.01), dt=0.02,
                T=4000.0 if quick else 8000.0, burn=2000.0 if quick else 4000.0, every=10)
    seeds = (1,) if quick else (1, 2, 3)
    runs = []
    for sd in seeds:
        runs.append(dict(base, tag='cold', kind='rest', seed=sd))
        for q in (0.2, 0.41, 0.82, 1.15):                 # k = 0.5, 2, 8, 16 for an isolated piece
            runs.append(dict(base, tag='free', kind='free', q=q, seed=sd))
        if not quick:
            for nu in (5.0, 50.0):                        # (gamma + gi)/(gamma + gi + nu) = 0.5, 0.09
                runs.append(dict(base, tag='collisional', kind='collisional', q=0.82, nu=nu, seed=sd))
    return runs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=3); ap.add_argument('--quick', action='store_true')
    ap.add_argument('--only', default=None); ap.add_argument('--set', default='main')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    runs = configs(args.quick, args.set)
    if args.only:
        runs = [r for r in runs if r['tag'] == args.only]
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(simulate, runs):
            res.append(r); c = r['cfg']
            print(f"[{time.monotonic() - t0:5.0f} s] {c['tag']:12s} q {c.get('q', 0):.4f} (k {k_nominal(c):.2f}) nu {c.get('nu', 0):4.1f} seed {c['seed']}: "
                  f"source output {r['src_out']:.4e} (bright {r['src_bright']:.3e}); pull {r['pull']:.4e} (quiet {r['pull_s']:.4e}, "
                  f"radiators {r['pull_B']:.3e}); lead {r['lead']:.3f}; sync {r['sync']:.3f}; energy residual "
                  f"{r['energy']['residual']:.1e}", flush=True)
    args.output.write_text(json.dumps(dict(experiment='round 16 step 3b: one kind of matter, every piece both sender and receiver',
                                           piece=PIECE, runs=res, seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print('wrote', args.output, f'({time.monotonic() - t0:.0f} s)')


if __name__ == '__main__':
    main()
