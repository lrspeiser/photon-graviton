"""Round 17, steps B and C: one kind of matter (round 16's, code/one_matter_v16.py) with the internal structures and the
time-reversal parity of the motion's mixing tested in code/rhythm_protect_v17.py, in the full energy-exact model.

Every piece (sources and receivers alike, same equations, same constants) holds its quiet oscillation s, its inversion
w and its finite store n as in round 16, and a list of internal fast families F_f (each a triplet, one mode per
direction) plus, optionally, a damped chamber R:
  bright family   radiates as a dipole into the shared wave through the single complex symmetric coupling M (its share of
                  the piece's dipole channel is sign_f F_f), loses `loss` inside, own frequency Delta_f
  dark family     the same without any coupling to the wave
  chamber         loses gR inside; the motion opens it from s; it passes on to one bright family by a fixed coupling gB
The wave sees each piece as (s, B_eff) with B_eff = sum over bright families of sign_f F_f; its forces are the round-16
ones on that vector. The motion's mixing into family f (weight c_f) has the parity of code/rhythm_protect_v17.py:
  even   dF_f/dt += (i/2)(-2 c_f delta) s,    ds/dt drive += -2 c_f delta . F_f    (round 16)
  odd    dF_f/dt += (i/2)(+2i c_f delta) s,   ds/dt drive += -2i c_f delta . F_f   (a velocity's parity)
both Hermitian: the energy the mixing moves between s's ensemble and a family is exactly accounted (the budget closes).
    ds/dt   = -(i/2) w (drive) + (G/2) w s - lambda s,      drive = (M z)_s + mixing terms
    dF_f/dt = (i/2) [sign_f (M z)_B (bright only) + mixing] - (loss_f + i Delta_f) F_f  [- i gB R, the chamber's family]
    dR/dt   = (i/2) [mixing]  - i gB F_target - gR R
    dw/dt   = W (1 - w) - g_par (1 + w) - 2 G |s|^2 - 2 Im(conj(s) drive),   dn/dt = -W (1 - w)/2
Energy: stores + ensembles ((1 + w)/2) + every family's and chamber's |F|^2 + internal output + radiated power, against
the start. Round 16's matter is structure 'single' with parity 'even'.

The runs (step B): cold, free (k = 2, 8, 16 for an isolated piece of that structure) and colliding sources, 3
arrangements, the same set-up as round 16 (48 sources in a ball of radius 3, 8 receivers of the same matter at radius 6);
and (step C) the pull on receivers at several distances, sources of several sizes, and receivers of several kinds.

    python code/one_matter_v17.py --set parity --output run-one-matter-v17/parity.json [--processes 3]
    python code/one_matter_v17.py --set best --output run-one-matter-v17/best.json
    python code/one_matter_v17.py --set distance --output run-one-matter-v17/distance.json
    python code/one_matter_v17.py --set mass --output run-one-matter-v17/mass.json
    python code/one_matter_summary_v16.py run-one-matter-v17/parity.json
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
from shared_wave_v16 import coupling, far_field, ball_min_sep, OMEGA   # noqa: E402
from one_matter_v16 import PIECE, sphere_points                        # noqa: E402
from rhythm_protect_v17 import STRUCTURES, q_for_k                     # noqa: E402

MIX = dict(even=(-2.0, -2.0), odd=(2.0j, -2.0j))       # (into the family, back onto s), times c_f delta


class Matter:
    def __init__(self, x, g0, piece, st, Delta):
        self.x = x; self.N = len(x)
        self.M, self.dM = coupling(x, g0, 1.0)
        self.Gam = self.M.imag
        self.p = piece; self.lam = (piece['W0'] + piece['g_par']) / 2
        self.st = st; self.Delta = Delta
        fams = st['families']; self.nf = len(fams)
        self.bright = np.array([f['kind'] == 'bright' for f in fams])
        self.sign = np.array([f.get('sign', 1.0) if f['kind'] == 'bright' else 0.0 for f in fams])
        self.c = np.array([f['c'] for f in fams]); self.loss = np.array([f['loss'] for f in fams])
        self.Dl = np.array([f['Delta'] for f in fams])
        self.valve = st.get('valve')
        self.kin, self.kback = MIX[st.get('parity', 'even')]

    def wave_vector(self, s, F):
        Beff = np.einsum('f,jfm->jm', self.sign, F)
        return np.concatenate([s[:, None], Beff], 1).ravel()

    def rhs(self, s, F, R, w, n, dl, n0):
        """s (N,), F (N, nf, 3), R (N, 3) or None, w, n (N,); dl (N, 3) the pieces' mixings."""
        p = self.p; N = self.N
        z = self.wave_vector(s, F)
        Mf = self.M @ z
        P = float(np.imag(np.vdot(z, Mf)))
        Mz = Mf.reshape(N, 4)
        dF_dot = np.einsum('jm,jfm->jf', dl, F)                           # delta . F_f per piece and family
        drive = Mz[:, 0] + self.kback * np.einsum('f,jf->j', self.c, dF_dot)
        dF = (0.5j * self.sign[None, :, None] * Mz[:, None, 1:] + 0.5j * self.kin * self.c[None, :, None] * dl[:, None, :] * s[:, None, None]
              - (self.loss + 1j * self.Dl)[None, :, None] * F - 1j * self.Delta[:, None, None] * F)
        dR = None
        if self.valve is not None:
            v = self.valve; t = v['into']
            drive = drive + self.kback * v['c'] * np.einsum('jm,jm->j', dl, R)
            dR = 0.5j * self.kin * v['c'] * dl * s[:, None] - 1j * v['gB'] * F[:, t, :] - (v['gR'] + 1j * self.Delta[:, None]) * R
            dF[:, t, :] += -1j * v['gB'] * R
        ds = -0.5j * w * drive + (0.5 * p['G'] * w - self.lam - 1j * self.Delta) * s
        W = p['W0'] * n / n0
        dw = W * (1 - w) - p['g_par'] * (1 + w) - 2 * p['G'] * np.abs(s) ** 2 - 2 * np.imag(np.conj(s) * drive)
        dn = -0.5 * W * (1 - w)
        sink = 0.5 * p['g_par'] * (1 + w) + p['G'] * np.abs(s) ** 2 + 2 * np.einsum('f,jfm->j', self.loss, np.abs(F) ** 2)
        if R is not None:
            sink = sink + 2 * self.valve['gR'] * np.sum(np.abs(R) ** 2, 1)
        return ds, dF, dR, dw, dn, float(np.sum(sink)), P


def energy(F, R, w, n):
    e = float(np.sum(0.5 * (1 + w)) + np.sum(np.abs(F) ** 2) + np.sum(n))
    return e + (float(np.sum(np.abs(R) ** 2)) if R is not None else 0.0)


def simulate(cfg):
    rng = np.random.default_rng(cfg['seed'])
    Ns, Nr = cfg['Ns'], cfg['Nr']
    g0 = cfg['gamma0']
    xs = ball_min_sep(Ns, cfg['Rb'], cfg.get('min_sep', 0.15), rng)
    xr = cfg['receivers'] if 'receivers' in cfg else sphere_points(Nr, cfg['rp'])
    xr = np.asarray(xr, float); Nr = len(xr)
    x = np.vstack([xs, xr]); N = Ns + Nr
    Delta = np.zeros(N)
    p = dict(PIECE, **cfg.get('piece', {}))
    st = STRUCTURES[cfg['structure']]
    mat = Matter(x, g0, p, st, Delta)
    w0 = 2 * mat.lam / p['G']
    amp0 = np.sqrt((p['W0'] * (1 - w0) - p['g_par'] * (1 + w0)) / (2 * p['G']))
    s = amp0 * np.exp(1j * rng.uniform(0, 2 * np.pi, N))
    F = np.zeros((N, mat.nf, 3), complex); R = np.zeros((N, 3), complex) if mat.valve is not None else None
    w = np.full(N, w0); n0 = cfg.get('n0', 1e5); n = np.full(N, n0)
    e_dir = rng.normal(size=(Ns, 3))
    dl = np.zeros((N, 3))
    kind, q, nu = cfg['kind'], cfg.get('q', 0.0), cfg.get('nu', 0.0)
    if kind in ('free', 'collisional', 'stopped'):
        dl[:Ns] = q * e_dir
    dt = cfg['dt']; nstep = int(round(cfg['T'] / dt)); nb = int(round(cfg['burn'] / dt))
    rho = np.exp(-nu * dt) if nu > 0 else 1.0; kick = q * np.sqrt(max(0.0, 1 - rho * rho))
    rc = xr - xs.mean(0) if cfg.get('center_on_sources') else xr
    rhat = rc / np.linalg.norm(rc, axis=1)[:, None]
    E_start = energy(F, R, w, n); sink_tot = 0.0; rad_tot = 0.0
    src = np.zeros(4 * N, bool); src[:4 * Ns] = True
    acc = dict(src_out=0.0, src_out_cold=0.0, src_bright=0.0, pull=0.0, pull_s=0.0, pull_B=0.0, lead=0.0, fed=0.0,
               amp_r=0.0, amp_s=0.0, w_r=0.0, w_s=0.0, E_r=0.0, sync=0.0, lead_src=0.0, fed_src=0.0, E_src=0.0, pull_s_src=0.0)
    cnt = 0; mom = []; pulls = np.zeros(Nr); pulls_src = np.zeros(Nr); s_prev = None; freq = np.zeros(N)
    src4 = np.repeat(np.arange(N) < Ns, 4)
    Mrs = mat.M[4 * Ns::4, :][:, src4]
    dMrs = mat.dM[:, 4 * Ns::4, :][:, :, src4]
    internal = dict(src_fam=0.0, src_chamber=0.0)
    every = cfg.get('every', 10)
    for stp in range(nstep):
        def f(s_, F_, R_, w_, n_):
            return mat.rhs(s_, F_, R_, w_, n_, dl, n0)
        k1 = f(s, F, R, w, n)
        def adv(h, k):
            return (s + h * k[0], F + h * k[1], (R + h * k[2]) if R is not None else None, w + h * k[3], n + h * k[4])
        k2 = f(*adv(0.5 * dt, k1)); k3 = f(*adv(0.5 * dt, k2)); k4 = f(*adv(dt, k3))
        comb = lambda i: (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6
        s = s + dt * comb(0); F = F + dt * comb(1)
        if R is not None:
            R = R + dt * comb(2)
        w = w + dt * comb(3); n = n + dt * comb(4)
        sink_tot += dt * comb(5); rad_tot += dt * comb(6)
        if kind == 'collisional':
            dl[:Ns] = rho * dl[:Ns] + kick * rng.normal(size=(Ns, 3))
        if kind == 'stopped' and stp == int(round(cfg['t_stop'] / dt)):
            dl[:] = 0.0                                              # the source is brought to rest
        if stp >= nb and stp % every == 0:
            z = mat.wave_vector(s, F); zz = z.reshape(N, 4)
            v = np.einsum('cab,b->ca', mat.dM, z)
            parts = (np.real(np.conj(z)[None, :] * v) / OMEGA).reshape(3, N, 4)
            Ftot = parts.sum(2).T; Fs_ = parts[:, :, 0].T
            pl = -(Ftot[Ns:] * rhat).sum(1); pls = -(Fs_[Ns:] * rhat).sum(1)
            acc['pull'] += pl.mean(); acc['pull_s'] += pls.mean(); acc['pull_B'] += (pl - pls).mean(); pulls += pl
            zsrc = z * src
            acc['src_out'] += float(np.real(np.conj(zsrc) @ (mat.Gam @ zsrc)))
            zb = zsrc.reshape(N, 4).copy(); zb[:, 0] = 0; zb = zb.ravel()
            acc['src_bright'] += float(np.real(np.conj(zb) @ (mat.Gam @ zb)))
            acc['src_out_cold'] += float(np.sum(2 * g0 * np.abs(s[:Ns]) ** 2))
            internal['src_fam'] += float(np.sum(2 * mat.loss[None, :, None] * np.abs(F[:Ns]) ** 2))
            if R is not None:
                internal['src_chamber'] += float(2 * mat.valve['gR'] * np.sum(np.abs(R[:Ns]) ** 2))
            Mz = (mat.M @ z).reshape(N, 4)
            Eext = Mz[Ns:, 0] - 2j * g0 * s[Ns:]
            acc['lead'] += float(np.mean(np.sin(np.angle(Eext) - np.angle(s[Ns:]))))
            acc['fed'] += float(np.mean(np.imag(np.conj(s[Ns:]) * Eext)))
            acc['E_r'] += float(np.mean(np.abs(Eext)))
            acc['amp_r'] += float(np.mean(np.abs(s[Ns:]))); acc['amp_s'] += float(np.mean(np.abs(s[:Ns])))
            acc['w_r'] += float(np.mean(w[Ns:])); acc['w_s'] += float(np.mean(w[:Ns]))
            acc['sync'] += float(np.abs(np.mean(np.exp(1j * np.angle(s[:Ns])))))
            Esrc = Mrs @ z[src4]
            acc['lead_src'] += float(np.mean(np.sin(np.angle(Esrc) - np.angle(s[Ns:]))))
            acc['fed_src'] += float(np.mean(np.imag(np.conj(s[Ns:]) * Esrc)))
            acc['E_src'] += float(np.mean(np.abs(Esrc)))
            gsrc = np.einsum('crb,b->rc', dMrs, z[src4])
            ps = -(np.real(np.conj(s[Ns:])[:, None] * gsrc) / OMEGA * rhat).sum(1)
            acc['pull_s_src'] += float(np.mean(ps)); pulls_src += ps
            if s_prev is not None:
                freq += np.angle(s * np.conj(s_prev))
            s_prev = s.copy()
            if cnt % 50 == 0 and cfg.get('momentum', True):
                Pf, pf = far_field(x, z, g0, 1.0, ndir=2000)
                mom.append(dict(force_sum=Ftot.sum(0).tolist(), momentum_out=pf.tolist(), power_far=Pf,
                                power_M=float(np.real(np.conj(z) @ (mat.Gam @ z)))))
            cnt += 1
    res = {k: v / cnt for k, v in acc.items()}
    res.update({k: v / cnt for k, v in internal.items()})
    span = (cnt - 1) * every * dt
    res['rhythm_sources'] = float(np.mean(freq[:Ns]) / span); res['rhythm_receivers'] = float(np.mean(freq[Ns:]) / span)
    res['rhythm_spread_sources'] = float(np.std(freq[:Ns]) / span); res['rhythm_spread_receivers'] = float(np.std(freq[Ns:]) / span)
    E_end = energy(F, R, w, n)
    res.update(cfg=cfg, pull_each=(pulls / cnt).tolist(), pull_s_src_each=(pulls_src / cnt).tolist(),
               receiver_distance=np.linalg.norm(rc, axis=1).tolist(),
               pull_se=float((pulls / cnt).std(ddof=1) / np.sqrt(Nr)) if Nr > 1 else None,
               energy=dict(start=E_start, end=E_end, internal_out=sink_tot, radiated=rad_tot,
                           residual=float((E_start - E_end - sink_tot - rad_tot) / (E_start - E_end))),
               store_used=float(1 - n.mean() / n0),
               momentum=[dict(m, balance=float(np.linalg.norm(np.array(m['force_sum']) + np.array(m['momentum_out']))
                                               / max(np.linalg.norm(m['momentum_out']), 1e-30))) for m in mom[:4]],
               bloch_max=float(np.max(4 * np.abs(s) ** 2 + w ** 2)))
    return res


def k_nominal(cfg):
    return cfg.get('k', 0.0)


BASE = dict(Ns=48, Nr=8, Rb=3.0, rp=6.0, gamma0=0.01, piece=dict(G=4.0, W0=1.6, g_par=0.01), dt=0.02,
            T=8000.0, burn=4000.0, every=10)


def shells(radii, per=4):
    """Receivers on several spheres around the source's centre: `per` directions per radius, turned from one radius to
    the next so that no two lie on one line of sight."""
    out = []
    for i, r in enumerate(radii):
        base = sphere_points(per, 1.0)
        a = 0.7 * i; ca, sa = np.cos(a), np.sin(a)
        rot = np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]]) @ np.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]])
        out.append(r * base @ rot.T)
    return np.vstack(out).tolist()


def case(structure, tag, k=0.0, nu=0.0, seed=1, **kw):
    q = q_for_k(STRUCTURES[structure], k, BASE['gamma0']) if k > 0 else 0.0
    kind = {'cold': 'rest', 'free': 'free', 'collisional': 'collisional', 'stopped': 'stopped'}[tag]
    dt = 0.01 if q > 5.0 else BASE['dt']                       # strong mixings (poor-scatterer pairs) need a finer step
    return dict(BASE, **dict(dict(structure=structure, tag=tag, kind=kind, k=k, q=q, nu=nu, seed=seed, dt=dt), **kw))


def configs(which):
    runs = []
    if which == 'parity':
        # round 16's matter with the velocity's parity: the decisive test again
        for sd in (1, 2, 3):
            for stn in ('single', 'single, odd'):
                runs.append(case(stn, 'cold', seed=sd))
                for k in (2.0, 8.0, 16.0):
                    runs.append(case(stn, 'free', k=k, seed=sd))
                for nu in (5.0, 50.0):
                    runs.append(case(stn, 'collisional', k=8.0, nu=nu, seed=sd))
        return runs
    if which == 'best':
        # the structures that kept warm sources in tune best in the reduced model, the same decisive test
        for sd in (1, 2, 3):
            for stn in ('paired, Delta = 4 (gamma + gi), odd', 'valve, gR = 20, gB = 5, odd'):
                runs.append(case(stn, 'cold', seed=sd))
                for k in (2.0, 8.0):
                    runs.append(case(stn, 'free', k=k, seed=sd))
                runs.append(case(stn, 'collisional', k=8.0, nu=50.0, seed=sd))
        return runs
    if which == 'distance':
        # step C: the pull at four distances (4 receivers each), cold, warm, colliding and stopped sources, 16,000 time
        # units (the far receivers re-time slowly); R(sigma, r) = (F(sigma, r)/F(0, r)) sqrt(I(0, r)/I(sigma, r))
        rec = shells((6.0, 9.0, 13.5, 20.0))
        for sd in (1, 2):
            for stn in ('single', 'paired, Delta = 4 (gamma + gi), odd'):
                kw = dict(receivers=rec, T=16000.0, burn=8000.0, momentum=False, seed=sd)
                runs.append(case(stn, 'cold', **kw))
                for k in (2.0, 8.0):
                    runs.append(case(stn, 'free', k=k, **kw))
                runs.append(case(stn, 'collisional', k=8.0, nu=50.0, **kw))
                runs.append(case(stn, 'stopped', k=8.0, t_stop=4000.0, **kw))
        return runs
    if which == 'mass':
        # step C: sources of 24, 48 and 96 pieces at the same density, receivers at the same distances
        rec = shells((12.0, 18.0))
        for sd in (1, 2):
            for Ns in (24, 48, 96):
                kw = dict(Ns=Ns, Rb=3.0 * (Ns / 48) ** (1 / 3), receivers=rec, T=8000.0, burn=4000.0, momentum=False, seed=sd)
                runs.append(case('single, odd', 'cold', **kw))
                runs.append(case('single, odd', 'free', k=8.0, **kw))
        return runs
    raise ValueError(which)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--processes', type=int, default=3); ap.add_argument('--set', default='parity')
    ap.add_argument('--only-structure', default=None)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    runs = configs(args.set)
    if args.only_structure:
        runs = [r for r in runs if r['structure'] == args.only_structure]
    res = []
    with Pool(args.processes) as pool:
        for r in pool.imap_unordered(simulate, runs):
            res.append(r); c = r['cfg']
            print(f"[{time.monotonic() - t0:5.0f} s] {c['structure']:>28} {c['tag']:12s} Ns {c['Ns']} k {c['k']:5.1f} (q {c['q']:.3f}) nu {c.get('nu', 0):4.1f} "
                  f"seed {c['seed']}: output {r['src_out']:.3e}, pull {r['pull']:+.3e} (quiet {r['pull_s']:+.3e}, from sources {r['pull_s_src']:+.3e}, "
                  f"radiators {r['pull_B']:+.3e}), keeping step {r['lead_src']:+.2f}, rhythm spread {r['rhythm_spread_sources']:.2e}, "
                  f"E res {r['energy']['residual']:.1e}", flush=True)
            args.output.write_text(json.dumps(dict(experiment=f'round 17: one kind of matter, set {args.set}', runs=res,
                                                   seconds=time.monotonic() - t0), indent=1, default=float) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
