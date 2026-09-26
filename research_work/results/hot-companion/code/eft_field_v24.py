"""Round 24: the Casimir-EFT candidate's local field equation on the collision maps (test A of
round24-casimir-eft.md, section 10).

The note's quasistatic scalar sector,
    div( (|grad phi| / a) grad phi ) = 4 pi G rho_phi,   rho_phi = rho_cold + (1 + k_eff) rho_free,
replaces the adopted law's S + g_hot magnitude and direction rule by one local equation. Here it is solved on the
suite's 3D collision grids (bullet_v4.kappa_map_v4's grid, densities and heat profiles; k_eff = the law's
k = 3 sigma^2/u^2 for stars and 0 for the gas, the note's two limits), in two ways:
  'qumond': g_phi = sqrt(a |Q|) Q/|Q|, Q the Newtonian field of rho_phi. Exact for round, flat or cylindrical
            sources; the adopted law's field form has the same structure;
  'aqual':  the exact solution, |g_phi| g_phi / a = Q + T with T divergence-free (the curl field), chosen so that
            g_phi is curl-free; found by iterating Helmholtz projections with isolated boundaries (_aqual).
The strong-field hold is not part of the note's equation. It enters as in the law, h = g_N + exp(-|g_N|/g_d) g_phi;
a law with lam -> infinity gives the note's equation as written. Lensing as for the law: Sigma_eff = -div h/(4 pi G),
projected along the line of sight.

Options (SETTINGS; MODES names the combinations used in round 24):
  source     'local':  today's gas and stars, with the heat of today's stars (the note's local equation);
             'memory': the law's memory, exactly as bullet_v4.kappa_map_v4 (the old companion rides with the
                       galaxies, a fresh sphere of radius u t around the stopped gas);
  magnitude  'vector': Q = F + g_hot, the heat's pull added as a vector (the note);
             'scalar': Q = F and |g_phi|^2 = a (|Q| + S), S the law's two-way scalar sum of the heat (the hot glow
                       heard from all around). As a field equation this is div[(|grad phi|/a - S/|grad phi|) grad phi]
                       = 4 pi G rho, from L = -|grad phi|^3/(12 pi G a) + S |grad phi|/(4 pi G) - rho phi (section 35);
             'excess': Q = F + g_hot (the note's source) and |g_phi|^2 = a (|Q| + S_ex), S_ex = S - |g_hot| >= 0 the part
                       of the hot glow's intensity that does not add up to a net flow: the note's equation plus one term,
                       L = -|grad phi|^3/(12 pi G a) + S_ex |grad phi|/(4 pi G) - rho_phi phi. Round sources: the law;
  direction  'unit': along Q; 'law': the law's (F + g_hot)/(|F| + |g_hot|) (qumond only);
  solver     'qumond' or 'aqual'.
MODES['law'] reproduces bullet_v4.kappa_map_v4 to float rounding (main() --check).

    python code/eft_field_v24.py --modes law,eft,eft_memory,eft_scalar --output run-eft-v24/bullet_modes.json
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import bullet_v3 as B                             # noqa: E402
import bullet_v4 as V                             # noqa: E402
import law as L                                   # noqa: E402

G = B.G
MODES = {
    'law':              dict(source='memory', magnitude='scalar', direction='law', solver='qumond'),
    'eft':              dict(source='local', magnitude='vector', direction='unit', solver='qumond'),
    'eft_aqual':        dict(source='local', magnitude='vector', direction='unit', solver='aqual'),
    'eft_memory':       dict(source='memory', magnitude='vector', direction='unit', solver='qumond'),
    'eft_scalar_local': dict(source='local', magnitude='scalar', direction='unit', solver='qumond'),
    'eft_scalar':       dict(source='memory', magnitude='scalar', direction='unit', solver='qumond'),
    'eft_scalar_aqual': dict(source='memory', magnitude='scalar', direction='unit', solver='aqual'),
    'eft_excess_local': dict(source='local', magnitude='excess', direction='unit', solver='qumond'),
    'eft_excess':       dict(source='memory', magnitude='excess', direction='unit', solver='qumond'),
    'eft_excess_aqual': dict(source='memory', magnitude='excess', direction='unit', solver='aqual'),
    'eft_excess_local_aqual': dict(source='local', magnitude='excess', direction='unit', solver='aqual'),
}
SETTINGS = dict(MODES['law'], aqual_iters=60, aqual_tol=1e-3, aqual_relax=0.3)   # relax 1.0 oscillates for the repair
                                                                                 # (run-eft-v24/bullet_aqual.json); 0.3 settles
REPORT = []          # one entry per map: the AQUAL iteration's history
_ORIG_MAP = V.kappa_map_v4


def install(mode='law'):
    """Make every collision map in the suite (bullet_v4.kappa_map_v4, read at call time by all collision codes) use
    this module's map in the given mode; 'law' restores the adopted law's own function."""
    if mode in (None, '', 'law'):
        V.kappa_map_v4 = _ORIG_MAP
        SETTINGS.update(MODES['law'])
        return
    SETTINGS.update(MODES[mode])
    V.kappa_map_v4 = kappa_map_eft


class _Conv:
    """bullet_v3.Conv with scipy's threaded FFT (for the AQUAL projections)."""
    def __init__(self, n, dx, fn, origin):
        import scipy.fft as sf
        self.sf = sf; m = 2 * n
        k = ((np.arange(m) - (np.arange(m) >= n) * m) * dx).astype(np.float32)
        r = np.sqrt(k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2, dtype=np.float32)
        with np.errstate(divide='ignore'):
            ker = np.where(r > 0, fn(np.maximum(r, 1e-6)), origin).astype(np.float32)
        del r
        self.n = n; self.K = sf.rfftn(ker, workers=-1); del ker

    def __call__(self, f):
        n = self.n; m = 2 * n
        pad = np.zeros((m, m, m), np.float32); pad[:n, :n, :n] = f
        return self.sf.irfftn(self.sf.rfftn(pad, workers=-1) * self.K, s=(m, m, m), workers=-1)[:n, :n, :n].astype(np.float32)


def _curl(v, dx):
    d = lambda f, ax: np.gradient(f, dx, axis=ax)
    return np.array([d(v[2], 1) - d(v[1], 2), d(v[0], 2) - d(v[2], 0), d(v[1], 0) - d(v[0], 1)], dtype=np.float32)


def _transverse(v, dx, inv_r):
    """The divergence-free part of v with isolated boundaries: curl A, A = (1/4 pi) int curl v / |x - x'| dV'."""
    w = _curl(v, dx)
    A = np.array([inv_r(w[i]) for i in range(3)], dtype=np.float32) * np.float32(dx ** 3 / (4 * np.pi))
    del w
    return _curl(A, dx)


def _norm(v):
    return np.sqrt(np.sum(v.astype(np.float64) ** 2, axis=0)).astype(np.float32) + np.float32(1e-30)


def _aqual(Q0, S, a, n, dx, iters, tol, relax):
    """Solve |g| g / a - S g/|g| = Q0 + T, curl g = 0, div T = 0 (S = 0: the note's equation). Returns the
    curl-free g and the history of (change in T, curl fraction of g) per iteration."""
    inv_r = _Conv(n, dx, lambda r: 1.0 / r, B.cube_average(lambda r: 1 / r) / dx)
    T = np.zeros_like(Q0); hist = []
    q_rms = float(np.sqrt(np.mean(Q0.astype(np.float64) ** 2)))

    def g_of(Q):
        Qm = _norm(Q)
        return np.sqrt(a * (Qm + S)) * Q / Qm
    for it in range(iters):
        g = g_of(Q0 + T)
        gT = _transverse(g, dx, inv_r)
        gL = g - gT
        gLm = _norm(gL)
        flux = (gLm ** 2 / a - S) * gL / gLm
        Tn = _transverse(flux, dx, inv_r)
        dT = float(np.sqrt(np.mean((Tn - T).astype(np.float64) ** 2))) / q_rms
        cf = float(np.sqrt(np.mean(gT.astype(np.float64) ** 2) / np.mean(g.astype(np.float64) ** 2)))
        T = T + np.float32(relax) * (Tn - T)
        hist.append(dict(it=it, dT=dT, curl_fraction=cf, T_over_Q=float(np.sqrt(np.mean(T.astype(np.float64) ** 2))) / q_rms))
        del g, gT, gL, gLm, flux, Tn
        if dT < tol:
            break
    g = g_of(Q0 + T)
    gL = g - _transverse(g, dx, inv_r)
    return gL, hist


def kappa_map_eft(current, ghost_gas, ghost_stars, pos, consts, n=192, dx=15.0, centre=(360., 50.),
                  fresh_kpc=30.0, heat=True, memory=True):
    """bullet_v4.kappa_map_v4's signature and grid, with the field from SETTINGS (module docstring)."""
    st = SETTINGS
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']; gd = lam * a
    x = (np.arange(n) - n / 2 + 0.5) * dx + centre[0]
    y = (np.arange(n) - n / 2 + 0.5) * dx + centre[1]
    z = (np.arange(n) - n / 2 + 0.5) * dx
    dV = dx ** 3
    rho_now = B.build_density(current, pos, x, y, z, 'gas') + B.build_density(current, pos, x, y, z, 'st')
    rho_gas_now = B.build_density(current, pos, x, y, z, 'gas')
    use_mem = st['source'] == 'memory'
    today = list(current.values())
    if use_mem:
        rho_ghost_gas = B.build_density(ghost_gas, pos, x, y, z, 'gas') if memory else rho_gas_now
    rho_ghost_st = np.zeros_like(rho_now); krho = np.zeros_like(rho_now)
    for comp, prof in ghost_stars:
        if not use_mem and not any(comp == c for c in today):
            continue                                  # local: only the stars that are there today
        rho_c = B.build_density({'c': comp}, pos, x, y, z, 'c')
        cx, cy = pos[comp['centre']]
        r3 = np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2 + z[None, None, :] ** 2)
        rho_ghost_st += rho_c
        if heat:
            krho += L.k_from_sig2(np.interp(r3, prof[0], prof[1]) ** 2, u).astype(np.float32) * rho_c
        del r3, rho_c
    inv_r = B.Conv(n, dx, lambda r: 1.0 / r, B.cube_average(lambda r: 1 / r) / dx)
    gN = -np.array(np.gradient(-G * inv_r(rho_now * dV), dx))
    if use_mem:
        F = -np.array(np.gradient(-G * inv_r((rho_ghost_gas + rho_ghost_st) * dV), dx))
    else:
        F = gN.copy()                                 # today's matter: the Newtonian field itself
    g_hot = -np.array(np.gradient(-G * inv_r(krho * dV), dx)) if heat else np.zeros_like(gN)
    del inv_r
    if use_mem and memory and fresh_kpc > 0:
        Kin = V.vec_conv(n, dx, rmax=fresh_kpc)
        F = F + G * dV * (V.apply_vec(Kin, rho_gas_now, n) - V.apply_vec(Kin, rho_ghost_gas, n))
        del Kin
    S = 0.0
    if heat and st['magnitude'] in ('scalar', 'excess'):
        inv_r2 = B.Conv(n, dx, lambda r: 1.0 / r ** 2, B.cube_average(lambda r: 1 / r ** 2) / dx ** 2)
        S = G * inv_r2(krho * dV)
        del inv_r2
    mag = np.sqrt(np.sum(gN ** 2, axis=0)) + 1e-30
    hold = np.exp(-mag / gd)
    Q = F if st['magnitude'] == 'scalar' else F + g_hot
    if st['magnitude'] == 'vector':
        S = 0.0
    elif st['magnitude'] == 'excess' and heat:
        S = np.maximum(S - np.sqrt(np.sum(g_hot ** 2, axis=0)), 0.0).astype(np.float32)
    if st['solver'] == 'aqual':
        g_phi, hist = _aqual(Q.astype(np.float32), S, a, n, dx, st['aqual_iters'], st['aqual_tol'], st['aqual_relax'])
        REPORT.append(dict(settings={k: st[k] for k in ('source', 'magnitude', 'solver')}, n=n, dx=dx, history=hist))
    else:
        Qm = np.sqrt(np.sum(Q ** 2, axis=0)) + 1e-30
        if st['direction'] == 'law':
            Fmag = np.sqrt(np.sum(F ** 2, axis=0)) + 1e-30
            hmag = np.sqrt(np.sum(g_hot ** 2, axis=0)) if heat else 0.0
            g_phi = np.sqrt(a * (Fmag + S)) * (F + g_hot) / (Fmag + hmag + 1e-30)
        else:
            g_phi = np.sqrt(a * (Qm + S)) * Q / Qm
    h = gN + hold * g_phi
    divh = sum(np.gradient(h[i], dx, axis=i) for i in range(3))
    Sig_eff = (-divh / (4 * np.pi * G)).sum(axis=2) * dx
    Sig_b = rho_now.sum(axis=2) * dx
    return x, y, Sig_eff, Sig_b


def hearing_for(mode):
    """The spherical parts (the stars' speeds from the Jeans equation, which set their heat) heard the same way as the
    map: Gauss's net flux when the heat is extra source mass ('vector'), the law's two-way sum when it is the scalar S."""
    return 'one_way_vector' if MODES[mode]['magnitude'] == 'vector' else 'two_way'


def bullet_rows(law, modes, n, dx, hearing='auto'):
    """The suite's Bullet case (bullet_static_v11.run_bullet, static distances) in each mode."""
    import bullet_static_v11 as BS
    import collisions_v10 as C10
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    rows = []
    geo0 = L.HOT_GEOMETRY
    for m in modes:
        install(m); t = time.monotonic(); REPORT.clear()
        L.HOT_GEOMETRY = hearing_for(m) if hearing == 'auto' else geo0
        try:
            r = BS.run_bullet(law, dict(f), n=n, dx=dx)
        finally:
            install('law'); geo = L.HOT_GEOMETRY; L.HOT_GEOMETRY = geo0
        ch = r['checks']
        row = dict(mode=m, settings=dict(MODES[m]), hearing=geo, seconds=time.monotonic() - t,
                   **{k: (v['value'] if isinstance(v, dict) else v) for k, v in ch.items()},
                   targets={k: v.get('target') for k, v in ch.items()}, aqual=list(REPORT))
        rows.append(row)
        print(f"{m:17s} kappa {ch['kappa_main']['value']:.3f}/{ch['kappa_sub']['value']:.3f}  gas {ch['gas_main']['value']:+.3f}/{ch['gas_sub']['value']:+.3f}  "
              f"peaks {ch['peak_main']['value']:.0f}/{ch['peak_sub']['value']:.0f} kpc ({ch['peak_main']['target']}, {ch['peak_sub']['target']})  "
              f"M250 {ch['m250_main']['value']:.3e} ({ch['m250_main']['target']}) / {ch['m250_sub']['value']:.3e} ({ch['m250_sub']['target']})  "
              f"[{row['seconds']:.0f} s]" + (f"  AQUAL: {len(REPORT[0]['history'])} its, curl {REPORT[0]['history'][-1]['curl_fraction']:.3f}, "
                                               f"T/Q {REPORT[0]['history'][-1]['T_over_Q']:.3f}" if REPORT else ''), flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--modes', type=str, default='law,eft,eft_memory,eft_scalar_local,eft_scalar')
    ap.add_argument('--law', type=str, default='round12')
    ap.add_argument('--n', type=int, default=192); ap.add_argument('--dx', type=float, default=15.0)
    ap.add_argument('--check', action='store_true', help='compare MODES["law"] with bullet_v4.kappa_map_v4 on the Bullet grid')
    ap.add_argument('--hearing', choices=['auto', 'law'], default='auto',
                    help="'auto': the spherical parts heard as the map's heat (hearing_for); 'law': as the --law says")
    ap.add_argument('--relax', type=float, default=None, help='AQUAL iteration relaxation (default SETTINGS: 0.3)')
    ap.add_argument('--iters', type=int, default=None, help='AQUAL iterations at most (default SETTINGS: 60)')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    from law_config import load_law
    import common as C
    law = load_law(args.law); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    out = dict(experiment='round 24: the Casimir-EFT local field equation on the Bullet Cluster (suite case, static distances)',
               law=law['name'], constants=dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'], u_kms=law['u_kms']), n=args.n, dx=args.dx)
    if args.check:
        import bullet_static_v11 as BS
        import collisions_v10 as C10
        f = C10.factors(0.296, 1.0, (70.0, 0.3))
        got = {}

        def grab(tag, fn):
            def w(*a, **k):
                res = fn(*a, **k); got[tag] = res[2]; return res
            return w
        V.kappa_map_v4 = grab('orig', _ORIG_MAP); BS.run_bullet(law, dict(f), n=args.n, dx=args.dx)
        SETTINGS.update(MODES['law']); V.kappa_map_v4 = grab('eft_law', kappa_map_eft); BS.run_bullet(law, dict(f), n=args.n, dx=args.dx)
        install('law')
        d = float(np.max(np.abs(got['orig'] - got['eft_law'])) / np.max(np.abs(got['orig'])))
        out['check_law_mode_max_rel_diff'] = d
        print(f'check: MODES["law"] against bullet_v4.kappa_map_v4, largest difference {d:.2e} of the peak surface density')
    if args.relax is not None: SETTINGS['aqual_relax'] = args.relax
    if args.iters is not None: SETTINGS['aqual_iters'] = args.iters
    out['rows'] = bullet_rows(law, args.modes.split(','), args.n, args.dx, args.hearing)
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
