"""Round 11: which static geometry do the gravity tests prefer?

Both variants keep the project's redshift law 1 + z = e^(alpha D) and its brightness law D_L = (1 + z) D
(energy loss and arrival-rate stretching); they differ in how an angle becomes a size:
    fixed   D_A = D            (the fixed-material transport branch; used since round 10)
    metric  D_A = D / (1 + z)  (the material-coasting geometry of five_candidate_tests/prior/
                                conformal_action_derivation_derivation.md, eq. 13; D_L = (1 + z)^2 D_A)
For each variant, with everything converted from the papers' own conventions:
  1. X-COP (stars deprojected): the constants refitted jointly (a and g_d on SPARC, u on X-COP);
  2. KiDS lensing (all, blue, red, disc, bulge, GAMA; the early/late gap) and Mistele's lensing speeds;
  3. the six SLACS lenses (light = matter; the stars their Einstein radii need);
  4. the far collisions MACS J0025, Abell 520 and El Gordo (published stars, on the Chabrier basis);
  5. the Bullet Cluster (round-5 case).

    python code/distance_variants_v11.py --output-dir run-distance-variants-v11
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                   # noqa: E402
import run as RUN                                 # noqa: E402
import run_v3 as R3                               # noqa: E402
import collisions_v10 as C10                      # noqa: E402
import collisions_v8 as V8                        # noqa: E402
import xcop_static_v11 as XS                      # noqa: E402
import kids_static_v11 as KS                      # noqa: E402
import bullet_static_v11 as BS                    # noqa: E402

GIRARDI = dict(P1=(811.0, 71.0, 278.0), P2=(749.0, 88.0, 186.0), P4=(579.0, 151.0, 523.0), P5=(668.0, 187.0, 570.0))
grade = lambda z: 'pass' if abs(z) <= 2 else ('close' if abs(z) <= 3 else 'fail')


def xcop_refit(law0):
    import common as C
    raw = json.loads((RUN.ROOT / 'research_work/results/path-memory/cl2-inputs-xcop-profiles.json').read_text())['clusters']
    cls = C.Context(tier='quick', verbose=False).xcop()
    dep, _ = XS.deproject_xcop(cls, raw)
    st = [XS.to_static(c, raw[c['name']]['header']['z']) for c in dep]
    gals = RUN.load_sparc()
    a, lam, u = law0['a_code'], law0['lam'], R3.fit_u3(st, law0['a_code'], law0['lam'])
    for _ in range(3):
        a, lam = RUN.fit_a_lam(gals, u, (a, lam))
        u = R3.fit_u3(st, a, lam)
    r = R3.resid(st, lambda c: R3.cluster_M3(c, a, u, lam))
    sc = RUN.sparc_score(gals, lambda g: RUN.galaxy_g(g, a, u, lam))[0]
    f = st[0]['factors']
    return dict(a_code=a, lam=lam, u_kms=u, xcop_rms=float(R3.rms(r)), sparc_rms=float(sc),
                size_A1644=f['size'], stars_A1644=f['stars'], gas_A1644=f['gas'])


def slacs(law, variant):
    import lenses_t35 as LT
    LT.A, LT.LAM, LT.U = law['a_code'], law['lam'], law['u_kms']
    orig, reader = LT.patched_readers('project')
    J = LT.M.J
    obs = {x['Name']: x for x in orig(J.INPUT_NAMES[3])}

    def reader2(rel):
        data = reader(rel)
        if variant == 'metric' and rel == J.INPUT_NAMES[4]:
            out = []
            for g in data:
                o = obs.get(g['Name'])
                if o is None or 'conditional_Dl_Mpc' not in g: out.append(g); continue
                g = dict(g)
                g['conditional_Dl_Mpc'] = g['conditional_Dl_Mpc'] / (1 + o['zFG'])
                g['conditional_Ds_Mpc'] = g['conditional_Ds_Mpc'] / (1 + o['zBG'])       # D_ls / D_s is unchanged
                out.append(g)
            return out
        return data
    J.read_json = reader2
    try:
        lenses = [LT.M.make_lens(n) for n in LT.M.LENSES]
    finally:
        J.read_json = orig
    rows = [LT.analyse(l, 'hot companion') for l in lenses]
    gaps = np.array([r['slip_gap_dex'] for r in rows]); dm = np.array([r['dm_lens'] for r in rows])
    return dict(light_equals_matter=float(gaps.mean()), se=float(gaps.std(ddof=1) / np.sqrt(len(gaps))),
                stars_needed_dex_vs_chabrier=float(dm.mean()), salpeter_factors=[float(10 ** (x - 0.25)) for x in dm])


def collisions(law):
    base = dict(macs0025=V8.macs0025(), abell520=V8.abell520(), el_gordo=V8.el_gordo())
    t_main = dict(macs0025=0.5, abell520=0.3, el_gordo=0.46)
    out = {}
    for name, spec in base.items():
        f = C10.factors(spec['z'], 1.0)
        st = C10.rescale(spec, f, star_extra=C10.chabrier_basis(name))
        sol = C10.solve(st, law, t_main[name], f['size'])
        m = C10.measure(name, st, sol, f['size'])
        rows = C10.compare(name, m, f, V8.observed_el_gordo(spec) if name == 'el_gordo' else None)
        d = {r['check']: float(r['z']) for r in rows}
        if name == 'abell520':
            zs = [(m['speeds'][k]['500kpc'] - o) / (hi if m['speeds'][k]['500kpc'] > o else lo) for k, (o, lo, hi) in GIRARDI.items()]
            d['speeds_rms_z'] = float(np.sqrt(np.mean(np.square(zs))))
        out[name] = d
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--variants', default='fixed,metric')
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law, with_constants
    law9 = load_law('round9')
    res = dict(experiment='round 11: the gravity tests under the two static geometries', variants={})
    for variant in args.variants.split(','):
        C10.VARIANT = variant
        print(f'== variant {variant}', flush=True)
        c = xcop_refit(law9)
        law = with_constants(law9, a_code=c['a_code'], lam=c['lam'], u_kms=c['u_kms'])
        print(f"   constants: a {law['a_SI']:.4e}, g_d {law['g_d_SI']:.4e}, u {law['u_kms']:.1f} km/s; SPARC {c['sparc_rms']:.2f} km/s,"
              f" X-COP rms {c['xcop_rms']:.3f} (A1644: size x{c['size_A1644']:.3f}, stars x{c['stars_A1644']:.3f})", flush=True)
        f25 = C10.factors(0.25, 0.75, KS.WMAP9)
        consts = dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'])
        k = KS.kids(consts, law['u_kms'], law['reach_kpc'], f25); mi = KS.mistele(consts, law['u_kms'], law['reach_kpc'], f25)
        print('   KiDS (z_l 0.25) median log10(obs/pred): ' + ', '.join(f"{s} {k[s]:+.3f}" for s in ('all', 'blue', 'red', 'disc', 'bulge', 'gama')) +
              f"; gap {k['gap_model']:.3f} (obs {k['gap_observed']:.3f}); Mistele rms z spirals {mi['LTG']['rms_z']:.2f}, ellipticals {mi['ETG']['rms_z']:.2f}", flush=True)
        s = slacs(law, variant)
        print(f"   SLACS: light = matter {s['light_equals_matter']:+.3f} +- {s['se']:.3f} dex; stars needed {s['stars_needed_dex_vs_chabrier']:+.3f} dex vs Chabrier"
              f" ({min(s['salpeter_factors']):.2f}-{max(s['salpeter_factors']):.2f} x Salpeter)", flush=True)
        col = collisions(law)
        for name, d in col.items():
            print(f"   {name}: " + ', '.join(f"{kk} {v:+.2f}" for kk, v in d.items()), flush=True)
        b = BS.run_bullet(law, dict(C10.factors(0.296, 1.0, (70.0, 0.3))))
        print('   Bullet: ' + ', '.join(f"{kk} {v['value']:.3g} ({('z %+.2f' % v['z']) if 'z' in v else ('ok' if v['ok'] else 'no')})" for kk, v in b['checks'].items()), flush=True)
        res['variants'][variant] = dict(constants=dict(a_SI=law['a_SI'], g_d_SI=law['g_d_SI'], u_kms=law['u_kms']), xcop=c,
                                        kids=k, mistele=mi, slacs=s, collisions=col, bullet=b['checks'],
                                        kids_factors=dict(size=f25['size'], stars=f25['stars'], lens_mass=f25['lens']))
    C10.VARIANT = 'fixed'
    res['seconds'] = time.monotonic() - t0
    (out / 'distance_variants_v11.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'distance_variants_v11.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()
