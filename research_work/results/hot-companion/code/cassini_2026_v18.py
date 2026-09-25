"""Round 18: Cassini's limit on a Galactic distortion of the Sun's field, updated.

Park, Hees, Famaey, Desmond & Durakovic (Phys. Rev. D, 28 July 2026; arXiv:2602.17884) re-estimated the quadrupole
Q2 with the DE440 ephemeris data, simultaneously with the other ephemeris parameters: Q2 = (1.6 +- 1.8) x 10^-27 s^-2
(1 sigma), 40% tighter than Hees et al. 2014's (3 +- 3) x 10^-27, which fixed our release length L in round 9.

For the adopted law (round-12 constants, the Galaxy's Newtonian pull at the Sun from 230 km/s at 8.2 kpc), this
computes Q2 and the wide-binary boosts for a range of release lengths L (the companion leaving an emitter is released
as R(r) = 1 - exp(-r/L)), with the same exact formula as the regression suite (regression/t_precision.py), and reports
how far the adopted L = 0.15 pc sits from the new measurement and which L the new measurement allows.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'regression'))
from law_config import load_law                                     # noqa: E402
from t_precision import sun_in_galaxy, galactic_pull, MSUN, AU      # noqa: E402

PC_AU = 206264.806
NEW = dict(value=1.6e-27, error=1.8e-27, source='Park, Hees, Famaey, Desmond & Durakovic 2026, Phys. Rev. D (28 July 2026), arXiv:2602.17884')
OLD = dict(value=3e-27, error=3e-27, source='Hees, Folkner, Jacobson & Park 2014, Phys. Rev. D 89, 102002')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists():
        raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    law = load_law()
    a, gd, L0 = law['a_SI'], law['g_d_SI'], law['release_length_au']
    ge, gobs = galactic_pull(law)
    print(f"adopted law {law['name']}: a = {a:.4e}, g_d = {gd:.4e} m/s^2, L = {L0:.0f} AU; Galactic g_N at the Sun {ge:.4e} (observed {gobs:.4e})")
    rows = []
    for L_pc in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0):
        Q2, boost = sun_in_galaxy(MSUN, ge, 0.0, 0.0, L_pc * PC_AU, a, gd)
        rows.append(dict(L_pc=L_pc, L_au=L_pc * PC_AU, Q2=Q2, z_new=(Q2 - NEW['value']) / NEW['error'], z_old=(Q2 - OLD['value']) / OLD['error'],
                         boost_3000=boost[3000], boost_7000=boost[7000], boost_20000=boost[20000]))
        r = rows[-1]
        print(f"L = {L_pc:5.2f} pc: Q2 = {Q2:.2e} s^-2 (new: {r['z_new']:+.2f} sigma; old: {r['z_old']:+.2f}); wide binaries +{100 * (r['boost_7000'] - 1):.1f}% at 7,000 AU, +{100 * (r['boost_20000'] - 1):.1f}% at 20,000 AU", flush=True)
    # the release lengths at which Q2 meets the new 1 sigma and 2 sigma bounds (Q2 falls with L)
    Ls = np.array([r['L_pc'] for r in rows]); Qs = np.array([r['Q2'] for r in rows])
    bounds = {}
    for nsig in (1, 2):
        lim = NEW['value'] + nsig * NEW['error']
        ok = Qs <= lim
        if ok.any() and not ok.all():
            i = int(np.argmax(ok))
            f = (Qs[i - 1] - lim) / (Qs[i - 1] - Qs[i])
            bounds[f'{nsig}sigma'] = dict(Q2_limit=lim, L_min_pc=float(Ls[i - 1] + f * (Ls[i] - Ls[i - 1])))
        else:
            bounds[f'{nsig}sigma'] = dict(Q2_limit=lim, L_min_pc=0.0 if ok.all() else None)
    Q0, b0 = sun_in_galaxy(MSUN, ge, 0.0, 0.0, L0, a, gd)
    res = dict(experiment='round 18: Cassini Q2 against the 2026 re-estimate', law=law['name'], g_e=ge, g_obs=gobs,
               measured_2026=NEW, measured_2014=OLD,
               adopted=dict(L_au=L0, L_pc=L0 / PC_AU, Q2=Q0, z_new=(Q0 - NEW['value']) / NEW['error'], z_old=(Q0 - OLD['value']) / OLD['error'],
                            boost_7000=b0[7000], boost_20000=b0[20000]),
               release_length_bounds=bounds, sweep=rows, seconds=time.monotonic() - t0)
    print(f"adopted L = {L0 / PC_AU:.3f} pc: Q2 = {Q0:.2e}, {res['adopted']['z_new']:+.2f} sigma from the 2026 value; "
          f"L for 1 sigma >= {bounds['1sigma']['L_min_pc']}, for 2 sigma >= {bounds['2sigma']['L_min_pc']} pc")
    (out / 'cassini_2026_v18.json').write_text(json.dumps(res, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
