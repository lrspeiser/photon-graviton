"""Galaxies: the 149 SPARC rotation curves (the data a and g_d were fitted on)."""
from __future__ import annotations
import numpy as np
import common as C
from checks import make, at_most, z_check

GROUP = 'galaxies'
MOND = dict(all=16.1328, train=15.8481, validation=19.8035, test=13.5164, bulge=30.3477)   # run-v3, MOND simple, km/s


def run(law, ctx):
    import run as R
    a, lam, u = law['a_code'], law['lam'], law['u_kms']
    gals = ctx.sparc()
    ours = lambda g: R.galaxy_g(g, a, u, lam)
    out = []
    ctx.log('SPARC scores')
    for split, key in ((None, 'all'), ('test', 'test'), ('validation', 'validation')):
        rms, msq = R.sparc_score(gals, ours, split)
        out.append(make(f'galaxies.sparc_rms_{key}', GROUP, f'SPARC typical speed miss ({key} galaxies)', rms,
                        crit=at_most(rms, MOND[key], 1.05 * MOND[key]), unit='km/s',
                        target=f'no worse than MOND ({MOND[key]:.2f}); within 5% is close',
                        refs='Lelli, McGaugh & Schombert 2016 (SPARC)'))
        if key == 'all':
            out.append(make('galaxies.sparc_msq', GROUP, 'SPARC mean squared log residual (the fitted statistic)', msq,
                            crit=('info', msq)))
    # median residual log10(g_obs/g_pred): overall and in the switch-off region
    K = C.K_SI
    A = np.vstack([np.c_[g['gN'] * K, (g['v'] ** 2 / g['r']) / ours(g)] for g in gals])
    lr, lg = np.log10(A[:, 1]), np.log10(A[:, 0])
    med = float(np.median(lr))
    out.append(make('galaxies.sparc_median', GROUP, 'SPARC median residual log10(observed/predicted pull)', med,
                    crit=z_check(med, 0.0, 0.02), unit='dex',
                    target='0 +- 0.02 dex (stellar mass-to-light zero point)'))
    m = (lg >= -10.25) & (lg < -9.5)
    mt = float(np.median(lr[m]))
    out.append(make('galaxies.sparc_median_switch', GROUP, 'SPARC median residual where the release switches (g_N 10^-10.25 to 10^-9.5)', mt,
                    crit=z_check(mt, 0.0, 0.02), unit='dex', target='0 +- 0.02 dex'))
    bf = np.array([np.max(g['vb2'] / np.maximum(g['vbar2'], 1e-9)) for g in gals])
    bulgy = [g for g, f in zip(gals, bf) if f > 0.5]
    rb = R.sparc_score(bulgy, ours)[0]
    out.append(make('galaxies.sparc_bulges', GROUP, f'bulge-dominated galaxies ({len(bulgy)}): typical speed miss', rb,
                    crit=at_most(rb, MOND['bulge'], 1.05 * MOND['bulge']), unit='km/s', target=f"no worse than MOND ({MOND['bulge']:.2f})"))
    return out
