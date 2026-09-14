"""Post-hoc diagnostic for RPG-1 T3 (not part of the declared protocol).

How much of the SLACS lens deficit comes from PF-1's distance and stellar-mass conventions? This recomputes
the Einstein radii with flat FLRW angular distances (H0 = 70, Omega_m = 0.3, the convention of the published
population masses) and the published masses, beside the declared PF-1 values.

    python geometry-diagnostic.py

Writes geometry-diagnostic.json next to this script.
"""
import json
import sys
from pathlib import Path
from astropy.cosmology import FlatLambdaCDM

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import baryons as B  # noqa: E402
import lensing as L  # noqa: E402


def main():
    cosmo = FlatLambdaCDM(H0=70, Om0=.3)
    obs, geo, light, masses = L.load()
    published = {}
    for r in json.loads(L.INPUTS['masses'].read_text()):
        published.setdefault((r['Name'], r['imf']), r['published_log10_stellar_mass'])
    rows = []
    for name in L.LENSES:
        o, g = obs[name], geo[name]
        zl, zs = o['zFG'], o['zBG']
        conventions = dict(
            pf1=(g['conditional_Dl_Mpc']*1000, g['conditional_Dls_over_Ds'], masses),
            flrw=(cosmo.angular_diameter_distance(zl).value*1000,
                  float(cosmo.angular_diameter_distance_z1z2(zl, zs)/cosmo.angular_diameter_distance(zs)), published))
        row = dict(name=name, z_lens=zl, z_source=zs)
        for key, (Dl, ratio, table) in conventions.items():
            r, frac, _ = L.light_profile(light[name]['components'], Dl)
            row[key] = dict(Dl_Mpc=Dl/1000, Dls_over_Ds=ratio)
            for imf in ('Chabrier', 'Salpeter'):
                M = 10**table[(name, imf)]
                row[key][imf] = dict(newtonian_ratio=L.einstein_radius(M, r, frac, Dl, ratio, None)[0]/o['bSIE'],
                                     rpg1_ratio=L.einstein_radius(M, r, frac, Dl, ratio, B.A_STAR)[0]/o['bSIE'],
                                     rpg1_mass_factor=L.mass_factor(M, r, frac, Dl, ratio, B.A_STAR, o['bSIE']))
        rows.append(row)
    out = dict(scope='Post-hoc diagnostic, not part of the RPG-1 protocol.', rows=rows)
    (HERE/'geometry-diagnostic.json').write_text(json.dumps(out, indent=1) + '\n', encoding='utf-8', newline='\n')
    for row in rows:
        print(row['name'], {k: {i: round(row[k][i]['rpg1_ratio'], 3) for i in ('Chabrier', 'Salpeter')} for k in ('pf1', 'flrw')},
              'FLRW mass factor', {i: round(row['flrw'][i]['rpg1_mass_factor'], 2) for i in ('Chabrier', 'Salpeter')})
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
