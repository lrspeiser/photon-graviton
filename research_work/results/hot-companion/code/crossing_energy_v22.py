"""Round 22, step 2a: the energy bill of the heat of crossing, and of the fast glow.

The law's companion carries l = a u / 2 watts per kilogram from cold matter (a = 2 l / u), and (1 + k) l from warm
matter. Round 16's heat of crossing gives the Bullet's stars k(t) = p (v^2 - 2 u v cos theta)/u^2 for the time they
spend inside the other system's flow; round 21 found that a glow spreading at about 600 km/s fits the Bullet only if
the collision puts about 3.3 times that power into it. Here the energy per kilogram of stars is integrated over the
crossing history (round 16's straight pass: receding at 3,900 km/s for 213 Myr, approaching at 3,000 km/s before),
shell by shell, and compared with:
  * the kinetic energy of the relative motion per kilogram of the smaller cluster (in the centre-of-mass frame);
  * the cold companion's emission over 13 Gyr (l x 13 Gyr);
  * what the settled heat of a cluster galaxy already costs over 10 Gyr (k = 3 sigma^2/u^2 at sigma = 1,000 km/s);
  * the stars' own rest-mass energy and light output (a solar-type 1 L_sun per M_sun, for scale).

    python code/crossing_energy_v22.py --output run-crossing-frame-v22/crossing_energy_v22.json
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import crossing_heat_v16 as CH                    # noqa: E402
import collisions_v10 as C10                      # noqa: E402
import crossing_frame_v22 as CF                   # noqa: E402

C_SI = 2.99792458e8
MYR_S = 3.15576e13
LSUN_PER_MSUN = 3.828e26 / 1.989e30               # W/kg


class _Captured(Exception):
    pass


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    from law_config import load_law
    import common as C
    law = load_law(); ctx = C.Context(verbose=False); C.apply_distances(law, ctx)
    f = C10.factors(0.296, 1.0, (70.0, 0.3))
    cap = {}

    def capture(cur, ghost_gas, ghost_stars, pos, consts, crossing=(), **k):
        cap['crossing'] = crossing; cap['cur'] = cur; cap['ghost_gas'] = ghost_gas
        raise _Captured

    orig_make, orig_map = CH.make_kernel, CH.kappa_map_cross
    CH.make_kernel, CH.kappa_map_cross = CF.FrameKernel, capture
    try:
        CH.run_bullet_cross(law, dict(f), cross=dict(v_out=3900.0, v_in=3000.0, b=150.0), n=64, dx=45.0)
    except _Captured:
        pass
    finally:
        CH.make_kernel, CH.kappa_map_cross = orig_make, orig_map
    a_si = float(law.get('a_SI', 6.297889390049439e-11)); u = float(law['u_kms'])
    ell = a_si * u * 1e3 / 2                                             # W/kg
    t = np.linspace(0.0, 1200.0, 24001); dt = t[1] - t[0]
    rows = []
    for comps, (r0, r1), kf in cap['crossing']:
        k = kf.k_of_t(t)
        E = float(np.sum(k) * dt * MYR_S * ell)                          # J/kg
        rows.append(dict(system=comps[0]['centre'], shell_kpc=[float(r0), float(r1)], k_max=float(k.max()),
                         k_time_integral_myr=float(np.sum(k) * dt), energy_J_per_kg=E,
                         peak_power_W_per_kg=float(k.max() * ell)))
    # the kinetic energy of the relative motion, per kg of the smaller cluster, in the centre-of-mass frame
    cur = cap['cur']
    M_main = sum(c['M'] for kname, c in cur.items() if kname.endswith('main') or kname.startswith('st_main'))
    M_sub = sum(c['M'] for kname, c in cur.items() if kname.endswith('sub') and not kname.startswith('st_main'))
    frac = M_main / (M_main + M_sub)
    ke = {f'{v:.0f} km/s': 0.5 * (frac * v * 1e3) ** 2 for v in (3000.0, 3900.0)}
    cold_13 = ell * 13e3 * MYR_S
    settled = 3 * (1000.0 / u) ** 2 * ell * 10e3 * MYR_S
    inner_sub = [r for r in rows if r['system'] == 'sub_bcg'][0]
    out = dict(experiment='round 22: the energy bill of the heat of crossing (Bullet, round 16 history)',
               ell_W_per_kg=ell, rows=rows, masses=dict(main=M_main, sub=M_sub, main_fraction=frac),
               kinetic_energy_sub_J_per_kg=ke, cold_companion_13Gyr_J_per_kg=cold_13,
               settled_cluster_heat_10Gyr_J_per_kg=settled, rest_energy_J_per_kg=C_SI ** 2,
               sunlike_light_W_per_kg=LSUN_PER_MSUN,
               summary=dict(inner_sub_energy=inner_sub['energy_J_per_kg'],
                            inner_sub_energy_x3p3=3.3 * inner_sub['energy_J_per_kg'],
                            over_kinetic_3900=inner_sub['energy_J_per_kg'] / ke['3900 km/s'],
                            x3p3_over_kinetic_3900=3.3 * inner_sub['energy_J_per_kg'] / ke['3900 km/s'],
                            over_cold_13Gyr=inner_sub['energy_J_per_kg'] / cold_13,
                            x3p3_over_settled_10Gyr=3.3 * inner_sub['energy_J_per_kg'] / settled,
                            x3p3_rest_fraction=3.3 * inner_sub['energy_J_per_kg'] / C_SI ** 2,
                            peak_power_over_sunlike=inner_sub['peak_power_W_per_kg'] / LSUN_PER_MSUN))
    for r in rows:
        print(f"{r['system']:8s} shell {r['shell_kpc'][0]:7.0f}-{r['shell_kpc'][1]:9.0f} kpc: k max {r['k_max']:6.1f}, "
              f"integral {r['k_time_integral_myr']:8.0f} Myr, energy {r['energy_J_per_kg']:.2e} J/kg, peak power "
              f"{r['peak_power_W_per_kg']:.2e} W/kg")
    s = out['summary']
    print(f"l = {ell:.3e} W/kg; smaller cluster's kinetic energy (CM frame): "
          + ', '.join(f'{k} {v:.2e} J/kg' for k, v in ke.items()))
    print(f"inner smaller-half stars: {s['inner_sub_energy']:.2e} J/kg = {s['over_kinetic_3900']:.1f} x kinetic (3,900 km/s), "
          f"{s['over_cold_13Gyr']:.1f} x the cold companion over 13 Gyr; x3.3: {s['inner_sub_energy_x3p3']:.2e} J/kg = "
          f"{s['x3p3_over_kinetic_3900']:.1f} x kinetic, {s['x3p3_over_settled_10Gyr']:.2f} x a cluster galaxy's settled heat over "
          f"10 Gyr, {s['x3p3_rest_fraction']:.1e} of rest energy; peak power (heat rule) {s['peak_power_over_sunlike']:.1f} x a Sun-like star's light per kg")
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()
