"""CG-0: does the clock field's gradient stay timelike across the transition to Newtonian gravity? See protocol.md.

    python timelike.py [--canonical] [--output-dir DIR]

X/X0 = 1 - (a_chi / c^2 alpha)^2, where a_chi is the acceleration supplied by the static clock-field perturbation:
    F1  the total acceleration (the model contract's C2 as written);
    F2  only the excess over Newtonian gravity, with the simple nu;
    F3  a limiting-gradient kinetic function: X stays positive but tends to zero where F1 would cross it, and
        the force saturates at c^2 alpha.
"""
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
sys.path.insert(0, str(RESULTS/'capture-to-orbit'))
sys.path.insert(0, str(RESULTS/'companion-extensions'))
import evidence_io  # noqa: E402
import inputs as I  # noqa: E402

C_MS, MPC_M, KPC_M = 299792458., 3.0856775814913673e22, 3.0856775814913673e19
ALPHA = 0.0002488993286382367                  # per Mpc: the repository's redshift rate (as in cr1.py)
C2ALPHA = C_MS**2*ALPHA/MPC_M                  # m/s^2, c times the rolling rate
A_STAR = 8.563335193921255e-11                 # m/s^2: the archived simple-mu a0 (RPG-1 baryons.py)
ACC = 1e6/KPC_M                                # (km/s)^2/kpc in m/s^2
GM_SUN, AU = 1.32712440018e20, 1.495978707e11
PLANETS_AU = dict(Mercury=.387, Venus=.723, Earth=1., Mars=1.524, Jupiter=5.203, Saturn=9.537,
                  Uranus=19.19, Neptune=30.07, Kuiper_40AU=40.)
SATURN_BOUND = 1e-14                           # m/s^2, Cassini bound quoted by Turyshev & Toth 2010
FLOOR = .1                                     # declared minimum X/X0


def nu(y):
    return .5 + np.sqrt(.25 + 1/y)


def x_over_x0(g_total, g_newton):
    """X/X0 under F1, F2 and F3 at points with total and Newtonian accelerations (m/s^2)."""
    g_newton = np.maximum(np.asarray(g_newton, float), 1e-30)
    f1 = 1 - (np.asarray(g_total, float)/C2ALPHA)**2
    f2 = 1 - (g_newton*(nu(g_newton/A_STAR) - 1)/C2ALPHA)**2
    return dict(F1=f1, F2=f2, F3=np.maximum(f1, 0.))


def sparc():
    ud, ub = I.UPSILON_DISK, I.UPSILON_BULGE
    g, gn, name = [], [], []
    for d in I.sparc_galaxies():
        R, V, Vg, Vd, Vb = (d['rotmod'][:, j] for j in (0, 1, 3, 4, 5))
        use = (R > 0) & (V > 0)
        g.append(V[use]**2/R[use]*ACC)
        gn.append((Vg*np.abs(Vg) + ud*Vd*np.abs(Vd) + ub*Vb*np.abs(Vb))[use]/R[use]*ACC)
        name += [d['name']]*int(use.sum())
    g, gn, name = np.concatenate(g), np.concatenate(gn), np.array(name)
    out = {}
    for key, x in x_over_x0(g, gn).items():
        bad = x <= 0
        out[key] = dict(radii=int(len(x)), fraction_below_floor=float(np.mean(x < FLOOR)), fraction_nonpositive=float(np.mean(bad)),
                        galaxies_with_nonpositive=int(len(set(name[bad]))), min=float(x.min()))
    out['max_observed_g_over_c2alpha'] = float(g.max()/C2ALPHA)
    return out


def milky_way():
    r = I.GRID[(I.GRID >= .5) & (I.GRID <= 25.)]
    M = np.interp(r, I.GRID, I.milky_way_receivers('I')['mass'])
    gn = I.G*M/r**2*ACC
    g = nu(gn/A_STAR)*gn
    x = x_over_x0(g, gn)
    f1 = x['F1']
    cross = float(np.interp(0., f1, r)) if f1.min() < 0 < f1.max() else None   # f1 rises outward
    return dict(model='baryon model I, spherically averaged, total g from the simple nu',
                radius_where_F1_X_is_zero_kpc=cross, g_over_c2alpha_at_half_kpc=float(g[0]/C2ALPHA),
                min={k: float(v.min()) for k, v in x.items()})


def slacs():
    data = json.loads((RESULTS/'radiation-polarized-gravity/rpg1-results.json').read_text(encoding='utf-8'))
    rows = []

    def walk(o, path):
        if isinstance(o, dict):
            if 'cases' in o and 'rpg1_constant_a' in o['cases']:
                y = o['cases']['rpg1_constant_a']['gN_at_bE_over_a_star']
                gn = y*A_STAR
                x = x_over_x0(nu(y)*gn, gn)
                rows.append(dict(entry=path, gN_over_a_star=y, **{k: float(v) for k, v in x.items()}))
            for k, v in o.items():
                walk(v, f'{path}/{k}')
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, f'{path}/{i}')
    walk(data, '')
    return dict(rows=rows, min={k: min(r[k] for r in rows) for k in ('F1', 'F2', 'F3')})


def solar_system():
    gn = {p: GM_SUN/(a*AU)**2 for p, a in PLANETS_AU.items()}
    rows = {p: {k: float(v) for k, v in x_over_x0(nu(g/A_STAR)*g, g).items()} for p, g in gn.items()}
    g_sat = gn['Saturn']
    anomaly = dict(F1=float(g_sat*(nu(g_sat/A_STAR) - 1)), F2=float(g_sat*(nu(g_sat/A_STAR) - 1)), F3=C2ALPHA)
    return dict(rows=rows, anomalous_acceleration_at_saturn=anomaly, bound=SATURN_BOUND)


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    s, mw, lens, sol = sparc(), milky_way(), slacs(), solar_system()
    verdict = {}
    for k in ('F1', 'F2', 'F3'):
        floor_ok = min(s[k]['min'], mw['min'][k], lens['min'][k], min(r[k] for r in sol['rows'].values())) >= FLOOR
        solar_ok = sol['anomalous_acceleration_at_saturn'][k] < SATURN_BOUND
        verdict[k] = dict(x_floor_met=bool(floor_ok), solar_system_bound_met=bool(solar_ok), passed=bool(floor_ok and solar_ok))
    result = dict(c2alpha_m_s2=C2ALPHA, a_star_m_s2=A_STAR, xi=A_STAR/C2ALPHA,
                  X_zero_at_g_over_a_star=C2ALPHA/A_STAR, X_over_X0_at_g_equal_a_star=1 - (A_STAR/C2ALPHA)**2,
                  sparc=s, milky_way=mw, slacs=lens, solar_system=sol, verdict=verdict,
                  any_formulation_passed=any(v['passed'] for v in verdict.values()))
    text = json.dumps(result, indent=1)
    print(json.dumps(dict(verdict=verdict, sparc={k: s[k] for k in ('F1', 'F2')}, milky_way=mw), indent=1))
    return evidence_io.finish(args, 'clock-gradient-cg0', text, HERE/'timelike-results.json')


if __name__ == '__main__':
    raise SystemExit(main())
