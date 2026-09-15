"""CF-1: gravitational focusing and capture of a slow companion population. A conditional transport-and-capture
diagnostic; see protocol.md.

    python cf1.py [--canonical] [--output-dir DIR]        (CF1_SMOKE=1 runs a reduced set)

Potentials come from ordinary matter plus a counted seed only: the Milky Way, six SLACS lens hosts, and Coma at
both ends of its bracketed gas-and-star inputs. An isotropic bath of trial speed u is specified at an exterior
boundary R_b.

The script first validates the no-capture transport (V1-V3). It then applies one declared interaction: elastic,
equal-mass companion-companion scattering with constant sigma/m. Results are per unit incident density; supply
requirements are reported separately. A full inferred Coma potential appears only as a separately labeled inverse
diagnostic. The run fails fast and logs its progress.
"""
import json
import math
import os
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
for _p in (RESULTS/'capture-to-orbit', RESULTS/'companion-extensions', RESULTS/'supported-reservoir', HERE):
    sys.path.insert(0, str(_p))
import evidence_io  # noqa: E402
import focus as F  # noqa: E402
import scatter as S  # noqa: E402
import inputs as I  # noqa: E402
import cr2  # noqa: E402  (lens hosts: Sersic stars, Auger masses, FLRW geometry)

G = F.G
SMOKE = os.environ.get('CF1_SMOKE') == '1'
TRIAL_U = (300., 1000., 3000.)          # km/s. 3,000 was chosen after seeing a desired contrast: a trial, not a prediction
F_SEED = .01                            # counted seed: 1% of the baryons, distributed like them, in every system
T_REF = 10.                             # Gyr, a reference time (assumption)
PER_GYR = 1.0227121650537077            # 1 (km/s)/kpc in 1/Gyr
KPC_M, MSUN_KG, MP_G = 3.0856775814913673e19, 1.98847e30, 1.67262192e-24
CM2_PER_G = .1*MSUN_KG/KPC_M**2         # 1 cm^2/g in kpc^2/Msun
RHO_MEAN = 3*(70/1e3)**2*.3/(8*math.pi*G)   # Msun/kpc^3, cosmic mean matter density for H0 = 70, Om = 0.3 (comparison unit)
EV_J, HBAR, C_MS = 1.602176634e-19, 1.054571817e-34, 299792458.
T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:7.1f}s] {msg}', flush=True)


def coma_baryons(ne0, mstar, r):
    """Beta-model gas (beta 0.75, r_c 296 kpc, mu_e 1.17) truncated at 3 Mpc, plus stars distributed like the gas."""
    beta, rc, rt = .75, 296., 3000.
    rho0 = 1.17*MP_G*ne0*(KPC_M*100)**3/(MSUN_KG*1000)       # g/cm^3 to Msun/kpc^3
    rho = np.where(r <= rt, rho0*(1 + (r/rc)**2)**(-1.5*beta), 0.)
    gas = cumulative_trapezoid(4*np.pi*r*r*rho, r, initial=0) + 4/3*np.pi*r[0]**3*rho[0]
    return gas*(1 + mstar/gas[-1]), float(gas[-1])


def nfw_total(r, M200=1e15, c=4.):
    """Inverse diagnostic only: an NFW total potential of literature scale, truncated at R200 so that its mass is the
    quoted M200 (continuing it to R_b = 30 Mpc would add 2.8e15 Msun)."""
    R200 = (3*M200/(4*math.pi*200*RHO_MEAN/.3))**(1/3)
    m = lambda x: np.log1p(x) - x/(1 + x)
    return M200*m(np.minimum(r, R200)/(R200/c))/m(c), R200


def systems():
    out = [dict(name='Milky Way (baryon model I)', r=I.GRID, M=I.milky_way_receivers('I')['mass'], R_b=300., kind='galaxy')]
    for name in (cr2.LENSES[:1] if SMOKE else cr2.LENSES):
        L = cr2.Lens(name, 'G1_flat_FLRW')
        out.append(dict(name=f'{name} (stars, Auger Chabrier, FLRW)', r=L.model.r, M=L.pop['Chabrier']*L.frac, R_b=300., kind='lens host'))
    r = np.geomspace(.1, 3e4, 6000)
    for ne0, ms in ((2.5e-3, .5e13), (4.5e-3, 2e13)):
        M, gas = coma_baryons(ne0, ms, r)
        out.append(dict(name=f'Coma (gas n_e0={ne0:g}, stars {ms:.1e})', r=r, M=M, R_b=3e4, kind='cluster', gas_mass=gas))
    M, R200 = nfw_total(r)
    out.append(dict(name='Coma inverse diagnostic (NFW 1e15, c=4)', r=r, M=M, R_b=3e4, kind='inverse diagnostic', R200=R200))
    return out


def seed_profile(pot):
    """The seed follows the baryons: density rho_b(r) and isotropic Jeans dispersion s(r) in the host potential."""
    r = pot.r
    rho = np.maximum(np.gradient(pot.M, r)/(4*np.pi*r*r), 0.)
    g = G*pot.M/r**2
    press = cumulative_trapezoid((rho*g)[::-1], r[::-1], initial=0)[::-1]
    s2 = np.where(rho > 0, -press/np.where(rho > 0, rho, 1.), 0.)
    return rho, np.sqrt(np.maximum(s2, 0.))


def kernel(pot, u, radii, rho, s, n):
    """Growth and heating kernels per unit (sigma/m) rho_inf, with the seed distributed like the baryons:
    K = int (n_in/n_inf) <|v_rel| net> rho_b dV / M_b  [km/s],   H = same with the bound energy change [(km/s)^3]."""
    rows = []
    for rr in radii:
        ve = math.sqrt(pot.v_esc2(rr))
        sr = float(np.interp(math.log(rr), pot.lr, s))
        o = S.outcomes(u, ve, sr, n=n)
        rows.append(dict(r_kpc=float(rr), v_esc=ve, s=sr, focus_density=math.sqrt(1 + ve*ve/(u*u)),
                         entry_factor=1 + ve*ve/(u*u), capture=o['capture'], eject=o['eject'], net=o['net_retained'],
                         rate_net=o['rate_weighted_net'], rate_heat=o['bound_energy_change']*o['mean_relative_speed']))
    r = np.array([q['r_kpc'] for q in rows])
    w = np.interp(np.log(r), pot.lr, rho)*4*np.pi*r**3          # rho dV per d ln r
    mb = float(np.trapezoid(w, np.log(r)))
    K = float(np.trapezoid(w*np.array([q['focus_density']*q['rate_net'] for q in rows]), np.log(r)))/mb
    H = float(np.trapezoid(w*np.array([q['focus_density']*q['rate_heat'] for q in rows]), np.log(r)))/mb
    bind = float(np.trapezoid(w*(np.interp(np.log(r), pot.lr, pot.phi) + 1.5*np.interp(np.log(r), pot.lr, s)**2), np.log(r)))/mb
    return dict(K_kms=K, H_kms3=H, seed_specific_energy=bind, rows=rows)


def growth(Ktab_u, Ktab_K, u, exposures):
    """Retained mass with feedback: dq/dt = X K(q) q with K(q) = sqrt(1+q) K0(u/sqrt(1+q)), q0 = F_SEED.
    X = (sigma/m) rho_inf per (km/s) in 1/Gyr units; integrated over T_REF."""
    lu = np.log(Ktab_u)

    def K(q):
        ue = u/math.sqrt(1 + q)
        return math.sqrt(1 + q)*float(np.interp(math.log(ue), lu, Ktab_K))
    out = []
    for X in exposures:
        def rhs(t, y):
            q = math.exp(y[0])
            return [X*K(q)*PER_GYR]
        def big(t, y):
            return y[0] - math.log(1e3)
        big.terminal = True
        sol = solve_ivp(rhs, (0, T_REF), [math.log(F_SEED)], rtol=1e-8, atol=1e-10, events=big)
        out.append(dict(exposure=X, q_final=math.exp(float(sol.y[0][-1])), runaway=bool(sol.status == 1),
                        t_end_Gyr=float(sol.t[-1])))
    return out


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    n_mc = 20000 if SMOKE else 60000
    result = dict(scope='Conditional transport-and-capture diagnostic; per unit incident density; potentials from '
                        'baryons plus a counted seed; 3,000 km/s is a trial parameter.', systems={})
    ok = True
    for sysd in systems():
        pot = F.Potential(sysd['r'], sysd['M'], sysd['R_b'], sysd['name'])
        rho, s = seed_profile(pot)
        Mb = float(pot.M[-1])
        r_half = float(np.interp(.5*Mb, pot.M, pot.r))
        vesc_inf = lambda rr: math.sqrt(pot.v_esc2(rr) + 2*G*Mb/pot.R_b)
        entry = dict(mass_within_R_b_Msun=Mb, r_half_kpc=r_half,
                     v_esc_centre_kms=math.sqrt(pot.v_esc2(pot.r[1])), v_esc_r_half_kms=math.sqrt(pot.v_esc2(r_half)),
                     v_esc_centre_to_infinity_kms=vesc_inf(pot.r[1]), kind=sysd['kind'])
        log(f"{sysd['name']}: M={Mb:.3g} Msun, r_half={r_half:.3g} kpc, v_esc centre {entry['v_esc_centre_kms']:.0f} km/s")
        # V1-V3: the no-capture transport
        # V1 compares the orbits with the exact effective-potential entry area; the analytic formula
        # pi R^2 (1 + v_esc^2/u^2) is recorded beside it and applies only when no barrier or boundary limits entry
        v1 = {f'{u:g}': abs(F.entry_by_orbits(pot, r_half, u)/F.entry_barrier(pot, r_half, u) - 1) for u in TRIAL_U}
        v1_formula = {f'{u:g}': F.entry_analytic(pot, r_half, u)/F.entry_barrier(pot, r_half, u) - 1 for u in TRIAL_U}
        v2 = {f'{f:g} r_half': abs(F.density_by_orbits(pot, f*r_half, TRIAL_U[1])/F.density_analytic(pot, f*r_half, TRIAL_U[1]) - 1)
              for f in (.1, .3, 1., 3., 10.) if f*r_half < pot.R_b}
        v3 = {f'{u:g}': F.no_capture_retained(pot, u, n=20 if SMOKE else 60) for u in TRIAL_U}
        passed = max(v1.values()) < 1e-3 and max(v2.values()) < 1e-3 and max(v3.values()) <= 1e-12
        ok &= passed
        entry['validation'] = dict(V1_entry_relative=v1, V1_formula_over_exact_minus_1=v1_formula, V2_density_relative=v2,
                                   V3_retained_fraction=v3, passed=bool(passed))
        log(f"  V1 max {max(v1.values()):.1e}  V2 max {max(v2.values()):.1e}  V3 max {max(v3.values()):.1e}  passed={passed}")
        radii = np.geomspace(max(pot.r[1], 1e-2*r_half), min(20*r_half, .9*pot.R_b), 16 if SMOKE else 32)
        # kernel table over effective bath speeds, for the trial speeds and for feedback
        u_tab = np.geomspace(30., 1e4, 8 if SMOKE else 22)
        tab = [kernel(pot, ue, radii, rho, s, n_mc // 4) for ue in u_tab]
        entry['kernel_table'] = dict(u_kms=u_tab.tolist(), K_kms=[t['K_kms'] for t in tab], H_kms3=[t['H_kms3'] for t in tab])
        entry['trials'] = {}
        for u in TRIAL_U:
            k = kernel(pot, u, radii, rho, s, n_mc)
            heat_per_retained = k['H_kms3']/k['K_kms'] if k['K_kms'] > 0 else None
            exps = np.geomspace(1e-6, 1e2, 25 if not SMOKE else 7)
            Ktab = np.array([t['K_kms'] for t in tab])
            g = growth(u_tab, Ktab, u, exps)
            reach = None
            j = next((i for i, q in enumerate(g) if q['q_final'] >= 1.), None)
            if j is not None:     # the exposure at which the retained mass reaches the baryonic mass, solved exactly
                lo = math.log(exps[j - 1]) if j > 0 else math.log(exps[0]) - 5
                reach = math.exp(brentq(lambda lx: math.log(growth(u_tab, Ktab, u, [math.exp(lx)])[0]['q_final']),
                                        lo, math.log(exps[j]), xtol=1e-6))
            entry['trials'][f'{u:g}'] = dict(
                K_kms=k['K_kms'], H_kms3=k['H_kms3'], seed_specific_energy=k['seed_specific_energy'],
                heat_per_net_retained_over_binding=(heat_per_retained/abs(k['seed_specific_energy']) if heat_per_retained else None),
                growth=g, exposure_for_retained_equal_to_baryons=reach,
                incident_density_for_that_at_1cm2g_over_mean=(reach/CM2_PER_G/RHO_MEAN if reach else None),
                radial=k['rows'])
            log(f"  u={u:g}: K={k['K_kms']:.4g} km/s, heating/binding per retained="
                f"{entry['trials'][f'{u:g}']['heat_per_net_retained_over_binding']}, exposure for q=1: {reach}")
        entry['wavelength_kpc'] = {f'{m:g} eV': {f'{u:g}': 2*math.pi*HBAR*C_MS**2/(m*EV_J)/(KPC_M*1e3)/u for u in TRIAL_U}
                                   for m in (1.34e-24, 1e-22)}
        result['systems'][sysd['name']] = entry
        if time.time() - T0 > 60*60:
            raise RuntimeError('wall-clock budget exhausted')
    gal = [v for v in result['systems'].values() if v['kind'] == 'galaxy']
    cl = [v for v in result['systems'].values() if v['kind'] == 'cluster']
    # signed kernels; a ratio is meaningful only when both systems retain (two eroding systems give a positive ratio)
    result['contrast'] = {}
    for u in TRIAL_U:
        kg = gal[0]['trials'][f'{u:g}']['K_kms']
        kc = [c['trials'][f'{u:g}']['K_kms'] for c in cl]
        result['contrast'][f'{u:g}'] = dict(milky_way_K_kms=kg, coma_K_kms=kc,
                                            coma_over_milky_way=[x/kg if (x > 0 and kg > 0) else None for x in kc])
    result['all_validation_passed'] = bool(ok)
    result['runtime_seconds'] = time.time() - T0
    text = json.dumps(result, indent=1, default=float)
    log(f'validation passed={ok}; contrast {result["contrast"]}')
    if SMOKE:
        return 0
    status = evidence_io.finish(args, 'gravitational-focusing-cf1', text, HERE/'cf1-results.json', ignore={'/runtime_seconds'})
    return status if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
