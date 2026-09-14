"""CR-1 driver: a driven collective reservoir from a counted seed. See protocol.md.

    python cr1.py [--output-dir DIR] [--canonical]

Validation (V1-V4); then, for each constituent mass, the Milky Way seed under
- channel A: no bound source;
- channel B1: stimulated growth into the occupied mode;
- channel B2: growth at the local phase;
- the required-rate diagnostic C: B1 and B2 scaled to reach the reference inventory in 10 Gyr.
Units: kpc, km/s, Msun, with time in kpc/(km/s) = 0.9778 Gyr.

Fails fast: a floating-point overflow, an invalid operation or a non-finite state raises at once. Every
stage logs timestamped progress, and a run that exceeds its wall-clock budget stops the driver.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent if (HERE.parent/'capture-to-orbit').exists() else \
    Path('C:/Users/henry/Documents/Codex/photon-graviton/research_work/results')
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(RESULTS/'capture-to-orbit'))
sys.path.insert(0, str(RESULTS/'companion-extensions'))
import gpp as P  # noqa: E402
import inputs as I  # noqa: E402
import evidence_io  # noqa: E402

G = I.G
C_KMS = 299792.458
C_MS = 299792458.
GYR = 0.9777922216807891                               # Gyr per kpc/(km/s)
ALPHA = 0.0002488993286382367                          # per Mpc
ALPHA_C = ALPHA*C_KMS/1000                             # per kpc/(km/s)
KPC_M = 3.0856775814913673e19
M_SUN_KG = 1.98840987e30
L_SUN_W = 3.828e26
EV_J = 1.602176634e-19
U_CMB = 4*5.670374419e-8*2.7255**4/C_MS                # J/m^3
RHO_CMB = U_CMB*KPC_M**3/C_MS**2/M_SUN_KG              # Msun/kpc^3 mass equivalent
MASSES_EV = (1.34e-24, 1e-22)
M_SEED = 1.0e9
T_END = 10/GYR
REF = dict(E=-0.05425640259, mu=-0.16276920779, r_half=3.92510134)   # archived companion-self-binding state
CHECKPOINTS = 5
RUNS = ('A_no_bound_source', 'B1_mode', 'B2_local', 'C_B1_required_rate', 'C_B2_required_rate')
RUN_BUDGET_S = 45*60
_T0 = time.time()


def log(msg):
    print(f'[{(time.time() - _T0)/60:7.2f} min] {msg}', flush=True)


def kappa(m_ev):
    """hbar/m in kpc km/s."""
    return 1.054571817e-34*C_MS**2/(m_ev*EV_J)/(KPC_M*1e3)


def validation():
    out = {}
    grid = P.Grid(60., 2047)
    s = P.GPP(grid)
    u = s.relax(grid.r*np.exp(-grid.r**2/20), 1.)
    e, rh = s.energies(u), P.half_mass_radius(grid, u)
    rel = dict(E=e['total']/REF['E'] - 1, mu=e['mu']/REF['mu'] - 1, r_half=rh/REF['r_half'] - 1)
    vir = (2*e['kinetic'] + e['self_gravity'])/abs(e['self_gravity'])
    out['V1'] = dict(relative_errors=rel, virial=vir, passed=bool(max(map(abs, rel.values())) < 1e-4 and abs(vir) < 1e-4))
    log(f"V1 {out['V1']}")
    t_dyn = 2*np.pi*rh/np.sqrt(.5/rh)
    dt, n = .05, int(100*t_dyn/.05)
    w, E0, drift = u.copy(), e['total'], 0.
    for i in range(n):
        w = s.step(w, dt)
        if i % 5000 == 0:
            drift = max(drift, abs(s.energies(w)['total']/E0 - 1))
    drift = max(drift, abs(s.energies(w)['total']/E0 - 1))
    change = float(np.max(np.abs(np.abs(w)**2 - np.abs(u)**2))/np.max(np.abs(u)**2))
    out['V2'] = dict(dynamical_time=t_dyn, steps=n, energy_drift=drift, density_change=change,
                     passed=bool(drift < 1e-8 and change < 1e-3))
    log(f"V2 {out['V2']}")
    g3 = P.Grid(200., 4095)
    s3 = P.GPP(g3, G=0., absorb_from=150., absorb_rate=1.)
    r = g3.r
    u3 = r*np.exp(-(r - 30)**2/8 + 2j*r)
    M0 = s3.mass(u3)
    for _ in range(15000):
        u3 = s3.step(u3, .02)
    inside = 4*np.pi*g3.dr*float(np.sum(np.abs(u3[r < 150])**2))/M0
    out['V3'] = dict(remaining=s3.mass(u3)/M0, reflected_inside=inside, absorbed=s3.absorbed/M0,
                     passed=bool(s3.absorbed/M0 > .999 and inside < 1e-3))
    log(f"V3 {out['V3']}")
    s4 = P.GPP(grid)
    w, rate, n = u.copy(), .1/(300*t_dyn), int(300*t_dyn/dt)
    for _ in range(n):
        w = s4.step(w, dt, ('rate', rate))
    e4 = s4.energies(w)
    E0M = REF['E']*e4['mass']**3
    merr = e4['mass']/(1 + rate*n*dt) - 1
    ledger = (e4['mass'] - 1 - s4.injected)/s4.injected
    exc = (e4['total'] - E0M)/abs(E0M)
    out['V4'] = dict(mass_error=merr, excess_over_ground=exc, ledger_error=ledger,
                     passed=bool(abs(merr) < 1e-8 and exc < 1e-3 and abs(ledger) < 1e-6))
    log(f"V4 {out['V4']}")
    return out, u, grid


def milky_way():
    """Variant I baryons (spherically averaged), their potential, and the photon mass-equivalent density."""
    r = I.GRID
    M = I.milky_way_receivers('I')['mass']
    _, gas, _ = I._milky_way_definitions()
    m_star = M - I._disk_enclosed(sum(S*1e6*np.exp(-hole/r - r/h) for S, h, hole in gas), r)
    C = cumulative_trapezoid(np.gradient(M, r)/r, r, initial=0.)
    phi = -G*(M/r + C[-1] - C)
    L_w = m_star[-1]/.5*L_SUN_W                          # L_bol = L[3.6], a declared upper bound
    dL = L_w*np.gradient(m_star, r)/m_star[-1]           # W per kpc
    x = np.sqrt(r[1:]*r[:-1])
    c_kpc_s = C_KMS/(KPC_M/1e3)
    u = np.array([np.trapezoid(dL/r*np.log((xi + r)/np.abs(xi - r)), r)/(8*np.pi*c_kpc_s*xi) for xi in x])   # J/kpc^3
    rho_star = u/C_MS**2/M_SUN_KG
    runs = [q for q in I.milky_way_runs() if q['baryons'] == 'I' and q['rd'] == 2.6 and q['lf'] == 1]
    return dict(r=r, M=M, phi=phi, L_bol_Lsun=float(L_w/L_SUN_W), x=x, rho_star=rho_star,
                starlight_over_cmb_at_8kpc=float(np.interp(8., x, rho_star)/RHO_CMB),
                reference_inventory=float(runs[0]['total_inventory']), eilers_R=runs[0]['R'], eilers_v=runs[0]['y'],
                baryon_v=runs[0]['vb'])


def rotation(mw, R_res, M_res):
    """Reservoir contribution added in quadrature to the archived baryon speeds at the 38 Eilers radii (reported)."""
    v_res = np.sqrt(G*np.interp(mw['eilers_R'], R_res, M_res)/mw['eilers_R'])
    total = np.sqrt(mw['baryon_v']**2 + v_res**2)
    return dict(reservoir_speed_range_kms=[float(v_res.min()), float(v_res.max())],
                total_RMSE_kms=float(np.sqrt(np.mean((total - mw['eilers_v'])**2))))


def calibrate_step(grid, k, phi_b, seed, vspan, window=.3/GYR, tol=1e-4):
    """Largest potential-phase fraction f (dt = f k/vspan) whose source-free, absorber-free evolution of the
    seed conserves energy to tol over the window. Strang splitting conserves a modified energy, whose
    offset from the true energy grows as f^2."""
    gs = P.GPP(grid, kappa=k, G=G, phi_ext=phi_b)
    E0 = gs.energies(seed)['total']
    for f in (.1, .05, .025, .0125, .00625, .003125):
        dt = f*k/vspan
        u, drift = seed.copy(), 0.
        for i in range(int(np.ceil(window/dt))):
            u = gs.step(u, dt)
            if i % 50 == 0:
                drift = max(drift, abs(gs.energies(u)['total']/E0 - 1))
        drift = max(drift, abs(gs.energies(u)['total']/E0 - 1))
        log(f'  step calibration: phase fraction {f} gives energy drift {drift:.2e} over {window*GYR:.2f} Gyr')
        if drift < tol:
            return f, drift
    raise RuntimeError(f'no time step conserved energy to {tol} over the calibration window')


def evolve(sys_, u, source, dr_stop, r_outer, relax, dt_of, label):
    """Real-time run to 10 Gyr with checkpoints and an adaptive step. Following the protocol's clause "until
    the solution leaves the resolved domain", it stops if r_half < dr_stop or if r_99 reaches the absorbing
    layer at r_outer. A non-finite state or an exhausted wall-clock budget raises."""
    rows, t, stop, t_start = [], 0., None, time.time()
    for tc in np.linspace(0., T_END, CHECKPOINTS + 1):
        k = 0
        dt = dt_of(u, source)
        while t < tc - 1e-12 and stop is None:
            h = min(dt, tc - t)
            u = sys_.step(u, h, source)
            t += h
            k += 1
            if k % 50 == 0:
                if not np.all(np.isfinite(u)):
                    raise FloatingPointError(f'{label}: non-finite state at {t*GYR:.4g} Gyr')
                if time.time() - t_start > RUN_BUDGET_S:
                    raise RuntimeError(f'{label}: wall-clock budget of {RUN_BUDGET_S/60:.0f} min exhausted at {t*GYR:.4g} Gyr')
                rh = P.half_mass_radius(sys_.grid, u)
                if rh < dr_stop:
                    stop = f'half-mass radius {rh:.4g} kpc below {dr_stop:.4g} kpc at {t*GYR:.4g} Gyr'
                    break
                Mc = P._cumulative(4*np.pi*np.abs(u)**2, sys_.grid.dr)
                r99 = float(np.interp(.99*Mc[-1], Mc, sys_.grid.r))
                if r99 > r_outer:
                    stop = f'r_99 {r99:.4g} kpc reached the absorbing layer at {r_outer:.4g} kpc at {t*GYR:.4g} Gyr'
                    break
                dt = dt_of(u, source)
        e = sys_.energies(u)
        g, bound = relax(u, e['mass'])
        row = dict(t_Gyr=t*GYR, mass_Msun=e['mass'], r_half_kpc=P.half_mass_radius(sys_.grid, u),
                   central_density=float(np.abs(u[0]/sys_.grid.r[0])**2), energy=e['total'],
                   excess_over_ground_over_abs=(e['total'] - g)/abs(g), ground_energy_from=bound,
                   injected_Msun=sys_.injected, absorbed_Msun=sys_.absorbed)
        rows.append(row)
        log(f"  {label} t={row['t_Gyr']:.2f} Gyr M={row['mass_Msun']:.4g} r_half={row['r_half_kpc']:.4g} "
            f"E_exc={row['excess_over_ground_over_abs']:.3g} absorbed={row['absorbed_Msun']:.3g}" + (f' STOP {stop}' if stop else ''))
        if stop:
            break
    return u, rows, stop


def run_mass(m_ev, mw, unit_state):
    k = kappa(m_ev)
    R, N, absorb, W0 = (800., 8191, 600., 30.) if m_ev < 1e-23 else (20., 4095, 15., 600.)
    grid = P.Grid(R, N)
    phi_b = np.interp(grid.r, mw['r'], mw['phi'])
    sys_ = P.GPP(grid, kappa=k, G=G, phi_ext=phi_b, absorb_from=absorb, absorb_rate=W0)
    span = float(phi_b.max() - phi_b.min())
    seed = sys_.relax(np.exp(-grid.r/5.)*grid.r, M_SEED, dts=tuple(k/span*f for f in (4., 1., .25, .06, .02, .005)))
    es = sys_.energies(seed)
    residual = sys_.stationary_residual(seed)
    log(f'm={m_ev:g} eV seed relaxed: r_half={P.half_mass_radius(grid, seed):.4g} kpc mu={es["mu"]:.6g} residual={residual:.2e}')
    mask0 = np.abs(seed)**2 > 1e-8*np.max(np.abs(seed)**2)
    fphase, fdrift = calibrate_step(grid, k, phi_b, seed, float(np.ptp((phi_b + sys_.self_potential(seed))[mask0])))
    q_star = ALPHA_C*np.interp(grid.r, mw['x'], mw['rho_star'])
    q = q_star + ALPHA_C*RHO_CMB                                            # Msun/kpc^3 per time unit
    rate0 = sys_.collection_rate(seed, q)
    rate0_star = sys_.collection_rate(seed, q_star)
    M_ref = mw['reference_inventory']
    K = (M_ref - M_SEED)/(rate0*T_END)

    def dt_of(u, source):
        """The calibrated fraction of a radian of potential phase across the occupied region, and at most 0.2%
        mass growth per step."""
        mask = np.abs(u)**2 > 1e-8*np.max(np.abs(u)**2)
        dt = fphase*k/max(float(np.ptp((phi_b + sys_.self_potential(u))[mask])), 1.)
        if source is not None:
            dt = min(dt, .002*sys_.mass(u)/max(source[2]*sys_.collection_rate(u, source[1]), 1e-300))
        return dt

    def ground(u, mass):
        """Ground-state energy at this mass, taken as the lower of two variational upper bounds:
        - the current density relaxed in imaginary time;
        - the self-gravitating soliton of this mass placed in the baryonic potential, evaluated on the unit
          soliton's own grid, so the bound holds even when that soliton is too compact for this grid.
        Returns the energy and which bound it came from."""
        gs = P.GPP(grid, kappa=k, G=G, phi_ext=phi_b)
        v = np.abs(u).astype(complex)
        v *= np.sqrt(mass/gs.mass(v))
        mask = np.abs(v)**2 > 1e-8*np.max(np.abs(v)**2)
        s = max(float(np.ptp(gs.potential(v)[mask])), 1.)
        relaxed = gs.energies(gs.relax(v, mass, dts=tuple(k/s*f for f in (.25, .06, .02))))['total']
        ell = k*k/(G*mass)
        w = 4*np.pi*np.abs(ug)**2
        soliton = REF['E']*G*G*mass**3/(k*k) + mass*float(np.sum(w*np.interp(gu.r*ell, mw['r'], mw['phi']))/np.sum(w))
        return (relaxed, 'relaxed') if relaxed <= soliton else (soliton, 'soliton bound')

    ug, gu = unit_state                                    # unit-mass self-gravitating ground state (kappa = G = M = 1)
    Mu = P._cumulative(4*np.pi*np.abs(ug)**2, gu.dr)
    ell = k*k/(G*M_ref)                                    # length unit of the reference-mass soliton
    out = dict(m_eV=m_ev, kappa_kpc_kms=k, grid=dict(R_kpc=R, N=N, dr_kpc=grid.dr, absorb_from_kpc=absorb),
               numerics=dict(seed_stationary_residual=residual, step_phase_fraction=fphase,
                             calibration_energy_drift=fdrift),
               seed=dict(mass_Msun=es['mass'], r_half_kpc=P.half_mass_radius(grid, seed), mu=es['mu'], energy=es['total'],
                         quanta=es['mass']*M_SUN_KG*C_MS**2/(m_ev*EV_J)),
               source=dict(collection_rate_Msun_per_Gyr=rate0/GYR, starlight_share=rate0_star/rate0,
                           collected_in_10Gyr_at_initial_rate_Msun=rate0*T_END, required_rate_factor_K=K,
                           reference_inventory_Msun=M_ref),
               analytic_endpoint=dict(kaup_limit_Msun=.633*k*C_KMS/G, reference_over_kaup=M_ref/(.633*k*C_KMS/G),
                                      self_gravitating_r_half_kpc=REF['r_half']*ell,
                                      schwarzschild_radius_kpc=2*G*M_ref/C_KMS**2,
                                      rotation=rotation(mw, gu.r*ell, M_ref*Mu/Mu[-1])))
    runs = {}
    for label, source in zip(RUNS, (None, ('mode', q, 1.), ('local', q, 1.), ('mode', q, K), ('local', q, K))):
        sys_.injected = sys_.absorbed = 0.
        t0 = time.time()
        log(f'm={m_ev:g} eV {label} start')
        u, rows, stop = evolve(sys_, seed.copy(), source, 10*grid.dr, absorb, ground, dt_of, f'm={m_ev:g} {label}')
        e = sys_.energies(u)
        Mr = P._cumulative(4*np.pi*np.abs(u)**2, grid.dr)
        runs[label] = dict(checkpoints=rows, stopped=stop, final_mass_Msun=e['mass'], growth_factor=e['mass']/es['mass'],
                           escaped_fraction=sys_.absorbed/(sys_.injected + es['mass']),
                           mass_ledger_error=(e['mass'] - es['mass'] - sys_.injected + sys_.absorbed)/max(es['mass'], sys_.injected),
                           rotation=rotation(mw, grid.r, Mr), seconds=time.time() - t0)
        log(f'm={m_ev:g} eV {label} done in {runs[label]["seconds"]/60:.1f} min, growth {runs[label]["growth_factor"]:.4g}')
    out['runs'] = runs
    return out


def main():
    np.seterr(over='raise', invalid='raise', divide='raise')
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    log('validation start')
    v, ug, gu = validation()
    mw = milky_way()
    log(f"Milky Way inputs ready: reference inventory {mw['reference_inventory']:.4g} Msun")
    result = dict(validation=v, validation_passed=bool(all(x['passed'] for x in v.values())),
                  milky_way=dict(L_bol_Lsun=mw['L_bol_Lsun'], starlight_over_cmb_at_8kpc=mw['starlight_over_cmb_at_8kpc'],
                                 reference_inventory_Msun=mw['reference_inventory']),
                  masses={f'{m:g}': run_mass(m, mw, (ug, gu)) for m in MASSES_EV},
                  protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),
                  source_sha256={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ('gpp.py', 'cr1.py')},
                  runtime_seconds=time.time() - t0)
    text = json.dumps(result, indent=1) + '\n'
    print(json.dumps(dict(validation_passed=result['validation_passed'],
                          masses={m: {kk: vv for kk, vv in d.items() if kk != 'runs'} for m, d in result['masses'].items()}),
                     indent=1))
    ignore = {'/runtime_seconds'} | {f'/masses/{m:g}/runs/{r}/seconds' for m in MASSES_EV for r in RUNS}
    status = evidence_io.finish(args, 'collective-reservoir', text, HERE/'cr1-results.json', ignore=ignore)
    return status if result['validation_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
