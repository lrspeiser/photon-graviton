"""RPG-1 driver: solver validation (V1), SPARC rotation (T1), Milky Way rotation and vertical force (T2),
lensing (T3), field energy and outer boundary (T4) and the PF-1 link (T5). See protocol.md.

    python rpg1.py [--output-dir DIR] [--canonical] [--workers N]

Regenerates rpg1-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive. Exits nonzero if a V1 tolerance fails or the regenerated
numbers differ from the archive.
"""
import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent if (HERE.parent/'capture-to-orbit').exists() else \
    Path('C:/Users/henry/Documents/Codex/photon-graviton/research_work/results')
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(RESULTS/'companion-extensions'))
import aqual as Q  # noqa: E402
import baryons as B  # noqa: E402
import evidence_io  # noqa: E402
import lensing as L  # noqa: E402

A = B.A_STAR
G = Q.G
C_MS = 299792458.
MPC_M = 3.0856775814913673e22
KPC_M = Q.KPC_M
L_SUN_W = 3.828e26
U_CMB = 4*5.670374419e-8*2.7255**4/C_MS                 # J/m^3
SN_C_ALPHA = 69.75912160673488                           # PF-1 T4 calibrated fit, km/s/Mpc
KZ_UNIT = 2*np.pi*G*1e6                                  # (km/s)^2/kpc per Msun/pc^2
SPLITS = ('train', 'validation', 'test')
_GALAXIES = {}


def _init():
    for d in B.I.sparc_galaxies():
        _GALAXIES[d['name']] = d


def solve_galaxy(task):
    name, hz_factor, refine, r_out_factor, full = task
    d = _GALAXIES[name]
    comp = B.sparc_components(d)
    M = B.sparc_total_mass(comp)
    grid, r_m = B.sparc_grid(d, comp, M, refine, r_out_factor)
    m, m_in = B.sparc_masses(grid, comp, hz_factor*comp['rd'])
    R = d['r']
    newton = Q.Solver(grid, m, m_in)
    newton.solve(0.)
    vN = newton.midplane_speed(R)
    gN = vN**2/R
    s = Q.Solver(grid, m, m_in)
    s.solve(A)
    out = dict(name=name, split=d['split'], R_kpc=R.tolist(), observed_kms=d['y'].tolist(),
               archived_algebraic_kms=d['mond_raw'].tolist(), newtonian_kms=vN.tolist(),
               algebraic_kms=np.sqrt(R*Q.nu(gN/A)*gN).tolist(), aqual_kms=s.midplane_speed(R).tolist(),
               iterations=len(s.history), converged=bool(s.converged), residual=s.residual,
               baryonic_mass_Msun=M, grid_mass_ratio=s.m_tot/M, r_M_kpc=float(r_m), grid=[grid.nr, grid.nth])
    if full:
        vf2 = float(np.sqrt(G*M*A))
        E, EN = s.energy_profile(), newton.energy_profile()
        lr = np.log(grid.rf)
        E10, E100 = np.interp(np.log([10*r_m, 100*r_m]), lr, E)
        phi_eq = s.phi.reshape(grid.nr, grid.nth)[:, -1]
        P10, P100 = np.interp(np.log([10*r_m, 100*r_m]), np.log(grid.rc), phi_eq)
        L_w = d['catalog']['L9']*1e9*L_SUN_W
        R_u = np.sqrt(L_w/(4*np.pi*C_MS*U_CMB))/KPC_M
        out.update(v_f_kms=vf2**.5, energy=dict(
            dE_dlnr_10_100_rM_over_third_Mvf2=float((E100 - E10)/np.log(10)/(M*vf2/3)),
            E_F_inside_10rM_over_Mvf2=float(E10/(M*vf2)),
            newtonian_gradient_energy_over_Mvf2=float(EN[-1]/(M*vf2)),
            potential_rise_10_100_rM_over_vf2_ln10=float((P100 - P10)/(vf2*np.log(10)))),
            radiation=dict(R_u_kpc=float(R_u), R_u_over_last_radius=float(R_u/R[-1]), R_u_over_r_M=float(R_u/r_m)))
    return out


def scores(rows, key):
    out = {}
    for split in SPLITS:
        dd = [r for r in rows if r['split'] == split]
        mse = [np.mean((np.array(r[key]) - np.array(r['observed_kms']))**2) for r in dd]
        le = [np.mean(np.log10(np.maximum(np.array(r[key]), 1e-3)/np.array(r['observed_kms']))**2) for r in dd]
        out[split] = dict(galaxies=len(dd), RMSE_kms=float(np.sqrt(np.mean(mse))), log_RMS=float(np.sqrt(np.mean(le))))
    return out


def difference_stats(rows, a, b):
    d = np.concatenate([np.array(r[a]) - np.array(r[b]) for r in rows])
    rel = np.concatenate([(np.array(r[a]) - np.array(r[b]))/np.array(r[b]) for r in rows])
    q = np.percentile(d, [5, 50, 95])
    return dict(points=len(d), p05_kms=float(q[0]), median_kms=float(q[1]), p95_kms=float(q[2]),
                max_abs_kms=float(np.max(np.abs(d))), median_relative=float(np.median(rel)),
                max_abs_relative=float(np.max(np.abs(rel))))


def validation():
    out = {}
    M, b = 1e10, 1.
    grid = Q.Grid(1e-3, 1e4, 400, 48)
    s = Q.Solver(grid, *Q.spherical_cell_masses(grid, lambda r: M*r**3/(r*r + b*b)**1.5))
    s.solve(A)
    gr = s.gradients(s.phi)[0]
    gN = G*M*grid.rf/(grid.rf**2 + b*b)**1.5
    ex = Q.nu(gN/A)*gN
    sel = (grid.rf >= .01*b) & (grid.rf <= 1000*b)
    err = float(max(np.max(np.abs(gr[sel, -1]/ex[sel] - 1)), np.max(np.abs(gr[sel, 0]/ex[sel] - 1))))
    out['plummer_aqual'] = dict(max_relative_error=err, iterations=len(s.history), residual=s.residual,
                                tolerance=1e-6, passed=err < 1e-6)
    grid = Q.Grid(1e-3, 3e4, 420, 80)
    R = np.geomspace(.3, 60, 40)
    Rz = np.linspace(4, 9, 11)
    z = np.full_like(Rz, 1.1)
    M, ak = 5e10, 3.
    thin = Q.thin_disk_cell_masses(grid, lambda x: M*(1 - ak/np.sqrt(x*x + ak*ak)))
    s = Q.Solver(grid, *thin)
    s.solve(A)
    d = np.sqrt(R*R + ak*ak)
    vex = np.sqrt(R*Q.nu(G*M/d**2/A)*G*M*R/d**3)
    dz = s.gradient_at(Rz, z)[1]
    dd = np.sqrt(Rz**2 + (z + ak)**2)
    kz = Q.nu(G*M/dd**2/A)*G*M*(z + ak)/dd**3
    ev, ek = float(np.max(np.abs(s.midplane_speed(R)/vex - 1))), float(np.max(np.abs(dz/kz - 1)))
    out['kuzmin_aqual_exact'] = dict(speed_max_relative_error=ev, kz_max_relative_error=ek, iterations=len(s.history),
                                     residual=s.residual, tolerance=5e-3, passed=ev < 5e-3 and ek < 5e-3)
    M, a_, b_ = 5e10, 3., .3

    def rho(x, zz):
        zb = np.sqrt(zz*zz + b_*b_)
        return b_*b_*M/(4*np.pi)*(a_*x*x + (a_ + 3*zb)*(a_ + zb)**2)/((x*x + (a_ + zb)**2)**2.5*zb**3)
    s = Q.Solver(grid, *Q.smooth_cell_masses(grid, rho))
    s.solve(0.)
    vex = np.sqrt(G*M*R*R/(R*R + (a_ + b_)**2)**1.5)
    zb = np.sqrt(z*z + b_*b_)
    kz = G*M*z*(a_ + zb)/(zb*(Rz*Rz + (a_ + zb)**2)**1.5)
    ev = float(np.max(np.abs(s.midplane_speed(R)/vex - 1)))
    ek = float(np.max(np.abs(s.gradient_at(Rz, z)[1]/kz - 1)))
    out['miyamoto_nagai_newtonian'] = dict(speed_max_relative_error=ev, kz_max_relative_error=ek,
                                           grid_mass_ratio=s.m_tot/M, tolerance=5e-3, passed=ev < 5e-3 and ek < 5e-3)
    return out


def convergence_selection(galaxies):
    """V1(d) galaxies from baryon properties only: gas-richest, median baryonic mass, most bulge-dominated."""
    props = []
    for d in galaxies:
        c = B.sparc_components(d)
        ms = B._cylinder_mass(c['sigma_star'], 1e4*c['rd'])
        mg = B._cylinder_mass(c['sigma_gas'], 1e4*max(c['rd'], c['h_gas'])) if c['sigma_gas'] else 0.
        mb = float(c['m_bulge'](np.array([1e6]))[0]) if c['m_bulge'] else 0.
        props.append((d['name'], mg/(ms + mg + mb), mb/(ms + mg + mb), ms + mg + mb))
    order = np.argsort([p[3] for p in props])
    return dict(gas_richest=max(props, key=lambda p: p[1])[0], median_mass=props[int(order[len(order)//2])][0],
                most_bulge_dominated=max(props, key=lambda p: p[2])[0])


def milky_way(variant, inputs):
    grid = Q.Grid(1e-2, 1000*max(25., np.sqrt(G*1.2e11/A)), int(np.ceil(np.log(1000*max(25., np.sqrt(G*1.2e11/A))/1e-2)
                                                                             /B.DLNR)), B.NTH)
    m, m_in, nominal = B.milky_way_masses(grid, variant)
    R = np.array([r['R_kpc'] for r in inputs['eilers']['rows']])
    v = np.array([r['vc_kms'] for r in inputs['eilers']['rows']])
    Rb = np.array([r['R_kpc'] for r in inputs['bovy']['rows']])
    kz_obs = np.array([r['Kz_over_2piG_Msun_pc2'] for r in inputs['bovy']['rows']])
    z = np.full_like(Rb, 1.1)
    newton = Q.Solver(grid, m, m_in)
    newton.solve(0.)
    s = Q.Solver(grid, m, m_in)
    s.solve(A)
    vN = newton.midplane_speed(R)
    gN = vN**2/R
    dRN, dzN = newton.gradient_at(Rb, z)
    preds = dict(newtonian=(vN, np.abs(dzN)/KZ_UNIT),
                 algebraic=(np.sqrt(R*Q.nu(gN/A)*gN), Q.nu(np.hypot(dRN, dzN)/A)*np.abs(dzN)/KZ_UNIT),
                 aqual=(s.midplane_speed(R), np.abs(s.gradient_at(Rb, z)[1])/KZ_UNIT))
    out = dict(variant=variant, grid_mass_Msun=s.m_tot, nominal_mass_Msun=nominal['stellar'] + nominal['gas'],
               iterations=len(s.history), residual=s.residual)
    for k, (vv, kk) in preds.items():
        out[k] = dict(rotation_RMSE_kms=float(np.sqrt(np.mean((vv - v)**2))), rotation_bias_kms=float(np.mean(vv - v)),
                      kz_RMS_Msun_pc2=float(np.sqrt(np.mean((kk - kz_obs)**2))), kz_bias_Msun_pc2=float(np.mean(kk - kz_obs)),
                      kz_chi2=float(np.sum(((kk - kz_obs)/np.array([r['Kz_error'] for r in inputs['bovy']['rows']]))**2)),
                      speeds_kms=vv.tolist(), kz_Msun_pc2=kk.tolist())
    E = s.energy_profile()
    r_m = np.sqrt(G*s.m_tot/A)
    vf2 = np.sqrt(G*s.m_tot*A)
    E10, E100 = np.interp(np.log([10*r_m, 100*r_m]), np.log(grid.rf), E)
    out['energy'] = dict(r_M_kpc=float(r_m), dE_dlnr_10_100_rM_over_third_Mvf2=float((E100 - E10)/np.log(10)/(s.m_tot*vf2/3)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path)
    ap.add_argument('--canonical', action='store_true')
    ap.add_argument('--workers', type=int, default=min(8, os.cpu_count() or 1))
    args = ap.parse_args()
    t0 = time.time()
    _init()
    galaxies = list(_GALAXIES.values())
    names = [d['name'] for d in galaxies]
    v1 = validation()
    pick = convergence_selection(galaxies)
    tasks = ([(n, .1, 1, 1., True) for n in names] + [(n, .2, 1, 1., False) for n in names]
             + [(n, .1, 2, 1., False) for n in pick.values()] + [(n, .1, 1, 10., False) for n in pick.values()])
    with ProcessPoolExecutor(args.workers, initializer=_init) as pool:
        solved = list(pool.map(solve_galaxy, tasks, chunksize=1))
    main_rows, thick = solved[:len(names)], solved[len(names):2*len(names)]
    refined, far = solved[2*len(names):2*len(names) + 3], solved[2*len(names) + 3:]
    by_name = {r['name']: r for r in main_rows}
    conv = {}
    for (label, n), rr, ff in zip(pick.items(), refined, far):
        base = np.array(by_name[n]['aqual_kms'])
        dres, dout = float(np.max(np.abs(np.array(rr['aqual_kms']) - base))), float(np.max(np.abs(np.array(ff['aqual_kms']) - base)))
        conv[label] = dict(galaxy=n, max_change_doubled_resolution_kms=dres, max_change_r_out_x10_kms=dout,
                           residual=by_name[n]['residual'],
                           passed=dres < .5 and dout < .1 and by_name[n]['residual'] < 1e-8)
    v1['sparc_convergence'] = conv
    v1_passed = all(v['passed'] for k, v in v1.items() if k != 'sparc_convergence') and all(c['passed'] for c in conv.values())
    t1 = dict(scores={k: scores(main_rows, k) for k in ('newtonian_kms', 'algebraic_kms', 'aqual_kms', 'archived_algebraic_kms')},
              thickness_0p2_Rd_scores={k: scores(thick, k) for k in ('newtonian_kms', 'algebraic_kms', 'aqual_kms')},
              aqual_minus_algebraic_same_baryons=difference_stats(main_rows, 'aqual_kms', 'algebraic_kms'),
              algebraic_minus_archived_algebraic=difference_stats(main_rows, 'algebraic_kms', 'archived_algebraic_kms'),
              thickness_0p2_minus_0p1_aqual=difference_stats([dict(a=t['aqual_kms'], b=m_['aqual_kms'])
                                                               for t, m_ in zip(thick, main_rows)], 'a', 'b'),
              all_converged=bool(all(r['converged'] for r in main_rows + thick)),
              max_residual=float(max(r['residual'] for r in main_rows + thick)),
              max_grid_mass_error=float(max(abs(r['grid_mass_ratio'] - 1) for r in main_rows)))
    inputs = json.loads((RESULTS/'milky-way-capture/inputs.json').read_text())
    t2 = {v: milky_way(v, inputs) for v in ('I', 'II')}
    lens_rows = L.run(A)
    summaries = []
    for imf in ('Chabrier', 'Salpeter'):
        for case in ('rpg1_constant_a', 'rpg1_linked_a', 'newtonian_baryons'):
            ratios = [row[imf]['cases'][case]['ratio_to_observed'] for row in lens_rows]
            summaries.append(dict(imf=imf, case=case, ratios=ratios,
                                  all_within_0p9_1p1=bool(all(.9 <= x <= 1.1 for x in ratios))))
    lens_consistent = {case: any(s['all_within_0p9_1p1'] for s in summaries if s['case'] == case)
                       for case in ('rpg1_constant_a', 'rpg1_linked_a')}
    t3 = dict(lenses=lens_rows, summaries=summaries)
    en = np.array([[r['energy'][k] for k in ('dE_dlnr_10_100_rM_over_third_Mvf2', 'E_F_inside_10rM_over_Mvf2',
                                              'newtonian_gradient_energy_over_Mvf2', 'potential_rise_10_100_rM_over_vf2_ln10')]
                   for r in main_rows])
    ru = np.array([[r['radiation']['R_u_over_last_radius'], r['radiation']['R_u_over_r_M']] for r in main_rows])
    t4 = dict(median_and_range={k: [float(np.median(en[:, i])), float(en[:, i].min()), float(en[:, i].max())]
                                for i, k in enumerate(('dE_dlnr_10_100_rM_over_third_Mvf2', 'E_F_inside_10rM_over_Mvf2',
                                                       'newtonian_gradient_energy_over_Mvf2',
                                                       'potential_rise_10_100_rM_over_vf2_ln10'))},
              radiation_radius=dict(median_R_u_over_last_radius=float(np.median(ru[:, 0])),
                                    galaxies_with_R_u_inside_last_radius=int(np.sum(ru[:, 0] < 1)),
                                    median_R_u_over_r_M=float(np.median(ru[:, 1])), u_CMB_J_m3=U_CMB,
                                    luminosity='L_bol = L[3.6] (declared upper bound; R_u scales as sqrt(L_bol/L[3.6]))'))
    c2a = C_MS**2*L.ALPHA0/MPC_M
    a_si = A/Q.SI_TO_CODE
    t5 = dict(c2_alpha_m_s2=c2a, xi_group_alpha=a_si/c2a, xi_supernova_alpha=a_si/(C_MS*SN_C_ALPHA*1e3/MPC_M),
              prediction='a*(z) = (1+z) a*(0) for a linear index',
              sparc_redshift_note='SPARC and the Milky Way lie at z < 0.03, so they cannot test the evolution',
              slacs_theta_linked_over_constant={
                  row['name']: {imf: (row[imf]['cases']['rpg1_linked_a']['theta_arcsec']
                                      / row[imf]['cases']['rpg1_constant_a']['theta_arcsec']
                                      if row[imf]['cases']['rpg1_constant_a']['theta_arcsec'] > 0 else None)
                                for imf in ('Chabrier', 'Salpeter')} for row in lens_rows})
    result = dict(V1=v1, V1_passed=bool(v1_passed), convergence_selection=pick, T1=t1, T2=t2, T3=t3,
                  T3_declared_lensing_consistent=lens_consistent, T4=t4, T5=t5, a_star_code_units=A,
                  protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),
                  source_sha256={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                                 for p in ('aqual.py', 'baryons.py', 'lensing.py', 'rpg1.py')},
                  galaxies=main_rows, runtime_seconds=time.time() - t0)
    text = json.dumps(result, indent=1) + '\n'
    t2_summary = {k: {m: {x: y for x, y in v[m].items() if x not in ('speeds_kms', 'kz_Msun_pc2')}
                      for m in ('newtonian', 'algebraic', 'aqual')} for k, v in t2.items()}
    print(json.dumps(dict(V1_passed=result['V1_passed'], T1=t1['scores'], T2=t2_summary, lens=lens_consistent), indent=1))
    status = evidence_io.finish(args, 'radiation-polarized-gravity', text, HERE/'rpg1-results.json',
                                ignore={'/runtime_seconds'})
    return status if v1_passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
