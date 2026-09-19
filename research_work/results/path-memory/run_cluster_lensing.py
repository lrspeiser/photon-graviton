"""Run CL-F1 without fitting lensing or modifying historical evidence.

python research_work/results/path-memory/run_cluster_lensing.py --demo-geometry
An explicit --geometry-json may instead contain lens_kpc, source_kpc,
lens_source_kpc and label. Without either flag export geometry-free quantities.
"""
import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import cumulative_trapezoid

import cl1 as C1
import cl2_sources as CS
import steady_field as SF
from cluster_lensing import BaryonShells, Geometry, WrittenClusterLens, solve_lens
from cluster_lensing_checks import checks, relative

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plain(x):
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, dict):
        return {k: plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [plain(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, (int, np.integer)):
        return int(x)
    return x


def build(name, extract, widths, amplitudes, n=1500, npj=4000):
    # Empty operator widths: reuse the frozen source construction without a
    # pressure fit or a redundant force-operator calculation.
    c = CS.build_cluster(name, extract['clusters'][name], [], n_grid=n,
                         frac_profile=CS.stellar_fraction_profile(extract))
    lens = WrittenClusterLens(BaryonShells(c['r'], c['M_b']), widths, amplitudes, npj)
    return c, lens


def coma_transfer(widths, amplitudes):
    kubo = json.loads(C1.KUBO.read_text(encoding='utf-8'))
    R = np.array([x['published_radius_h_inverse_Mpc'] for x in kubo['rows']])*1000/C1.H_COMA
    y = np.array([x['shear_t'] for x in kubo['rows']])
    e = np.array([x['plotted_sigma_t'] for x in kubo['rows']])
    output = []
    for ne0, ms in C1.COMA_BRACKET:
        model = C1.coma_model(ne0, ms)
        r, rho = model['r'], model['rho']
        mass = cumulative_trapezoid(4*np.pi*r*r*rho, r, initial=0) + 4*np.pi/3*r[0]**3*rho[0]
        lens = WrittenClusterLens(BaryonShells(r, mass), widths, amplitudes)
        p = lens.profile(R)
        variants = {'baryons': (p['baryon_projected_mass_Msun']/(np.pi*R*R)
                                - p['baryon_surface_density_Msun_kpc2'], p['baryon_surface_density_Msun_kpc2']),
                    'frozen_pressure_response': (p['delta_sigma_Msun_kpc2'], p['sigma_Msun_kpc2'])}
        fits = {}
        for key, (ds, sigma) in variants.items():
            amp = max(float(np.sum(ds*y/e**2)/np.sum((ds/e)**2)), 0.)
            pred = amp*ds
            reduced = pred/(1-amp*sigma)
            fits[key] = dict(shape_chi2=float(np.sum(((pred-y)/e)**2)),
                             fitted_inverse_sigma_crit=amp, fitted_weak_shear=pred,
                             delta_sigma_Msun_kpc2=ds,
                             reduced_minus_weak_max_sigma=float(np.max(np.abs(reduced-pred)/e)))
        output.append(dict(ne0_cm3=ne0, stellar_mass_Msun=ms, total_baryon_mass_Msun=mass[-1],
                           radius_kpc=R, observed_shear=y, plotted_errors=e, fits=fits))
    return dict(status='Exposed-data shape test only; one fitted geometry nuisance per bracket. No absolute lensing validation.',
                brackets=output)


def plot_results(out, clusters, coma):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 4, figsize=(13, 8), constrained_layout=True)
    for ax, (name, row) in zip(axes.flat, clusters.items()):
        p = {k: np.asarray(v, float) for k, v in row['profile'].items()}
        R = p['radius_kpc']
        baryon = p['baryon_projected_mass_Msun']/(np.pi*R*R) - p['baryon_surface_density_Msun_kpc2']
        ax.plot(R, baryon, color='#718096', label='Ordinary matter')
        ax.plot(R, p['delta_sigma_Msun_kpc2'], color='#bf5a28', label='Frozen written field + matter')
        ax.axvline(row['source_outer_radius_kpc'], color='#bfc5ce', ls=':', lw=1)
        ax.set(xscale='log', title=name, xlim=(10, 15000))
        if np.all(baryon > 0) and np.all(p['delta_sigma_Msun_kpc2'] > 0):
            ax.set_yscale('log')
            ax.set_ylim(1e5, 2e9)
        else:
            ax.set_yscale('symlog', linthresh=1e6)
        ax.grid(alpha=.15)
    axes[0, 0].legend(fontsize=7)
    fig.supxlabel('Projected radius [kpc] · dotted line: finite source boundary')
    fig.supylabel('DeltaSigma [solar masses / kpc²] · geometry independent')
    fig.suptitle('Conditional cluster-lensing predictions — pressure fit frozen, no dark-matter source', fontsize=13)
    fig.savefig(out/'cluster-profiles.png', dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    first = coma['brackets'][0]
    ax.errorbar(np.asarray(first['radius_kpc'])/1000, first['observed_shear'], yerr=first['plotted_errors'],
                fmt='o', color='#222f40', label='Existing Coma figure-reconstructed bins')
    for j, row in enumerate(coma['brackets']):
        for key, color in [('baryons', '#718096'), ('frozen_pressure_response', '#bf5a28')]:
            fit = row['fits'][key]
            ax.plot(np.asarray(row['radius_kpc'])/1000, fit['fitted_weak_shear'], color=color, ls='-' if j == 0 else '--',
                    label=f'{key.replace("_", " ")} · bracket {j+1} · chi²={fit["shape_chi2"]:.2f}')
    ax.set(xscale='log', xlabel='Projected radius [Mpc]', ylabel='Tangential shear',
           title='Coma transfer: shape test with one geometry nuisance per curve')
    ax.legend(fontsize=8)
    ax.grid(alpha=.15)
    fig.savefig(out/'coma-transfer.png', dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path)
    geo = ap.add_mutually_exclusive_group()
    geo.add_argument('--demo-geometry', action='store_true')
    geo.add_argument('--geometry-json', type=Path)
    ap.add_argument('--cluster', action='append', help='Select an X-COP cluster; repeat for several.')
    args = ap.parse_args()
    geometry = None
    if args.demo_geometry:
        geometry = Geometry(200000., 1000000., 800000., 'Fictional static Euclidean demonstration; NOT the actual cluster distances')
    elif args.geometry_json:
        geometry = Geometry(**json.loads(args.geometry_json.read_text(encoding='utf-8')))
    out = args.output_dir or ROOT/'research_work/generated'/datetime.now(timezone.utc).strftime('cluster-lensing-%Y%m%dT%H%M%SZ')
    if out.exists() and any(out.iterdir()):
        ap.error('Output directory must be new or empty; previous evidence is never overwritten.')
    out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    pinned = [HERE/'cl2-results.json', CS.XCOP, C1.KUBO, HERE/'cl1-results.json',
              HERE/'protocol-cluster-lensing-forward.md', HERE/'cluster_lensing.py',
              HERE/'cluster_lensing_checks.py', Path(__file__)]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in pinned}
    archive = json.loads((HERE/'cl2-results.json').read_text(encoding='utf-8'))
    spectrum = archive['E1_spectrum']['clusters_alone']['spectrum']
    widths = [s['w_kpc'] for s in spectrum]
    amplitudes = [s['Lambda'] for s in spectrum]
    extract = CS.load_xcop()
    names = args.cluster or list(extract['clusters'])
    if any(n not in extract['clusters'] for n in names):
        ap.error('Unknown cluster; choices: ' + ', '.join(extract['clusters']))
    validation = checks()
    print('Analytic and independent optics gates:', validation['numerical_verification_passed'], flush=True)
    if not validation['numerical_verification_passed']:
        (out/'failed-gates.json').write_text(json.dumps(plain(validation), indent=2, allow_nan=False), encoding='utf-8', newline='\n')
        return 1
    clusters, refinement = {}, {}
    for name in names:
        c, lens = build(name, extract, widths, amplitudes)
        R = np.geomspace(1., max(5*c['r_out'], 8*max(widths)), 360)
        p = lens.profile(R) if geometry is None else lens.observables(R, geometry)
        row = dict(source_outer_radius_kpc=c['r_out'], total_baryon_mass_Msun=float(c['M_b'][-1]),
                   stellar_source=c['star_source'], below_10_kpc='Extrapolated and not resolution-certified',
                   profile=p, lens_solution=solve_lens(lens, geometry, R, sampled=p) if geometry else None)
        clusters[name] = row
        with (out/f'{name}.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerow(p.keys())
            writer.writerows(zip(*(plain(v) for v in p.values())))
        print(f'{name}: exported {len(R)} radii', flush=True)

    # The declared real-source gate is always run even for a selected cluster.
    for name in ('A1795', 'A2319'):
        c, lo = build(name, extract, widths, amplitudes)
        _, hi = build(name, extract, widths, amplitudes, n=3000, npj=8000)
        R = np.geomspace(10., max(5*c['r_out'], 8*max(widths)), 150)
        a, b = lo.profile(R), hi.profile(R)
        errors = {key: relative(a[key], b[key]) for key in ('alpha_rad', 'delta_sigma_Msun_kpc2')}
        refinement[name] = dict(errors=errors, passed=all(e < 3e-3 for e in errors.values()))
        print(f'{name} refinement: {refinement[name]}', flush=True)
    numerical = validation['numerical_verification_passed'] and all(x['passed'] for x in refinement.values())
    coma = coma_transfer(widths, amplitudes)
    unchanged = all(sha(ROOT/p) == digest for p, digest in hashes.items())
    data = dict(experiment='CL-F1', base_commit='37dab0641ddd0af99b762d4a23cf24e8efba0e79',
                input_and_code_sha256=hashes, inputs_preserved=unchanged, frozen_pressure_spectrum=spectrum,
                fit_parameters_added_for_cluster_lensing=0,
                geometry=asdict(geometry) if geometry else None,
                coupling='L1 is conditional and un-derived; ordinary matter plus the existing written field.',
                numerical_verification=validation, real_source_refinement=refinement,
                numerical_verification_passed=bool(numerical and unchanged),
                observational_status='Not validated: no matched X-COP shear/source-geometry dataset; Coma is shape-only. CL-2 cross-scale failures remain.',
                clusters=clusters, coma_transfer=coma, elapsed_seconds=time.monotonic()-start)
    (out/'results.json').write_text(json.dumps(plain(data), indent=2, allow_nan=False), encoding='utf-8', newline='\n')
    plot_results(out, clusters, coma)
    lines = ['# Cluster lensing: executable, conditional, and still testable', '',
             f'Base `37dab06`; frozen pressure response; {len(clusters)} X-COP sources; no dark-matter source.', '',
             f'Numerical gates: **{"PASS" if numerical and unchanged else "FAIL"}**. Archive/input integrity: **{unchanged}**.', '',
             'The solver calculates physical deflection and DeltaSigma, then (when distances are supplied) convergence, shear, reduced shear, signed magnification, critical curves and point-source image positions. CSV files contain the full profiles.', '',
             'The five footprint amplitudes come only from the archived pressure fit (8.55 chi² per pressure point). No lensing fit changes them. Negative effective field density is retained; it is a property of the potential, not added matter.', '',
             f'Geometry: {geometry.label if geometry else "none; only geometry-free quantities exported"}.', '']
    if geometry:
        lines += [f'D_l={geometry.lens_kpc/1000:g} Mpc, D_s={geometry.source_kpc/1000:g} Mpc, D_ls={geometry.lens_source_kpc/1000:g} Mpc. Source offset: 5 arcsec.', '',
                  '| Cluster | Tangential critical radii (arcsec) | Images found |', '|---|---:|---:|']
        for name, row in clusters.items():
            sol = row['lens_solution']
            rings = ', '.join(f'{r["theta_arcsec"]:.3f}' for r in sol['critical_curves']['lambda_t']) or 'none in scan'
            lines.append(f'| {name} | {rings} | {len(sol["images"])} |')
        lines += ['', 'These angular results are conditional on the stated geometry. The search finds sign-bracketed roots; tangent roots and images outside the recorded radial range are not certified.', '']
    lines += ['![Geometry-free cluster profiles](cluster-profiles.png)', '',
              '## Coma: frozen-response transfer', '',
              '| Ordinary-matter bracket | Baryons-only shape chi² | Written-field shape chi² |', '|---|---:|---:|']
    for j, row in enumerate(coma['brackets']):
        f = row['fits']
        lines.append(f'| {j+1} | {f["baryons"]["shape_chi2"]:.3f} | {f["frozen_pressure_response"]["shape_chi2"]:.3f} |')
    lines += ['', 'Six previously exposed, figure-reconstructed bins. Each curve fits one nonnegative inverse-critical-density nuisance; this tests shape only and does not establish absolute bending. No untouched observational holdout is available.', '',
              '![Coma shape transfer](coma-transfer.png)', '', '## What remains unresolved', '',
              'The light rule has not been derived from a photon interaction. Photon supply, field-energy accounting and formation at these strengths remain open. CL-2 still fails the joint galaxy/lens/cluster test. X-COP optical profiles are predictions with finite source boundaries, spherical symmetry and no matched shear catalog. Below 10 kpc the input source is extrapolated. A cluster merger needs a nonspherical, evolving source calculation.', '',
              'The existing finite footprints lose their additional deflection beyond the ordinary source plus their widths. More normalization cannot repair that outer-shape limitation. A creative successor must change how the photon-written field propagates or how light samples it, and must earn its energy budget rather than introduce an independent invisible mass profile.', '',
              '## Reproduce', '', '```sh',
              'python research_work/results/path-memory/cluster_lensing_checks.py',
              'python research_work/results/path-memory/run_cluster_lensing.py --demo-geometry',
              '```', '', 'For adopted geometry, replace `--demo-geometry` with `--geometry-json path.json`; use `--cluster A1795` for a single cluster. Without a geometry flag, physical bend and surface-density diagnostics are exported without assigning angular predictions.']
    (out/'report.md').write_text('\n'.join(lines)+'\n', encoding='utf-8', newline='\n')
    print(f'Report: {out / "report.md"}', flush=True)
    return 0 if numerical and unchanged else 1


if __name__ == '__main__':
    raise SystemExit(main())
