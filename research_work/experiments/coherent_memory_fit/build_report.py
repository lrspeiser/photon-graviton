"""Build the human-readable CMF reports from preserved evidence."""
import argparse,json
from pathlib import Path
import laws as L
D=L.D;HERE=L.HERE
def main(update_status=False):
    e1=HERE/"evidence-v1";e2=HERE/"evidence-gradient-v1"
    selection=D.read(e2/"selection-before-test.json");selected=selection["selected"]
    pred=D.read(e2/"predictions.json");refs=D.read(e1/"references.json")
    paired=D.read(e2/"paired-bootstrap.json");light=D.read(e2/"lensing.json")
    audit=D.read(HERE/"audit.json");chart=D.read(HERE/"observed-comparison-data.json")
    g=chart["galaxy"]
    rows=[]
    for name,scores in [("Ordinary matter",refs["scores"]["ordinary"]),("Matched simple MOND",refs["scores"]["MOND"]),
                         ("Previous IH-1 shared law",refs["scores"]["IH1"]),("CMF shared law",pred["joint"]),
                         ("CMF galaxy-only law",pred["galaxy"])]:
        vs=" / ".join(f"{scores[s]['galaxy_rmse']:.2f}" for s in ("train","validation","test"))
        ps=" / ".join(f"{scores[s]['cluster_chi2_per_point']:.2f}" for s in ("train","validation","test"))
        rows.append(f"| {name} | {vs} | {ps} |")
    table="\n".join(rows)
    primary=[r for r in light["models"] if r["primary"] and r["model"] in ("joint","galaxy")]
    lens_table="\n".join(f"| {r['model']} | {r['ne0']:.4f} | {r['beta_unbounded']:.3f} | {r['shape_only']['chi2']:.3f} | {r['bounded']['chi2']:.3f} |" for r in primary)
    p=paired["joint"];q=paired["galaxy"]
    params="\n".join(f"| {n} | {v:.12g} |" for n,v in zip(selected["joint"]["basis"],selected["joint"]["theta"]))
    report=f"""# CMF-1 / CMF-2: the shared solution is still unresolved

We executed two declared campaigns, **132 fitting configurations and 396 optimization starts**, without dark matter, expansion or changing adopted distances. These are 11 structural basis/memory choices crossed with regularization and objective weights, not 396 independent theories. All optimizers reported termination success; none of that establishes physical validity.

The selected shared law improves the validation galaxy score but performs worse than the matched MOND benchmark on the test-labelled galaxies. The best galaxy-only memory law has a small point-score advantage over MOND, with an uncertainty interval spanning no improvement and a severe cluster-pressure failure. Adding ordinary-source profile slopes did not replace either selected candidate.

[Requested observation comparison](observed-comparison.png) | [PDF](observed-comparison.pdf) | [SVG](observed-comparison.svg) | [Exact plotted data](observed-comparison-data.json) | [Independent audit](audit.json)

## The requested chart

![One shared candidate versus observations and Newtonian baselines](observed-comparison.png)

NGC 3198 was chosen by name before inspecting its fit as a familiar extended-disk example, not because it fit well. It belongs to the historically exposed training partition. At its outermost measured radius **{g['radius_kpc'][-1]:.2f} kpc**, observed circular rotation is **{g['observed_kms'][-1]:.1f} +/- {g['error_kms'][-1]:.1f} km/s**, the shared candidate predicts **{g['candidate_kms'][-1]:.1f} km/s**, and ordinary-matter Newtonian gravity predicts **{g['newtonian_kms'][-1]:.1f} km/s**. The measured rotation field is a proxy for the circular orbit a star would follow; these data, often obtained from gas, do not track one individual star. Pressure/asymmetric-drift corrections and noncircular stellar motions are not solved here.

Coma shows tangential shear, the distortion of background images, rather than a directly measured trajectory of one photon. The blue/orange bands span the two fixed ordinary-source brackets, not statistical confidence bands. Curves use fixed source efficiency beta=1, the maximum in the adopted static thin-lens geometry; finite source distances lower this efficiency. The observed points and errors are reconstructed from Kubo et al.'s published figure. No lens normalization was fitted for the chart.

For the requested Newtonian light comparison, a ballistic test ray moving at c has half the standard weak-field deflection. Orange shows that explicitly defined approximation. Dotted gray also shows standard light bending from ordinary matter alone, so the stronger conventional baseline is visible. Neither includes dark matter.

## Matched observation scores

| Model | Galaxy train / validation / test RMSE, km/s | Cluster pressure train / validation / test chi2 per point |
|---|---|---|
{table}

The recalculated simple-MOND scale is **{refs['MOND_a0_SI']:.12g} m/s^2**, fitted only to the same training objects using the same equal-object squared-velocity objective. This replaces an unmatched historical comparison using a0=8.563335193921255e-11 and a 3,150-row reduction. The new baseline both uses all 3,152 rows and refits its coefficient; the score change must not be attributed only to adding two rows. MOND is credited as a reference, not presented as this project's mechanism.

All partitions were historically exposed. The second campaign explicitly used the first campaign's findings to motivate its new feature. Each stage saved its selection before its own test/Coma evaluation, but that ordering does not make reused observations blind. Cluster scores condition on one fitted nonnegative boundary pressure per object, including test-labelled clusters. Only quoted diagonal pressure errors are available.

The paired object bootstrap gives shared-minus-MOND test RMSE **{p['delta_rmse_kms']:.3f} km/s**, with conditional 95% interval **[{p['percentile95'][0]:.3f}, {p['percentile95'][1]:.3f}]**. The galaxy-only difference is **{q['delta_rmse_kms']:.3f} km/s**, interval **[{q['percentile95'][0]:.3f}, {q['percentile95'][1]:.3f}]**. These 10,000-resample intervals are conditional on the exposed sample and fixed selected models; they do not correct for selection, missing baryonic uncertainties or repeated experimentation.

The shared choice has {selected['joint']['n_parameters']} coefficients; the galaxy-only choice has {selected['galaxy']['n_parameters']}, versus one fitted acceleration scale for this simple-MOND reference. Source and input nuisance assumptions also matter. A small raw-score advantage is not a simplicity or evidence advantage. No active dark-matter comparison was run, so superiority to dark-matter models is not established.

## The selected equation

The primary shared selection is **{selected['joint']['id']}**. Its selected memory length is zero: the best shared candidate from this family does not require outward profile memory. Define ordinary acceleration g_b in m/s^2, physical radius r in kpc, source-mass proxy M in solar masses:

    x = tanh[-ln(max(g_b,1e-30)/1e-10)/4]
    y = tanh[ln(r/10)/4]
    z = tanh[ln(M/1e11)/4]
    P = sum theta_j F_j
    g = g_b exp[8 tanh(P/8) / (1+(g_b/1e-7)^2)]
    v_circular = sqrt(r_physical * g)

| Basis F_j | theta_j |
|---|---:|
{params}

The hypothetical radial-memory versions solve ell dH/dln(r)=x-H with their declared inner boundary; h=H-x enters the response. This is radial profile memory, not a time step, vector graviton alignment or a microscopic hitchhiking derivation. CMF-2 added only derivatives of the ordinary-source acceleration, never observed velocity or pressure. Coefficients of either sign allow enhancement and suppression while keeping total inward acceleration positive. The exponential cap and high-acceleration release are declared modeling choices.

## Removing the arbitrary lens normalization

The same total force enters the stipulated weak-field integral:

    alpha(b) = (4/c^2) integral_0^infinity g(sqrt(b^2+u^2)) b/sqrt(b^2+u^2) du
    gamma_t = (D_l beta/2) [alpha/b - d alpha/db]
    beta = 1-D_l/D_s, 0 <= beta <= 1

The fixed adopted Coma distance is 100 Mpc. We preserve the previous h=0.7 radial conversion; their cross-catalog angular consistency has not been recalibrated. Source photometric redshifts are not converted into distances through expansion. The source model is truncated at 3 Mpc; the extra force has a declared 9 Mpc reach followed by an inverse-square tail. The small feature near the gas cutoff in the plot comes from this assumed truncation.

| Candidate | Gas central density, cm^-3 | Preferred unbounded beta | Free-shape chi2 | Physically bounded chi2 |
|---|---:|---:|---:|---:|
{lens_table}

For the shared candidate, the preferred amplitude needs beta greater than one. Moving sources farther away cannot supply that best-fit normalization in this geometry. However, the six shear bins have broad errors: the bounded high-source case is not by itself a decisive statistical exclusion. Missing covariance, source weights, bin edges and independent geometry still prevent an absolute lensing-validation claim. The galaxy-only candidate can reach an admissible beta, but fails the pressure transfer; it is not a joint solution.

## Verification and preserved history

Both source implementations were committed before their fits: CMF-1 protocol f9a5928, source 9726084; CMF-2 protocol 7121a21, source 115f3fa. A concurrent remote merge brought the older cluster branch onto main at 3a920b0. We preserved it in merge e5960c3 before the second implementation; source/input hashes verify the present calculations' dependencies.

All **33 preliminary controls** pass. The first and second campaigns pass **14 and 16** selected pressure/optical refinement checks respectively. Source-grid and inner-boundary sensitivities are reported separately from pass/fail checks; the second campaign also checks the galaxy-only memory law rather than inferring its robustness from the zero-memory joint winner.

The audit checks **{audit['evidence_files']} evidence hashes**, **{audit['pinned_source_input_checks']} inputs/sources against their recorded Git commits**, and **149 unchanged distances**. An independently written selected polynomial reproduces velocities to {audit['independent_velocity_error_kms']:.3g} km/s and pressure to {audit['independent_pressure_error']:.3g} in the archived units. Independent ODE memory agrees within {audit['independent_memory_error']:.3g}; adaptive optical quadrature agrees within {audit['independent_optics_relative_error']:.3g} relative. The adaptive integrator emitted {len(audit['adaptive_integration_warnings'])} roundoff/tolerance warnings near piecewise interpolation knots; these are saved in audit.json. The comparison passes the audit's 1e-4 threshold, not a claim that the integrator attained every tighter internal tolerance.

The first chart rendering had an overlapping annotation/footer; visual inspection led to spacing corrections. Fitting evidence and plotted values were unchanged. Run archives refuse overwrite. No failed candidate or negative shear bin was removed.

## What remains before there is a solution

The work narrows the mathematical requirements; it does not solve the requested universe. Shared predictions must improve galaxy transfer without losing cluster pressure and lens strength. A physical transport law must derive the ordinary-source field, its vector following/swirl behavior, and matter/light coupling. Source fuel, field energy, attachment/release energy, recoil and boundary flow must balance. The empirical positive force law does not supply that accounting.

Fresh matched cluster pressure, ordinary-matter and calibrated lensing measurements are needed, along with independent galaxy validation. New geometry cannot be tuned to repair this result. Increasing polynomial flexibility has already been tried here and did not yield a shared win.

## Attribution

[SPARC data](https://arxiv.org/abs/1606.09251), [X-COP pressure](https://arxiv.org/abs/1805.00042), [Kubo et al. Coma shear](https://arxiv.org/abs/0709.0506), and [Alabi et al. fixed Coma distance convention](https://academic.oup.com/mnras/article/496/3/3182/5859958) supply observations or stated distance inputs. The latter paper's expansion model is not adopted.

[MOND](https://adsabs.harvard.edu/pdf/1983ApJ...270..365M), [standard line-of-sight lensing](https://arxiv.org/abs/astro-ph/9912508), polynomial regression, saturation, relaxation, ridge penalties and bootstrap are established work. [Earlier symbolic acceleration-law searches](https://doi.org/10.1093/mnras/stad597) also exist. This project owns its implementation and declared candidate combinations, not these established tools, and asserts no historical novelty.
"""
    (HERE/"report.md").write_text(report,encoding="utf8",newline="\n")
    summary=dict(selected={k:{kk:f[kk] for kk in ("id","family","ell","theta","basis","n_parameters","weight","ridge")} for k,f in selected.items()},
        scores={**{k:{s:{kk:v[s][kk] for kk in ("galaxy_rmse","cluster_chi2_per_point")} for s in ("train","validation","test")} for k,v in pred.items()},
                "MOND":{s:{kk:refs["scores"]["MOND"][s][kk] for kk in ("galaxy_rmse","cluster_chi2_per_point")} for s in ("train","validation","test")}},
        paired_bootstrap=paired,fits=132,optimization_starts=396,structural_basis_memory_choices=11,
        complete_solution=False,beats_MOND=False,beats_dark_matter_established=False,
        dark_matter=False,expansion=False,distances_adjusted=False,energy_funding_derived=False)
    L.save(HERE/"summary.json",summary)
    readme="""# Shared memory fits and observed comparison

Start with [the chart](observed-comparison.png), [report](report.md), [exact summary](summary.json), and [audit](audit.json).
The selected shared formula remains unsuccessful as a full galaxy/cluster solution.

## Reproduction
Python 3.13.5; NumPy 2.2.6; SciPy 1.16.1; Matplotlib 3.10.5.
Set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. From the repository root:

    python -B research_work/experiments/coherent_memory_fit/audit.py
    python -B research_work/experiments/coherent_memory_fit/make_comparison.py
    python -B research_work/experiments/coherent_memory_fit/build_report.py

These verify evidence or rebuild derived artifacts, not the fitted archives.
For a fresh campaign, use an isolated checkout at the source commits:
- CMF-1: protocol f9a5928, source 9726084; run.py creates evidence-v1.
- CMF-2: protocol 7121a21, source 115f3fa; gradient_run.py uses the preserved first archive and creates evidence-gradient-v1.
Both drivers refuse to overwrite their archive directories. Keep Git history for provenance checks.

The 132 configurations comprise 84 plus 48 fits, with three starts each.
All partitions are exposed; no fresh validation, full energy conservation, microscopic attachment or historical originality is claimed.
No dark-matter model or expanding universe is active. No catalog distance is adjusted.
"""
    (HERE/"README.md").write_text(readme,encoding="utf8",newline="\n")
    if update_status:
        block=("**CMF-1 / CMF-2 executed, 19 September 2026.** [Observed lensing/rotation chart](research_work/experiments/coherent_memory_fit/observed-comparison.png), "
               "[report and equations](research_work/experiments/coherent_memory_fit/report.md). "
               "132 fit configurations, 396 starts, 11 basis/memory choices; all exposed observations, fixed distances, no dark matter or expansion. "
               "The matched simple-MOND reference is now refitted on identical 3,152 rows and object weights: validation/test RMSE 26.14/16.09 km/s. "
               "The selected shared candidate gives 23.99/21.72 km/s and cluster-pressure validation/test chi2 per point 9.83/5.06. "
               "Galaxy-only memory reaches 15.65 km/s but its paired improvement interval includes zero and cluster-pressure test score is 195.64. "
               "The shared Coma model prefers unphysical source efficiency beta=2.94/1.72; constraining beta<=1 yields chi2 12.77/7.59 for the two ordinary-source brackets. "
               "These are conditional six-bin calculations, not an absolute lensing validation or decisive exclusion. "
               "At NGC 3198's 44.08 kpc outer data point: observed 149+/-3, candidate 119.3, Newtonian ordinary matter 64.7 km/s. "
               "33 controls, 30 selected refinements, 21 evidence hashes and 27 pinned source/input checks pass; adaptive quadrature warnings are explicitly archived. "
               "Concurrent cluster merge 3a920b0 was preserved in e5960c3. Full shared solution, vector field generation, microscopic light coupling and energy funding remain open.\n\n")
        for name,prefix in (("CURRENT-STATUS.md",block),("CHANGELOG.md","## 2026-09-19: matched source-memory search and observed charts\n\n"+block)):
            path=L.ROOT/name;raw=path.read_bytes()
            if b"**CMF-1 / CMF-2 executed" in raw:continue
            pos=raw.index(b"\n")+1
            path.write_bytes(raw[:pos]+b"\r\n"+prefix.replace("\n","\r\n").encode("utf8")+raw[pos:])
    print("Report, summary and reproduction instructions generated.")
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--update-status",action="store_true")
    main(parser.parse_args().update_status)
