# SM-1: a screened square-root field, with reciprocal dynamics

**21 September 2026. Decision: retain as a concrete effective candidate, not a superior observational theory.**

SM-1 adds an explicit high-acceleration completion to the project's square-root law, derives its total-potential field equation and energy, and implements reciprocal matter/field evolution. The declared numerical checks pass. The exposed SPARC algebraic screen is mixed, not an improvement overall. Actual synthetic disk solutions depart from the spherical shortcut by up to 30%, so that screen is not the final disk prediction.

## Provenance and reproduction

- Inspected baseline: `28f154242aa6e6ae010a97af8124d01cc63fed3b`.
- Declaration: `fbc02535abfe707c2b347e7ddc3f50d38d63d244`, [protocol.md](protocol.md).
- Executed implementation: `c43a5c647603e26e4fe577a27e9bde75eb4d8478`, [model.py](model.py), [run.py](run.py).
- Completed workflow: [35647197394](https://github.com/lrspeiser/photon-graviton/actions/runs/35647197394), conclusion **success**. This is the new SM-1 workflow, not a rerun of all historical repository suites.
- Durable compact evidence: [results-summary.json](results-summary.json). Full 149-galaxy JSON: Actions artifact `10660875299`, also downloaded into the originating conversation as `sm1-results-c43a5c6.zip`. The summary is not a substitute for that complete per-galaxy archive; Actions retention is finite.
- Original artifact SHA-256: `a3bf8dd3c415957554d0ac2a5569b2cb197382519b5c9357a12ae0ee60240ae2`.
- Extracted full JSON SHA-256: `c44b013153d80ba165ac6fe4708e9fd1cf8a43a282c880c8c112974049bcfff0`.

From the repository root, using NumPy 2.3.5 and SciPy 1.17.0:

```sh
python research_work/results/screened-memory/run.py --sparc --disks --output-dir sm1-replay-001
```

Use a fresh output directory. The runner refuses to overwrite prior evidence. The complete replay ran on CPython 3.13.15; local algebra and dynamics also ran on 3.13.5. Source and input hashes are recorded in the results.

## What was already known, and what is added

[PM-1](../path-memory/report.md) had already fitted a local square-root acceleration law. [PM-2A](../path-memory/report-pm2a.md) had already constructed both total-potential and independently sourced auxiliary-field completions. Neither is a discovery of SM-1. The microscopic [Phase Junction program](../../../phase_junction_network/README.md) remains separate and retains its own open constraints, matter-loop and electromagnetic issues.

The additions here are a particular **decaying high-acceleration excess**, its invertible constitutive relation, its consistent gradient functional and reciprocal time evolution, a frozen-data screen, and actual nonspherical field solutions. These are an effective construction. The chosen interpolation is not uniquely derived from a finite junction algebra, and the new scalar is not identified with the microscopic companion by these tests.

## 1. Repair the screening implementation, not the entire mechanism

For the review's isolated, spherical, independently sourced companion, write

\[
 g_N=\mu_\chi(g_\chi/a_0)g_\chi.
\]

Local ellipticity requires

\[
 \frac{d g_N}{d g_\chi}=\mu_\chi+y\mu_\chi'>0,
 \qquad y=g_\chi/a_0.
\]

Its inverse therefore has an increasing companion force. Choosing a higher power of the companion gradient makes the **fractional** correction shrink but does not make its **absolute** force decrease. This is a restriction on that particular separately sourced scalar with those boundary assumptions, not a no-go theorem for environmental, nonlocal or memory gravity, nor a proof that every increasing companion force is observationally excluded.

SM-1 instead screens the **total-potential completion** already allowed by PM-2A. It does not silently reuse the independent-companion equation with a different label.

## 2. One frozen equation

Define

\[
 x=\frac{g_N}{a_0},\qquad z=\frac{g}{a_0},\qquad
 h(x)=x+\frac{\sqrt{x}}{1+x^2},\qquad
 a_0=6.54\times10^{-11}\ {\rm m\,s^{-2}}.
\]

The spherical prediction is

\[
 \boxed{g=g_N+\frac{\sqrt{a_0g_N}}{1+(g_N/a_0)^2}.}
\]

The acceleration scale is the **rounded previous PM-1 fitted scale**. The crossover and exponent are fixed phenomenological choices, not new fitted parameters in this run and not microscopic predictions. No galaxy receives its own halo, transition or acceleration scale.

At low acceleration the excess approaches the original square root. At high acceleration,

\[
 g-g_N\sim\frac{a_0^{5/2}}{g_N^{3/2}},
\]

so the absolute excess vanishes. The familiar isolated outer scaling \(v_\infty^4=GMa_0\) follows from the chosen asymptote; it is built in, not an independent discovery.

The source-based, nonspherical completion is

\[
 \mu(z)=\frac{h^{-1}(z)}{z},\qquad
 \boxed{\nabla\cdot[\mu(|\nabla U|/a_0)\nabla U]=4\pi G\rho_b.}
\]

The measured force is \(-\nabla U\). One may define an excess potential \(U-\Phi_N\), but it does not then obey the review's independently sourced companion equation.

### Invertibility and gradient energy

For every \(x>0\),

\[
 h'(x)=1+\frac{1-3x^2}{2\sqrt{x}(1+x^2)^2}
 >1-\frac{75}{128}\left(\frac35\right)^{3/4}
 =0.6005478789\ldots>0.
\]

The bound follows by bounding the negative term with \(3x^{3/2}/[2(1+x^2)^2]\), whose maximum occurs at \(x^2=3/5\). Hence

\[
 \mu>0,\qquad \mu+z\mu'=1/h'(x)>0.
\]

A consistent gradient functional is

\[
 W(z)=\int_0^z\mu(v)v\,dv,\qquad
 E_{\rm grad}=\int\frac{a_0^2}{4\pi G}W(|\nabla U|/a_0)\,d^3x.
\]

The isolated low-acceleration point-mass asymptote has a logarithmic potential and a logarithmically divergent gradient energy at infinite radius. Energy statements therefore need a finite domain or a declared environmental boundary and its flux; no finite whole-universe energy has been established here.

These results establish convexity of the scalar gradient term at nonzero gradient. They do not establish boundedness of the complete gravitational energy including attractive matter coupling, or nonlinear stability of the full matter/field system. The principal symbol is degenerate at zero gradient.

### Two further exact consequences

Derived after the primary run, not retroactively declared as preregistered gates:

**The isolated spherical extra acceleration has a ceiling**

\[
 \Delta g_{\max}=\frac{3^{3/4}}4 a_0
 =3.7269940\times10^{-11}\ {\rm m\,s^{-2}},
 \qquad g_N/a_0=1/\sqrt3.
\]

**Every circular test orbit around an isolated point mass is radially stable in the prescribed static potential.** Since \(x=GM/(a_0r^2)\),

\[
 \kappa^2=\frac{3g}{r}+\frac{dg}{dr}
 =\frac{a_0}{r}\left[x+
 \frac{2\sqrt{x}(1+3x^2)}{(1+x^2)^2}\right]>0.
\]

The identity was also checked numerically over 24 decades in \(x\), with relative discrepancy \(8.9\times10^{-16}\). This proves neither circularization of arbitrary trajectories nor stability once the source and field co-evolve.

## 3. Reciprocal evolution and the energy budget

The declared preferred-frame effective extension is

\[
 \frac{U_{tt}+\Gamma U_t}{c_*^2}
 -\nabla\cdot[\mu(|\nabla U|/a_0)\nabla U]
 =-4\pi G\rho_b,
 \qquad m\ddot{\mathbf q}=-m\nabla U.
\]

For \(\Gamma=0\), its energy contains

\[
 H=\sum_i\frac{p_i^2}{2m_i}+\sum_i m_iU(\mathbf q_i)
 +\int\frac{1}{4\pi G}
 \left[\frac{U_t^2}{2c_*^2}+a_0^2W\right]d^3x.
\]

The discrete code uses extended particles and differentiates the same normalized source interaction for both field sourcing and particle reaction. There is no one-way force or inserted initial field excitation. For damping, the code additionally records

\[
 \dot Q=\int\frac{\Gamma U_t^2}{4\pi Gc_*^2}\,d^3x\ge0.
\]

This closes an **energy ledger**, not a derived microscopic bath, bath momentum budget, photon-conversion mechanism or long-lived cluster deposit.

At a nonzero constant background gradient, the scalar's squared principal speeds are proportional to \(\mu\) and \(\mu+z\mu'\), both positive. The chosen toy value \(c_*=0.7\) lies below the reference unit-speed bound; it is not a derivation of a physical shared light/gravity cone. Scalar radiation and preferred-frame observational constraints remain untested.

### Measured implementation checks

| Check | Measured maximum error |
|---|---:|
| Constitutive inverse, 28 decades of input | \(7.8\times10^{-15}\) |
| Independent finite-difference gradient-energy derivative | \(6.0\times10^{-11}\) |
| Full discrete field Hamiltonian derivative | \(7.6\times10^{-11}\) |
| Full particle interaction derivative | \(6.8\times10^{-10}\) |

The declared one-dimensional empty-field run uses two initially resting extended particles, a fixed finite box and duration 2. Conservative normalized energy errors at timesteps 0.004, 0.002 and 0.001 are **3.70e-6, 9.40e-7 and 2.36e-7**. The damped energy-plus-heat errors are **2.04e-6, 5.10e-7 and 1.28e-7**. Heat is nondecreasing.

**There is no resolved particle motion in that primary duration-2 run.** The field develops, but both positions remain at +/-1. The positive field energy is offset by negative matter-field interaction energy; no positive energy is created without an offset.

A clearly labeled post-primary duration-8 exploration shows the positions moving to about +/-0.35515 in the conservative 161-node box. Halving the timestep changes the positive position by 2.4e-6; doubling the grid gives +/-0.34864. The longer run is not spatially certified, does not isolate wall effects with single-source controls, and is not a stable-orbit or physical-front test. It is retained as evidence of actual coupled motion with an energy budget, not as an astrophysical result.

## 4. Exposed SPARC screen: mixed, not an overall improvement

The runner uses the frozen 149 galaxy names and 89/29/31 split, all **3,152 positive-radius raw rows**, and the archived mass-to-light ratios 0.5 for disks and 0.7 for bulges. Signed gas contributions are retained; two negative total baryon squared-speed rows are clipped to zero consistently across models and counted. No observed velocity enters any predictor.

The table reports the **arithmetic mean of each galaxy's velocity RMSE, in km/s**, not a percentage error, likelihood, historical headline metric or uncertainty-normalized chi-square. All models use the same rows and weights within this comparison.

| Model | Train: 89 | Validation: 29 | Archived test: 31 |
|---|---:|---:|---:|
| Baryons only | 45.93 | 48.64 | 41.68 |
| Unscreened square root, same rounded \(a_0\) | 16.17 | 21.85 | 14.78 |
| **SM-1 screened law** | **17.35** | **20.32** | **15.33** |
| Archived simple-MOND comparison | 15.78 | 19.79 | 13.50 |

Screening changes the mean error by **+7.30%, -7.02% and +3.68%**, respectively. Across all 149 galaxies it changes **16.98 to 17.50 km/s**. The simple-MOND comparison uses its previously selected scale, 8.563e-11 m/s^2; neither scale is fitted again here.

This is exposed-data reuse, not a new blind test. Mean, median and pooled scores are all available in the full JSON; no metric was selected to hide the worsened means. No paired significance claim is made. In particular, **SM-1 is not promoted as a better galaxy fit**.

## 5. Actual disk solutions: the spherical shortcut is not the same theory

The same constitutive law was solved on thick exponential sources with radial scale 1, vertical scale 0.15, and compactness \(GM/(a_0R_d^2)=0.1\) or 10. The inner unresolved mass uses a spherical integral rather than the inherited thin-cylinder approximation. Newtonian, unscreened and screened forces share the same source on each grid.

All **14 finite-volume solves converge**. At the declared seven radii:

| Source compactness | Largest screened grid change | Outer-boundary change | Largest screened PDE departure from algebraic radial mapping |
|---|---:|---:|---:|
| 0.1 | 0.4504% | 0.00255% | -29.90% |
| 10 | 0.4943% | 0.00207% | -11.94% |

The grid tests double radial and angular resolution; the boundary control doubles outer radius without degrading logarithmic radial resolution. Both sources pass the fixed 2% refinement and 0.5% boundary gates.

The largest shortcut departures occur at \(r=0.3R_d\). For the low-compactness disk they are about -14.08% at \(R_d\), -5.78% at \(2R_d\), and +1.30% at \(8R_d\). These are synthetic model predictions, not observed galaxy discrepancies.

**Consequence:** the SPARC table is only an algebraic screening result. It cannot be relabeled as the fit of the actual nonspherical field equation. A 149-source field comparison must solve each frozen density and recompute its Newtonian baseline; substituting force-equivalent mass or merely applying the spherical formula would repeat PM-1's earlier interpretation error.

## 6. What remains genuinely unresolved

At Earth's orbital radius, the isolated point-Sun excess changes from **6.23e-7 to 7.57e-23 m/s^2**. This demonstrates the intended monopole suppression, **not a Solar-System pass**. The Galactic external field can generate a quadrupole even when the local monopole correction is negligible.

The updated Cassini/DE440 analysis by Park et al., revised 28 July 2026, reports \(Q_2=(1.6\pm1.8)\times10^{-27}\,s^{-2}\). It is a relevant next constraint, not a number this experiment has matched. Its published tensions apply to the formulations and interpolation families analyzed there and are not a substitute for evaluating this exact law and declared environment.

Still open: that external-field calculation; a full same-density SPARC field fit; relativistic matter/light coupling and lensing; scalar/binary radiation; a tensor-sector integration; a microscopic derivation of the interpolation and acceleration scale; and memory capable of predicting merger offsets and deposits without separately assigned halos. No cluster, cosmology, photon-conversion, fundamental-gravity or universal-stability success is claimed.

## Next decisive work

Keep SM-1 frozen as a benchmark. The useful next comparison is the **actual same-density disk field plus the Solar external-field quadrupole under this same law**, not a search over interpolation coefficients until one score improves. In parallel, only a sourced, reciprocal memory mechanism with a physical energy/transport interpretation can turn the broader companion hypothesis into a distinct theory rather than another static acceleration prescription.

## Prior art and external constraints

- Bekenstein and Milgrom (1984), *Does the missing mass problem signal the breakdown of Newtonian gravity?* Aquadratic variational gravity predates this project. The present total-potential construction belongs to that broad framework; it is not claimed as a new foundational theory.
- Desmond, Hees and Famaey (2024), [arXiv:2401.04796](https://arxiv.org/abs/2401.04796), joint radial-acceleration/Solar-quadrupole constraints.
- Park et al. (2026), [arXiv:2602.17884v2](https://arxiv.org/abs/2602.17884v2), updated Cassini radio-tracking constraint. External references provide context, not evidence that SM-1 passes their tests.
