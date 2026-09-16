# PM-1: path-memory gravity against the data already in the repository

Declared before execution, 15 September 2026. Baseline: `main` at a09827a.

**Why it is run.** The owner's proposal after reading the walking-droplet work: treat matter and the gravitational field as one system with memory, so the field at time t depends on where matter has been. They set out six candidate families (A–F) with toy-model calculations, and asked for them to be tested against the data already on this machine rather than against new analogies. The repository holds the two datasets that can decide them: SPARC's 149 galaxies and 3,150 radii with a frozen baryon model and a frozen train/validation/test split, and the Milky Way's 38 Eilers bins with two baryon fiducials. Both already carry archived comparison scores, so a new law can be judged the way this repository judges every other one — **universal against universal, at equal freedom**.

**The owner's own framing is kept.** These are mathematical hypotheses. A candidate that reproduces a flat curve but misses the baryonic mass–speed relation is not sufficient, and an equation that permits an orbit is not the same as an orbit that survives perturbation.

## What is fixed

- **SPARC.** The archived 149 galaxies through `capture-to-orbit/inputs.py`: the rotmod curves, the catalogue (R_d, L₃.₆, M_HI), and the enclosed baryonic mass M_b(<r) that `sparc_receivers` builds from the stellar disk (Υ = 0.5), the bulge (Υ = 0.7) and the helium-corrected gas (1.33 M_HI). No mass-to-light ratio is refitted; [[sparc-stellar-ML-resolved]] settled that Υ\* = 0.5 is a global nuisance.
- **The split.** The archived 89/29/31 train/validation/test division, used as it stands. These data are exposed; this is not a blind test.
- **The Milky Way.** The 38 Eilers bins and baseline I, as every other experiment here uses them.
- **The comparison scores**, from CURRENT-STATUS: baryons alone 52.57 / 58.22 / 47.77 (SPARC train/validation/test) and 52.57 / 62.34 (Milky Way I / II); simple MOND with a fitted a₀ 19.89 / 26.88 / 16.40 and 9.53 / 12.04; the exact-third reference 29.03 / 32.50 / 23.59 and 6.76 / 10.41.
- **Scoring.** Equal-galaxy RMSE in km/s, as the frozen benchmarks use, with one universal constant fitted on the training galaxies only and then applied unchanged.

## The candidates, and what each run decides

### A, causal memory of density — the control

Ψ_A is an exponential memory of ρ at retarded time. For a stationary source it returns αΦ_N, so the force law is unchanged and only G is rescaled. **Run:** verify numerically on one SPARC galaxy that the memory average of a static density reproduces Φ_N to 10⁻¹⁰, and that a rigidly rotating axisymmetric density gives the identical field. **Decides:** that delay alone explains nothing, which is the null every other candidate must beat.

### B, orbit-generated waves — the candidate with data traction

A single remembered wavelength gives an on-ring force A k J₀(kR)J₁(kR), whose envelope falls as 1/R but alternates in sign. A band with equal potential weight per logarithmic wavenumber gives, exactly,

    g_band(R) = (C/2R)·[J₀(k_min R)² − J₀(k_max R)²] → C/2R  when  k_min R ≪ 1 ≪ k_max R,

and, for one fixed ring, a logarithmic potential whose force outside the ring is C/r. Summing rings makes the coefficients add, so the force at r comes from the material inside r.

**Two sourcing rules are declared, because the composition rule is the physics at issue:**
- **B-linear:** each ring contributes C_i ∝ m_i, so g_mem(r) = κ·M_b(<r)/r. This must fail the mass–speed relation (v_f⁴ ∝ M², not M), and it is run to show that.
- **B-root:** the amplitude saturates as the owner's Ċ = sS − bC² suggests, so g_mem(r) = k·√(M_b(<r))/r with **one universal k for all 149 galaxies and the Milky Way**. Then v² = √(M_b) k outside the matter, which is the mass–speed scaling the data require.

**Runs.** For each rule, the total acceleration is g_N + g_mem with g_N from the same baryon model, v(r) = √(r·g). One constant is fitted on the 89 training galaxies, then applied unchanged to validation, test and the Milky Way. The banded force is used in its exact J₀² form, not only its 1/R limit, so the cut-offs k_min and k_max are explicit; the wide-band limit is reported as the special case.

**Also run, as the honest control:** the single-wavelength law, to show the alternating bands it produces on a real rotation curve.

**Declared readings:** the fitted k and its implied a-scale k² = G a\*, the three SPARC RMSEs and the two Milky Way ones, the residual trend with radius, and the baryonic mass–speed slope the fitted law predicts against the data's own.

### C, a receiver that changes its response — the instability

Equilibrium gives g = ½[g_N + √(g_N² + 4a\*g_N)], the familiar simple-ν form, which this repository has already scored at 22.16 / 29.42 / 16.90 in RPG-1. The new content is the owner's stability result. **Run:** reproduce their cubic Tλ³ + (1+q₀)λ² + Tλ + (1+3q₀) = 0 and its root 0.1533 ± 1.2858i at q₀ = 0.9, T = 1; map the growth rate over q₀ ∈ (0,1) and T ∈ [10⁻³, 10³], including the T → 0 and T → ∞ limits, where it must go marginal; and evaluate q₀ at the actual outer radii of the SPARC galaxies, since q₀ → 1 there is what makes the e-folding time physically fatal. **Decides:** whether any memory time is stable, and at which real radii the instability bites.

### E, energy deposited as a gravitating state — the budget

v_u² = GL/c³ from the catalogue's L₃.₆, per galaxy, against the observed v_f. **Decides:** the shortfall factor, with its own number rather than an illustrative one.

### D and F, declared but not run

D's guided field needs H/f = 2√(GM/a\*), tens of kiloparsecs, and SPARC measures no such thickness; assigning each galaxy the H that reproduces its speed would be circular. F needs a trajectory-history integrator and a separate light prescription. Both are named here so that their absence is a decision, not an oversight.

## Gates

- **G1, the shape.** A candidate passes only if one universal constant beats baryons alone on all three SPARC splits and on both Milky Way fiducials.
- **G2, the mass–speed relation.** The same fitted law must reproduce the baryonic mass–speed slope within the scatter of the data it is scored on; a candidate that passes G1 and fails G2 is reported as failing.
- **G3, universality.** No per-galaxy parameter. The only freedom is the single constant, fitted on the training split.
- **G4, stability.** For any candidate whose mechanism is a feedback with delay, the linearised orbit must not have a growing mode at the radii where the mechanism is needed.

Passing G1–G3 would make B a law worth deriving, not a law demonstrated: the spectrum, the confinement, the phase coherence and the energy accounting would all still be owed, and the protocol says so in advance.

## Scope: what this does not do

- It does not fit mass-to-light ratios, distances or inclinations, and it does not use SPARC's optional dark-matter fits.
- It does not derive the spectrum, the cut-offs, the phase coherence or the energy supply of B's modes. It tests the force law that a declared spectrum implies.
- It does not test lensing, the Solar System, free-fall universality or gravitational-wave timing, all of which any surviving candidate must face.
- The data are exposed. Nothing here is a blind test, and the split is reused from earlier work.
- The owner's own package (`path_memory_gravity_lab.zip`) could not be read from this session, so this is an independent implementation of their equations from the description, not a check of their code.

## Corrections, after the owner's review of the first run

Declared before the corrected run. The first run's fits are unchanged by corrections 1, 4, 5 and 6, which are reporting and proof; corrections 2 and 3 add calculations and replace one that was wrong.

**Correction 1, the acceleration scale is a re-expression of the fitted amplitude, not a second prediction.** With g_mem = β√M_b(<r)/r and the monopole comparison field g_mono = GM_b(<r)/r², identically g_mem = √(a\*·g_mono) with a\* = β²/G. Fitting β *is* fitting a\*, and the asymptote v_f⁴ = β²M_b is built into the √M sourcing choice rather than discovered. What the fit establishes is narrower and still worth stating: one normalization serves 149 galaxies with no per-galaxy freedom, its value lands at 0.55 of MOND's canonical scale, and the finite-radius slope and the individual curve shapes remain genuine tests.

**Correction 2, independent ring saturation does not give the law that was fitted.** With Ċ_i = s·m_i − b·C_i² each ring equilibrates at C_i ∝ √m_i, so a sum of rings gives Σ√m_i, not √(Σm_i): splitting the same matter into N equal pieces multiplies the field by √N. Measured here on the six best-sampled galaxies, the per-ring rule exceeds the cumulative rule by 5.5–8.0 at the outermost radius, and the excess grows with the sampling. The cumulative rule is therefore a fitting law, not a consequence of the declared mechanism. Two things follow, both run:
- **The collective control**, with no extra parameter: one shared mode saturating as C_total = β√M_total with its weight distributed in proportion to mass, so g = β·M_b(<r)/(r·√M_total). It shares B-root's outer normalization and differs inside by √(M(<r)/M_total).
- **A derivation that does work.** The spherical limit of the nonlinear field equation ∇·[(|∇ψ|/a\*)∇ψ] = 4πGρ_b gives g_ψ = √(G a\* M_b(<r))/r, which is exactly the fitted form. The square root comes from a collective nonlinear field, not from assigning a square-root charge to each ring. That is a benchmark to compare against, not a claim of path memory, and in a disk it will not equal the enclosed-mass prescription.

**Correction 3, the finite-band row used the on-ring expression at the observation radius.** (C/2R)[J₀(q_min R)² − J₀(q_max R)²] is exact only for r = R_s. The field of a ring at R_s observed at r is C∫J₁(qr)J₀(qR_s)dq. At r = 2, R_s = 1 over [0.01, 100] that is 0.49877 against the on-ring expression's 0.24989, a factor of two, reproduced independently here. The banded calculation is replaced by the source-integrated sum over rings, and the wide-band limit is verified directly rather than assumed: the ring integral is 1/r outside the ring and near zero inside (0.4988 at r = 2, 0.2000 at r = 5, 0.011 at r = 0.5, −0.005 at r = 0.2, for a ring at R_s = 1).

**Correction 4, the luminosity is band-specific.** SPARC's catalogue L₃.₆ is a 3.6 μm luminosity in that band's solar units, and multiplying it by the bolometric solar constant assumes a spectral conversion. Candidate E's deficit is relabelled as conditional on that stated conversion, and the band ratio is reported beside it. A bolometric correction of order a few does not move a deficit of 10⁹, but the number must say what it is.

**Correction 5, the radius count.** 3,152 counts rotmod rows with R > 0; the frozen comparison arrays hold 3,150. Both are reported, with the mask that produces each.

**Correction 6, candidate C's instability is analytic, and it does not transfer.** Routh–Hurwitz on Tλ³ + (1+q₀)λ² + Tλ + (1+3q₀) requires (1+q₀)T > T(1+3q₀), that is q₀ < 0. So every q₀ > 0 with finite T > 0 has a growing mode — stronger than the grid, which only sampled. The growth rate vanishes in both limits (≈ q₀T/(1+q₀)² as T → 0, ≈ q₀/T as T → ∞), so the e-folding time depends on the actual memory time and "unstable" does not mean "disrupts in one orbit at every memory time". And C's failure says nothing about B-root: in B-root's frozen potential Φ = −GM/r + K ln(r/r₀) with K = β√M, circular orbits are radially stable, κ² = GM/r³ + 2K/r² > 0. Whether they survive when matter and field evolve together is PM-2's question, not PM-1's.

## Correction 7, after the owner's review of the landed run: the mass is force-equivalent

The code builds `M_force = R·v_bar²/G = R²·g_N/G`, the mass a spherical Newtonian source would need to produce the baryon model's radial force. That is an enclosed mass only in spherical symmetry, and the SPARC inputs are disks. Substituting it into the fitted rule gives, identically,

    g = g_N + β√(M_force)/R = g_N + √(a\*·g_N),   a\* = β²/G,

verified here to 4×10⁻¹⁶ across all 149 galaxies. **So PM-1 fitted a local acceleration law — a pointwise function of the local Newtonian field — not a force sourced by the matter inside r.** The broadband ring construction motivates the latter; the code implements the former; they coincide only in spherical symmetry. No observed speed leaked into the predictor: `v_bar` is the ordinary-matter model's own prediction. The scores stand as scores of a one-parameter empirical acceleration law, and the archived run is preserved under that label.

Three consequences are carried forward rather than patched over:
- **The mass–speed regression is against force-equivalent mass.** It uses `M_force` at the outermost radius, so the reported slopes relate observed speed to that quantity, not to an independently integrated stellar-plus-gas mass. The advertised baryonic mass–speed test has not yet been done.
- **The collective control's total is not a total source mass.** `M_force` at the last sampled radius depends on where the rotation curve stops. A physical total cannot be defined by where someone stopped observing, and the invariance test is explicit in PM-2A: removing outer observation rows must not change the predicted force at the remaining inner radii.
- **The √8 demonstration is a mathematical warning, not a physical annulus test.** It differences `M_force` and clips negative increments. A source test must take annular masses from a stated nonnegative density with mass conservation checked.

## Files and reproduction

`pm1.py` runs every candidate and writes `pm1-results.json`; `checks.py` is a suite job that reruns the cheap parts and the regression anchor. About a minute on one core.
