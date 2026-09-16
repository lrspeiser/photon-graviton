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

## Files and reproduction

`pm1.py` runs every candidate and writes `pm1-results.json`; `checks.py` is a suite job that reruns the cheap parts and the regression anchor. About a minute on one core.
