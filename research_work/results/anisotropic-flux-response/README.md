# JR-10: the companion's stress state, not its amount, sets lensing versus kinematics

**22 September 2026. An explicitly empirical search, requested as such, followed by the
separation of what is exact general relativity from what is a fitted law.
Historically examined objects; not a blind test.**

## Why this direction

JR-9 left two hard facts. The companion's amplitude is pinned by resolved stellar
motions measured to 1–2 percent, so making it stronger is blocked. And the energy
ledger forbids paying for the companion's mass out of starlight by about 4.9×10⁴ per
lens — consistent with RC-1's 10⁵–10⁸ against the microwave background and earlier
10³–10⁶ comparisons.

So the companion cannot be static dust whose mass we have to fund. This experiment
tests the alternative that costs no extra energy at all: **the same stress-energy,
arranged anisotropically.**

## The part that is exact, with no fitting anywhere

In linearised GR with `ds² = −(1+2Φ)dt² + (1−2Ψ)dx²`:

```
∇²Φ = 4πG (ρ + T^k_k)      slow stars feel Φ
∇²Ψ = 4πG  ρ               light feels Φ + Ψ
```

Writing `w_r = p_r/ρ` and `w_t = p_t/ρ` and `τ = 1 + w_r + 2w_t` for the Tolman
active-mass factor:

```
γ_χ  ≡  Ψ_χ / Φ_χ  =  1 / τ
```

and the companion enters the deflection at weight `(1 + γ_χ)/2` — which is exactly
the R11 form, now with a derived meaning instead of a free coefficient.

| Stress state | w_r | w_t | τ | γ_χ |
|---|---:|---:|---:|---:|
| static dust (what R10 assumed) | 0 | 0 | 1 | 1 |
| relativistic gas, p = ρ/3 | ⅓ | ⅓ | 2 | ½ |
| free-streaming radial null flux | 1 | 0 | 2 | ½ |
| radial tension | −½ | 0 | ½ | 2 |
| **static radial field gradient** | **1** | **−1** | **0** | **infinite** |

The last row is exact and is the one the data has been asking for since JR-8:
`T^k_k = −ρ` cancels `ρ` in the Tolman active mass, so **a static radial field
gradient exerts no force at all on slow stars while its energy density still deflects
light.** Lensing without kinematics, paid for out of stress rather than new energy.

## Stage 1: what each lens actually demands

`code/stress.py` → [`run-v1/required-stress.json`](run-v1/required-stress.json)

Under the JR-9 cross-prediction protocol — stellar-mass offset and anisotropy fixed by
the resolved kinematics alone, universal R10 coefficients frozen — solve for the single
stress ratio that reproduces each catalog Einstein angle.

| Lens | z | Angle error at τ=1 | Required τ | Required γ_χ | w (if w_t=0) | In band |
|---|---:|---:|---:|---:|---:|:---:|
| J1204+0358 | 0.1644 | +5.68% | 1.2180 | 0.821 | +0.218 | yes |
| J1402+6321 | 0.2046 | +4.63% | 1.1423 | 0.875 | +0.142 | yes |
| J0037-0942 | 0.1955 | +4.29% | 1.1290 | 0.886 | +0.129 | yes |
| J1621+3931 | 0.2449 | −11.69% | 0.7371 | 1.357 | −0.263 | yes |
| J1112+0826 | 0.2730 | −14.61% | 0.6594 | 1.517 | −0.341 | yes |
| J1630+4520 | 0.2479 | −15.81% | 0.6495 | 1.540 | −0.350 | yes |

**Every lens is satisfiable with a modest, physically allowed stress anisotropy.** The
required τ span 0.65 to 1.22 — well inside the band the dominant energy condition
permits — and the implied stress ratios are of order ±0.3, not exotic values. The
over-bent three want mild pressure, the under-bent three want mild tension. No extra
energy, no extra mass, no per-lens gravity multiplier.

Six numbers from six angles is a measurement of the requirement, not a test. The test
is whether one law of an independently measured property reproduces them.

## Stage 2: only redshift predicts it

`code/law_search.py` → [`run-v2-law/law-search.json`](run-v2-law/law-search.json)

Every candidate is scored by **leave-one-out prediction of the Einstein angle**: fit
the law on five lenses, predict the sixth, never having used its angle. Predictors
built from the observed Einstein angle (the Einstein radius in kpc, θ_E/R_e) are
excluded as circular — τ is solved from that same angle — and reported separately.

| Predictor | Form | Params | In-sample | **Leave-one-out** |
|---|---|---:|---:|---:|
| *static dust, τ=1* | *baseline* | *0* | *10.58%* | ***10.58%*** |
| lens redshift | exponential | 2 | 3.36% | **4.89%** |
| lens redshift | linear | 2 | 3.66% | 5.66% |
| source redshift | exponential | 2 | 8.06% | 9.98% |
| stellar mass, R_e, compactness, surface density, aperture σ, axis ratio, Sérsic n, … | any | 1–2 | 9.94% | **11.73%** |

**No structural property of the galaxy predicts its stress state.** Mass, size,
compactness, velocity dispersion, surface density and light-profile shape all collapse
to the constant fit, which scores *worse* than assuming static dust. Only redshift
carries information.

## Stage 3: the locked formula, bounded by construction and permutation-tested

`code/bounded_law.py` → [`run-v3-bounded/bounded-law.json`](run-v3-bounded/bounded-law.json)

Stage 2's best form extrapolated to τ = 4.16 at z = 0 — outside the dominant energy
condition. A law that leaves the physical band is not the law. So reparameterise using
the two exact endpoints of a radial stream:

```
w_r = 1,   w_t = −f,   τ = 2(1 − f),   f ∈ [0,1]
```

`f` is the fraction of transverse tension developed — how far the stream has gone from
free-streaming (f=0, τ=2) to fully anchored (f=1, τ=0). Fitting `f` through a logistic
link keeps it in [0,1], so **τ stays in (0,2) and γ_χ stays above ½ for any input, at
any redshift, permanently. No fitted value can violate the energy conditions.**

### The locked formula

```
f(z) = sigmoid( 4.3550 + 2.7552 · ln z )
τ(z) = 2 (1 − f)
γ_χ  = 1 / τ
```

| | |
|---|---:|
| In-sample angle RMS | **3.35%** |
| **Leave-one-out angle RMS** | **4.66%** |
| Static-dust baseline (0 parameters) | 10.58% |
| Permutation test, 2000 shuffles of the whole search grid | null median 8.04%, 5th pct 5.37% |
| **Look-elsewhere corrected p** | **0.0165** |

Leave-one-out, per lens, each predicting its own angle having never seen it:

| Lens | τ predicted | γ_χ | Angle error |
|---|---:|---:|---:|
| J0037-0942 | 1.0527 | 0.950 | +2.41% |
| J1112+0826 | 0.6059 | 1.651 | +3.79% |
| J1204+0358 | 1.4282 | 0.700 | −3.83% |
| J1402+6321 | 0.9764 | 1.024 | +5.52% |
| J1621+3931 | 0.7741 | 1.292 | −2.13% |
| J1630+4520 | 0.7852 | 1.274 | −7.79% |

The law crosses the static-dust value τ = 1 at **z = 0.2058**, matching the split
JR-9B found independently at z ≈ 0.22 before this search existed. Its z → 0 limit is
**τ = 2 exactly** — the free-streaming radial null flux state.

## Stage 4: what this derives, and what it plainly does not

`code/derive_stress.py` → [`run-v4-derivation/derivation-and-bounds.json`](run-v4-derivation/derivation-and-bounds.json)

### Stress engineering cannot solve the energy problem

The dominant energy condition bounds each w to [−1,1], so τ ≤ 4. The energy density a
fixed stellar force needs scales as 1/τ, so **the largest saving any physically allowed
stress state can produce is a factor of four.** Against a measured shortfall of
4.9×10⁴, that leaves 1.2×10⁴.

**Anisotropic stress fixes the direction of the response, not its cost.** The energy
arrow is untouched by anything in JR-10 and remains the deepest open problem in the
portfolio. This result is worth having because it closes a direction: no rearrangement
of stress, however clever, bridges four orders of magnitude.

### The driving variable is not identified

| Pair | Pearson | Rank |
|---|---:|---:|
| lens redshift vs lens distance | **+1.000** | **+1.000** |
| lens redshift vs Einstein radius (kpc) | +0.860 | +0.886 |
| Einstein radius vs lens distance | +0.861 | +0.886 |

Under this project's distance convention, lens distance is a deterministic function of
redshift, so **"the law depends on redshift" and "the law depends on lens distance" are
literally the same statement in these data.** The physical Einstein radius is collinear
too. Six systems spanning z = 0.164–0.273 cannot separate them.

This matters for the next step. Deriving a redshift mechanism now would mean deriving a
dependence on a variable the data has not identified. **The next acquisition, not the
next derivation, is what advances this.** What breaks the degeneracy: lenses at matched
redshift with different Einstein radii, or matched Einstein radii at different
redshifts. Neither exists in this six-system set.

### The discipline is thinner than it looks

Because γ_χ only changes lensing, no V_rms prediction moves and not one of the 149
SPARC rotation curves is affected — which also means **SPARC cannot constrain this
law.** The discipline comes only from the energy-condition band (which every required
value satisfies), the leave-one-out score, and the permutation test. It does not come
from a second observable. That is a genuine weakness and it is why the p = 0.0165 is
suggestive rather than decisive.

### The local prediction, with its caveat

The law sends τ → 2 as z → 0, so around nearby galaxies the companion should deflect
light at (1 + ½)/2 = **¾ of the rate a static-dust companion of the same stellar force
would produce** — a testable statement about local weak lensing at fixed measured
rotation, which no lens in the fitted set could have produced.

The caveat is real: the z → 0 limit follows from the log-linear form, which the six
systems cannot distinguish from the linear form inside the observed range (LOO 4.66%
versus 5.17%). The linear form does **not** give τ = 2 at z = 0. This prediction is
conditional on a form choice made after seeing the fits, and should be labelled that
way wherever it is quoted.

## Provenance and honesty about ordering

This package was **not** preregistered. It was an explicitly requested empirical search
— find the formula that works, then derive it — and it is labelled as such. What
protects it from being a free fit is recorded rather than assumed: the search grid is
declared in the code (16 predictors × 3–4 forms, 63 candidates scored), the circular
predictors are excluded by name, scoring is leave-one-out throughout, and the
look-elsewhere correction is a permutation test over that same grid rather than a
hand-waved allowance.

Everything in "the part that is exact" is textbook linearised GR and was written down
before the search. Everything in stages 2 and 3 is a fit.

## Reproduction

```
cd code
python stress.py        --output-dir ../run-v1            # ~5 s
python law_search.py    --output-dir ../run-v2-law
python bounded_law.py   --output-dir ../run-v3-bounded     # ~45 s incl. 2000 permutations
python derive_stress.py --output-dir ../run-v4-derivation
```

`code/model.py` is a symlink to the JR-9 loader, which reproduces the published R10
lens predictions to zero relative difference before anything else runs. Requires NumPy
and SciPy; no network access.
