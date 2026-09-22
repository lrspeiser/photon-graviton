# T0.2, relativistic part: where the slip can live

**22 September 2026.** The field equations in [THEORY.md](../../../THEORY.md) are
non-relativistic. This asks what a relativistic theory must look like to produce a
gravitational slip η ≠ 1 of the size measured, while surviving GW170817. Every identity
below is checked symbolically with SymPy. Every number comes from `run-v1/results.json`.

## 1. Slip is radial stress: an identity

Take any static, weak-field metric in any theory,
`ds² = −(1+2Φ)dt² + (1−2Ψ)(dr² + r²dΩ²)`. The radial Einstein equation then gives, exactly,

```
M_Φ(<R) − M_Ψ(<R)  =  4π R³ p_r(R) / c²
```

* M_Φ = R²Φ′/G is the mass that sets orbits.
* M_Ψ is the mass read from the spatial curvature.
* p_r(R) is the radial stress, at R, of whatever beyond ordinary matter sources the metric.

The script checks four things symbolically, and all four return `True`:
* the linearised G^t_t;
* the linearised G^r_r;
* the static conservation law, dp_r/dr + 2(p_r − p_t)/r = 0;
* the mass identity itself.

**The identity is not ours.** It is the radial Einstein equation. Faber & Visser
(MNRAS 372, 136, 2006) proposed combining rotation curves and lensing to measure the
pressure of dark matter this way. What is new is reading it for extra gravity without
dark matter: **inside R, lensing and dynamics can disagree only if the extra gravity carries a
radial stress at R.**
* η > 1, where light is bent more than stars are pulled, needs radial *tension*.
* η < 1 needs radial *pressure*.

The toy galaxy, 10¹¹ M☉ Hernquist under QUMOND with η = 1.34:

| Radius | g_N/a₀ | radial stress / extra energy density | sideways stress / extra energy density |
|---:|---:|---:|---:|
| 1 kpc | 7.44 | −0.130 | −0.062 |
| 10 kpc | 0.70 | −0.166 | −0.044 |
| 100 kpc | 0.011 | −0.235 | −0.009 |
| 300 kpc | 0.001 | −0.247 | −0.003 |

Far out, the stress tends to `1/η − 1 = −0.254` along the radius and to zero sideways, so it
points purely along the radius. The identity holds at every radius to 4 × 10⁻¹⁶.

## 2. The four slots

| Where the extra gravity lives | Light feels it at | η | GW170817 | Our data |
|---|---:|---:|---|---|
| Scalar, uniform stretch of spacetime (conformal): relativistic AQUAL | 0% | −1 | passes | **excluded**: we measure 91–127% |
| Separate metric for matter and light (disformal): TeVeS | 100% | 1/(2k−1), with k = 1 in TeVeS | **fails** on both speed and Shapiro delay | — |
| One shared metric, no stress: AeST | 100% | 1 | passes | misses the measured ±0.4 |
| **AeST + conformal coupling κ of matter to its scalar** | 1/(1+κ) | (1−κ)/(1+κ) | **passes** | **fits** (table below) |

Sources checked against the papers themselves, not summaries:

* **Light is blind to a conformal factor.** Skordis & Złośnik (PRL 127, 161302, 2021;
  arXiv:2007.00082) say so directly: "null geodesics are unaltered by conformal
  transformations."
* **TeVeS used a disformal metric to force equal potentials.** Stated in the same paper.
* **AeST has Ψ = Φ in its weak-field quasistatic limit.** Also stated there: "(5) leads to Ψ = Φ"
  and "Since Ψ = Φ, (6) leads to the right lensing." The earlier "as far as we can
  establish" hedge is removed.
* **Separate-metric theories are dead.** Boran, Desai, Kahya & Woodard (PRD 97, 041501, 2018)
  show GW170817 rules out theories that put matter on a different metric from
  gravitational waves: |γ_GW − γ_EM| < 9.8 × 10⁻⁸.

**Why the last slot passes GW170817.** A conformal factor does not move null cones, so light and
gravitational waves still share the same paths. They also share the same speed and the same
Shapiro delay. Both feel only Φ + Ψ, and in this construction that combination is exactly
what AeST gives. The slip lives entirely in Φ − Ψ, which only slow-moving matter feels.

**Within AeST-type theories this is the natural slot:**
* A nonminimal G₄(ϕ)R coupling is the same slot seen in another frame.
* Aether (vector) terms give no slip at this order: Einstein-aether has γ = 1 (Foster & Jacobson 2006).
* Couplings of type G₄(X) or G₅ change the speed of gravitational waves.
* Beyond-Horndeski theories offer other routes, not explored here.

## 3. What each measurement requires

| Systems | η | Light / matter | κ | Radial stress / mean extra density | Far-field w_r |
|---|---:|---:|---:|---:|---:|
| 3 over-bent SLACS lenses | 0.82 – 0.88 | 0.91 – 0.94 | +0.064 – +0.099 | +0.045 – +0.073 | +0.14 – +0.22 |
| 11 X-COP clusters | 1.27 – 1.42 | 1.14 – 1.21 | −0.174 – −0.119 | −0.099 – −0.071 | −0.30 – −0.21 |
| 3 under-bent SLACS lenses | 1.36 – 1.54 | 1.18 – 1.27 | −0.213 – −0.153 | −0.117 – −0.088 | −0.35 – −0.27 |

In words:
* In clusters and the under-bent lenses, stars and gas feel the AeST part minus a 12–21%
  outward uniform-stretch push, while light feels the whole AeST part.
* In the over-bent lenses the uniform-stretch part pulls inward instead, by 6–10%.

## 4. Checks

* **Solar system.** At the Cassini ray, 1.6 R☉, the extra gravity is 8.6 × 10⁻¹² of the
  Newtonian deflection. Slip moves γ by at most **3.7 × 10⁻¹²**, against a bound of
  2.3 × 10⁻⁵.
* **The one published direct cluster slip measurement.** Pizzuti et al.
  (arXiv:1602.03385) used galaxy kinematics against lensing in the relaxed cluster MACS J1206.2-0847 and
  found η(r₂₀₀) = 1.01 +0.31/−0.28 for the *total* potential. Our cluster slip predicts
  η_total = 1 + f(η − 1) = **1.21–1.37** for f = 0.79–0.87. That is 0.7–1.2σ above their centre:
  consistent, at the edge, and not a confirmation. Galaxy orbits do not feel gas pressure, so
  this route escapes the slip-versus-non-thermal-pressure degeneracy.

## 5. What this does not do

1. **Set κ.** A single κ cannot give η > 1 in clusters and η < 1 in three lenses.
   * Either κ depends on something not yet identified, or the over-bent lenses carry a
     systematic. The stellar-mass convention, roadmap T3.5, moves lens results by 2.4×.
   * One public lead: Brouwer et al. (A&A 650, A113, 2021) find early- and late-type galaxies
     of equal stellar mass on different weak-lensing acceleration relations, at ≥6σ. A
     universal modification cannot produce that. A type-dependent κ could. Not yet tested.
2. **Background drift.** If the scalar's background value drifts in time, particle masses
   drift at κ times that rate, and lunar laser ranging bounds such drifts near 10⁻¹³ per
   year. A static background avoids the issue.
3. **Prior art.** Both ingredients are published: AeST, and the conformally coupled scalar
   of Bekenstein & Milgrom's 1984 relativistic AQUAL. We have not found the combination
   proposed as a source of slip. A full literature check belongs to roadmap T3.3.
4. **The simple ν in the solar system.** This is separate from slip. With the simple interpolation
   function the extra pull tends to a₀ at high acceleration rather than switching off. Hees et al.
   (MNRAS 455, 449, 2016) find that Cassini rules out several popular MOND transition functions
   and leaves others viable. Whether ours survives is roadmap T2.3.

## Reproduction

```
cd code
python slip_slots.py --output-dir ../run-v1
```

Requires NumPy and SymPy. No input data and no network access are needed. The symbolic
Einstein-tensor step takes about a minute.
