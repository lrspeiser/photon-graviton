# Nonthermal escape and occupation of the returning modes

13 September 2026. Optically thin emission diagnostic before reabsorption feedback.

## Result

Simply allowing low-energy relaxation waves to escape does not make their coupled modes approximately empty in the tested galaxy-scale examples. At the earlier orbital heat budget and very small quantum energies, the free-streaming population is large even for long emission histories. Therefore the protected-state calculation cannot justify negligible upward excitation merely by saying the radiation is nonthermal and outgoing.

This does not determine the actual optical depth or reject every escaping channel. The result assumes the outgoing waves occupy the same modes that couple to the reverse transition. Their spectral, angular, polarization and spatial overlap must be derived. Absorption can reduce the outgoing population, but doing so returns energy to the receiving material; it cannot be counted as successful cooling without tracking that return.

## Geometry and mode count

**Assumptions:** uniform emission throughout a sphere, isotropic initial directions, free propagation at c, two bosonic polarizations, no absorption or scattering, and a constant occupation n across a stipulated logarithmic energy band. This is not a thermal distribution. It is also not the actual three-dimensional source distribution of the Milky Way.

**Geometrical result:** the mean distance to the boundary from uniformly sampled volume and direction is 3R/4. Thus, for a steady luminosity L, mean radiation energy inside the sphere is L*3R/(4c), and volume-averaged energy density is

\[
\bar u=\frac{L\,t_{\rm esc}}V,\qquad
t_{\rm esc}=\frac{3R}{4c},\quad V=\frac{4\pi R^3}3.
\]

The code independently integrates the path length over position and angle. For a source active less than the largest crossing time, it replaces each path residence time by the smaller of that time and the source duration. The finite-fill correction prevents applying a steady field before it exists.

**Known massless-mode density**, applied to a nonthermal constant occupation rather than a Planck distribution:

\[
\bar u=\frac{g n}{8\pi^2\hbar^3c^3}
(\varepsilon_{\rm hi}^4-\varepsilon_{\rm lo}^4).
\]

Here epsilon_hi/lo=epsilon_0*exp(+/-b/2), with logarithmic band width b. The density-of-states factor is standard quantum statistics; see [Tong's quantum-gas notes](https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S3.html). The spectral shape, shared coupling and two-polarization choice are hypotheses.

In the weak incoherent bosonic coupling used by the protected-state model, upward and downward factors are n and n+1. Their ratio n/(n+1) therefore requires n<=0.001001 to keep upward/downward excitation at or below the earlier design ratio 0.001. This uses occupation of the coupled modes, not a temperature assignment. A volume-averaged occupation only diagnoses the stated homogeneous approximation; it is not a measured local transition rate.

## Numerical example

Use the previously calculated orbital heat Q=4.471e50 J, an illustrative bright gap E_a=1e-8 eV, and release quantum centered on epsilon_0=E_a/2. In a 120-kpc sphere, mean free residence is about 293,541 years. If this heat is emitted uniformly over one billion years, L=1.417e34 W.

| Logarithmic emission-band width | Mean occupation after free streaming | Required emission duration for n<=0.001001 |
|---|---:|---:|
| 0.01 | 4.67e16 | 4.67e28 years |
| 1 | 2.58e14 | 2.57e26 years |

Both examples have upward/downward factors extremely close to unity if the same states couple as assumed. This is not proof of immediate reabsorption: the absolute coupling and spatial overlap remain missing. It is a failure of the negligible-occupation assumption for this emission model.

The broader band distributes energy among more available modes but is an optimistic alternative to a sharply defined transition. A literal single narrow gap cannot automatically emit across that whole band. The required duration follows from Q divided by the maximum luminosity compatible with the occupation allowance, not from an assumed universe age. Giving the universe more time would require adopting and testing that actual source history; no conventional age or horizon is imposed.

## Scope and checks

The output covers both ordinary-matter baselines, two excitation energies, the prior orbital heat and the separate conditional all-deposit-formation heat, three container radii and two bandwidths: 48 configurations. Each includes source durations 1e6,1e9,1e12 years as diagnostic histories, giving 144 emitted-population examples. These histories are not measurements or adopted formation ages. The all-deposit case retains the explicit assumption that the inventory energy comes from protected excitations, rather than a separately existing seed.

No-absorption residence times are verified by angular/volume integration; finite source duration is retained where required. Energy in flight is source power times the effective residence. The mode-count integral then fixes occupation. These are source-population calculations, not a completed radiative-transfer solution or a fit to images, spectra or luminosities.

Uniform free escape is not a universal lower bound on occupation at every receiving site. Actual source geometry and selection rules can reduce coupling in some locations or directions, while trapping can increase residence. A coupled transport calculation must determine those effects together. Coherent emission would also invalidate the simple incoherent rate treatment rather than automatically inherit its protection result.

## Decision

Do not use direct escape alone to justify the empty-bath approximation at these energy scales and release powers. The next physically discriminating step is emission/absorption transport with the receiving-state populations included. It should calculate whether occupation suppresses net cooling, redirects emission, or empties protected states, using common couplings. A different high-energy or weakly coupled outlet is another hypothesis requiring an energy-conserving transition, not an arbitrary new sink.

The paper retains the conditional storage possibilities but now records both thermal capacity and nonthermal occupation constraints. The microscopic state spectrum, source supply, protected lifetime, redshift timing and gravitational predictions remain unresolved; the exact-third empirical reference is unchanged.

Reproduce with `python research_work/results/companion-extensions/escape-occupation.py`. [Source](escape-occupation.py), [complete results](escape-occupation-results.json).
