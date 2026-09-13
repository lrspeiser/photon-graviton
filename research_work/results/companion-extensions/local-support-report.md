# Local support audit: the ideal phase interpretation needs revision

13 September 2026. Goal progress on settling, support and companion states. Exact-third capture and the shared fitted phase parameters are frozen; no stellar-speed parameter is refitted.

## Finding

The successful inner rotation fit is a useful density profile, but it is not a demonstrated ideal-Bose condensate equilibrium. Computing the pressure needed at each radius exposes a large mismatch with the transported compact fraction. This strengthens the distinction already made between a phase-inspired selector and a physical phase solution.

Within 0.1-30 kpc, 90.4%/82.3% of the settled companion mass under matter baselines I/II lies where required pressure exceeds the ideal condensed-phase ceiling at the fitted effective mass 17.78 eV/c squared. Above that ceiling an ideal gas could be in a hotter noncondensed state; the result does not mean that no gas can support the profile.

At 5 kpc the pressure is below the ceiling, but only a small condensed fraction is compatible with it:

| Quantity at 5 kpc | Baseline I | Baseline II |
|---|---:|---:|
| Assigned compact fraction of local density | 81.9% | 85.4% |
| Condensed fraction implied by ideal pressure support | 2.53% | 8.31% |
| Extra pressure needed if assigned fraction is retained | 94.0% | 95.3% |
| Required isotropic velocity scale | 151.0 km/s | 143.9 km/s |

Local fractions differ from the earlier whole-inventory mobile fractions: inward transport concentrates the selected component. A compact transport label need not remain a thermodynamic condensed state. If the label is only a history tag, the present mismatch does not invalidate transport, but the condensation explanation for its stopping radius is lost.

## Equations and provenance

For the stipulated spherical reservoir and an angularly averaged ordinary-matter force, isotropic support requires

\[
\frac{dP}{dr}=-\rho g,\qquad P(r)=\int_r^{r_{\max}}\rho(x)g(x)\,dx,
\quad g=g_{b,\rm monopole}+GM_d(<r)/r^2.
\]

These are known hydrostatic/isotropic Jeans relations, with zero imposed outer pressure. They are applied to our proposed density, not derived as a new law. The [Jeans derivation](https://galaxiesbook.org/chapters/I-04.-Equilibria-of-Collisionless-Stellar-Systems_4-The-Jeans-equations.html) explains the connection between force and dispersion. Spherical averaging cannot establish three-dimensional equilibrium in a flattened galaxy. A Jeans solution also does not guarantee a positive distribution function or stability.

For a homogeneous, noninteracting Bose gas below its critical temperature, thermal pressure is

\[
P=\frac{k_BT}{\lambda_T^3}\zeta(5/2),\qquad
P_c=\frac{\zeta(5/2)}{\zeta(3/2)}n k_BT_c,
\]
\[
\frac{P}{P_c}=(T/T_c)^{5/2},\qquad
f_{\rm cond}=1-(P/P_c)^{3/5}\quad(P\leq P_c).
\]

These are established ideal-gas formulas; see [Tong's quantum-gas treatment](https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S3.html). Using them for companions is hypothetical. Below condensation the classical equality P=n kB T is not the correct ideal-Bose equation of state. Earlier T=m sigma squared/kB was explicitly an assumed kinetic selector, and must not be upgraded into a thermodynamically supported condensate without this check. Interactions and finite-size wave-gradient stresses are excluded by the ideal homogeneous approximation.

The settled density is reconstructed as

\[
\rho_{\rm new}(r)=[1-f_0(r)]\rho_0(r)+s^{-3}f_0(r/s)\rho_0(r/s).
\]

The factor s^-3 is the known volume Jacobian for the proposed contraction. It preserves the original capture inventory. No new density normalization was fitted.

If the transported compact fraction f_assigned is interpreted as the present condensed fraction, it supplies thermal pressure Pc(1-f_assigned)^(5/3). The remainder P_required minus this pressure must be supplied by some other stress or the profile must change. This is a pressure accounting diagnostic, not a derived additional force.

## Creative alternatives now have a quantified target

1. Interacting bound states. Repulsive interactions can contribute pressure beyond ideal thermal pressure. The required radial residual is calculable; a shared interaction law must reproduce it and be solved self-consistently, rather than assigning a free pressure at every radius.
2. Orbital support. Compact material might be a population of bound objects or excitations with angular momentum rather than a cold condensate. This retains a possible settling history but needs a physically nonnegative orbital distribution and an explanation of angular-momentum transport.
3. Separate components. A cold bound population could coexist with a hot or wave-supported reservoir. Its cross-component coupling must transfer support consistently; pressure in one uncoupled population does not automatically hold up another.
4. Different effective mass or dispersion. At the fixed 5-kpc profiles, matching the assigned fraction by changing only effective mass would require about 6.21/5.64 eV/c squared. These are local inverse estimates, not new globally fitted masses. Changing mass changes the original selector and density, so a full rerun is required. A radius-by-radius mass choice would not explain the system.
5. Wave gradients or nonequilibrium states. A collective excitation may not obey the homogeneous ideal-gas model. Derive its stresses and formation history explicitly. Merely calling the reservoir a wave does not provide pressure or guarantee stability.

The next useful experiment is to test a small family of shared interaction-pressure laws against the required residual, then solve their equilibrium profiles. This can determine whether an interaction offers a physical stopping rule. The propagation/time and nuclear/diffuse emission tracks remain open; none is canceled by this local support result.

## Verification and scope

The [protocol](local-support-protocol.md) preceded the calculation. [Executable](local-support.py) and `local-support-results.json` preserve four profiles, radial samples and refinements. An independent Plummer analytic-pressure integral agrees within 1e-5 relative tolerance. Density integration reproduces capture mass within 4e-6. The maximum local derivative residual in the reporting region is 0.20%; sharper selector transitions dominate it. Doubling resolution changes sampled pressure ratios by less than 0.0016%, and halving the distant boundary changes them by less than 0.0006%. Thresholded mass fractions change by less than 0.027 percentage points. Results are numerically robust within the stated approximation.

No new observational agreement is claimed. This is evidence that the simple equilibrium-condensate interpretation needs revision, while the prior empirical rotation fit and energy-sign ledger remain what they were. A global virial balance is insufficient to establish local phase support. The full research goal remains active.
