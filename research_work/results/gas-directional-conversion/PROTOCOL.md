# JR-5: gas-mediated directional conversion

21 September 2026, America/Los_Angeles. Baseline a5805452b97ef0856d364de90e394572a8fd2146. Exploratory construction: toy simulations precede observational fitting. Earlier results and raw measurements are preserved.

## Theory first

Matter/radiation may exchange energy with a propagating companion and a bound companion reservoir. Gas may aid or disrupt that exchange in directions selected by its actual spatial distribution. This experiment adds a gas-dependent conversion/return mechanism, not a new independently fitted halo or lensing multiplier. It does not narrow the broad research portfolio to this implementation.

## Explicit transport mechanism

Along a ray, photon-like energy L and traveling companion C have the same trial speed v. A bound reservoir B captures C at rate kc and releases it back at kr. Their energy-density-per-ray equations are

    L_t + v L_r = -k_plus L + k_minus C
    C_t + v C_r = k_plus L - k_minus C - kc C + kr B
    B_t = kc C - kr B.

The sum has only source/boundary flux. All rates are nonnegative; the full finite-volume evolution retains recoil/thermal channels only as unresolved microscopic physics, not implicitly solved. Equal trial travel speeds and phenomenological transition rates are model postulates. No photon energy is counted again as an independently generated gravitational reservoir.

In steady state B=(kc/kr)C and the traveling companion flux fraction f obeys

    df/dr = a(r,n)*(1-f) - b(r,n)*f,
    a=alpha*rho_g/(1e6 Msun/kpc^2),
    b=beta*(rho_g/rho_ref)*rho_g/(1e6 Msun/kpc^2),
    rho_ref=1e7 Msun/kpc^3.

alpha,beta are nonnegative universal coefficients. Linear forward exchange and density-quadratic reverse exchange are hypotheses motivated by encounter versus correlated-encounter sensitivity, NOT derived quantum cross sections. A comparison has both rates linear. Exact constant-cell transfer is f_next=f_eq+(f_previous-f_eq)*exp[-(a+b)dr]. It preserves 0<=f<=1 and the sum of photon/companion traveling flux. The supplied central outgoing mixture f(0)=1/2 represents a fixed emitter conversion assumption. It is not a claim gas-free conversion is microscopically solved.

The gas-free capture/residence kernel is fixed to the earlier R10 spherical reservoir, so B_new(r,n)=2*f(r,n)*B_R10(r). This follows conditionally from holding kc/kr and the incoming total flux fixed. It retains R10's empirically inferred source normalization and spatial capture kernel; it does NOT derive the old envelope, source luminosity or energy supply from scratch. In the stationary equal-speed model the sum of traveling energies is unchanged by gas, while stored energy can change because capture/release delays differ from free passage. Gravitational work and full matter backreaction remain outside this kinetic toy.

## Toy campaign

Use equal-total-mass gas distributions: sphere, thin exponential disk, thick disk and axisymmetric gas ring. Central and oblique rays see different columns/densities. Sweep forward and reverse coefficients over 0,0.03,0.3,3; retain controls and all combinations. Predict the companion's equatorial and polar radial forces and face-on/edge-on lens deflection from ONE three-dimensional axisymmetric mass distribution using even Legendre multipoles, not independent directional force multipliers. Compare a monopole-only approximation as an attribution control. Gas-free, forward-only, reverse-only, rotation-of-observer and local energy controls are required. Run an explicit source-fed initially empty three-channel finite-volume evolution and compare its stationary flux fractions with the corresponding finite-volume steady solution; no finite-age astrophysical interpretation is asserted.

## Observational transfer

Keep all 149 archived SPARC galaxies, the original 89/29/31 roles, observed speeds/errors, distances, inclinations, stellar normalization, gas mass and universal R10 parameters. Estimate one exponential gas radial scale per galaxy by matching its published gas-only force component at its catalog gas mass; this uses no observed total rotation residual. Its vertical exponential height is provisionally 0.1 gas radial scales. These are inferred axisymmetric gas profiles, not resolved HI maps; mismatches to the gas-force input and scale-bound fits must be reported. Sensitivities use heights 0.05 and 0.2 without fitting individual heights. Real directional asymmetry cannot be validated from the one-dimensional rotation archive.

Compare null, forward-only, reverse-only, two linear rates, and linear-forward/quadratic-reverse rates. Fit at most two universal nonnegative conversion coefficients (0..30) on 89 training galaxies; select on the 29-role mean equal-galaxy fractional-squared error including the null. Report 31 comparison results only after selection. All data were previously examined, so no fresh blindness or confirmation claim. Also report ordinary velocity RMSE, median fractional RMS, number below 10/20 percent and quoted-error chi-square. No uncertainty floor or source edits. Preserve unfavored models and optimizer bound hits.

The same energy-partition law and angular density-to-force calculation used in the toys supplies the observational screening. Local source/capture assumptions and the axisymmetric gas reconstruction are separate limitations. Lensing toys use the same potential, but the six archived SLACS lenses have no usable gas maps in this input package: do not set their gas to zero and call unchanged predictions validation. Clusters also require their actual gas distribution and are not claimed tested without it.

## Numerical reliability and output

Double radial/angular resolution and multipole order at frozen coefficients; distinguish a numerical uncertainty from an observational discrepancy. Use independent direct volume force or lens integration for selected toy controls. Output code, input hashes, toy sweep, source-only gas reconstructions, all per-galaxy predictions, fitting history, declared selection and failures. A promising fit is an effective source-interaction candidate, not proof of a microscopic conversion law, adequate fuel, planetary recovery or gas/star-formation feedback.
