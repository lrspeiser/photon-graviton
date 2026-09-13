# Can a continuum of receiving states remove the color dependence?

13 September 2026. Analytic restriction on the existing conditional point-response interaction.

## Outcome

No: replacing its single gap with any positive, incident-energy-independent distribution of gaps cannot give constant fractional redshift. Wherever its forward fractional loss is nonzero and finite, its logarithmic energy slope is at least three. A twofold increase in photon energy produces at least an eightfold increase in the fractional loss rate through the same receiving environment.

This is a restriction on a particular interaction and approximation, not on every hypothetical photon-companion process. It does not use an expanding universe, dark matter, a universe age or a fitted astronomical linewidth limit. It shows that tuning only this interaction's receiving spectrum cannot meet our common-redshift target.

## Assumptions and provenance

The [existing derivation](../pair-production-balance/derivation.md) uses a weak-coupling, heavy-store, point-electric-dipole interaction. Forward creation of a mode pair with gap Delta gives sigma(E,Delta) proportional to E*(E-Delta)^3 above threshold. This rate was derived for our stipulated operator, not measured in nature. Its transition-rate method is the known golden-rule calculation; see [Tong, section 3.6.1](https://davidtong.org/pdfs/teaching/quantum-field-theory/qft3.pdf).

Let W(Delta)>=0 absorb mode density, receiving-store density and squared coupling, with constants chosen so alpha below has units of inverse path length. W can depend on gap but must not depend on incident photon energy. This calculation includes only forward creation into fixed available states, with no inverse transfer, changing occupation, spatial form factor or coherent response. We assume the integrals exist.

**Known energy-loss accounting applied to that proposed rate:** multiply each event rate by its fractional photon loss Delta/E:

\[
\alpha(E)=C\int_0^E W(\Delta)\,\Delta\,(E-\Delta)^3\,d\Delta,\qquad C>0.
\]

**Analytic consequence derived here; no claim of literature novelty:** differentiating gives

\[
E\alpha'(E)-3\alpha(E)
=3C\int_0^E W(\Delta)\,\Delta^2\,(E-\Delta)^2\,d\Delta\geq0.
\]

The moving-upper-bound term vanishes because the outgoing phase-space factor is zero at Delta=E. Hence

\[
\frac{d\ln\alpha}{d\ln E}\geq3,\qquad
\frac{\alpha(E_2)}{\alpha(E_1)}\geq\left(\frac{E_2}{E_1}\right)^3
\quad(E_2>E_1,\ \alpha(E_1)>0).
\]

Equivalently, alpha/E^3 is nondecreasing. This version extends directly to positive discrete sums or mixed spectra. A zero-gap component carries no energy in this forward formula; a divergent spectrum invalidates finite-rate assumptions rather than supplying a finite achromatic solution.

## Numerical examples and checks

The executable uses 12 illustrative spectra W proportional to Delta^q with q=-2,-1,0,1, lower gap 1e-8 eV and upper gaps 0.001,1,100 eV. It evaluates seven photon energies from 0.01 to 10 eV. All 84 calculations obey the inequality; slopes range from approximately 3.00003 to 6.00000. Fractional-loss ratios between 10 eV and 1 eV range from about 1,000 to 1,000,000. The rates are normalized at 1 eV to compare shapes, not calibrated to observations.

Doubling logarithmic quadrature from 128 to 256 nodes changes rates by less than 1e-10 fractionally. An independent polynomial antiderivative checks the q=0 case to 1e-12 relative accuracy. These computations illustrate the analytic result; the sampled power laws do not establish its generality—the positive-integral argument does.

The energy range is a diagnostic grid, not a claim that we fitted actual radio-to-gamma data. The analytic bound applies across any range where the stipulated interaction and approximations are valid. If the approximation breaks down, its replacement must supply the rate rather than extrapolating this expression.

## What would need to change

1. **Spatial or energy-dependent interaction.** If a multiplicative squared-amplitude factor F(E)^2 were added, flattening the mean would require d ln(F^2)/d ln E=-d ln alpha/d ln E, at most -3. In a small-gap limit this corresponds approximately to F proportional to E^(-3/2). This is an inverse requirement, not a derived form factor. The output records the exact compensating factor relative to 1 eV; multiplying by it would guarantee a flat curve by construction and would not validate physics. A real spatial calculation must also predict the angular distribution and test image preservation.
2. **Occupation and inverse transfer.** Net loss can involve terms of both signs, so the positivity argument no longer applies. The earlier production/loss balance already shows finite stationary occupation for a conditional mode. A proposed cancellation must follow from the common rates and illumination history, preserve stored energy accounting and work across colors; it cannot use separate fitted subtractions for each wavelength.
3. **Coherent evolving field.** A field that changes optical phase continuously falls outside this incoherent forward-event model. It must derive its own evolution, source/receiver clocks, energy exchange and event-duration behavior. It is not licensed to retain this forward rate while claiming its spectral bound has vanished.

The previously written scale-invariant fractional-jump kernel is therefore a target requiring different microscopic structure, not something obtained just by increasing the number of bound modes in the point interaction.

## Research decision

Do not spend further fits on positive gap spectra alone as a repair to this direct-redshift interaction. Preserve the calculation for production into a reservoir or as an ingredient of a changing field, where its role and observations differ. Before comparing a replacement to redshift data, derive the energy and angular dependence of an explicit spatial or coherent interaction. The one-third deposition reference and separate saturation work remain unchanged.

Reproduce with `python research_work/results/companion-extensions/continuum-color.py`. [Source](continuum-color.py) and [all numerical results](continuum-color-results.json).
