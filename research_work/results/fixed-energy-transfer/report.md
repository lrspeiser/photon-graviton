# Many small transfers: an explicit mean-redshift process

2026-09-10. This is a proposed statistical energy-transfer mechanism, not a first-principles action or a validated unified theory. It advances candidate discrimination and spectral predictions using an already exposed radio observation. [Protocol](protocol.md), [calculation](run.py), [results and hashes](results.json).

**Main result:** a photon can lose the same small energy amount on each event, with event rate proportional to its remaining energy, and reproduce the existing exponential mean-redshift rule. Unlike stationary whole-mode conversion, this process lowers the energies of surviving photons. It also necessarily broadens their energy distribution. The chosen unchanged-speed implementation still fails to produce event stretching.

## Proposed rule and known solution

**New candidate postulate, without a novelty claim:** restrict this idealized model to an energy ladder E=m*epsilon. Each downward step m->m-1 transfers epsilon to companions. Its probability per unit path is alpha*m, where the same previously calibrated alpha remains fixed. The proportional rate is an assumption, not something derived from gravity or time. More energetic photons take more steps and therefore lose the same fraction of energy on average.

**Known linear-death-process mathematics, conditional on that postulate:**

\[
\frac{dP_m}{ds}=\alpha[(m+1)P_{m+1}-mP_m],\qquad
m_{\rm out}\sim\mathrm{Binomial}(m_0,q),\quad q=e^{-\alpha D}.
\]

Consequently,

\[
\langle E_{\rm out}\rangle=E_0q,\qquad
\mathrm{Var}(E_{\rm out})=\epsilon E_0q(1-q),\qquad
\frac{\mathrm{Var}(E_{\rm out})}{\langle E_{\rm out}\rangle^2}
=\frac{\epsilon}{\langle E_{\rm out}\rangle}(1-q).
\]

The mean frequency falls by q, giving the desired mean relation S=1/q=exp(alpha D). This is an achromatic mean relation with frequency-dependent stochastic width, not a proof that every measured line centroid follows that equation exactly after flux weighting and instrumental reduction.

Each history balances energy exactly: m*epsilon remains in the photon and (m0-m)*epsilon has been transferred. The destination is a labeled companion reservoir; no gravitational field, capture dynamics or momentum exchange has been derived. The binomial distribution is over photon-energy states, not the number of photons in the original beam. At m=0 the photon is exhausted: P0=(1-q)^m0. For huge m0 and these transfer fractions this boundary probability is negligible, but the model is not exactly photon-number-preserving for arbitrary lengths. A physical continuous-energy and low-energy boundary completion is still required.

## Real radio lines constrain the step size

The previously archived methanol measurement uses a 48.3724558-GHz rest-frequency line, absorber redshift 0.885801 and fitted Gaussian FWHM 13.85 km/s. The underlying study cautions about different line profiles and uses similar-frequency transitions for its robust comparison. We reuse the already inspected summary; this is not a new spectrum or independent rate calibration. [Kanekar et al.](https://arxiv.org/abs/1412.7757)

For this diagnostic, attribute S=1.885801 to the transfer and allocate a stated part of the observed variance to conversion. No distance is inferred from redshift. In the large-energy-ladder limit, the conditional Gaussian-equivalent conversion width is

\[
\sigma_v^2\simeq c^2\frac{\epsilon/h}{\nu_{\rm obs}}(1-1/S).
\]

Known Gaussian width conversion gives FWHM=2*sqrt(2 ln2)*sigma_v. Giving conversion the entire observed variance yields epsilon/h <= **21.0 Hz**, or epsilon <= **8.69e-14 eV**. Allocating 10% of the variance gives 2.10 Hz and 8.69e-15 eV instead. These are illustrative noise-budget ceilings, not posterior confidence bounds; gas motions, covering factors, instrument effects and covariance have not been fitted. The astrophysical application uses the large-m limit; rounding energy to an integer ladder would be far below the resolved width here, but that does not derive a physical quantization rule.

| Assumed step energy divided by h | Predicted added FWHM, lower-frequency line | Higher-frequency line |
|---|---:|---:|
| 1 Hz | 3.02 km/s | 2.70 km/s |
| 10 Hz | 9.55 km/s | 8.54 km/s |
| 100 Hz | 30.21 km/s | 27.01 km/s |
| 1000 Hz | 95.53 km/s | 85.40 km/s |

Large steps conflict with the narrow-line allowance in this particular model. Smaller steps are not excluded by this summary alone. The predicted difference in widths is a future raw-profile test; the published tied-width summary is not two independent width measurements. A full line likelihood could change these allowances.

The Hz unit here is **energy divided by Planck's constant**. If every transfer produced one massless graviton, it would be that graviton's frequency, but neither this identity nor its numerical scale establishes a coherent gravitational wave detectable by LIGO. At the 10-Hz illustrative step, the radio photon undergoes about 2.27 billion mean transfers over the specified redshift factor. The physics must explain those events and their coherence, not merely assign their energy.

## What the same candidate predicts elsewhere

At 100 million light-years, the fixed alpha gives the unchanged target z=0.0076605. With the 10-Hz step, a 48.37-GHz input has about 0.888 km/s of added width, while a 6e14-Hz optical input has about 0.00798 km/s. These are synthetic predictions, not local-observation passes. The process predicts greater fractional broadening for lower-frequency photons even though their mean redshift is the same.

With separately assumed straight, unchanged-speed transport and ordinary clocks, event separation remains 1 times the emitted separation. No additional photon/companion delay occurs if companions are emitted along the ray and also travel at c; these kinematics are assumptions and do not establish a nonzero physical emission amplitude or momentum-complete interaction. The energy flux factor is q rather than q^2 because this implementation does not stretch arrival times. It therefore remains incomplete as an explanation of supernova timing and the shared-transport brightness relation.

This fixed-energy step is different from the previously tested fixed-*fractional*-energy jumps. Their noise formulas and bounds cannot be interchanged. Neither bound applies automatically to deterministic coherent conversion.

## Verification and next scientific requirement

Nine finite-state matrix evolutions agree with exact binomial probabilities within 1.5e-15 and with their moments within 9e-14. Probability normalization and per-history energy bookkeeping pass. These verify the statistical model, not optical phase coherence, field stress, graviton identity or observations across the spectrum.

A useful candidate cause must derive the rate alpha*m, the transfer energy or its distribution, angular/momentum behavior and clock/event response together. Adding an unrelated timing multiplier would not complete it. This branch can reproduce the mean redshift and constrain spectral noise, but cannot yet pass the complete stage-1/3 requirements. Independent physical void inference and frozen full-model holdouts remain outstanding; all four goal stages remain incomplete.
