# Can one light-bending response repair the stellar-motion model?

**Result:** under the retained six-galaxy spherical models and inner-stellar-motion mass posteriors, no common extra-sector lensing multiplier matches all six published lens-angle summaries within the displayed central mass ranges. A universal reduction to time-only bending also does not repair them. This is an exposed-data diagnostic, not a joint image likelihood or a general exclusion of the companion hypothesis.

## Shared force, separate measurable metric response

Keep the already fitted dynamical potential Phi=Phi_b+Phi_c unchanged. As a diagnostic weak-field completion, use spatial potential Psi=Phi_b+gamma_c Phi_c. Ordinary matter retains the equal-potential prescription. Slow stellar motion responds to Phi, while light bending responds to Phi+Psi. This is established weak-field metric mathematics; see [Will's review](https://arxiv.org/abs/1403.7377). Assigning an independently adjustable constant gamma_c to the extra component is an optional phenomenological postulate, not a derived deposited-energy field equation.

Relative to the previous equal-potential extra-force lens calculation,

    f=(1+gamma_c)/2
    reduced_deflection(theta)=B_b(theta)+f B_c(theta).

B_b and B_c include the existing source/lens distance ratio. f=1 is the old equal-potential completion; f=0.5 corresponds to an extra component affecting only the temporal potential. Altering f leaves the stipulated stellar model unchanged at this weak-field order. It is not permission to fit gravity and light independently in a completed theory: the diagnostic asks what one shared field completion would have to supply.

At a published angle theta_cat, the spherical lens equation demands

    f_required = [theta_cat-B_b(theta_cat)]/B_c(theta_cat).

This is known linear algebra applied to the two computed deflection components. It does not assert that the Einstein radius itself scales linearly with f. We evaluate the deflections at the catalog angle rather than incorrectly multiplying the previous predicted Einstein radius.

## Comparison with existing observations

Use all six previously exposed resolved-kinematics systems, the frozen empirical force law, light components, cutoff and conditional distance geometry. Stellar mass and orbit anisotropy were inferred only from inner velocity bins; neither the outer bin nor the lens angle was used in that inference. Here the published lens angles are explicitly used to calculate the required response, so these are inferred requirements, not blind predictions.

| Galaxy | Required f at median stellar mass | Range from central 95% stellar mass interval |
|---|---:|---:|
| J0037-0942 | 0.875 | 0.751-1.000 |
| J1112+0826 | 1.941 | 1.770-2.113 |
| J1204+0358 | 0.756 | 0.573-0.934 |
| J1402+6321 | 0.800 | 0.691-0.910 |
| J1621+3931 | 1.632 | 1.370-1.887 |
| J1630+4520 | 1.989 | 1.841-2.139 |

The intersection is empty: the largest lower endpoint is 1.841 while the smallest upper endpoint is 0.910. The low-response and high-response groups require opposite changes to the current f=1 completion. In particular f=0.5 falls below all these ranges.

These intervals propagate only the existing stellar mass posterior, under the uniform-log-mass and uniform-anisotropy prior, holding the lens summary exact for the diagnostic. They are **not full confidence intervals on f**, nor simultaneous confidence bands. Lens image/model uncertainty, nonsphericity, line-of-sight structure, spatial mass-to-light gradients, more general stellar orbits and geometry uncertainty are omitted. The catalog SIE angle is a fitted nonspherical model summary, not an exact observed spherical ring radius. No rejection significance follows from the empty intersection alone.

At median masses, the f=0.5 lens-equation residuals are negative in all six systems: the reduced deflection at the catalog angle is too small. This is a statement about that completion with those masses, not a direct measurement of a universal gravity parameter. The extra force is still an empirical template rather than a calculated companion deposition profile.

## Implication for the full theory

The answer is not to assign the six listed factors to the six galaxies. That would remove the predictive test. One physical field equation must explain both the dynamical potential and light-bending response, together with the ordinary-matter distribution and its uncertainties. A different mass structure can change both inferred motion and lensing; its stellar fit must be recalculated rather than kept artificially fixed.

The newly derived cosmological beam areas also remain separate from the conditional distance geometry used here. They have not yet been derived for these lenses, and the six-system diagnostic is not a test of a complete global lapse field. Existing source and distance assumptions are retained to isolate the response question.

**Decision:** a single added constant light-bending multiplier is not an adequate repair within the current approximations. Proceed toward a specified spatial source/metric model, rather than changing a global normalization or fitting one per object. response.py and response-results.json preserve deflection components, mass quantiles, required responses and fixed-completion residuals. No final holdout was opened. All six objectives remain open.
