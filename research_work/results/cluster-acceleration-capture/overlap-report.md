# Capture when two wells overlap

## Result and limitation

The acceleration-squared candidate is not additive across receivers. Two identical Plummer wells separated by 2a have a combined weak-capture area 1.660 times the sum of their isolated weak-capture areas. When exactly superposed, the factor is 2: total acceleration doubles, so local opacity quadruples instead of merely doubling. This is a static optically thin comparison, not a collision simulation or arbitrary-strength absorption result.

At the midpoint between separated identical wells, however, the gravitational accelerations cancel and the candidate predicts zero local capture. The potential remains deep. This exposes the difference between a rule based on net acceleration and the user's broader idea of capture by deep gravitational wells. The present candidate must not be described as a well-depth rule.

| Center separation / a | Combined weak area / sum of isolated weak areas |
|---|---:|
| 0 | 2.000 |
| 0.5 | 1.963 |
| 1 | 1.871 |
| 2 | 1.660 |
| 5 | 1.325 |
| 10 | 1.168 |
| 20 | 1.085 |

## Established field identity and candidate consequence

For kappa=chi|g|^2 and optical depth much less than one along relevant rays, effective area is approximately integral kappa dV. Thus

    sigma_pair=chi integral |g1+g2|^2 dV
              =sigma1+sigma2+2 chi integral g1.g2 dV.

The cross term can be locally negative where pulls oppose one another, but its full integral is positive for these two positive-mass wells. Known integration by parts and Poisson's equation give

    integral g1.g2 dV=-4 pi G integral rho1 Phi2 dV.

The right-hand side is positive for Phi2 referenced to zero at infinity. This is a Newtonian mathematical identity, not a claim that negative gravitational binding energy is an additional source of capturable photon energy.

For one Plummer well,

    integral |g|^2 dV=3 pi^2 G^2 M^2/(4a).

In units G=M=a=1, the cross integral is evaluated by averaging the second Plummer potential over angular directions around the first:

    <1/sqrt(1+|r-D|^2)> =
    2/[sqrt(1+(r+D)^2)+sqrt(1+(r-D)^2)].

That average is integrated against the first well's known density. Independent direct volume integration of g1.g2 agrees within 1.40e-14 relative for seven separations; doubling angular nodes changes results by at most 1.25e-11. The co-located limit gives exactly twice the sum of isolated areas. At large separation the cross contribution falls as 1/D, so interaction becomes small but is not identically absent.

No chi value or incoming energy density is fitted. The calculation supplies the coefficient of area in the weak-absorption limit. At finite optical depth, full attenuation saturates and requires ray integration through the combined field; the tabulated factors must not be extrapolated to opaque receivers.

## Consequences for the broader model

The earlier beta=sum(n_i sigma_i) competition approximation presumes independently characterized receivers. A nonlinear capture law based on the summed field can require overlap corrections or a continuous opacity calculation. Simply adding isolated collecting areas during a merger is inconsistent with this candidate.

The local cancellation also provides a useful discriminator between possible capture rules. A potential-depth or tidal/curvature-based rule can behave differently from net acceleration; each requires a precise physical definition and a separate test. A uniform external acceleration raises additional frame-dependence issues for an acceleration-based rule, reinforcing that this is not yet a generally covariant microscopic interaction.

Next assess a physically defined gravity observable and retain combined-field transport when comparing separated and interacting systems. Source motion, deposited support, field backreaction, finite supply and observed collision/lensing structure remain unmodeled in this static test. All six objectives remain open. No final observational holdouts were opened.

Formula status: field superposition, Poisson's equation and the integral identities are established mathematics. The acceleration-squared opacity remains a provisional hypothesis; no novelty or observational confirmation is claimed.

Run `python research_work/results/cluster-acceleration-capture/overlap.py`.
