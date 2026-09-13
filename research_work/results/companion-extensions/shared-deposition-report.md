# One shared deposition scaling: no acceptable transfer yet

13 September 2026. Six exposed fitted halos; 40 existing stellar-motion bins. This bounded experiment attempts to compress the required distributions into one population-based formula. It does not fit a microscopic capture rate or calculate new incoming companion flux.

## Proposed rule and provenance

    rho_d(r) = A / [4 pi rs^3 x (1+x)^2], x=r/rs
    A = K eta Mpop (W/Wstar)^p
    rs = b Re (W/Wstar)^q
    W = G Mpop/Re, Wstar=(300 km/s)^2

The density shape is the known NFW profile. W is a characteristic binding-scale proxy, not the actual resolved potential. Mpop is the archived population mass proxy, Re the half-light radius, and eta the original exact-third retention mapping. The power laws and density mathematics are known; their shared companion deposition interpretation is a project hypothesis. K,b,p,q are common to all galaxies. No independently adjustable galaxy halo amplitude or radius remains in the candidate rule. Archived population masses and eta have their existing calibration uncertainty; they are not newly measured luminosities.

Under stationary capture/loss bookkeeping this would specify the required product kappa*u_c*tau = c*rho_d. It does not determine those three factors separately. Permanent retention instead needs an integrated illumination history. Neither source accounting, physical outer truncation, stable support nor microscopic interactions are solved. The formal NFW tail has divergent total mass, so this is an effective profile diagnostic, not a completed finite-energy deposition law.

## Test design

Three nested families have two shared parameters (p=q=0), three (free p,q=0), or four (free p,q). Fit inverse enclosed-mass targets at the observed annular midpoint radii, using equal-galaxy squared asinh(Mhalo/Mpop) differences. This includes the zero-halo target without a log-zero workaround. These target masses are inferred from earlier fits, not direct mass observations. All bounds, starts and objective were declared in shared-deposition-protocol.md before fitting.

Each family is fitted to all six systems for description and then six times to five systems, predicting the omitted sixth. Replay motions and lensing with original halo-fit stellar masses, gradients, orbital parameters and geometry frozen. Do not recalibrate stellar mass to restore the lens fit. Thus a failed replay tests this frozen closure; a later joint stellar/orbital refit would be a different experiment. Previously exposed data and weakly constrained/boundary halo targets prevent blind-validation or general impossibility claims.

## Results

Lower motion chi-squared and lens fractional RMS are better. Lens RMS here is an unweighted fractional discrepancy, not a lens likelihood.

| Shared parameters | All-six motion chi-squared | Omitted-system motion chi-squared | Omitted-system lens RMS |
|---|---:|---:|---:|
| 2 | 9301.66 | 12667.19 | 41.08% |
| 3 | 6497.66 | 8820.79 | 41.36% |
| 4 | 6222.34 | 25235.99 | 65.85% |

The separate fitted halos give motion chi-squared 49.713 and lens closure by calibration. That flexible baseline is not an equally parameterized competitor, but establishes the predictions the shared rule was intended to recover. All-six lens RMS is 29.99%, 30.64% and 31.02% respectively. None of these candidates is acceptable as a replacement. The three-parameter family has lower omitted motion error, but still large discrepancies and slightly worse lens RMS than the two-parameter family. Do not select it as a successful theory.

Descriptive all-six coefficients:

| Family | K | b | p | q |
|---|---:|---:|---:|---:|
| 2 | 0.586024 | 0.211610 | 0 | 0 |
| 3 | 0.055139 | 0.142849 | 3.286341 | 0 |
| 4 | 0.710867 | 2.893030 | -0.071577 | -4.000000 |

The fourth parameter reaches its declared lower bound. Allowing it reduces descriptive motion error but worsens omitted-system transfer sharply. No expanded bounds or galaxy exclusions were tried. Omitted-fit coefficients differ from these descriptive coefficients and are fully saved in JSON.

## What we learned

Identical imposed halo density still reproduces halo gravity, as the preceding independent equivalence check established. This experiment fails at predicting which density to impose from a simple shared mass/size scaling. It does not identify a gravity-solver bug, nor prove that no universal capture law exists. The zero-halo case and widely varying, partly unconstrained target scales also make this a difficult and uncertain inverse target set.

Keep the original exact-third reference; do not adopt these coefficients. A next physically distinct candidate could calculate how incident companion supply and interception vary with environment, rather than assigning the same supply from population mass and size alone. First match ordinary-matter calibration and propagate target-profile uncertainty: a unique fitted halo shape must not be treated as a uniquely measured deposition requirement. Current results do not justify merely adding more exponents to this power law.

Run shared-deposition.py to reproduce all 21 multistart fits and the frozen observable replay. Results include per-system predictions, omitted fits, optimizer success counts, boundary flags and source hashes. No production gravity functions or prior results were changed. The v1.5 PDF predates this supplement.
