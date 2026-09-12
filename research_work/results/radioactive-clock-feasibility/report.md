# Can radioactive decay provide an independent supernova clock?

Radioactive-tail measurements are a possible source-physics constraint for the joint-light test, not an established solution. They offer a route beyond simply assuming distant explosions have the same fitted width as nearby ones. However, observed luminosity is not automatically equal to the radioactive heating rate.

## What the literature supports

[Jeffery (1999)](https://arxiv.org/abs/astro-ph/9907015) develops approximate radioactive energy-deposition models for late supernova light curves, including time-varying gamma-ray escape. [Dimitriadis et al. (2017)](https://arxiv.org/abs/1701.07267) find that the late optical/near-infrared light curve of SN 2011fe requires additional energy escape or redistribution beyond a simple full-trapping model and includes multiple radioactive chains. These are source-physics results, not expanding-universe assumptions we must import.

The possibility requires an explicit hypothesis that source-frame nuclear lifetimes retain their laboratory meaning in the emitting environment. Our companion propagation postulates do not by themselves establish that condition.

## Conditional source model

**Known radioactive-decay law:** a single dominant isotope supplies heating proportional to exp(-t_e/tau), where tau is its source-frame mean life, not its half-life. The calculation uses 111.3 days as an illustrative Co-56-like benchmark; it is not a new precision measurement of a nuclear constant.

**Illustrative deposition approximation from the established class of radioactive-tail models:**

\[
d(t_e)=f_p+(1-f_p)\left[1-e^{-(t_0/t_e)^2}\right],\qquad
L(t_e)=L_0e^{-t_e/\tau}d(t_e).
\]

Here t_0 controls gamma-ray escape and f_p is a retained charged-particle energy fraction. The numerical examples choose t_0=35 or 70 days and f_p=0.035; these are declared examples, not inferred values for the DES sample. The ejecta move outward within a static-universe model. No cosmological expansion is used in this calculation.

**Our conditional propagation mapping and known radiometry:** let A=S^b be the event stretch. For constant propagation conditions during the event,

\[
F(t_o)=\frac{L_0}{4\pi D^2 S^{1+b}}e^{-t_o/(A\tau)}d(t_o/A).
\]

The flux normalization disappears from the logarithmic decline rate:

\[
k_o=-\frac{d\ln F}{dt_o}
=\frac1A\left[\frac1\tau-\frac{d\ln d}{dt_e}\right].
\]

These equations are conditional deductions combining known decay/deposition mathematics with our proposed event mapping. They are not a new microscopic explanation for redshift.

## Why the simple clock reading can be wrong

If deposition decreases, the supernova fades faster than the radioactive clock alone. Ignoring escape would estimate A_naive=1/(tau*k_o), which is less than the actual A. Across the 36 declared numerical cases, A_naive/A ranges from 0.386 to 0.912. Thus even perfect decline-rate measurements can give a badly biased stretch estimate if the emitting source is modeled incorrectly.

The analytic slope agrees with finite differences to relative error below 4.2e-11. This checks the calculation, not the suitability of the deposition approximation for an actual event.

With a single isotope, prompt bolometric reradiation and nonincreasing deposition, one obtains the conditional lower bound

\[
A\geq\frac{1}{\tau k_o}.
\]

It is not a general supernova bound. Additional isotopes, delayed release of stored energy, light echoes, changing observed-band fractions or increasing deposition can invalidate its assumptions. In a single optical filter, changes of spectral energy distribution add another logarithmic derivative; the bolometric formula cannot simply be applied to that band.

## Route to a useful test

Use observer-time multi-band late light curves and spectra to constrain deposition, spectral redistribution and extra power sources. Retain the explosion-epoch uncertainty and fit the same source/transport model jointly, rather than declaring an optical decline rate to be a nuclear lifetime. Nearby sources can test the source model; a separated distant sample is needed to test cumulative propagation. If the permitted source model leaves arbitrary time-dependent deposition or luminosity, the earlier identifiability problem returns.

This identifies a physically motivated constraint worth developing; it does not yet supply an independent calibrated clock for our dataset. It also avoids requiring an ultimate explanation of time before testing an effective model. The existing timing-calibration batch is unchanged, and none of the six scientific demonstrations is complete.

```powershell
python research_work/results/radioactive-clock-feasibility/run.py
```
