# Companion capture scale and stored density released in the matched comparison

13 September 2026. Six fitted local parameters per galaxy, matching the free-NFW diagnostic. This is an inverse profile fit, not a physical supply prediction.

## Outcome

With capture scale and stored-density normalization fitted separately for each galaxy, the companion-shaped model gives total motion chi-squared 26.49, compared with 49.71 for the six-parameter NFW fit and 111.27-113.56 for the transferred companion law with four stellar parameters. All versions use the same forty motion bins, six systems, lens calibration and stellar model freedoms at their stated parameter counts.

The lower descriptive residual demonstrates flexibility of this profile family. It does not establish the exact-third law, a radiation source, sufficient supply, a positive orbital distribution function or stability. Free normalization absorbs the one-third factor. The required parameter shifts are substantial: two compact fits need density normalizations roughly 500 and 2400 times the reference, while J1621+3931 reaches the maximum capture radius and implies about 8.42e16 solar masses of effective stored mass within the numerical range.

Thus the earlier fixed companion law is not a demonstrated solution even though a freely adjusted companion-like density can fit the motions. No parameter set from this inverse diagnostic replaces the shared reference.

## Formula and provenance

Keep the reference opacity coefficient k0 fixed, and recompute the angular incident field at each trial capture scale ac:

    kappa(r) = k0/[1+(r/ac)^2]^2,
    J(r;ac,k0) = one-half integral over incoming directions of exp(-optical depth),
    rho_d(r) = D J(r;ac,k0)/[1+(r/ac)^2]^2.

Attenuation, angular integration and the gravitational integrals are known mathematics. The opacity shape and its deposited-density interpretation remain project hypotheses. Changing ac changes both geometry and optical depth; this experiment does not merely stretch a precomputed profile while leaving attenuation fixed.

The reference D is A*eta(X), with corrected A=2 C0 and the existing training-only amplitude adjustment. Here D is freed by fitting the deposited component's fraction f of required bending at the catalogue lens radius. Unit-density deflection d_unit(ac) gives D=f*alpha_required/d_unit(ac), while the remaining fraction calibrates stellar mass. These are lens-calibration identities, not a derivation of deposited energy from photons.

Because D is free, changing eta(X) can be canceled by changing the source normalization. The fit consequently cannot measure or validate the one-third exponent. Both historical population-proxy labels now share the same physical fitting problem; their reference values are retained only to quantify departures from the old law.

## Matched freedoms and data use

Fit capture scale ac/Re in [0.1,100] and bending fraction f in [0,0.95], plus the same h, beta0, beta_infinity and ra/Re stellar parameters as the free-NFW calculation. This is six local parameters per galaxy, 36 across forty motion bins. Both models consume catalogue lens angles to determine stellar normalization. All motion bins, including outer bins, are fitted.

The local parameter counts match, but physical assumptions differ. The companion family retains a fixed galaxy-trained k0; the NFW family has its own stipulated density shape. This is not a formal complexity-adjusted or blind predictive comparison. Necessary central and sampled orbital conditions remain enforced, with the same caveat that passing them does not prove a nonnegative distribution function or stable equilibrium.

## Fitting outcomes

| Galaxy | Free NFW chi-squared | Free companion chi-squared | Conditional outer residual | ac/Re | Deposited bending fraction |
|---|---:|---:|---:|---:|---:|
| J0037-0942 | 9.265 | 4.615 | 0.041 | 1.9626 | 0.1391 |
| J1112+0826 | 17.641 | 9.060 | 1.095 | 2.5315 | 0.4803 |
| J1204+0358 | 13.945 | 8.169 | 1.957 | 0.1000 | 0.5774 |
| J1402+6321 | 4.494 | 0.420 | -0.254 | 2.0778 | 0.4079 |
| J1621+3931 | 2.820 | 2.762 | -0.153 | 100.0000 | 0.0847 |
| J1630+4520 | 1.548 | 1.468 | 0.228 | 0.1000 | 0.2617 |

The total splits into inner chi-squared 21.32 and conditional outer residual-square sum 5.17. The improved outer residuals are not predictions because they were used in fitting. J1402+6321 can now be described with a small motion residual, but that outcome relies on free local density and scale rather than a successfully transferred source law.

J1204+0358 and J1630+4520 reach ac/Re=0.1; J1621+3931 reaches ac/Re=100. Stellar-gradient and orbital bounds also remain active. A boundary fit is not an independently measured physical scale, and local optimization does not establish a global optimum.

## Departures that the physical mechanism would need to explain

Reference ratios below use the Chabrier proxy. The complete output includes both proxy baselines; neither is selected as preferred evidence.

| Galaxy | Fitted ac / reference ac | Fitted D / reference D | Stored effective mass inside numerical grid (solar masses) |
|---|---:|---:|---:|
| J0037-0942 | 1.189 | 1.6932 | 4.394e+12 |
| J1112+0826 | 1.533 | 4.9028 | 7.6224e+12 |
| J1204+0358 | 0.06057 | 2444.6 | 9.3023e+10 |
| J1402+6321 | 1.259 | 5.5381 | 1.1092e+13 |
| J1621+3931 | 60.57 | 1.2309 | 8.4243e+16 |
| J1630+4520 | 0.06057 | 513.9 | 1.0935e+11 |

The two compact configurations need density normalizations about 2352-2445 and 493-514 times their population-dependent references. These ratios are not total-mass multipliers: changing the capture scale also changes the spatial integral and attenuation. Their fitted stored masses are about 9.30e10 and 1.09e11 solar masses, respectively.

J1621+3931 is a different degeneracy: its density normalization remains near the reference, but its capture scale grows by about 60.6 times to roughly 877 kpc. The implied stored mass inside the inherited radial grid is about 8.42e16 solar masses. The motion fit alone cannot justify assigning such an extended reservoir to this galaxy. It creates a major source-inventory and large-radius gravitational/lensing requirement that must be independently tested.

## Finite inventory and energy bookkeeping

The mass inside the numerical grid is a computed integral, not an observed inventory. The grid reaches at least 490 capture scales for these fits, and much farther for the compact cases. Since 0<=J<=1, the remaining positive mass beyond radius R is bounded by

    M_tail <= 4*pi*D*ac^4/R.

This follows by bounding the density by D*ac^4/r^4 and integrating 4*pi*r^2*rho. It is known integration applied to the proposed profile, not a new capture law. It controls the exterior tail only; numerical uncertainty of the interior integral and physical validity of the profile remain separate.

For J1621+3931 the tail bound is 1.30% of the computed interior mass, giving an approximate interior-plus-tail upper inventory of 8.53e16 solar masses under the adopted profile. Under the reference closure E_d=M_d*c^2, the interior inventory corresponds to about 1.51e64 joules. This is an energy-equivalent requirement of that closure, not a claim that gravitational field energy has a unique local density. No universe age or size is assumed, and no source history supplying it has been demonstrated.

Other fitted inventories and their energy equivalents are recorded in free-companion-inventory-results.json. A low chi-squared does not supply those energies or establish their retention.

## Numerical method and checks

The initial search tabulates unit-density enclosed mass and bending at 257 logarithmic capture scales and interpolates them logarithmically. The best three successful candidates are polished using direct profile evaluation. There are 22 starts per galaxy: both reference companion cases, the stellar-only fit, the free-NFW optimum as a seed, and eighteen grid starts. Twenty-one or twenty-two converge, with two or three successful direct polishing runs. Final predictions and scores use direct profiles.

Reference companion scores are reproduced to within 1.09e-7; the zero-density control also passes its reproduction check. The final lens equation closes within 2.23e-16 fractionally and the covariance decomposition within 2.67e-15. These algebraic closure precisions are not the physical accuracy of the model. The initial interpolation differs at final parameters by up to 2.90e-6 in bending and 3.38e-4 in unit enclosed mass; final direct evaluations remove reliance on those interpolation errors.

Adaptive integration issued a roundoff warning during tabulation at the requested 1e-9 relative tolerance. A separate refinement calculation therefore checked the fitted profiles directly. Increasing the incoming-direction quadrature from 96 to 192 and 384 changes fitted-parameter motion predictions by at most 1.53e-7 km/s. Independent fixed deflection quadrature over unsplit angles differed by up to 1.13e-6 for the most extended case; splitting intervals at capture-scale crossings reduced all comparisons with adaptive quadrature to within 6.87e-9. These checks support the reported score precision, not every requested integration tolerance. The stellar projection/radial grid is inherited, not independently refined in this calculation.

All final central and refined 8193-point slope constraints pass. This is still a necessary-condition check, not distribution-function inversion or stability. The broad numerical radial extent also does not constitute observed coverage at those radii.

## Assessment

At matched local fitting freedom, the companion density family is capable of lower descriptive residuals than the bounded NFW fit tested here. That is more informative than comparing a frozen companion profile with freely fitted halos. It is not a derived or independently confirmed solution: the free parameters expose large and diverse departures from the shared galaxy law and, in one case, an enormous extended reservoir.

The useful next distinction is between profiles supported by plausible supply and independently allowed large-radius gravity, and profiles selected only because they fit a small inner/lens dataset. The model must predict these density and scale requirements from sources, transport and retention, or explain why another physically specified branch applies. The exact-third reference remains unchanged, and no success is claimed for redshift, event timing, stable capture or energy supply on the basis of this stellar fit.

## Reproduction

Run free-companion.py, free-companion-refinement.py and free-companion-inventory.py in research_work/results/companion-extensions/. The outputs preserve all historical fits and record parameters, every motion bin, reference ratios, optimizer outcomes, numerical checks and input hashes. The pre-execution scope and choices are in free-companion-protocol.md.
