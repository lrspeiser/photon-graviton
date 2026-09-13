# M87 finite-emitter calculation

## Result and significance

Replacing the central point emitter with an observationally informed extended light distribution changes M87's total conditional cluster deposition from 33,687,363 to 33,684,549 solar luminosities: a decrease of 0.00835%. Source extent resolves the point-source singularity but does not materially increase this emitter's energy supply. This is a sensitivity calculation, not an observed deposit, accumulated reservoir, gravitational field or lensing fit.

The spherical approximation gives a finite central volumetric deposition power of 9.21662e9 Lsun/Mpc^3. Independent analytic and numerical convolution values differ by 6.22e-15 relative. This numerical agreement is not an observational uncertainty. The full projected cluster map has NOT been replaced: other emitters remain unresolved, and this calculation proves boundedness rather than measuring a lensing peak.

## Observational input versus assumptions

[Simon, Cappellari and Hartke, Table 2](https://arxiv.org/html/2303.18229v2#S4.T2), DOI 10.1093/mnras/stad3309, supplies eleven Gaussian components of deconvolved F850LP stellar surface brightness. Their photometry combines Hubble and SDSS imaging and removes AGN/jet contamination. Only the photometric table is imported; their dark halo and dynamical mass estimates are not inputs.

The original tabulated values are preserved in profile.json. The catalog distance 16.74942788 Mpc and bolometric luminosity 1.504e11 Lsun come from the unchanged DustPedia THEMIS input. Input hashes are saved in results.json. The following are explicit modeling choices, not measurements:

- Circularize each elliptical Gaussian with sigma_c=sigma*sqrt(q), preserving its integrated light; assume a spherical three-dimensional Gaussian whose projection gives that circularized profile. This loses ellipticity and does not uniquely infer depth from an image.
- Assign bolometric luminosity in proportion to F850LP Gaussian light. This assumes no spatial color variation and treats the pilot luminosity as following the stellar shape; excluded AGN emission needs separate accounting if relevant to the bolometric input.
- Retain the stipulated 1 Mpc receiver, constant alpha=0.0002488993265191759/Mpc, kappa=10/Mpc, stationary illumination, lossless companions and permanent capture. No age, mass target or lensing profile is fitted.

The resulting circularized projected half-light radius is 5.691 kpc. Component scales come from the published light fit rather than tuning to improve the deposit map. Finite profiles for the other galaxies, deprojection uncertainty and wavelength-dependent emission remain outstanding.

## Equations and checks

All Gaussian projection, convolution, geometry and integration equations below are established mathematics. Applying them to photon-to-companion conversion is conditional on this project's postulated transport; none is claimed to be a unique fundamental law.

Weights are w_j proportional to I_j sigma_j^2 q_j, normalized to unity. The spherical volume luminosity profile is

    j(r) = L sum_j w_j exp[-r^2/(2 s_j^2)] / (2 pi s_j^2)^(3/2).

For an internal infinitesimal emitter and separation s within the uniform receiver, the previously derived power-density kernel is

    K(s) = kappa alpha [exp(-alpha s)-exp(-kappa s)]
           / [4 pi (kappa-alpha) s^2].

The extended source gives the convolution q(x)=integral j(y) K(|x-y|) d^3y. At the center the Gaussian integral yields

    q(0) = L kappa alpha/[4 pi (kappa-alpha)]
           * sum_j w_j/s_j^2
           * [erfcx(alpha s_j/sqrt(2))-erfcx(kappa s_j/sqrt(2))].

Here erfcx(x)=exp(x^2) erfc(x). Direct radial quadrature independently verifies this expression. The analytic center formula extends the uniform internal-source kernel over the Gaussian tails. Numerical source integration ends at 12 component sigmas, all inside the receiver; omitted luminosity fraction is 5.19e-31. Thus source light outside the receiver is negligible in this approximation, not an imposed physical cutoff in M87.

Because companion energy fraction C(s)<=alpha s, K(s)<=kappa alpha/(4 pi s). Convolving 1/s with a normalized spherical Gaussian gives erf(r/[sqrt(2)s_j])/r, bounded everywhere by sqrt(2/pi)/s_j. Consequently the extended emitter's density and its projection through the finite receiver are bounded. The resulting conservative projected bound, 1.87207e10 Lsun/Mpc^2, is NOT a predicted central value or observed surface density.

For total retained power, integrate the existing internal-source angular capture fraction over each Gaussian's Maxwell radial distribution. Independent 64/128 angular-node calculations agree at reported precision. Positive energy fractions and the center bound are asserted in run.py. No additional light or energy is introduced by extending the source.

## Next physical requirements

Evaluate the resolved projection and measured apertures, extend source profiles across relevant galaxies, and propagate source/color/deprojection uncertainty. Source histories and a supported reservoir with a specified gravitational response are still needed to predict both lensing and motion. This result does not establish cluster sufficiency or shortage. All six objectives remain open; no final holdouts were opened.

Run `python research_work/results/cluster-m87-finite-emitter/run.py`.
