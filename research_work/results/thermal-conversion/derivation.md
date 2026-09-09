# Photon-energy conversion and the microwave spectrum

## Why this test is needed

The redshift rule must say what happens to a bath of radiation, not only to individual starlight photons. A blackbody has a linked temperature, energy density and photon number. Lowering each photon's energy while leaving the number unchanged need not produce another blackbody. Without this check, we could claim the observed microwave spectrum is preserved by a rule that predicts the wrong normalization.

This pass distinguishes two mechanisms already present in the project: constant-speed partial photon-energy conversion, and a changing propagation field that changes the electromagnetic density of modes. It uses an initially thermal bath as an explicit conditional comparison, not a Big-Bang premise or a derived origin of the microwave background. A different source history, escape, thermalization or photon emission/absorption requires a separate calculation.

## 1. Fixed volume, fixed light speed, conserved photon number

Let N_nu be photon number density per frequency. Suppose every photon frequency decreases by the factor q, where 0<q<=1, while physical volume, standard light speed and total photon number remain fixed. Conservation under the frequency-variable transformation gives

    N_final(nu) = N_initial(nu/q)/q.

Multiplying by photon energy h nu gives

    u_final(nu) = u_initial(nu/q).

For an initial Planck spectrum at T_i, the standard frequency-cubed prefactor then yields

    u_final(nu) = q^-3 u_Planck(nu, q T_i).

The color temperature decreases to q T_i, but the normalization is greater than a standard blackbody at that temperature by q^-3. Integrating the spectrum confirms that photon number remains unchanged and total energy decreases by q. A standard blackbody at q T_i would instead have q^3 times as many photons and q^4 times as much energy in the same volume.

The code independently integrates these number and energy spectra at three q values. This is a homogeneous-bath calculation; it does not import an expansion dilution term or substitute luminosity-distance factors for a physical-volume balance.

## 2. How additional physics could change the result

One mathematical way to recover a Planck normalization is to retain only the fraction q^3 of photons, independently of frequency, after applying the energy scaling. The final energy fraction is then q^4. This is an extra photon-removal process, not a consequence of partial energy conversion.

For continuous energy drift dE/dt=-h E in a fixed volume, Planck-preserving evolution at temperature drift dlnT/dt=-h requires a photon-removal rate 3h per photon in this simple homogeneous comparison:

    dN/dt = -3h N
    dU/dt = -4h U
    energy into the partial-redshift channel per time = h U
    energy carried by removed photons per time = 3h U.

All removed energy needs a destination. The relative integrated channel energies are 1:3 for this continuously concurrent prescription; they are not generally the same as first shifting a finite photon ensemble and only afterward removing photons. A mechanism that removes photons also changes transmitted brightness and must pass the other observation requirements. No such removal rate is adopted here.

Other possibilities include escape, emission/replenishment, an explicitly calculated thermalizing receiver, a nonthermal initial spectrum, or a different mode structure. Each changes the prediction and must include its energy and photon ledger. Redshifting arbitrary radiation alone does not create a thermal spectrum. We have not derived the initial state or the origin of the observed microwave background.

## 3. Why the earlier propagation-field result is different

In the homogeneous nondispersive index model, omega=c k/n and the number of modes per physical volume per frequency scales as n^3 nu^2. Keeping wavevectors and mode occupations fixed while n changes preserves a thermal occupation law with n T constant. Its equilibrium densities scale as

    photon number proportional to n^3 T^3
    energy density proportional to n^3 T^4.

If q=n_initial/n_final, the same number of photons can therefore have final energy q U_initial while retaining the index model's thermal spectrum. That does not contradict section 1: its electromagnetic mode density changed, whereas section 1 holds the standard mode density fixed.

The earlier [field-ray](../../../fixed_atom_time_field/fixed_atom_time_field.md) and [matched-wave](../matched-wave/derivation.md) results retain this conditional distinction. They still need consistent atomic standards and a detector response. A modified mode density is not a free correction that can be attached to the constant-speed conversion model while leaving all its other laws unchanged.

## 4. Observation product and declared comparison

We use the recovered `time_revision/data/firas.txt`, all 43 numeric channels. [NASA LAMBDA's product description](https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html) states that its monopole column is reconstructed from a 2.725 K blackbody plus published residuals; it also provides residual uncertainties and a modeled Galaxy spectrum. We fit the residuals, not the rounded reconstructed intensity column.

The [protocol](protocol.json) was written before this pass's fit. The product was used in earlier research, so this is exploratory reuse, not new blind validation. Input and protocol hashes are recorded in the results.

For each fixed q, fit a final color temperature and a signed nuisance coefficient g multiplying the released Galaxy residual template. The predicted residual in kJy/sr is

    R_model(nu) = 1000 [q^-3 B_nu(T_color)-B_nu(2.725 K)] + g G_nu,

where B_nu is expressed in MJy/sr. Frequency in cm^-1 is converted to Hz using 100 c. The score sums squared residuals divided by the published diagonal variances. T_color is fitted independently: the initial thermal temperature is not fixed. The amplitude q^-3 is not allowed to float independently in the conversion hypotheses.

We also fit a descriptive positive free amplitude, temperature and the same Galaxy template. There is no supplied full covariance, calibration/gain prior, instrument model or full foreground likelihood. Therefore the scores below are diagnostic comparisons, not official FIRAS exclusion limits, confidence levels or inferred cosmic ages. In particular, a free calibration gain would be degenerate with a uniform spectral amplitude unless constrained independently.

## 5. Results

| Photon energy survival q | Required Planck amplitude q^-3 | Fitted color temperature, K | Diagonal chi-square |
|---|---:|---:|---:|
| 1 | 1 | 2.725001 | 45.02 |
| 0.9999 | 1.000300 | 2.724791 | 86.00 |
| 0.999 | 1.003006 | 2.722905 | 4339.94 |
| 0.99 | 1.030610 | 2.704052 | 438288.66 |
| 0.9 | 1.371742 | 2.515693 | 51713898.88 |
| 0.5 | 8 | 1.605381 | 3170943771.30 |

Each fixed-q case uses two fitted nuisance parameters. The free-amplitude comparison returns amplitude 1.0000073 and chi-square 44.99. This does not measure q independently of the stated history and calibration assumptions. At large losses, the fitted Galaxy-template coefficients become very large and should not be regarded as physically established foreground corrections; allowing them did not rescue those diagnostic fits.

Two different numerical optimizers agree on the profiled scores. Successful numerical convergence means the specified comparison was executed, not that the physical assumptions or likelihood are complete.

## What this adds to the theory requirements

A constant-speed, photon-number-preserving conversion law cannot claim generic Planck-spectrum preservation in a fixed volume. Under the specified initially thermal, unreplenished history it generates a normalization mismatch that is substantial in the recovered FIRAS residual comparison.

This is not a theorem excluding every conversion cosmology or every fictional background history. The missing ingredients are the actual microwave-frequency conversion rate, initial/source spectrum, escape and photon-number processes, any thermalization mechanism and the energy received by it, and a complete observational calibration model. Angular correlations and structure are not tested here. All original 20 tasks and 32 observational areas remain in scope.

Run `python -X utf8 research_work/results/thermal-conversion/check_thermal_conversion.py`. Saved evidence is in `thermal-conversion-results.json`. The original FIRAS file is preserved unchanged.
