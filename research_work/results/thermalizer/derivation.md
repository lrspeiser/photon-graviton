# A heated thermal absorber/emitter during photon-energy conversion

## Why this calculation is needed

The preceding inverse-source calculation specified the emission needed to hold a chosen background unchanged. A real emitter cannot choose its spectrum freely. Here we replace that chosen source with thermal emission tied to absorption at a material temperature. This checks whether an ordinary thermalizing medium could maintain the background, and what energy and transparency costs would follow. Without it, saying that matter reprocesses light into microwaves leaves the central spectrum and heating questions unanswered.

This is a specified macroscopic comparison, not an adopted cosmic material or a first-principles photon-to-companion vertex. No expansion, Big Bang, dark matter or assumed thermal companion population enters it. Uniform material temperature, gray opacity and an external heat source are explicit assumptions. A real material's composition, density, opacity and heating source remain to be supplied.

## 1. Tie emission to absorption

Take constant fractional photon-frequency loss d nu/dt=-H_c nu and a homogeneous, isotropic, fixed-volume bath with standard fixed light speed. The material is in local thermal equilibrium at T_d. Its emission and effective absorption satisfy Kirchhoff's law, with source function B_nu(T_d). This is the standard thermal radiative-transfer relation reviewed in [NRAO's radiative-transfer lecture](https://www.aoc.nrao.edu/~smyers/courses/astro12/L10.html) and [Dullemond's continuum transfer treatment](https://researchmgt.monash.edu/ws/portalfiles/portal/311391053/311390078_oa.pdf).

For frequency-independent rate kappa, the number spectrum obeys

    partial_t N_nu = H_c partial_nu(nu N_nu)
                    + kappa [N_P(nu,T_d)-N_nu].

kappa is the effective thermal absorption/relaxation rate per time; stimulated emission, when relevant, is incorporated consistently into the usual net opacity. It is not a separately derived microscopic collision count. There is no conservative elastic scattering term or boundary escape in this comparison. A prescribed material thermostat is not a free energy source: its required heating is computed below.

Set a=kappa/H_c>0. With the high-frequency boundary contribution vanishing, the stationary characteristic solution is

    N_nu = a integral_0^infinity exp[(1-a)s] N_P(nu exp(s),T_d) ds.

For the energy spectrum or isotropic intensity this becomes

    I_nu = a integral_0^infinity exp(-a s) B_(nu exp(s))(T_d) ds
         = integral_0^infinity exp(-y) B_(nu exp(y/a))(T_d) dy.

This is a positive spectrum supplied by the thermal medium. It is a weighted combination of photons emitted at different times and shifted by different amounts, rather than an assumed blackbody. At a tending to infinity it approaches B_nu(T_d).

## 2. Exact photon and energy moments

Integrating the kinetic equation or the characteristic solution gives

    N_total = N_P,total(T_d),
    U_gamma = a/(a+1) U_P(T_d).

Thus the bath has the photon count of a blackbody at T_d but less energy. It cannot be exactly a standard blackbody at any temperature for finite a: equality of photon count would require that temperature to equal T_d, contradicting the energy moment. A free normalization is a different observational model, not standard Planck normalization.

For large a, on a bounded frequency interval,

    I_nu = B_nu + (nu/a) partial_nu B_nu + O(a^-2)
         = (1+3/a) B_nu[T_d(1-1/a)] + O(a^-2).

The leading distortion therefore contains both a temperature shift and an amplitude change. Instrumental gain uncertainty matters when comparing small distortions with observations. The expansion is not uniform out to arbitrarily large frequency; the numerical fits use the full integral.

## 3. Heating and receiving sectors

The material's gross thermal emission power density is kappa U_P(T_d), while its absorbed power density is kappa U_gamma. To maintain T_d it needs net external heating

    P_heat = kappa [U_P(T_d)-U_gamma] = H_c U_gamma.

The same power leaves the photon bath through conversion. At stationarity this is energy flowing from the heating source, through matter and light, into the receiving sector. Absorption/re-emission can recycle a much larger power internally; that recycled power is not a new net supply.

This is different from treating every removed photon as permanently lost to another reservoir in the preceding inverse-source example. If absorbed energy is re-emitted, it cannot also remain permanently deposited. If converted energy is captured into permanent gravitational deposits, the deposits continue to gain energy and the heating source must pay for that gain. Neither the thermal medium nor a source-independent perpetual steady universe has been derived.

## 4. Exploratory comparison with the recovered microwave spectrum

The protocol was written before this fit. It reuses the 43 already-exposed channels of `time_revision/data/firas.txt`, whose [NASA product documentation](https://lambda.gsfc.nasa.gov/product/cobe/firas_monopole_spect.html) explains the residual, uncertainty and Galaxy-template columns. We fit the residuals relative to the reference 2.725 K Planck spectrum, not the rounded reconstructed monopole column. For each fixed a, we fit T_d and a signed coefficient of the supplied Galaxy template.

Only diagonal published errors are used. The full covariance, instrument gain/calibration prior and foreground likelihood are absent. There is no fresh held-out data, confidence interval, official exclusion limit or inferred universal rate. A large fitted Galaxy coefficient does not establish a physical foreground.

| kappa/H_c | Fitted material T_d (K) | Diagonal chi-square |
| ---: | ---: | ---: |
| No conversion, Planck comparator | 2.72500054 | 45.019 |
| 10 | 2.79577691 | 39,705,581 |
| 100 | 2.73140502 | 429,993 |
| 1,000 | 2.72563128 | 4,332.026 |
| 10,000 | 2.72506351 | 85.994 |
| 100,000 | 2.72500684 | 45.239 |

Rapid thermalization approaches the comparator in this simplified comparison. It is not correct to report 100,000 as a measured lower bound: this is a preset coarse grid with incomplete observational uncertainties and strong material assumptions.

## 5. Erasing directional information is a separate cost

With a fixed homogeneous thermal source, an incoming directional perturbation is attenuated relative to conversion-only transport by exp(-kappa t). If its frequency survival from conversion is q=exp(-H_c t), this additional attenuation is q^a. Equivalently, for the conversion frequency ratio S=1/q, it is S^-a.

At large a the background photons are repeatedly replaced before appreciable conversion accumulates. An observed angular pattern must then be generated by the actual spatially varying thermal medium, heating or another specified process; preserving an arbitrary incoming pattern cannot be assumed. No angular correlation calculation has been performed here.

The gray assumption covers the modeled radiation band. It does not establish optical opacity. Extending the same large absorption rate to visible starlight would also erase distant beams, so any proposed material must demonstrate the required frequency dependence against both microwave and optical observations rather than extrapolating this fit silently.

## 6. Verification and next physical inputs

The script checks the characteristic integral against adaptive integration and doubled quadrature order, independently checks the differential transport equation, and integrates both moments. It then repeats the temperature/template minimization using two optimizers. These checks establish the stated calculation, not a successful universe model.

Next supply a material opacity and temperature distribution, a source-funded heating history, and the companion production/capture law. The same model must predict spectrum, directional structure and transparency while accounting for material, photon and deposited energy. Do not replace these physical inputs with the source spectrum reconstructed from the desired answer.
