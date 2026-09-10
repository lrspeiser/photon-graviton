# The lens source fraction does not transfer unchanged to disk galaxies

Applying the lens-calibrated source fraction f=0.6 universally fails a necessary outer-gravity bound for **147 of 149** archived SPARC disk galaxies under the existing baryonic baseline. Even after increasing stellar masses by the common lens-benchmark factor 1.645, **134 of 149** fail.

This is a strong conditional failure of the **constant source fraction plus spherical Newtonian-response prescription**. It is not a rejection of every possible companion theory, an assessment of the total photon-energy supply, or a fresh blind test of SPARC. It shows that improving central lens observables did not establish a model that explains outer disk gravity.

## The bound in plain language

A spherical positive-mass source cannot pull a star harder than it would if all of its mass were inside the star's orbit. Moving the wave around or changing its particle mass cannot exceed that limit while total source mass and the gravitational response stay fixed.

The **known Newtonian shell relation** gives

\[
g_c(R)=\frac{GM_c(<R)}{R^2}\le\frac{GM_{c,\mathrm{total}}}{R^2}.
\]

The lens pilot assumes M_c,total=f M_star. Combining the source with the archived ordinary-matter radial field gives the most generous possible speed bound

\[
v_{\max}^2(R)=v_b^2(R)+\frac{GfM_{\star,\max}}R.
\]

Conversely a **conditional necessary source-fraction bound**, derived from that known relation, is

\[
f_{\min}(R)=\frac{R\,[v_{\mathrm{obs}}^2(R)-v_b^2(R)]_+}{G M_{\star,\max}}.
\]

The positive part [ ]_+ sets a negative extra-force requirement to zero; such a point passes this necessary mass bound without establishing a successful model. If f_min>0.6, no radial rearrangement of a spherical source with the stipulated total mass can produce the measured central-value speed. A finite supported wave can provide less force than the bound, so passing it is not sufficient for a wave solution.

The **ordinary galaxy remains a disk plus any bulge and gas in the archived mass model**. Only the added companion source is assumed spherical here. We do not use a spherical enclosed-mass formula for the disk's own gravity.

## Inputs and generous mass allowance

We use the same 149 eligible galaxies and 3,150 radial rows as the prior SPARC audit: quality flag <=2, inclination >=30 degrees, positive disk scale and at least five valid rotation points. All systems and their historical roles are retained. They were already analyzed under earlier empirical models, so this is a transfer check on exposed data, not a newly blind 149-system validation.

The archived ordinary squared-speed decomposition is

\[
v_b^2=v_{\rm gas}|v_{\rm gas}|+0.5v_{\rm disk}|v_{\rm disk}|+0.7v_{\rm bulge}|v_{\rm bulge}|.
\]

Signed contributions are retained, including any outward-directed component. The disk and bulge coefficients are adopted mass-to-light ratios, not new laws. The gas contribution is included as provided. The data and mass-model context are described by [Lelli, McGaugh and Schombert, SPARC (2016)](https://arxiv.org/abs/1606.09251).

Rather than infer a bulge/disk luminosity split solely to obtain the source total, we deliberately allow

\[
M_{\star,\max}=0.7\,L_{3.6}
\]

in the corresponding solar units. Applying the larger of the adopted stellar mass-to-light ratios to **all** the light bounds the actual 0.5-disk/0.7-bulge mass from above under those assumptions. It therefore gives the source more mass than a disk-dominated galaxy would normally receive at f=0.6. The catalog luminosities are in 10^9 solar luminosities; that factor is restored explicitly.

The second scenario multiplies both the stellar force contributions and this stellar-mass allowance by lambda=1.6453621586, while keeping gas unchanged. This is the previously calibrated lens mass-correction sensitivity, not an independently validated SPARC population correction. Lens and SPARC population calibrations are not identical; a unified population treatment remains necessary. Both scenarios are conditional on their declared mass conventions.

## Results at the outermost valid measured radius

The outermost point is used for the headline count so galaxies with many radial measurements do not receive more weight.

| Ordinary stellar-mass scenario | Galaxies | Fail f=0.6 bound | Still fail after lowering observed speed by one quoted speed error | Median lower bound on f |
|---|---:|---:|---:|---:|
| Archived 0.5 disk / 0.7 bulge baseline | 149 | **147** | 143 | **5.58** |
| Stellar masses increased by 1.645 | 149 | **134** | 126 | **3.04** |

These are necessary lower bounds on total source fraction, not fitted halo masses or measured conversion efficiencies. The actual stellar mass can be below our generous allowance, and an extended source has some mass outside the observed radius, so a successful spherical model may require still more total source mass. No source fraction is changed in this test.

Lowering the speed by one quoted error is only a speed-error sensitivity. It does not propagate correlated inclination, distance, luminosity, gas, population or noncircular-motion uncertainty and is not a statistical rejection probability.

The problem is especially severe in faint galaxies:

| Total 3.6-micron luminosity | Count | Baseline failures | Median minimum f, baseline | Failures after stellar-mass increase |
|---|---:|---:|---:|---:|
| Below 10^9 solar luminosities | 37 | 37 | 17.48 | 37 |
| 10^9 to below 10^10 | 44 | 44 | 6.99 | 43 |
| At least 10^10 | 68 | 66 | 2.67 | 54 |

This trend is a diagnostic for a future formation model, not a fitted new luminosity law. Inferring a different f for each galaxy from its required gravity would describe the observations without explaining the deposits.

Forty galaxies have catalog distance methods identified as Cepheid or tip-of-the-red-giant-branch distances. All 40 fail under the original baseline and 38 under the increased-stellar-mass scenario. Thus the main failure is not confined to galaxies assigned Hubble-flow distances. Other published distances are retained as allowed fictional-universe inputs with their method flags visible; their original assumptions are not silently removed.

## Verification and limits

The script reconstructs every archived baryonic rotation value from the original component table; the largest discrepancy over 3,150 points is 5.7e-14 km/s. An independent units check agrees with Astropy's gravitational constant in kpc, km/s and solar masses. The necessary-fraction calculation agrees with direct comparison against the maximum-speed bound to roundoff, and every galaxy is retained in both scenarios.

The local catalog is read as whitespace-separated records in the documented column order. Its stored spacing does not conform to the stated fixed-width byte layout, so a strict Astropy MRT reader rejected it; no catalog value was altered to force that reader to work. Records are checked against all 175 names/metadata rows and the exact previously used 149-system rotation cohort. The file hashes and all 298 galaxy/scenario records are saved.

The bound assumes a spherical nonnegative source and the same Newtonian gravitational response used in the nonrelativistic wave pilot. It is not valid as a universal force bound for arbitrary nonspherical deposits, modified metric couplings, substantial relativistic stresses or sources outside this isolated spherical description. Those alternatives need their own dynamical and lensing equations. The finite wave particle mass m does not enter this generous bound; changing m alone cannot fix a failure while f and these assumptions stay fixed.

This is **not the deferred photon-supply budget**. We did not integrate stellar emission, choose an age or calculate conversion efficiency. We tested whether a source normalization already adopted for the lens fit can satisfy a necessary gravitational requirement elsewhere.

## Consequence for the overall concept

The calibrated f=0.6 cannot presently be promoted to a universal accumulation rule. It was an empirical choice that helped the lens sample, not a first-principles prediction. Under the same spherical response, a successful shared theory needs a physically derived dependence of accumulation on environment or history that can explain much larger source-to-stellar ratios in disk galaxies, particularly faint ones. Alternatively, changing gravity or source geometry changes the theory and must be derived and tested consistently.

This is a bigger priority than further fine-tuning lens orbital preferences. We need to connect photon transfer, companion travel, capture, storage/support and the resulting field normalization to one another. Otherwise a locally useful lens fit and a separately useful rotation fit remain unrelated descriptions. The redshift and event-timing mechanism, conservation, source stability, broad observations and fresh tests after development are still open requirements.

Reproduce with `python research_work/results/wave-disk-transfer-bound/run.py` and `verify.py`. The runner requires NumPy; verification additionally uses Astropy. `galaxy-bounds.json` contains every bound, speed, distance-method flag and historical split; `results.json` contains grouped results. No new parameter is fitted and no original data or calibrated coefficient is overwritten.
