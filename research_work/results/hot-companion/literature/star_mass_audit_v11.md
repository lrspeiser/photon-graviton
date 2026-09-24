# How published cluster stellar masses were computed: a literature audit

Written 2026-09-24 for the hot-companion project (round 11; results README section 21.2). The goal is to put every stellar mass the project uses on one basis.

**Labels.** **STATED** means the paper says it, with a pointer to where. **INFERRED** means my own arithmetic or reading. "Not stated" means I looked and the paper does not say it.

**Versions read.** I read the LaTeX source of the latest arXiv version of each paper, unless noted. Section, table and equation numbers follow that version. For X-COP I also read the published A&A PDF and its Corrigendum.

The papers were read from their arXiv sources (LaTeX); the X-COP stellar-mass release files are the ones linked from the X-COP data page.

---

## 0. Corrections to the brief

1. **arXiv:2012.03904 is not the X-COP stellar-mass paper.** It is McConnachie & Venn (2020), a Research Note on Gaia EDR3 proper motions of dwarf galaxies.
   - The X-COP paper is **Ghizzardi et al. 2021, arXiv:2007.01084**, published as **A&A 646, A92 (2021)**.
   - It has a **Corrigendum, A&A 652, C3 (2021)**. The Corrigendum only fixes a unit-conversion error in the *expected* iron yield (Maoz & Graur values). **The stellar masses are unchanged.**
   - The published Table 3 (stellar masses) is identical to arXiv v1.
2. **arXiv:1001.3129 is not Treu et al. 2010.** It is a mathematics paper, "Lasry-Lions regularization and a Lemma of Ilmanen". Treu et al. 2010, "The initial mass function of early-type galaxies", is **arXiv:0911.3392**.
3. **"Auger et al. 2010" is two papers.** The IMF result is in **arXiv:1007.2409** (Auger, Treu, Gavazzi et al. 2010b). The SLACS X correlation paper is arXiv:1007.2880.
4. **Sifón et al. 2013 (arXiv:1201.0991) has no stellar masses.** I searched the source for "stellar", "IMF" and "M/L" and found nothing. It is a dynamical-mass paper.
5. **The key later El Gordo stellar mass is Hilton et al. 2013, arXiv:1301.0780.** This is the "Hilton et al., in prep." that Menanteau et al. 2012 cite for the ACT IRAC stellar-mass programme.

---

## 1. X-COP stellar masses (Ghizzardi et al. 2021, arXiv:2007.01084)

Section numbers below are those of the published paper: Sect. 4 "Optical data" and Sect. 4.1 "Stellar mass profiles". The arXiv v1 text is the same.

### STATED

**Sample**
- Seven X-COP clusters are also in MENeaCS: A1795, A85, A644, A2319, ZW1215, A2029 and A2142 (Table 3). The optical analysis uses the data set of van der Burg et al. 2015 ("vdB15", **arXiv:1412.2137**).

**Photometry and bands**
- g and r from CFHT MegaCam, plus u and i from the INT Wide Field Camera. vdB15 Sect. 2.1 adds archival MegaCam u for 7 clusters and i for 2.
- Sources are detected in r.
- Aperture fluxes use Gaussian weights on PSF-homogenised stacks.
- Typical 5σ limits in u, g, r, i: 24.3, 24.8, 24.2 and 23.3 mag.

**M/L method: SED fitting, not a colour–M/L relation**
- ugri aperture fluxes are fitted with Bruzual & Charlot (2003) libraries.
- Star formation history: SFR ∝ e^(−t/τ), with τ from 10 Myr to 10 Gyr.
- Chabrier (2003) IMF, solar metallicity, Calzetti et al. (2000) dust law.
- "Using the total flux measured in the r-band, luminosities are converted into stellar masses." So the fitted M/L is applied to the total r flux.
- vdB15 Sect. 3.1 names the fitting code as FAST (Kriek et al. 2009) and says it gives M/L in the r band.
- **The age grid is not stated** in either Ghizzardi or vdB15.

**Distance and membership**
- "Initially each galaxy is assumed to be part of the cluster (to set the luminosity distance)."
- Fore- and background galaxies are then subtracted statistically using the COSMOS ugri catalogue (Muzzin et al. 2013).
- vdB15 Sect. 3.1 adds the details:
  - an EAZY photometric-redshift cut at z < 0.3 is applied to both the cluster fields and the COSMOS reference;
  - the reference is an r-selected COSMOS catalogue covering 1.62 deg²;
  - the cluster distance modulus is assigned to every galaxy.
- There is no galaxy-by-galaxy membership.

**Which galaxies**
- All galaxies with M* > 10⁹ M☉ (Sect. 4.1, footnote 4). The footnote says this captures more than 98% of the stellar mass of a Schechter function with α ≈ −1 and M* = 10^10.8 M☉ (from van der Burg et al. 2018).
- **"the BCG is included in our mass computation, while the Intra-Cluster Light (ICL) is not."** (Sect. 4.1)
- Profiles are centred on the X-ray centroids, not on the BCGs as in vdB15. The two centres differ by at most 2% of R500.

**Projection**
- The cumulative profiles in Fig. 16 are "derived by integrating within a projected radius".
- To get the mass inside a sphere of radius R500, they assume a gNFW galaxy distribution with c = 0.72 and α = 1.64 (vdB15). They find that 75% of the line-of-sight mass lies inside the R500 sphere and "multiply the stellar mass estimates by a factor 0.75". These spherical values go into Table 3.
- In vdB15 (Sect. 5, Table 2), c = 0.72 and α = 1.64 are the gNFW fit to the *number density* of galaxies with M* > 10^10 M☉. The fit to the stellar-mass density is c = 0.64 and α = 1.63.

**Cosmology**
- Flat ΛCDM with H0 = 70 km s⁻¹ Mpc⁻¹, Ωm = 0.3 and ΩΛ = 0.7 (end of Sect. 1). vdB15 uses the same.

**Values: Table 3, M_star,500 (spherical), in 10¹² M☉**

| Cluster | M_star,500 |
|---|---|
| A1795 | 3.02 ± 0.28 |
| A85 | 2.10 ± 0.33 |
| A644 | 3.70 ± 0.38 |
| A2319 | 5.11 ± 0.46 |
| ZW1215 | 3.34 ± 0.39 |
| A2029 | 6.51 ± 0.47 |
| A2142 | 6.97 ± 0.51 |

**Stellar fraction (Sect. 4.2)**
- "All our clusters are within" M_star,500/M500 = 0.4% to 1%, "except for A85 which is slightly below".
- Table 4 gives the fits of M_star,500 against M_gas,500 and M500. The slopes are 1.11 (+0.60/−0.52) and 1.39 (+0.68/−0.75), with about 21–22% scatter.

**Systematic uncertainties they quote**
- **IMF.** A Salpeter-diet IMF "provides values that exceed by a factor ∼2 those obtained using a Chabrier IMF" (Sect. 4.2).
- **Comparison with other samples.** Chiu et al. 2018 and Lin et al. 2012 (the latter rescaled from Kroupa to Chabrier by ×0.76) are higher than X-COP. At M_gas,500 = [6.5, 8.8, 15] × 10¹³ M☉ the excesses are [90%, 58%, 12%] for Chiu and [78%, 56%, 24%] for Lin. Their conclusion: "assuming a systematic discrepancy of 50%–60% … seems reasonable".
- **Table 7 budget.** The stellar masses may be underestimated by up to 60%. The ICL could add up to 50%; they take the DES stack of Zhang et al. (2019) conservatively, and note that the SDSS stack of Zibetti et al. (2005) gives about 10%. The mass-loss factor r_o carries 20%.
- **Radial extent (Sect. 5.3).** M(<R200)/M(<R500) is about 1.8 for stars against 1.5 for gas.
- **Statistical errors.** 100 bootstrap perturbations of the fluxes per galaxy, plus Poisson noise and cosmic variance (Moster et al. 2011). Fig. 16 shows the statistical and cosmic-variance errors separately.

### INFERRED

- **Stellar-to-gas ratio at R500 (spherical, from Table 3): 0.023 to 0.049, median 0.044, mean 0.040.** A85 is lowest at 0.023; A2029 and A2142 are highest at 0.049 and 0.047.
- **M*/M500: 0.37% to 0.78%, median 0.65%.**
- **The released profiles are projected, not spherical.** The release file `<cluster>_mstar.fits` has two tables:
  - `MSTAR_RAW`, in R/R500, with statistical and cosmic-variance errors given separately;
  - `MSTAR_SMOOTHED`, in kpc, with total errors.

  The header says only "Stellar mass profile used in Ghizzardi et al. 2020". But 0.75 × (release value at R/R500 = 1) reproduces Table 3 to 0.4–4% for six clusters. A644 is the exception: 0.75 × release gives 4.14 against the tabulated 3.70.
- **The project mixes projected stars with spherical gas.** In `research_work/results/hot-companion/code/run.py`, `load_xcop()` interpolates the release `Mstar` on a *spherical* radius grid. It also forms the median star/gas ratio `FRAC` from projected stars divided by spherical gas. At R500 that ratio is **0.035 to 0.073 (median 0.061)**, against 0.023 to 0.049 on the paper's spherical basis. This matches the "X-COP 0.035–0.073" quoted in `THEORY.md`.

  So the project's X-COP stellar masses at R500 are too high by about 1/0.75 ≈ 1.33. The excess grows inward. For the vdB15 gNFW (my integration, satellites only, with R500/R200 between 0.62 and 0.68), the ratio of spherical to cylindrical mass is:

  | Radius | Sphere / cylinder |
  |---|---|
  | R500 | 0.71–0.80 (0.75 for a line-of-sight cut-off near 1.2–1.4 R200) |
  | 0.5 R500 | 0.62–0.67 |
  | 0.3 R500 | 0.58–0.62 |
  | 0.1 R500 | 0.53–0.54 |

  The paper does not say how far along the line of sight it integrates. Near the centre the BCG dominates; its projected and spherical masses are nearly equal, so the true inner correction is closer to 1 than this table.
- **How the masses scale with the distance law.** At fixed fitted M/L, M* ∝ D_L², and each radius in kpc ∝ D_A. At z < 0.1 the SED-fitted M/L barely depends on distance.

---

## 2. MACS J0025.4−1222 (Bradač et al. 2008, arXiv:0806.2320 v2)

### STATED

**The numbers (Table 3)**
- The table is titled "2-D projected enclosed mass within a radius of 300 kpc centered on BCG1 and BCG3 and within 500 kpc centered on the gas peak".
- Galaxies (stars): **SE–BCG1 0.027 ± 0.008 × 10¹⁴ M☉; NW–BCG3 0.019 ± 0.006 × 10¹⁴ M☉.** Around the gas peak within 500 kpc: 0.05 ± 0.01 × 10¹⁴ M☉.
- The table note says "stellar mass from F814W data".

**Method (Sect. 6)**
- Members come from colour–colour cuts, using DEIMOS spectroscopic redshifts where available. The cuts (Sect. 4.2) are 0.0 < F555W − F450W < 1.4 and 1.8 < F555W − F814W < 3.0. They "primarily select red cluster members".
- Observed I-band (ACS F814W) MAG_AUTO magnitudes are converted to rest-frame K luminosities using:
  - M_K,☉ = 3.28;
  - Galactic extinction A_I = 0.054;
  - a K-correction of 2.464 from a Coleman et al. (1980) elliptical template, extended with GISSEL (Bruzual & Charlot 1993);
  - "(we assume zero evolutionary correction)".
- The authors add: "our template is likely slightly redder than the cluster members (possibly giving a larger stellar mass)".
- M/L: "we follow Drory et al. (2004) and assume … M*/L_K = 0.74 ± 0.30 (see also Bell et al. 2003)."
- The paper also gives:
  - total M/L_K(< 500 kpc) = 130 (+40/−80) and a stellar fraction of 1.0 (+0.7/−0.4)%;
  - M/L_K(< 300 kpc) = 70 (+40/−50) for SE and 100 (+40/−60) for NW;
  - "Both subclusters have a stellar-to-total mass ratio of 1 per cent within 300 kpc".
- Cosmology: ΛCDM with Ωm = 0.3, ΩΛ = 0.7, H0 = 70, giving 6.61 kpc/arcsec at z = 0.586.

**Not stated**
- The IMF.
- A magnitude limit for the stellar-mass members.
- Any ICL treatment.
- Whether the BCGs are in the stellar sums. The apertures are centred on them, but the weak-lensing catalogue of Sect. 4.2 removes objects brighter than the BCGs.
- Where the ±0.30 on M/L_K comes from.
- Whether Table 3 was measured on the smoothed luminosity map (Fig. 4 is smoothed with an 80 kpc FWHM Gaussian).

### INFERRED

- **The 0.74 is Drory's value for massive field galaxies at z = 0.5.** Drory et al. 2004 (astro-ph/0403041), Table 1, lists a mean M/L_K of **0.74** (and M/L_B of 2.58) at **z = 0.5 for M > 10¹¹ h⁻¹ M☉**. The exact match makes this the very likely source.
- **Drory's M/L rests on a Salpeter IMF.** Drory Sect. 3 uses a Salpeter IMF (0.1–100 M☉), Maraston (1998) models, solar metallicity, τ-models with 28 ages from 0.001 to 14 Gyr, and A_V from 0 to 3.
  - So Bradač's stellar masses are effectively **Salpeter-based**, for field galaxies, from SED fits.
  - To put them on a Chabrier basis, divide by about 1.7–1.8, i.e. 0.24–0.25 dex (Hilton 2013 and Treu 2010 respectively).
  - Bradač also cites Bell et al. 2003, which uses a "diet Salpeter" IMF. Their stated source, though, is Drory.
- **The K-band light is probably overestimated.** Observed F814W at z = 0.586 samples rest ≈ 5100 Å. The conversion to rest K uses a non-evolving local elliptical colour, which overestimates L_K if the members are bluer, as the authors themselves note. The Drory M/L_K already refers to z = 0.5 populations.
- **Implied K luminosities:** 0.027 × 10¹⁴ / 0.74 ≈ 3.6 × 10¹² and 0.019 × 10¹⁴ / 0.74 ≈ 2.6 × 10¹² L_K,☉.
- **ICL is almost certainly excluded,** since MAG_AUTO galaxy photometry does not capture it.

---

## 3. El Gordo, ACT-CL J0102−4915 (z = 0.87)

### 3a. Menanteau et al. 2012 (arXiv:1109.0953 v3), Sect. 3.1.4 "Stellar Mass from SED Fitting" and Sect. 3.2

**STATED**

- **Photometry:** optical griz (SOAR, plus FORS2 RIz, co-added into an "über" set; SExtractor MAG_AUTO) and Spitzer/IRAC 3.6 and 4.5 μm.
  - The IRAC magnitudes are 4″ apertures with aperture corrections of −0.35 and −0.37 mag.
  - IRAC 80% completeness is about 22.6 AB in both channels.
  - Optical detection is on the über i-band at a 1.5σ threshold.
- **Members:**
  - 89 spectroscopic members, "augmented" to 411 in total using photometric redshifts. The photo-zs are from BPZ, recalibrated on 517 redshifts in seven clusters, with δz ≈ 0.04.
  - Galaxies within |Δz| ≤ 0.06 of the cluster mean are selected.
  - No magnitude limit is stated for this selection.
- **SED grid:**
  - Bruzual & Charlot 2003, **solar metallicity only**.
  - Exponentially declining SFHs with **20 values of τ from 0.1 to 20 Gyr**.
  - **53 ages from 0.001 to 7.0 Gyr.**
  - **Chabrier 2003 IMF.**
  - Calzetti 2000 dust, with E(B−V) from 0 to 0.5 in steps of 0.02.
  - Fits are by χ² with an analytic normalisation. Errors come from Monte Carlo.
- **The age grid is stated as 0.001–7.0 Gyr.** It is not "the age of the universe at z = 0.87".
- **On the fitted parameters:** "we do not expect the values of these parameters [age, τ, E(B−V)] to be very robust. However, stellar mass is well constrained by the IRAC photometry…" **The paper does not mention fits piling up at the maximum age.**
- **The factor-of-two statement is confirmed, verbatim:** "uncertainties in the IMF and the modeling of thermally pulsating AGB stars … is likely to lead to the stellar mass estimates only being accurate to within a factor of two; these systematic uncertainties are not taken into account in quoted errors on M*."
- **Total:** "The total stellar mass for all cluster members within r200 = 2111 h70⁻¹ kpc is **M*_200 = 1.31 ± 0.26 × 10¹³ M☉**. This suggests a ratio of stellar mass to total mass of < 1% within r200."
- **Components (Sect. 3.2):** the galaxies are split along a line between two stated coordinates. Each component, within its own r200, has **7.5 ± 1.4 × 10¹² M☉ (NW) and 5.6 ± 1.3 × 10¹² M☉ (SE)**.
- The BCG is a spectroscopic member at z = 0.87014. It is about 1.5 mag brighter than any other galaxy, blue, E+A+[O II], with R_e = 10.4 kpc (Sect. 3.4).
- **Not stated:** whether the BCG is explicitly in the sum (though it is a member), any ICL treatment, and whether r200 is a projected or a spherical radius.
- **Cosmology:** flat ΛCDM with Ωm = 0.27, ΩΛ = 0.73, H0 = 70 h70. That gives D_A = 1617 Mpc and 470 kpc/arcmin.

**INFERRED**

- **The 7.0 Gyr cap is slightly older than the universe in their own cosmology.** The age of the universe at z = 0.87 for their parameters is 6.60 Gyr (my integration). If ages above 6.6 Gyr were allowed, the cap is not an age-of-universe prior. There is no information on how many fits reach the edge. In a static universe with no age limit, a fit that wanted older ages would be capped at 7 Gyr, which biases M/L low. It is unknown whether this happens.
- **The r200 sum may be incomplete in the outer parts.** r200 corresponds to a 4.5′ radius. The FORS2 field is 6.8′ × 6.8′, a half-width of about 3.4′ or 1.6 Mpc. Photo-z membership needs the optical data, so the sum probably misses members beyond about 1.6 Mpc. Hilton et al. 2013 (Sect. 3.1) note that for most ACT clusters "the optical data do not provide coverage out to R500".
- **Stellar fraction:** 1.31 × 10¹³ / M200a (2.16 × 10¹⁵) ≈ 0.6%.

### 3b. Hilton et al. 2013 (arXiv:1301.0780 v2): ACT IRAC stellar content, including El Gordo

**STATED**

- **Values:**
  - **El Gordo M*_500 = 18.8 (+3.5/−4.1) × 10¹² M☉**, from observed-frame m_3.6 = 14.22 ± 0.09 (Table 4, Sect. 5).
  - **BCG M* = 16.2 (+2.8/−3.4) × 10¹¹ M☉** (Table 3).
  - R500 = 1.1 Mpc, and dynamical M500 = 9.8 ± 2.3 × 10¹⁴ M☉ (Table 1).
- **Method:**
  1. Background-subtracted 3.6 μm flux within R500, with the Extended Groth Strip as the background field.
  2. Summed to m* + 2, then extrapolated to m* + 5 assuming α = −0.8. The extrapolation adds less than 10%.
  3. The BCG is added.
  4. A constant M/L from a single τ = 0.1 Gyr burst (BC03, solar metallicity, **Salpeter IMF**, formed at z_f = 3) converts light to mass.
  5. The result is multiplied by 0.73 to deproject, assuming an NFW profile with c = 2.8 integrated to 3 R500.
  6. "Note that our measurement does not include any contribution from the intracluster light (ICL)."
- **Error bars:** they reflect z_f = 2 and z_f = 5 only. "Systematic errors due to the choice of stellar population model and/or IMF are neglected."
- **IMF conversion:** "if we adopted a Chabrier rather than Salpeter IMF, our stellar mass estimates would be **0.24 dex lower**". Models with more TP-AGB light would also lower the masses (Sect. 4.2 and Sect. 6).
- **Cosmology:** Ωm = 0.3, ΩΛ = 0.7, H0 = 70.

**INFERRED**

- On a Chabrier basis (−0.24 dex): M*_500 ≈ **10.8 (+2.0/−2.4) × 10¹² M☉** and BCG ≈ 9.3 × 10¹¹ M☉.
- f*_500 ≈ 1.9% (Salpeter), or about 1.1% (Chabrier).
- The z_f = 3 burst is about 4.2 Gyr old at z = 0.87 in their cosmology (3.1 Gyr for z_f = 2, 5.2 Gyr for z_f = 5). This M/L is tied to the ΛCDM timeline.
- Menanteau's 1.31 × 10¹³ M☉ (Chabrier, within r200, probably projected) and Hilton's 1.08 × 10¹³ M☉ (Chabrier-equivalent, spherical, within R500) agree within the stated factor of two.

### 3c. Other El Gordo papers checked

- **Jee et al. 2014 (arXiv:1309.5097), Sect. 4.5 and Table 3:** light only, no stellar mass.
  - Rest-frame B from B_rest = −0.565 (F775W − F850LP) + 1.38 − DM, using Kinney et al. 1996 templates.
  - Within r = 386 h70⁻¹ kpc apertures: L_B = 3.22 × 10¹² (NW) and 3.40 × 10¹² h70⁻² L_B,☉ (SE).
  - Total lensing M/L_B = 130 ± 10 (NW) and 88 ± 7 (SE).
- **Kim et al. 2021 (arXiv:2106.00031):** no stellar-mass estimate. It only cites Menanteau's stellar-mass density map.
- **Sifón et al. 2013 (arXiv:1201.0991):** none.
- I also searched the local texts of arXiv:1304.0455, 1310.6786, 1405.2617, 1412.1826, 1511.02578, 1607.04641, 1608.05413, 1711.08438, 1905.00025, 2012.03950, 2209.02718, 2308.00744, 2309.10374 and 2503.18613. None gives an independent cluster-total stellar mass. Caminha et al. 2023 and Diego et al. use total mass-to-light scalings of galaxies in their lens models.

---

## 4. Abell 520 (z = 0.201)

**STATED: none of these papers gives a stellar M/L or stellar masses for member galaxies or clumps.** Every "M/L" below is a *total* (lensing) mass-to-light ratio. Some are computed after subtracting the gas mass.

- **Mahdavi et al. 2007 (arXiv:0706.3048), Sect. 2.1 and Table 2.**
  - Rest-frame **B** from red-sequence galaxies in CFHT g′ and r′: r < 22, and g′ − r′ no more than 0.25 mag bluer than the red edge of the sequence.
  - Table 2 ("Masses and Mass-to-light Ratios", 150 kpc apertures) gives L_B for each peak.
  - For the "Cluster" aperture (710 kpc): L_B = 21.57 × 10¹¹ h70⁻² L_B,☉ and M/L_B = 232 ± 25.
- **Jee et al. 2012 (arXiv:1202.6368), Sect. 3.2 and Table 1.**
  - Rest-frame **B** from CFHT g and g − r via Kinney et al. 1996 templates.
  - Members follow the M07 definition.
  - Within 710 kpc, L_B = 1.7 × 10¹² L_B,☉ and M/L_B = 263 ± 28.
  - Table 1 gives L_B and total M/L_B for P1–P6 (P3: 588 ± 56).
- **Clowe et al. 2012 (arXiv:1209.2143), Sect. 3 (Table 1, smoothed maps; Table 2, aperture densitometry) and Sect. 4.**
  - The light is **observed-frame F814W** ("roughly a restframe R passband"). The table header's unit "L_z☉" does not match this.
  - Members are colour-selected with HST and Magellan colours.
  - Interlopers contribute about 3 × 10⁹ L☉ per aperture, less than 2%.
  - M/L is computed after subtracting the gas mass. The weighted mean for structures 1, 2, 4 and 5 is 108 ± 24.
  - Cosmology: Ωm = 0.27, ΩΛ = 0.73, H0 = 70.
- **Jee et al. 2014 (arXiv:1401.3356), Sect. 5.1 and Table 2.**
  - Gives **L_B, L_R and L_F814W** and the corresponding total M/L after subtracting the gas upper limit.
  - R_rest = F814W − 0.083 (F606W − F814W) + 0.42 − DM, from Kinney templates.
  - Members are those whose F606W − F814W and F435W − F606W colours match the spectroscopic members.
  - For P1, P2, P4 and P5 the mean is M/L_R ≈ 114 and M/L_B ≈ 131.
  - The paper states that Clowe et al. 2012 applied no K-correction.
- **Girardi et al. 2008 (arXiv:0809.3139):** no luminosity or M/L estimate. It is a galaxy-dynamics paper.
- **Okabe & Umetsu 2008 (astro-ph/0702649), Sect. 4.2 and Table 8:** gives the cluster total M/L in the **i′** band for A520. I did not transcribe the value.

**INFERRED:** to get A520 stellar masses, the project must supply its own stellar M/L for red-sequence galaxies at z = 0.2. The published B, R and F814W luminosities are the inputs.

---

## 5. Bullet Cluster (1E 0657−56, z = 0.296)

### STATED

**Clowe et al. 2006 (astro-ph/0608407), Sect. 3 and Table 2 ("Component Masses", 100 kpc apertures)**
- Stellar masses: main BCG 0.54 ± 0.08, main plasma 0.23 ± 0.02, subcluster BCG 0.58 ± 0.09, subcluster plasma 0.12 ± 0.01, all in 10¹² M☉.
- Method, verbatim: "Stellar masses are calculated from the **I-band luminosity of all galaxies equal in brightness or fainter than the component BCG**. The luminosities were converted into mass assuming (Kauffmann et al. 2003) **M/L_I = 2**. The assumed mass-to-light ratio is highly uncertain (**can vary between 0.5 and 3**)…"
- The quoted errors cover only the luminosity measurement.
- There is no colour selection, so "these measurements are an **upper limit** on the stellar mass" because they include non-members.
- Cosmology: Ωm = 0.3, λ = 0.7, H0 = 70, giving 4.413 kpc/arcsec.
- **Not stated:** the IMF, whether "I-band" means observed or rest-frame I, and which data set supplied the I-band photometry.

**Bradač et al. 2006 (astro-ph/0608408)**
- **No stellar M/L is assumed and no stellar mass is computed.**
- The 20 brightest members (F606W, colour-selected with Magellan data) enter the lens model as non-singular isothermal spheres, scaled by σ ∝ L^(1/4) and r_c ∝ L^(1/2). This is total galaxy mass.
- Cosmology: Ωm = 0.3, ΩΛ = 0.7, H0 = 70.

**Also checked: Paraficz et al. 2016 (arXiv:1209.0384)**
- Galaxy haloes are scaled with the Faber–Jackson and fundamental-plane relations. These are total masses, with no stellar M/L.

### INFERRED

- **Kauffmann et al. 2003 (astro-ph/0204055) uses a Kroupa (2001) IMF** (Sect. 4, "A Library of Star Formation Histories"). They state that Salpeter would give about twice the mass. So Clowe's M/L_I = 2 is nominally on a Kroupa basis.
- Kauffmann's M/L values are in the z band. Their colour bins top out at log(M/L)_z > 0.2, i.e. above about 1.6. An M/L_I of 2 therefore sits at the heavy end for old populations on that basis.

---

## 6. Fundamental-plane (FP) evolution of M/L for cluster early types

### STATED rates (rest-frame B)

- **van Dokkum & van der Marel 2007 (astro-ph/0609587).** Abstract and Sect. 5.3: **d log(M/L_B)/dz = −0.555 ± 0.042** for cluster galaxies with M ≳ 10¹¹ M☉.
  - 16 clusters at 0.02 ≤ z ≤ 1.28; the Coma zero point is at z = 0.024.
  - The slope moves only between −0.54 and −0.56 for reasonable FP coefficients.
  - The 90% upper limit on intrinsic cluster-to-cluster scatter is 0.057 dex.
  - Per-cluster offsets are in Table 3, for example:

    | Cluster | z | Δlog(M/L_B) |
    |---|---|---|
    | A2218 | 0.176 | +0.009 ± 0.037 |
    | A2390 | 0.228 | −0.035 ± 0.065 |
    | CL1358+62 | 0.327 | −0.162 ± 0.029 |
    | CL0016+16 | 0.546 | −0.281 ± 0.032 |
    | MS2053−04 | 0.583 | −0.287 ± 0.056 |
    | MS1054−03 | 0.831 | −0.427 ± 0.040 |
    | RXJ0152−13 | 0.837 | −0.449 ± 0.048 |
    | RXJ1226+33 | 0.892 | −0.558 ± 0.056 |

  - Field galaxies (same table) evolve slightly faster.
  - Cosmology: Ωm = 0.3, ΩΛ = 0.7, H0 = 71.
- **Holden et al. 2010 (arXiv:1009.4479)**, MS 1054−03 at z = 0.831 against Coma:
  - At fixed σ (virial estimator, Sect. 4.2): **Δlog(M/L_B) = −0.50 ± 0.03, or −0.60 ± 0.04 per unit z** (abstract).
  - From the FP offset (Sect. 4.1.3): Δlog(M/L_B) = −0.44 ± 0.03.
  - Cosmology: Ωm = 0.27, ΩΛ = 0.73, H0 = 71.
- **Saglia et al. 2010 (arXiv:1009.0645 v2, including the 2016 erratum), EDisCS**, abstract:
  - Clusters: **Δlog(M/L_B) = (−0.54 ± 0.01) z = (−1.61 ± 0.01) log(1+z)**, independent of cluster velocity dispersion.
  - Field: (−0.76 ± 0.01) z.
  - The evolution is milder after an incompleteness correction or once size and σ evolution are included.
- **Treu et al. 2005 (astro-ph/0503164):**
  - Field spheroidals in GOODS-N: ⟨d log(M/L_B)/dz⟩ = −0.72 (+0.07/−0.05) ± 0.04, and the rate depends on mass.
  - For cluster E/S0s they adopt **−0.46 ± 0.04** from van Dokkum & Stanford 2003 (their Eq. for the cluster fiducial; van Dokkum et al. 1998 give −0.49 ± 0.05).
- **Other bands:** none of these FP papers gives a cluster FP rate outside B.
  - A model-based field comparison: Drory et al. 2004, Table 1, SED fits with Salpeter IMF and Maraston 1998 models, M > 10¹¹ h⁻¹ M☉. M/L_K = 0.74, 0.68, 0.62, 0.55 and M/L_B = 2.58, 2.27, 1.83, 1.23 at z = 0.5, 0.7, 0.9, 1.1.
  - INFERRED from that table: K evolves about 2.5 times more slowly than B.

### INFERRED: implied fading from each linear rate

Values are Δlog(M/L_B) relative to z = 0, with the ratio (M/L_B)(z)/(M/L_B)(0) in parentheses. vDvdM07 fit a free intercept (14 degrees of freedom for 16 points), so slope × z is an approximation; the Table 3 points above are the direct measurements.

| z | vDvdM07 −0.555 | Holden −0.60 | Saglia −0.54 | Saglia −1.61 log(1+z) | Cluster fiducial −0.46 |
|---|---|---|---|---|---|
| 0.20 | −0.111 ± 0.008 (0.77) | −0.120 (0.76) | −0.108 (0.78) | −0.127 (0.75) | −0.092 (0.81) |
| 0.30 | −0.167 ± 0.013 (0.68) | −0.180 (0.66) | −0.162 (0.69) | −0.183 (0.66) | −0.138 (0.73) |
| 0.59 | −0.327 ± 0.025 (0.47) | −0.354 (0.44) | −0.319 (0.48) | −0.324 (0.47) | −0.271 (0.54) |
| 0.87 | −0.483 ± 0.037 (0.33) | −0.522 (0.30) | −0.470 (0.34) | −0.438 (0.37) | −0.400 (0.40) |

These are dynamical M/L values inside R_e, interpreted as fading of the stellar populations. They are **not** stellar M/L normalisations.

### How the FP M/L depends on the assumed distances

**STATED**
- vDvdM07, Sect. 5.1: they fit log r_e = a log σ + b log I_e + c with a = 1.20 and b = −0.83 (Jørgensen et al. 1996, B band). Then **Δlog(M/L) = (c_z − c₀)/b**, with c_{z,i} = log r_e − a log σ − b log I_e.
- r_e is converted from arcseconds to kpc with the adopted cosmology, i.e. through D_A. log I_e comes from the rest-frame B surface brightness μ_e "corrected for (1+z)⁴ cosmological surface brightness dimming" (Sect. 4).
- The analysis "depends on the value of Ω but not on the Hubble constant" (Sect. 3.2.2).
- The virial estimator is M/L_B = 5σ²/(2πG r_e ⟨I_e⟩) (Holden et al. 2010, Sect. 4.2).

**INFERRED: conversion recipe**
- Physically, r_e = θ_e D_A and I_e ∝ F/θ_e² × (D_L/D_A)². The "(1+z)⁴ correction" equals (D_L/D_A)² only when the Etherington relation D_L = (1+z)² D_A holds, as in any FRW metric.
- For a new distance law, define δA(z) = log₁₀[D_A,new(z)/D_A,ΛCDM(z)] and δL(z) = log₁₀[D_L,new(z)/D_L,ΛCDM(z)]. Measure both relative to the local reference at z₀ ≈ 0.024, i.e. ΔδA = δA(z) − δA(z₀) and ΔδL = δL(z) − δL(z₀). This cancels any overall change of scale, such as H0.
- **Virial estimator (M/L ∝ D_A/D_L²):** Δlog(M/L)_new = Δlog(M/L)_published + ΔδA − 2 ΔδL.
- **FP-offset estimator (b = −0.83):** Δlog(M/L)_new = Δlog(M/L)_published + (2 + 1/b) ΔδA − 2 ΔδL = published **+ 0.795 ΔδA − 2 ΔδL**.
  - Check: when Etherington holds (ΔδL = ΔδA), this reduces to (1/b) ΔδA = −1.20 ΔδA, as expected for an Ω-only dependence.
- Both recipes assume the FP coefficients and the local zero point stay fixed.
- Stellar masses from photometry and SED fits scale as M* ∝ D_L² at fixed M/L, and physical apertures scale with D_A. Rest-frame colours and 4000 Å breaks depend only on z, not on distance.

### Do red-sequence colours and 4000 Å breaks get bluer or weaker with z? STATED: yes

- **Holden et al. 2010, Sect. 4.3 (cluster, most cosmology-independent).**
  - Red-sequence E and S0 galaxies with log σ > 2.2, measured at **fixed σ**, in "redshifted" U and V bandpasses. These are rest-frame colours with no distance involved.
  - **Δ(U − V)_z = −0.24 ± 0.02 mag** between Coma (z = 0.023) and MS 1054−03 (z = 0.831). INFERRED: about −0.30 mag per unit z.
- **Bell et al. 2004 (astro-ph/0303394), COMBO-17, all environments, Sect. 4.**
  - The red-sequence colour at M_V − 5 log h = −20 follows **⟨U − V⟩ = 1.40 − 0.31 z** including local points, or 1.48 − 0.40 z from COMBO-17 alone.
  - The red sequence "become[s] redder by ∼0.3 mag" between z = 1.1 and z = 0.2.
  - This is at fixed absolute magnitude, so it depends weakly on distance through M_V; the colour–magnitude slope is only −0.08 mag/mag.
- **Moresco et al. 2012 (arXiv:1201.3609), Table 2, massive passive early types (mostly field).**
  - The median D4000_n in the high-mass bins falls with z: 1.992 (z = 0.16), 1.947 (0.22), 1.875 (0.53), 1.854 (0.67), 1.771 (0.83), 1.71 (1.02) and 1.575 (1.24).
  - Caveats: the mass bins differ between subsamples, and the stellar masses behind the bin selection assume ΛCDM distances.

---

## 7. IMF of massive ellipticals (σ ≈ 250–300 km/s)

Conversions I apply, all as stated in the cited papers:
- Salpeter to Chabrier: −0.24 dex (Hilton 2013) or −0.25 dex (Treu 2010; Auger 2010b table note).
- Kroupa is about 0.06 dex heavier than Chabrier (Auger 2010b, Sect. 4.1).

### STATED

- **Treu et al. 2010 (arXiv:0911.3392): strong lensing plus dynamics, SLACS, 56 lenses, NFW haloes.**
  - Relative to Chabrier: **⟨log α⟩ = 0.25 ± 0.03 ± 0.02**. My arithmetic: 10^0.25 = **1.78 times the Chabrier stellar mass**.
  - Relative to Salpeter: 0.00 ± 0.03 ± 0.02. Other stellar-population choices give 0.03–0.06 relative to Salpeter and 0.27 relative to Chabrier (Table 2).
  - A "tentative" trend (Sect. 3.2): log α_Salp = (1.31 ± 0.16) log σ* − 3.14 ± 0.01.
- **Auger et al. 2010b (arXiv:1007.2409): lensing, dynamics and weak lensing, 53 early types.**
  - "the data clearly prefer a Salpeter-like IMF over a lighter IMF (e.g., Chabrier or Kroupa), irrespective of the choice of DM halo."
  - Free-IMF fits (Table 1; Salpeter corresponds to log α = 0, Chabrier to about −0.25): log α = +0.03 ± 0.03 for an NFW halo (η = 0.08 ± 0.04), −0.08 ± 0.02 for Gnedin contraction and −0.11 ± 0.02 for Blumenthal contraction.
- **Cappellari et al. 2012 (arXiv:1202.3308, Nature Letter; main text and Fig. 2): dynamics of the ATLAS3D galaxies.**
  - The IMF normalisation varies "by up to a factor of three in mass". "For increasing (M/L)_stars the normalization of the inferred IMF varies from the one of Kroupa/Chabrier up to an IMF more massive than Salpeter."
  - "A Salpeter normalization at larger (M/L)_stars is consistent on average with results from strong gravitational lensing … (σ ≳ 200 km/s)."
  - The Letter plots the trend against (M/L)_stars, colour-coded by σ_e. **It gives no explicit fit against σ.**
- **Cappellari et al. 2013, ATLAS3D XX (arXiv:1208.3523), Sect. 4.3** (the σ fit from the same team):
  - **log[(M/L)_stars/(M/L)_Salp] = −0.12 ± 0.01 + (0.35 ± 0.06) log(σ_e / 130 km/s)**, with about 0.10 dex observed scatter.
  - This implies a transition from Kroupa to Salpeter over σ_e ≈ 90–290 km/s.
  - The fit uses galaxies with M_JAM < 2 × 10¹¹ M☉. The authors note their slope is about 3.5 times shallower than Treu et al. 2010's trend.
- **Conroy & van Dokkum 2012 (arXiv:1205.6473): absorption-line spectroscopy, 38 early types plus the M31 bulge.**
  - The IMF becomes more bottom-heavy with σ and [Mg/Fe]. At the highest σ it is "more bottom-heavy than even a Salpeter IMF".
  - They quote the IMF as (M/L_K)/(M/L_K)_MW, where MW means Kroupa 2001. In these units "a Salpeter IMF has a value of ≈1.6" (Sect. 4.1, Fig. 5). Table 2 lists per-galaxy σ, M/L_K and (M/L_K)_MW, with typical errors of 7%.

### INFERRED: my evaluation of those relations and tables

| Source | σ = 250 km/s | σ = 300 km/s |
|---|---|---|
| Treu 2010 trend | α_Salp ≈ 1.00 → ≈ **1.8 × Chabrier** | α_Salp ≈ 1.27 → ≈ **2.3 × Chabrier** |
| Cappellari 2013 fit | ≈ 0.95 × Salpeter → ≈ **1.7 × Chabrier** | ≈ 1.02 × Salpeter → ≈ **1.8 × Chabrier** |
| Auger 2010b (not a function of σ) | NFW: **1.9 × Chabrier**; contracted haloes: 1.4–1.5 × Chabrier | same |

**Conroy & van Dokkum 2012, Table 2, by σ bin:**

| σ bin (km/s) | Galaxies | Mean (median) ratio to Kroupa | ≈ ratio to Chabrier |
|---|---|---|---|
| 250–300 | 5 | **1.62 (1.66)** | ≈ **1.85** |
| 200–250 | 14 | 1.51 | ≈ 1.7 |
| M87, σ = 385 | 1 | 1.90 | – |

**Consensus:** at σ ≈ 250–300 km/s, the stellar mass is about **1.7–1.9 times** the Chabrier value (≈ Salpeter), with a plausible range of about 1.4–2.3.

---

## Key numbers

The IMF column gives the basis of each number as published. "Chabrier-equiv" is my conversion.

| Object | Source (arXiv; location) | Stellar mass | Aperture | IMF | Light and M/L | BCG | ICL | Cosmology |
|---|---|---|---|---|---|---|---|---|
| X-COP, 7 clusters | 2007.01084; Table 3 | M*,500 = 2.10–6.97 × 10¹² M☉ (spherical; 0.75 × projected) | R500 sphere | Chabrier | ugri SED (FAST, BC03, τ-models, solar Z, Calzetti); M/L_r on total r | incl. | excl. | 70 / 0.3 / 0.7 |
| X-COP (derived) | inferred from Table 3 | M*/Mgas = 0.023–0.049 (median 0.044); M*/M500 = 0.37–0.78% (stated 0.4–1%) | R500 | Chabrier | – | incl. | excl. | – |
| X-COP release file | FITS `_mstar` | **projected**: 1.33 × the spherical values at R500 | cylinder | Chabrier | – | incl. | excl. | 70 / 0.3 / 0.7 |
| MACS J0025 | 0806.2320; Table 3, Sect. 6 | 0.027 ± 0.008 (SE), 0.019 ± 0.006 (NW) × 10¹⁴; 0.05 ± 0.01 × 10¹⁴ (500 kpc) | 300 kpc projected | not stated; M/L_K = 0.74 from Drory04 (Salpeter, Maraston98), inferred | F814W → rest K, local E template, no evolution | not stated | not incl. (inferred) | 70 / 0.3 / 0.7 |
| El Gordo | 1109.0953; Sect. 3.1.4, 3.2 | 1.31 ± 0.26 × 10¹³ (all); NW 7.5 ± 1.4, SE 5.6 ± 1.3 × 10¹² | r200 = 2111 kpc | Chabrier | griz + IRAC SED; BC03; solar Z; τ 0.1–20 Gyr; ages ≤ 7.0 Gyr; ×2 systematic | member (not explicit) | not mentioned | 70 / 0.27 / 0.73 |
| El Gordo | 1301.0780; Tables 3, 4 | M*,500 = 18.8 (+3.5/−4.1) × 10¹² (Salpeter) → ≈ 10.8 × 10¹² Chabrier-equiv; BCG 16.2 × 10¹¹ | R500 = 1.1 Mpc sphere (×0.73) | Salpeter | IRAC 3.6 μm; z_f = 3 burst M/L | incl. | excl. | 70 / 0.3 / 0.7 |
| Bullet | astro-ph/0608407; Table 2 | 0.54, 0.23, 0.58, 0.12 × 10¹² (upper limits) | 100 kpc | nominally Kroupa (via Kauffmann03), inferred | I band, M/L_I = 2 (range 0.5–3) | incl. | not incl. | 70 / 0.3 / 0.7 |
| A520 | 0706.3048, 1202.6368, 1209.2143, 1401.3356 | none; total M/L only (B, R, F814W, i′) | 150 kpc | – | – | – | – | – |
| FP M/L_B rate | astro-ph/0609587 | −0.555 ± 0.042 per unit z → ×0.77, 0.68, 0.47, 0.33 at z = 0.2, 0.3, 0.59, 0.87 | – | – | – | – | – | – |
| FP M/L_B rate | 1009.4479 / 1009.0645 | −0.60 ± 0.04 (fixed σ) / −0.54 ± 0.01 per unit z | – | – | – | – | – | – |
| Colour change | 1009.4479 | Δ(U − V) = −0.24 ± 0.02 mag at z = 0.83 (fixed σ) | – | – | – | – | – | – |
| IMF at σ 250–300 | 0911.3392, 1007.2409, 1208.3523, 1205.6473 | ≈ 1.7–1.9 × Chabrier (range ≈ 1.4–2.3) | – | – | – | – | – | – |

Cosmology is given as H0 / Ωm / ΩΛ.

---

## Not verified, and open issues

- **Ages and truncation.** Neither Ghizzardi nor vdB15 states the FAST age grid. The X-COP deprojection factor of 0.75 has no stated line-of-sight truncation.
- **Release file.** The X-COP release header does not say "projected". That conclusion rests on the 0.75 match. A644 is the one mismatch (12%), which may mean a different profile version.
- **The 0.74 identification.** Bradač's 0.74 is matched to Drory Table 1 by value and redshift; Bradač does not name the bin. The source of the ±0.30 is not stated.
- **Clowe's I band.** Whether it is rest-frame or observed I, and which telescope supplied it, is not stated.
- **Okabe & Umetsu.** The A520 i′ M/L in their Table 8 was not transcribed.
- **IMF conversions vary between papers.**
  - Kroupa to Chabrier: ×0.76 in Ghizzardi (citing Chiu et al. 2018); −0.06 dex in Auger 2010b; about −0.11 dex implied by Hilton's +0.13 dex for Kroupa to Salpeter.
  - Kroupa to Salpeter: about ×1.6 in K-band M/L (Conroy & van Dokkum 2012) and ×2 in Kauffmann 2003.
  - Use a single convention when rescaling.
- **Fading versus normalisation.** The FP rates measure how dynamical M/L fades with redshift. They are not stellar-M/L normalisations. Size and σ evolution and progenitor bias change them at the 0.05–0.1 dex level; see Saglia et al. 2010 and Holden et al. 2010 Sect. 4.2.
