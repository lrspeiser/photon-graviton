# Fit a photon-removal revision, then transfer it to more distant supernovae

**Outcome:** one constant removal coefficient improves its training sample but does not improve the frozen central predictions across the more distant sample on the same covariance basis. Its prediction bands widen when parameter uncertainty is propagated. This is a concrete observed-data comparison, not another synthetic receiver check. It does not establish a successful joint redshift/timing/brightness mechanism.

## Proposed revision and known formulas

Fix the previously fitted nearby rate alpha=0.0002488993286382367/Mpc and the timing requirement S_t=S_E=S. Retain static Euclidean geometric dilution. Add the optional whole-photon survival probability

    S=exp(alpha D),  P=exp(-epsilon alpha D)=S^(-epsilon), epsilon>=0
    F_bol=L_bol P/[4 pi D^2 S^2]
    effective brightness distance = D S^(1+epsilon/2).

Fractional attenuation, logarithmic magnitudes and flux accounting are known mathematics. Interpreting epsilon as a new companion-producing process is an unproved hypothesis; no novelty is claimed for an opacity law. **This revision abandons exact photon-number conservation.** It is not silently substituted for the existing baseline. Achromatic whole-photon removal does not itself shift surviving photon frequencies or arrival stages; the original time/redshift mechanism is still needed.

For fixed redshift the additional magnitude is 2.5 epsilon log10(S). If removed photon energy also feeds companions, that energy needs an explicit receiving current and must be added to the energy ledger. Photon survival times per-photon energy survival gives total surviving beam energy fraction S^(-1-epsilon), before arrival-rate dilution. Event stretching must not be counted a second time in integrated energy.

## Observations, calibration, and declared split

Use the already cached, previously exposed Pantheon+SH0ES release and its full STAT+SYS covariance. The [official column definitions](https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/README) identify m_b_corr as a standardized magnitude, CEPH_DIST as a Cepheid host modulus, and explain the shared covariance, including repeated supernovae and calibrator systematics. We use none of MU_SH0ES as independent distance evidence. No dark halo or expansion-derived distance law is adopted in the prediction.

These products nevertheless retain SALT2, selection, host/color, intrinsic-scatter and peculiar-velocity assumptions. Applying the ideal flux scaling to these standardized B-band products is conditional, not an independent bolometric measurement or raw-light-curve fit. Published corrections, including velocity errors derived with a fiducial cosmology, have not been reconstructed under our hypothesis. This limitation cannot be removed by renaming the data model independent.

The protocol was written before running this revision, after the historical data and results were already known. Thus the following is a **retrospective transfer test**, not untouched validation:

- Absolute magnitude: 77 Cepheid-host light-curve rows, 43 distinct supernovae.
- Opacity fit: 466 noncalibrators with 0.1<=zHD<0.3.
- Frozen transfer comparison: 494 noncalibrators with zHD>=0.3.

For calibrators, treat published Cepheid moduli as stipulated geometric distances and include our own predicted attenuation at those distances. Their local measured redshifts are not all assigned to conversion. For the more distant rows, use the inherited conditional convention D=ln(1+zHD)/alpha and S=1+zHEL. This infers distance from our formula; these high-redshift rows are not independent measured-distance checks. Frame and standardization assumptions remain material.

## Fit and per-group results

Generalized least squares gives epsilon=0.48119 with formal conditional standard error 0.10281, equivalent removal rate 0.000119768/Mpc. The fitted absolute magnitude is -19.26740. The model would retain about 71.6% of photons at z=1, conditional on the proposed removal law. These are fitted quantities, not derived physical constants.

Positive residuals mean observed standardized magnitudes are fainter than predicted. Negative values mean the revised model dims them too much.

| Redshift group | Rows | Baseline mean residual, mag | Frozen revised mean residual, mag |
|---|---:|---:|---:|
| 0.1-0.3, fitted | 466 | +0.1249 | +0.0336 |
| 0.3-0.6 | 365 | +0.1869 | +0.0234 |
| 0.6-1.0 | 104 | +0.2048 | -0.0717 |
| 1.0-3.0 | 25 | +0.3497 | -0.0917 |

These are correlated GLS bin means, not independent measurements or separate calibration targets. No parameter was retuned on these latter groups.

Training chi-square decreases from 468.030 to 446.126 after fitting the one coefficient. On all 494 transfer rows, using the same baseline residual covariance for both central predictions, chi-square changes from **424.482 to 432.554**. There is no central-prediction improvement by that comparison, despite smaller mean offsets in some bins. Correlations and within-bin residuals matter.

Propagating fitted-coefficient uncertainty and its cross covariance gives a revised transfer score of 403.728. That uses a different residual covariance and must not be advertised as a direct point-fit improvement or a calibrated model-selection likelihood. We do not assign an extreme significance, add scores from overlapping tests, or claim that this globally excludes the revision. The conditional uncertainty and data-reduction limitations remain substantial.

## Reproducible statistical calculation

For calibrator weights w=C_cc^(-1)1/[1^T C_cc^(-1)1], define K_c=(2.5/ln10) alpha D_c. Then

    M(epsilon)=w^T(m_c-mu_c-2K_c)-epsilon w^T K_c
    m_pred,h(epsilon)=m_pred,h(0)+epsilon X_h
    X_h=2.5 log10(1+zHEL_h)-w^T K_c.

The residual covariance V includes absolute-magnitude calibration variance and both calibrator/evaluation cross terms, as in the earlier analysis. On the declared training rows,

    a=V_tt^(-1)X_t/[X_t^T V_tt^(-1)X_t]
    epsilon_hat=a^T r_t(0),  v=a^T V_tt a
    V_prediction=V_hh+v X_h X_h^T
                 -(V_ht a)X_h^T-X_h(V_ht a)^T.

The fitted value is interior to epsilon>=0; no boundary-normal approximation was needed. Distances, alpha and the published reduction remain fixed, so the formal error does not include all physical uncertainty. opacity.py writes all 960 predictions and input hashes. Historical brightness artifacts are unchanged.

**Decision:** do not adopt constant whole-photon removal as a demonstrated brightness repair. It can improve a calibration range and still overshoot more distant observations. A different physical rule must earn its behavior rather than use per-group adjustments. The unresolved local time/energy mechanism, gravity deposition and final untouched predictions remain open alongside this brightness result.
