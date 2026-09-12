# Diagnosis of the first flagged faint-source timing fits

The ongoing 160-case calibration has produced boundary fits. This separate diagnosis freezes the first two flagged cases and profiles their intrinsic duration scatter. It changes neither the running estimator nor its thresholds or retained results.

## What failed

Both optimizers converged from all three starts. In both cases the fitted intrinsic log-width scatter sigma reached the imposed lower limit 0.03, making `interior_solution` false. No numerical exception or missing event caused these flags. The injected population scatter is 0.1.

| Artificial sample | Injected timing b | Fitted b | Fitted sigma | Reason flagged |
| --- | ---: | ---: | ---: | --- |
| Plain asymmetric Gaussian, SNR 5, seed 902 | 1 | 1.03003 | 0.03 | Scatter floor |
| Plain asymmetric Gaussian, SNR 5, seed 903 | 0 | -0.04309 | 0.03 | Scatter floor |

In plain language, these faint samples prefer less natural variation among explosions than the estimator allows. That can arise from sampling, uncertain light curves or model assumptions. This calculation identifies the active constraint; it does not establish which of those is the ultimate cause.

## Executed scatter profile

For each fixed sigma from 0.03 to 0.6, we refit the mean width and timing exponent using three starts. The known profile-likelihood statistic is twice the loss of log likelihood relative to the saved best fit. This is a statistical diagnostic, not a new physical equation.

| Sample | Profile statistic at sigma=0.04 | At injected sigma=0.1 | Refit b at sigma=0.1 |
| --- | ---: | ---: | ---: |
| Stretch, seed 902 | 0.05272 | 2.96124 | 1.02112 |
| No stretch, seed 903 | 0.03074 | 3.23994 | -0.04430 |

The likelihood changes very little between sigma=0.03 and 0.04, and the timing estimate remains similar when scatter is fixed at its injected value. The latter is a diagnostic using known artificial truth, not a correction available for real supernovae. The statistic at sigma=0.1 is not assigned a calibrated confidence level here, particularly with a boundary and possible model mismatch.

At the original timing fits, the nominal 95% likelihood-ratio regions would include the injected b in both cases if the validity screen were ignored. We do not remove that screen after seeing the answers.

## Interpretation of the batch counts

The frozen calibration counts invalid/boundary fits as uncovered. Its reported rate therefore measures the fraction of trials delivering an accepted interval that contains the truth. It combines failures to deliver an accepted interval with intervals that miss the truth. It must not be described simply as evidence that two otherwise valid intervals excluded the answer.

At least one cell has already failed the declared zero-invalid-fits screening condition. Later successful cases cannot undo that failure. The full batch continues as specified to establish frequencies across source shapes, brightness levels and timing rules; these first cases are not grounds for stopping selectively or changing thresholds.

## What a repair would require

A defensible revision may allow zero intrinsic scatter explicitly and calibrate timing uncertainty in the presence of that boundary, or replace source/shape assumptions if broader tests identify bias. Merely lowering the numerical floor is insufficient: when scatter becomes smaller than the width-grid spacing, the discrete population quadrature can itself become inaccurate. A continuous narrow-scatter limit or independently checked integration would be needed.

Neither fixing scatter to the injected 0.1 nor discarding faint cases is a demonstrated repair. Any revised estimator needs a separately frozen repeat of the controls and coverage experiment. Original failures must remain reported. No observed supernova flux has been analyzed here, and the joint redshift/timing/brightness demonstration remains incomplete.

The result file retains every fixed-scatter fit, the input array hashes and a snapshot of the two selected cases. The running batch's checkpoint is not modified.

```powershell
python research_work/results/timing-boundary-diagnosis/run.py
```
