# Timing uncertainty: nominal profile-likelihood intervals

This is an uncertainty-method calculation on the six previously exposed artificial controls and one numerical refinement. No new observations or injected samples are used. The preceding control experiment showed timing separation; this step asks how wide the inferred uncertainty ranges are.

## What the intervals mean

For each proposed timing exponent b, the calculation refits the average intrinsic duration and intrinsic scatter. The known profile-likelihood statistic is q(b)=2[log L(best fit)-log L(best fit at fixed b)]. We use known chi-square-one thresholds of approximately 1 and 3.84 for nominal 68.27% and 95% ranges. These thresholds rely on asymptotic regularity assumptions; shape mismatch, selection and small samples can make actual coverage differ.

In plain terms, the range contains timing rules that fit nearly as well as the best answer after allowing the nuisance parameters to adjust. Calling a range nominal 95% does not establish that it will contain the truth in 95% of realistic repeated surveys. That requires a separate calibration.

| Artificial case | Injected b | Fitted b | Nominal 95% region(s) | Contains injected b? |
| --- | ---: | ---: | --- | --- |
| b0-seed701 | 0 | 0.1125 | [-0.0608, 0.2856] | True |
| b0-seed702 | 0 | -0.1126 | [-0.3270, 0.0990] | True |
| b0-seed703 | 0 | -0.0061 | [-0.1723, 0.1595] | True |
| b1-seed701 | 1 | 1.1097 | [0.9270, 1.2934] | True |
| b1-seed702 | 1 | 0.9099 | [0.7035, 1.1143] | True |
| b1-seed703 | 1 | 1.0221 | [0.8466, 1.1985] | True |
| b0-seed702-refined | 0 | -0.1113 | [-0.3248, 0.0997] | True |

The nominal 68.27% regions contain the injected exponent in 3/6 exposed base cases. These six cases share three paired seeds, so this fraction is not a reliable independent estimate of coverage.

The nominal 95.00% regions contain the injected exponent in 6/6 exposed base cases. These six cases share three paired seeds, so this fraction is not a reliable independent estimate of coverage.

## Numerical verification and limits

At every profile evaluation, at least one of two nuisance-parameter starts converged. The likelihood at the saved unrestricted best b agrees within 1e-6; no sampled profile point improves that optimum beyond tolerance. All reported crossing residuals are below 1e-5. Interval roots use Brent refinement after a 101-point scan including the fitted exponent and 0/1 controls. A finite scan does not prove that no narrow disconnected interval was missed. All profile evaluations and boundary truncation flags are retained.

The saved fixed-b likelihood ratios compare two timing laws only within this declared statistical model. They are not likelihood ratios for a complete expanding-universe theory versus a complete companion theory. Neither source-population evolution nor survey selection has been marginalized here.

## Consequence for the research program

The implementation now supplies candidate uncertainty intervals for repeated-injection calibration. A frozen larger ensemble, multiple source shapes and noise levels, and checks of the residual positive shift are still needed before applying those intervals scientifically. This step cannot establish physical redshift causation, calibrated brightness, or completion of the joint light demonstration.

## Reproduce

```powershell
python research_work/results/timing-profile-uncertainty/run.py
python research_work/results/timing-profile-uncertainty/summarize.py
```
