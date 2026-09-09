# Observer-time estimator: a failed cadence feasibility gate

## Why test artificial events first?

The proposed photon conversion must eventually predict event durations, not only spectral shifts. Before confronting the real supernova light curves, we need an estimator that does not manufacture or selectively remove duration trends. This pass tests a width estimator using actual DES observation times and error patterns but artificial fluxes with known timing behavior.

The frozen feasibility gate failed. That is a limitation of this first estimator and selection procedure, not evidence against time dilation or against the photon-companion theory. No real-flux timing exponent was fitted, and the failed gate was not relaxed.

## Frozen design and inputs

`protocol.json` was written before running the injections. Input files were checked against the pinned timing-foundation audit. Candidate events were labelled real and spectroscopically Ia in the release. The pass selected at most one band per event whose central wavelength divided by 1+z was within 5% of 450 nm. It required at least eight observations within 80 observer days of the header peak estimate, with at least two before and three after it and positive finite flux errors.

This yielded **98 cadence slots**, with heliocentric redshifts from about 0.0176 to 0.83. The real FLUXCAL column was never used. Artificial flux amplitudes had a declared nominal amplitude/error ratio of 20, while keeping the relative errors and uneven observation times. Header peak dates set only the observer-time origin; the artificial true peak also received an independent offset between -5 and +5 days.

Two artificial shapes were used: an asymmetric Gaussian and the same curve with an added shoulder. Three paired random seeds supplied independent intrinsic log-width scatter and Gaussian measurement noise. Each shape was tested with propagation exponent b=0 and b=1, for twelve runs total. Intrinsic widths had no imposed redshift evolution in these recovery tests.

Actual filter throughput, real brightness selection, classification completeness and source spectral evolution were not validated by the central-wavelength shortcut. This was a preliminary cadence/shape test, not the final inference protocol and not a blind data set.

## Estimator and outcome definition

Each artificial light curve was fit in observer days with a five-parameter exponential/logistic pulse of the Bazin family. This is a phenomenological shape, not an explosion model or a physical conversion law; a primary application describes the function in [Dai et al.](https://academic.oup.com/mnras/article/477/3/4142/4978464). The injected curves deliberately have a different functional form, so this is not just recovery of the fitter's own formula.

The optimizer used three initial rise times and fixed observer-time parameter bounds. No times were divided by 1+z in the estimator. Width was the fitted baseline-subtracted full width at half maximum (FWHM), accepted only if both crossings lay inside the observed time coverage and parameters were away from their bounds. Each run regressed log width on log(1+z), including a free intercept.

The declared requirements were at least 30 eligible events, at least 80% acceptance in every run, and median recovered exponent within 0.15 of the injected value for each shape/truth group. These are feasibility thresholds chosen before the run, not observational tolerances or confidence limits.

## Results

| Artificial shape | Injected b | Accepted events across the three runs | Median recovered b |
| --- | ---: | --- | ---: |
| Asymmetric Gaussian | 0 | 91, 92, 95 of 98 | -0.0947 |
| Asymmetric Gaussian | 1 | 82, 84, 80 of 98 | 1.1431 |
| Gaussian plus shoulder | 0 | 88, 90, 91 of 98 | -0.0732 |
| Gaussian plus shoulder | 1 | 71, 77, 70 of 98 | 1.1282 |

The group medians met the bias threshold, but all three stretched-shoulder runs failed the 80% acceptance requirement. Individual recovered b values also varied widely: the stretched-shoulder runs returned approximately 1.5235, 1.0110 and 1.1282. Three realizations do not establish statistical coverage or a calibrated uncertainty.

For those failing runs, incomplete fitted half-maximum coverage rejected 17, 15 and 20 events respectively; parameter boundaries rejected another 10, 6 and 8. The surviving sample is thus sensitive to duration and shape. This is why fitting only accepted widths would require an explicit treatment of selection before interpreting a redshift slope.

An exploratory diagnostic of the first failing run found acceptance of 5/11, 59/73 and 7/14 in redshift intervals [0,0.2), [0.2,0.5) and [0.5,1). These bins were inspected after the failure to diagnose it; they are not predeclared inference bins or significance tests. The endpoint samples are small and band selection is uneven.

## Intrinsic evolution remains exactly degenerate with propagation in this test

The artificial time scale is

    scale_i = intrinsic_width_i (1+z_i)^(b_propagation+e_intrinsic).

Assigning a given exponent to propagation or to intrinsic duration evolution produces exactly the same noiseless artificial light curve when their sum is the same. The script checks this equality. Widths alone cannot decide which physical explanation generated that trend. A physical source model, external source information or a justified restricted evolution assumption is required; an arbitrary redshift-dependent source duration cannot count as a prediction.

## What must change next

Retain the original gate and all failed trials as history. A revised estimator should account for partially observed light curves rather than silently discarding their width information, and should accommodate shape variation. A joint likelihood on observed fluxes or an explicitly censored width likelihood is a candidate method to test. It must be checked with the same b=0/b=1 controls plus a wider set of noise, peak-offset and shape scenarios under a newly recorded protocol.

The next run must also address wavelength matching and duration-dependent selection before real-data interpretation. Merely keeping the group median near its injected value is insufficient while acceptance varies strongly with injected stretching and shape. No tuning of the physical photon interaction was performed to repair this estimator.

## Reproduction and verification

First acquire the pinned inputs with the timing-foundation command. Then, from the repository root:

```sh
python -m pip install -r research_work/results/timing-injection/requirements.txt
python research_work/results/timing-injection/run_injections.py
```

The script writes to `research_work/generated/timing-injection` by default; `PHOTON_GRAVITON_RESULTS` selects another output base. Saved canonical results preserve every trial's acceptance status and fitted width. A failed scientific feasibility gate is reported in JSON and is not converted into a claim that the estimator passed because the program exited normally.

This optional data-dependent test is separate from the 31-job offline physics suite. It requires the verified DES cache and Astropy. The full theory, timing inference and independent observational validation remain unfinished.
