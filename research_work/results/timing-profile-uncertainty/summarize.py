"""Report approximate intervals without claiming calibrated coverage."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
lines=['# Timing uncertainty: nominal profile-likelihood intervals', '',
    'This is an uncertainty-method calculation on the six previously exposed artificial controls and one numerical refinement. No new observations or injected samples are used. The preceding control experiment showed timing separation; this step asks how wide the inferred uncertainty ranges are.', '',
    '## What the intervals mean', '',
    'For each proposed timing exponent b, the calculation refits the average intrinsic duration and intrinsic scatter. The known profile-likelihood statistic is q(b)=2[log L(best fit)-log L(best fit at fixed b)]. We use known chi-square-one thresholds of approximately 1 and 3.84 for nominal 68.27% and 95% ranges. These thresholds rely on asymptotic regularity assumptions; shape mismatch, selection and small samples can make actual coverage differ.', '',
    'In plain terms, the range contains timing rules that fit nearly as well as the best answer after allowing the nuisance parameters to adjust. Calling a range nominal 95% does not establish that it will contain the truth in 95% of realistic repeated surveys. That requires a separate calibration.', '',
    '| Artificial case | Injected b | Fitted b | Nominal 95% region(s) | Contains injected b? |', '| --- | ---: | ---: | --- | --- |']
for case in r['cases']:
    interval=next(v for v in case['intervals'] if v['nominal_level']==.95)
    regions='; '.join(f"[{v['low']:.4f}, {v['high']:.4f}]"+(' (truncated)' if v['lower_truncated'] or v['upper_truncated'] else '') for v in interval['regions'])
    lines.append(f"| {case['label']} | {case['injected_b']} | {case['fitted_b']:.4f} | {regions} | {interval['contains_injected_truth']} |")
base=[v for v in r['cases'] if not v['label'].endswith('-refined')]
for level in [.6826894921370859,.95]:
    count=sum(next(i for i in v['intervals'] if i['nominal_level']==level)['contains_injected_truth'] for v in base)
    lines+=['',f'The nominal {100*level:.2f}% regions contain the injected exponent in {count}/6 exposed base cases. These six cases share three paired seeds, so this fraction is not a reliable independent estimate of coverage.']
lines+=['', '## Numerical verification and limits', '',
    'At every profile evaluation, at least one of two nuisance-parameter starts converged. The likelihood at the saved unrestricted best b agrees within 1e-6; no sampled profile point improves that optimum beyond tolerance. All reported crossing residuals are below 1e-5. Interval roots use Brent refinement after a 101-point scan including the fitted exponent and 0/1 controls. A finite scan does not prove that no narrow disconnected interval was missed. All profile evaluations and boundary truncation flags are retained.', '',
    'The saved fixed-b likelihood ratios compare two timing laws only within this declared statistical model. They are not likelihood ratios for a complete expanding-universe theory versus a complete companion theory. Neither source-population evolution nor survey selection has been marginalized here.', '',
    '## Consequence for the research program', '',
    'The implementation now supplies candidate uncertainty intervals for repeated-injection calibration. A frozen larger ensemble, multiple source shapes and noise levels, and checks of the residual positive shift are still needed before applying those intervals scientifically. This step cannot establish physical redshift causation, calibrated brightness, or completion of the joint light demonstration.', '',
    '## Reproduce', '', '```powershell', 'python research_work/results/timing-profile-uncertainty/run.py', 'python research_work/results/timing-profile-uncertainty/summarize.py', '```','']
(HERE/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
print('Wrote nominal-interval report for '+str(len(r['cases']))+' cases.')
