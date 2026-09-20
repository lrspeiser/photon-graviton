# SE-T: total scalar-source accounting

Completed states: 18/21. Checks pass: True.

| Campaign/case | Matter source | Other fields | X/Y waves | Total coupling source | Wave fraction |
|---|---:|---:|---:|---:|---:|
| evidence-v1/no-emission | 6.3299013 | 0.02334689 | 0 | 6.3532482 | 0% |
| evidence-v1/no-excitation | 6.3242883 | 0.023307391 | 0 | 6.3475957 | 0% |
| evidence-v1/emission-only | 6.3284291 | 0.023349672 | 0.0022235933 | 6.3540023 | 0.034995% |
| evidence-v1/exchange-0 | 6.3284829 | 0.023313652 | -0.0077077366 | 6.3440889 | -0.12149% |
| evidence-v1/exchange-50 | 6.3284622 | 0.023321197 | -0.0039466418 | 6.3478368 | -0.062173% |
| evidence-v1/exchange-200 | 6.3284338 | 0.023335996 | 0.0017998152 | 6.3535696 | 0.028328% |
| evidence-v1/exchange-negative | 6.3285127 | 0.023304034 | -0.011826464 | 6.3399902 | -0.18654% |
| evidence-v1/sign-mirror | 6.3284338 | 0.023335996 | 0.0017998152 | 6.3535696 | 0.028328% |
| evidence-v1/time-refinement | 6.3284338 | 0.023335998 | 0.0017998146 | 6.3535696 | 0.028328% |
| evidence-v1/space-refinement | 6.3290212 | 0.023002146 | 0.0016032895 | 6.3536266 | 0.025234% |
| evidence-v1/rotation | 6.3286686 | 0.023247689 | 0.0018043335 | 6.3537206 | 0.028398% |
| emitter-v1/emitter-mixed-chi0 | 6.3286562 | 0.023279689 | -0.020370321 | 6.3315656 | -0.32173% |
| emitter-v1/emitter-mixed-chi50 | 6.3285606 | 0.023311447 | 0.00011009762 | 6.3519821 | 0.0017333% |
| emitter-v1/emitter-mixed-chi200 | 6.3284317 | 0.023372331 | 0.029619003 | 6.3814231 | 0.46414% |
| emitter-v1/emitter-Y-chi0 | 6.3286189 | 0.023355922 | 0.0014556288 | 6.3534305 | 0.022911% |
| emitter-v1/emitter-Y-chi50 | 6.3285335 | 0.023382663 | 0.021244154 | 6.3731603 | 0.33334% |
| emitter-v1/emitter-Y-chi200 | 6.3284193 | 0.023432688 | 0.049054909 | 6.4009069 | 0.76637% |
| emitter-v1/no-excitation | 6.3242883 | 0.023307391 | 0 | 6.3475957 | 0% |

These integrated scalar coupling derivatives are not mass, force or lensing.
The phi restoring term is recorded separately in JSON and is not included in
the table total. Field geometry and time evolution must still be solved.
A large wave-sector ratio can coexist with a small share of the total source.
The full energy derivative is independently reconstructed with full-grid
source weights at two perturbation sizes, and checked against the integrated
implemented scalar equation. All raw-state hashes and errors are preserved.
Partial diagnostic success is not campaign completion or observational success.
This evaluates the candidate Hamiltonian, using established differentiation;
older gravity formulas are not acceptance targets.
