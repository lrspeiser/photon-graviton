# Repeated assembly: mean time and recycled energy

13 September 2026. Extension of the constant-rate ladder; not a galaxy fit.

## Result

A small probability of completing one attempt does not imply that a packet can never form. Allowing repeated attempts gives a finite mean completion time for every finite ladder and positive forward rate. However, the required interaction rate and circulated energy become enormous when reverse transfer dominates. Even the unbiased example processes about 200,000 times its net assembly energy. Recycled energy must not be confused with fresh stellar energy, but losses from recycling require replacement.

## Explicit hypothesis and derivation

Keep the previous ladder with states j=0,...,K, energy j delta, upward rate a and downward rate b. At zero only upward transfer occurs; K absorbs the process and immediately converts half the packet energy into protected storage and half into outgoing radiation. The protected endpoint is stipulated. No protected decay during assembly, finite pump depletion, feedback from returned energy, or delay of terminal relaxation is included.

**Known backward-equation method, applied to our hypothetical ladder:** write t_j for the mean first time to K. Then

\[
-1=a(t_{j+1}-t_j)+b(t_{j-1}-t_j),\quad
-1=a(t_1-t_0),\quad t_K=0.
\]

Setting r=b/a and summing the differences gives

\[
F\equiv at_0=\sum_{j=0}^{K-1}(K-j)r^j
=\frac{K-(K+1)r+r^{K+1}}{(1-r)^2},
\]

with F=K(K+1)/2 when r=1. This is standard birth-death first-passage mathematics, not a new fundamental interaction law. The finite sum independently verifies the closed expression in the small cases.

The expected number of upward transfers is F, because rate a applies in every unfinished state. The number of downward transfers is F-K: net progress is exactly K. Consequently

\[
E_{forward}=F\delta,\quad E_{returned}=(F-K)\delta,\quad
E_{net}=K\delta=E_{protected}+E_{release}.
\]

Define gross/net traffic O=F/K. An externally lost fraction ell of returned energy would add replacement cost ell(O-1) times the net assembly energy, if the ideal rates were maintained. Requiring that cost below the net input gives ell<1/(O-1). This is a loss-budget condition, not proof that those rates can survive imperfect recycling. Returned radiation also needs its own propagation and occupation model.

## Physical benchmark

We import the 24 packet designs from collective-assembly-results.json and evaluate four reverse/forward ratios, giving 96 cases. We separately choose a mean first-fill duration of 10^12 years to display rates and powers. This is neither a cosmic age nor a protection lifetime. Changing it rescales the required forward rate and powers inversely, while O is unchanged.

For K=399,951, delta=1e-8 eV (baseline I, narrow-band packet design):

| Reverse/forward | Gross/net energy traffic | Required forward rate per unfinished site (1/s) | Returned-energy loss fraction ceiling |
|---:|---:|---:|---:|
| 0.9 | approximately 10 | 1.27e-13 | approximately 0.111 |
| 1 | 199,976 | 2.53e-9 | approximately 5.00e-6 |
| 1.000001 | approximately 229,530 | 2.91e-9 | approximately 4.36e-6 |
| 1.030301 | approximately 10^5182.46 | approximately 10^5168.56 | approximately 10^(-5182.46) |

The last row is an effectively unusable constant-rate design at this benchmark, not a proof against every collective mechanism. The ratios are sensitivity assumptions, not measured or derived assembly constants.

If the whole fitted inventory M=2.09261e11 solar masses is built into half-energy protected packets, total net assembly energy is 2Mc^2. Dividing by the illustrative mean duration gives 2.37014e39 W; half ends as stored energy, half as terminal release. This is an inventory/time scale, not a predicted luminosity history or an identified supply. At r=1 the corresponding gross forward traffic scale is about 4.74e44 W, much of it ideally recycled. Completion times are distributed: a mean first-fill time does not guarantee that all sites finish by that time.

## Interpretation and next physical decision

The favorable r=0.9 case shows why a driven assembly mechanism is worth testing: repeated return need not dominate the energy budget if accumulation is biased forward. The bias cannot be declared free. It must arise from a source population, selection rules, or a specified interaction, with its power and reverse channels included. State-dependent rates, coherent transfer and metastability are alternatives outside this constant-rate model.

Neither these energy scales nor recycled traffic establish an available external photon supply. The incoming field may change as it is captured or returned. The next calculation should link a proposed bias to that field and track its depletion and replenishment. The exact one-third empirical reference remains unchanged; the one-half terminal split is specific to this diagnostic.

## Verification

Twelve independent backward-equation matrix solves and finite sums check the analytic mean-time formula to relative tolerance 1e-8. Logarithmic evaluation retains extremely large finite results without silently setting completion probability to zero. All 96 outputs distinguish mean assembly duration, protection-lifetime design, gross traffic and net input. This is mathematical verification of an explicit toy process, not observational validation.

Files: `assembly-recycling.py`, `assembly-recycling-results.json`. Predecessor: [collective release and assembly](collective-assembly-report.md).
