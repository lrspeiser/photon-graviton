# Collective release and reversible packet assembly

13 September 2026. Conditional extension; no new observational fit.

## What changes

Combining many small incoming transfers into a larger excitation could ease the earlier outgoing-mode bottleneck, provided relaxation releases a larger quantum. Mere grouping does not help: for fixed released energy, unchanged outgoing quantum energy, bandwidth and escape geometry, the number of emitted quanta and their mode occupation are unchanged.

We test a distinct hypothesis: K transfers of energy delta build one bright excitation K delta. Half is retained in a protected state and half exits as one collective relaxation quantum. This assumes a permitted collective transition and its reverse; it is not derived from the earlier single-transfer interaction.

## Conditional scaling

The preceding finite-site result has lifetime bound proportional to M/(R^2 E_a^4) at fixed fractional bandwidth, protected fraction and relaxation fraction. Replacing E_a by K delta gives

\[
L_K=L_1/K^4,\qquad
K_{min}=\left\lceil(L_1/T)^{1/4}\right\rceil.
\]

This is an algebraic consequence of that hypothesis, not a universal storage law: site number falls as 1/K and outgoing mode number rises as K^3. It retains the homogeneous occupation and escape assumptions. T is a selected lifetime benchmark, not the age of the universe. Matching a necessary bound does not guarantee actual stability or adequate pumping. The protected decay here removes an entire aggregate; fractional leakage would require a different model.

Energy bookkeeping for a completed packet is K delta incoming, K delta/2 retained, K delta/2 released. Formation rates, momentum recoil and reverse collective excitation remain to be specified. A seed's pre-existing mass is not newly supplied energy.

For the fitted baseline-I inventory in a 120 kpc sphere, 90% protected occupation and T=10^12 years:

| Small transfer (eV) | Log bandwidth | Contributions needed | Aggregate bright energy (eV) |
|---:|---:|---:|---:|
| 1e-08 | 0.01 | 399,951 | 0.00399951 |
| 1e-08 | 1 | 115,872 | 0.00115872 |
| 2e-10 | 0.01 | 19,997,532 | 0.0039995064 |
| 2e-10 | 1 | 5,793,555 | 0.001158711 |

At fixed other inputs, the required aggregate energy is independent of the chosen small transfer energy except for integer rounding. Smaller individual transfers require more assembly steps. The broad-band case remains favorable mode counting, not proof that a sharp transition couples to the whole band.

## Can assembly finish?

As a separate toy model, let a partly assembled packet contain j units, with constant upward rate a and downward rate b. States 0 and K end the attempt; start at j=1. A downward step returns delta to the supplying channel rather than destroying it. At K the collective relaxation is assumed to occur immediately. Its actual competition with reverse transitions is omitted, making this an optimistic assembly endpoint.

**Known mathematics:** the first-passage probability satisfies

\[
a(P_{j+1}-P_j)+b(P_{j-1}-P_j)=0,\quad P_0=0,\quad P_K=1,
\]
\[
P_1=\frac{1-r}{1-r^K}\quad(r=b/a\ne1),\qquad P_1=1/K\quad(r=1).
\]

This is the classical gambler's-ruin result; see [Cambridge probability notes, random walk and gambler's ruin](https://www.statslab.cam.ac.uk/~rrw1/prob/prob-weber.pdf). Applying it to companion assembly is our hypothesis. These constant rates do not automatically follow from the previous finite-site or bosonic occupation models.

For K=399,951 (the narrow-band 1e-8 eV example):

| Reverse/forward rate | Completion probability per seeded attempt |
|---:|---:|
| 0.9 | approximately 0.1 |
| 1 | approximately 2.50e-6 |
| 1.000001 | approximately 2.03e-6 |
| 1.030301 | approximately 10^(-5186.53) |

The last ratio is borrowed only as an illustrative sensitivity value from the earlier reciprocal-rate example. We have not derived that ratio for assembly. A tiny upward bias can change the result dramatically; its physical origin and energy source must be identified. The eventual success probability is not a formation rate, elapsed time or energy efficiency. Repeated attempts, recycling, seed counts and absolute rates matter. With infinitely many attempts a nonzero probability need not mean permanent failure.

## What this establishes and what comes next

A higher-energy release channel provides a concrete way to reduce mode crowding without making each photon transfer large. It moves the unresolved work to the nonlinear accumulation and relaxation mechanism. Simply assigning a large packet energy does not supply it.

The next discriminating task is a driven ladder with state-dependent rates or a metastable collective transition, derived from an interaction or explicit energy landscape. Track reverse transfers, pump power and the release spectrum together. A retained fraction of one-half in this diagnostic is not the empirical one-third retention exponent; that reference is unchanged. No result here improves a galaxy or lensing fit by itself.

## Verification

Twenty-four aggregate designs cover two inventories, two small-transfer energies, two bandwidths and three free lifetime benchmarks (10^9, 10^12, 10^15 years). Each checks the smallest integer K meeting the conditional bound, with 96 assembly probabilities. Twelve independent finite Markov-chain linear solves agree with the analytic first-passage expression to relative tolerance 1e-10. Aggregate energy invariance under changing the small transfer is checked to the integer-rounding allowance. These validate the calculations, not the proposed microscopic physics.

Executable: `collective-assembly.py`; output: `collective-assembly-results.json`. Source inventory and bound: [site/mode budget](site-mode-budget-report.md).
