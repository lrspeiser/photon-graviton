# Collective coherence: normalization, preparation and reverse flow

13 September 2026. Audit of a possible compensation for angular-rate suppression.

## Finding

An initially empty collection of N sites does not receive an automatic N-squared enhancement for its first inelastic excitation. A normalized collective final state gives a rate factor N. At fixed total constituent density, grouping the sites therefore does not automatically restore the rate lost by narrowing the angular kernel. Stronger collective factors are possible for prepared, partly excited symmetric states, but the reverse transition is also enhanced and the preparation energy must be supplied.

This keeps collective mechanisms open while ruling out a normalization shortcut. It does not derive a working companion interaction, a stable deposit or the exact-third reference.

## Known physics and our application

Collective radiative states and enhanced transitions are established quantum-optical concepts: [Dicke, Coherence in Spontaneous Radiation Processes (1954)](https://link.aps.org/doi/10.1103/PhysRev.93.99). The paper does not establish our companion-gravity hypothesis. We apply ideal symmetric two-level-state algebra as a candidate replacement for the constant-matrix-element loading ladder.

Assume N identical two-level constituents with excitation energy delta, phase-matched coupling and no dephasing. Let |i> denote the state with only constituent i excited. Acting on the empty domain gives

\[
S_+|0\rangle=\sum_i|i\rangle=\sqrt N\,|W\rangle,
\quad |W\rangle=N^{-1/2}\sum_i|i\rangle.
\]

Thus the squared matrix element is N, not N^2. The latter would follow from failing to normalize |W>. Summing all orthogonal single-excitation final states also gives N for equal constituent couplings. This statement concerns inelastic excitation; coherent elastic scattering into the same unchanged material state is a different process.

For normalized symmetric states with m excitations, known collective-spin algebra gives

\[
|\langle m+1|S_+|m\rangle|^2=(N-m)(m+1),
\]
\[
|\langle m-1|S_-|m\rangle|^2=m(N-m+1).
\]

At fixed total number of constituents, the upward rate factor per constituent is (N-m)(m+1)/N. It equals one for m=0 and is about N/4 near half filling. That state already contains approximately N delta/2 of excitation energy. These factors describe ideal symmetric states, not arbitrary incoherent mixtures with the same mean excitation.

## Reciprocity remains

Multiplying the earlier paired-mode rates by these matrix factors gives

\[
a_m=A(N-m)(m+1),\qquad
b_m=B m(N-m+1),\quad r=B/A.
\]

Across a particular link, a_m/b_(m+1)=1/r: the collective factor cancels. For a closed finite ladder with maintained field coefficients and no protected export, stationary probabilities therefore obey

\[
P_m=Z^{-1}r^{-m},\quad m=0,...,N.
\]

This is the symmetric manifold's stationary distribution under the stipulated rates; it is not the binomial distribution of N independently thermalized sites. A physical realization would need to preserve that manifold.

| N | r | Mean excitations | Mean upward factor per constituent | Stationary net accumulation |
|---:|---:|---:|---:|---:|
| 10,000 | 0.9 | 9991 | 8.9838 | zero |
| 10,000 | 1 | 5000 | 1667 | zero |
| 10,000 | 1.030301 | 33.0022 | 33.7778 | zero |

Large gross transition rates do not establish large net capture. Initial filling and finite-time transients are not excluded by stationary balance. Sustained export to another protected reservoir could keep a domain partly filled and maintain current, but its transition, reverse channel and capacity must be included. An absorbing endpoint cannot be introduced and then described as a consequence of coherence alone.

## Connection to the angular-rate problem

The previous Gaussian branch lost rate through a shrinking angular cone and longitudinal mismatch. Increasing domain size alone does not compensate that loss for an empty receiving domain at fixed constituent density. A prepared collectively enhanced state might offset some angular suppression under additional assumptions about domain geometry and population, but it does not by itself remove the Gaussian momentum-mismatch factor or derive phase matching.

The next useful model would explicitly couple loading, collective occupation and protected export, with the same matrix factors applied to both directions. It must show where the preparation energy comes from, what stops dephasing, and whether the resulting angular kernel and total rate satisfy the earlier constraints. The lifetime of the prepared state cannot be inferred from its large transition rate; enhanced release can also make it short-lived.

## Verification

Twenty-four full tensor-product number-state checks for N=2,4,6,8 reproduce both collective matrix factors. Twelve larger example states record upward/downward factors and preparation energy. Nine stationary distributions for N=10,100,10000 and three rate ratios verify detailed balance and zero net current. Three explicit normalized bright-state checks verify the first-excitation N scaling. This is established algebra tested in our application, not evidence for gravitons or a new observational fit.

Files: `collective-normalization.py`, `collective-normalization-results.json`. Related: [Gaussian angular-rate closure](kernel-rate-closure-report.md), [mode-driven loading](mode-driven-bias-report.md).
