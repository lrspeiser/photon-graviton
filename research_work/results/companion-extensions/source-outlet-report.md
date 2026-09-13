# Self-consistent source and outlet during collective loading

13 September 2026. Steady driven-loader hypothesis, not a complete deposit theory.

## Result

An externally replenished input field and an escaping output field can maintain loading without prescribing their occupations independently. Slow escape lets the output population build up and reduces throughput. The result still assumes an instantaneous terminal export of each completed packet into protected storage plus relaxation radiation. That terminal process, its inverse and capacity are not derived here. Continuous export/reset is a new loader architecture, distinct from a site permanently holding one protected excitation.

## Equations and assumptions

Reuse the reciprocal paired-mode rates a=Gamma n_H(1+n_L) and b=Gamma n_L(1+n_H), with r=b/a. A ladder of K steps resets to zero after completion. For fixed field populations, its exact steady cycle rate is a/F, where

\[
F=\sum_{j=0}^{K-1}(K-j)r^j.
\]

Thus the net quantum-transfer rate per loader is J=K a/F. The first-passage and renewal mathematics are known; this driven-loader application and terminal sink are hypotheses. This is a stationary closure, not a proof of dynamical stability.

Let m_H and m_L be mode counts per loader, tau_H the source relaxation time and tau_L the outlet escape time. An external reservoir maintains source occupation n_s; the output boundary has zero occupation. Assume uniform populations within each mode set:

\[
\dot n_H=(n_s-n_H)/\tau_H-J/m_H,\qquad
\dot n_L=J/m_L-n_L/\tau_L.
\]

At stationarity define A=Gamma tau_H/m_H, B=Gamma tau_L/m_L and j=J/(Gamma n_s). Then

\[
n_H=n_s(1-Aj),\quad n_L=n_sBj,
\]
\[
j=\frac{K n_H(1+n_L)}{n_s F(K,r)}.
\]

The calculation solves this scalar consistency equation, retaining the finite ladder boundary effects. In a long, forward-biased ladder where boundary effects are small, J approximately equals a-b=Gamma(n_H-n_L), giving the useful approximation

\[
j\simeq\frac1{1+A+B}.
\]

This approximation is not used for the exact roots and need not hold for a nearly balanced or reverse-biased short ladder. Some tested roots have r>1 because the imposed absorbing terminal sink still allows rare completed packets; r>1 is not equivalent to zero loading.

## Numerical examples

Thirty-six cases cover K=20 and 399,951; source occupations 1e-12, 1 and 1000; source loads A=0.1 and 10; and outlet loads B=0.1, 10 and 1000. These are dimensionless sensitivity choices, not measured galaxy parameters.

For K=399,951 and n_s=1e-12:

| A | B | Reverse/forward ratio | J/(Gamma n_s) | Source depletion fraction |
|---:|---:|---:|---:|---:|
| 0.1 | 0.1 | 0.090909 | 0.833334 | 0.083333 |
| 0.1 | 10 | 0.909093 | 0.0900903 | 0.009009 |
| 0.1 | 1000 | 0.999004 | 0.000998904 | 0.000099890 |
| 10 | 0.1 | 0.090909 | 0.0900901 | 0.900901 |
| 10 | 10 | 0.909093 | 0.0476191 | 0.476191 |
| 10 | 1000 | 0.999004 | 0.000989122 | 0.009891 |

Large A permits severe source depletion when loading is efficient. Large B makes the two coupled occupations nearly equal and throttles loading. Neither Gamma, the mode counts nor the physical escape geometry has been derived, so these are not formation times or galaxy predictions.

## Energy and source requirement

Each net step removes E_H from the incoming mode set, adds E_L=E_H-delta to the outgoing set and stores delta in the unfinished ladder. In stationary cycling, half the assembled energy is exported to protected storage and half as terminal relaxation radiation:

\[
P_{in}=E_HJ,\quad P_{out}=E_LJ,\quad
P_{protected}=P_{relax}=\delta J/2.
\]

The ledger closes exactly. Internal forward and reverse traffic is additional to these net boundary flows and is already counted in the ladder rates. Store recoil and a momentum-carrying spatial interaction remain absent.

For E_H=2 eV and delta=1e-8 eV, P_in/(P_protected+P_relax)=2e8. Most processed light therefore exits carrying almost all its original energy; it is not heat. If the prior illustrative inventory/time scale of 2.37014e39 W net assembly were supplied by this one-step optical route, its corresponding processed input would be approximately 4.74e47 W. This is a conditional required luminosity, not a measured supply, an extra emitted luminosity, or a claim about a fixed cosmic age. Reusing photons in successive loading stages could change the incident-light requirement, but would require modeling those stages and their changing spectrum. If the incoming carriers are lower-energy companions instead, their source spectrum must replace this optical example.

## Implications

The previous favorable ratio r=0.9 can arise from source and outlet dynamics in a stipulated open system. This removes the need to choose n_H and n_L independently in that example. It does not derive the external reservoir, outlet geometry or protected sink. An empty boundary is an assumption about a physically distinct mode set, not permission to erase returning waves.

The next physical test is whether the required mode separation and outlet rate can exist without excessive angular scattering or an unobserved spectral component. The very large processed-power ratio must also be compared with the available source field. The traveling-companion connection and exact one-third reference remain unresolved and unchanged, respectively.

## Verification

All 36 roots satisfy the finite-ladder consistency tolerance and the stationary energy ledger. Eighteen K=20 cases are independently checked by solving the stationary cyclic Markov generator; the maximum relative current difference is below 8e-12. A near-unity-r series avoids cancellation in the mean-time expression. This verifies the stationary calculation, not stability, absolute supply or physical realization.

Files: `source-outlet.py`, `source-outlet-results.json`. Predecessor: [occupation-driven bias](mode-driven-bias-report.md).
