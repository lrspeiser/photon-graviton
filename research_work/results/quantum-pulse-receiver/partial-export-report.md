# Fixed partial export preserves coherence but has no stable finite reserve

The full-discharge problem has a limited repair: export a fixed energy d, rather than emptying the receiver. This translates the receiver's coherent state instead of replacing it with a sharp-energy state. But if d equals the mean photon input, the receiver's energy variance grows with use and its finite reserve eventually encounters the lower energy boundary. The simplest constant-release prescription therefore cannot be promoted as a steady physical mechanism.

## Candidate formula and provenance

Retain the stipulated receiver spectrum P_R=E_R-U, with E_R>=U. Propose, on states E_R>=U+d,

    |E_R, E_R-U>_R |vacuum>_c
      -> |E_R-d, E_R-U-d>_R |d,d>_c.

This is a project postulate using known energy/momentum-preserving translation mathematics. Both totals are unchanged. The outgoing companion state is identical for every receiver component, so it does not record the receiver's energy. Inner products and coherence are preserved on the permitted subspace. This is an asymptotic map with no derived interaction rate, locality, or graviton identity. It cannot act this way on energies below U+d: an explicit different boundary operation is required there.

Let X=E_R-U be the available reserve before a photon encounter and L=(1-1/S)E the photon energy transferred in the existing converter. Conversion followed by this export gives

    X_next = X + L - d.

In the ideal interior, the receiver state is a mixture of translated versions of its initial coherent state after outgoing photons are traced. For fresh independent packets with the same energy probabilities, its diagonal energy distribution obeys the known random-walk recursion. Before boundaries matter,

    mean(X_N) = mean(X_0) + N [mean(L)-d]
    Var(X_N) = Var(X_0) + N Var(L).

Setting d=mean(L) balances mean input and output but does not bound the fluctuations. The earlier translation-overlap result preserves single-use coherence in the interior; it does not guarantee support remains in that interior. These probability identities and the recurrence of centered finite-variance one-dimensional lattice walks are known mathematics, not new physics; see [Lalley's random-walk notes](https://www.stat.uchicago.edu/~lalley/Courses/312/RW.pdf).

For the nondegenerate finite symmetric increment distribution used here, a finite initial reserve eventually crosses the lower boundary with probability one in the unlimited-use idealization. This is a statement about this chosen mathematical process, not an assumed age or size of the universe, nor a general exclusion of finite-duration operation.

## Calculation of the allowed operating branch

Use the preceding packet example without fitting any astronomical parameters: S=2; nine photon energies 9.6 through 10.4 with the existing Gaussian-shaped weights; mean transferred energy 5; transfer variance 0.00959016. The initial reserve has mean 5 and standard deviation 0.5. Set d=5. The increments are -0.2 through +0.2 in steps of 0.05.

partial-export.py propagates these discrete probabilities exactly by convolution and removes weight when the prescribed operation would require X<0. Removed weight is an absorbing **failure flag**, not a physical negative-energy state or a repair rule. The upper numerical boundary has negligible escaped weight (3.1e-22 by the last row); neither loss is renormalized away.

| Number of uses | Weight that has not reached the failure boundary | Lower-boundary failure weight |
|---:|---:|---:|
| 100 | 0.999996 | 0.000004 |
| 1,000 | 0.894700 | 0.105300 |
| 10,000 | 0.395681 | 0.604319 |

These are dimensionless mechanism calculations, not lifetimes, distances, observed probabilities or a constraint on the universe's age. They implement an absorbing-history diagnostic; they do not give the unconditional quantum evolution of an unspecified boundary completion. The unconstrained reference variance reaches 96.1516 after 10,000 uses, despite its unchanged mean reserve of 5. Survivor distributions are biased by failure and must not be assigned that unconstrained mean.

## Consequences for a continuously operating model

A larger reserve postpones boundary encounters but does not create a stationary finite distribution for the centered walk. Releasing less than the mean input introduces secular accumulation; releasing more drains the reserve. A zero-variance energy input can avoid this diffusion in the ideal interior but does not represent general broadband packets or varying stellar spectra.

Adaptive release or a coherent exchange with a traveling field could change these conclusions, but requires its own quantum operation. Measuring the exact transferred energy and exporting it generally records which photon energy was present, threatening the coherence used by this particular pulse-stretching construction. That is not a proof that all coherent feedback is impossible. Nor is arbitrary receiver motion or replenishment included here.

**Decision:** fixed partial export is a valid finite-domain bookkeeping map, not a stable completion of the cosmic model. The next physical mechanism must specify state-dependent energy flow and preserve the needed arrival/phase transformation together. We should not fit a galaxy-dependent reserve or silently reset a receiver to make it succeed. No observational calibration or held-out prediction has been added; all six objectives remain open.
