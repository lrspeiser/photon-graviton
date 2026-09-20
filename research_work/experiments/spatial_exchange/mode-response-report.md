# SE-1M: conversion can change response per unit energy in this candidate

20 September 2026. This is a local calculation from the SE-1 Hamiltonian,
not an observational fit or a complete galaxy/cluster solution.

## Result

All 288 declared fixtures pass. The maximum scaled finite-difference error is
4.02070144e-6 against the declared 1e-5 threshold. A separate eigensolver and
direct matrix-derivative audit passes 580 checks; maximum independent scaled
response error is 1.13686838e-13. Every tested mode has positive energy and
group speed no greater than the candidate's local front speed.

The massless branch has source coefficient S=2 throughout. Of 144 massive
branch fixtures, 72 have S>2, 51 have S<2, and 21 have S=2. Of the suppressed
cases, 44 have negative S. Across all cases S ranges from -108.16297473 to
243.91816818. A negative coefficient is a reversed scalar source at positive
mode energy, not negative mode energy. Neither sign establishes stability.

The first execution failed when serializing NumPy integer counts. Its failure
record and source commit remain preserved. The second execution changes only
serialization and output directory; no equations or thresholds changed.

## What changed relative to the earlier conversion experiment

RC-1 showed reversible population exchange but no enhanced environmental force
for its identical-response symmetric channels. Here the companion's restoring
energy and mixing depend on the field itself. Differentiating that interaction
adds a reciprocal scalar source. With

    kappa = k0 tanh(U/U0)
    nu = mu^2 exp(2 chi U) (1+kappa^2)

the massive local mode gives

    S = 2 + [chi + kappa kappa_U/(1+kappa^2)] nu/(K^2+nu).

S is the derivative of the wave Hamiltonian at fixed canonical wave state,
divided by its energy, averaged over phase. It is not a lensing multiplier.
The interaction can amplify, suppress or reverse the source without assigning
the mode additional energy by hand. This is a concrete candidate change worth
testing; it does not demonstrate an affordable astrophysical field.

## Spatial evidence available at this checkpoint

Six of eleven SE-1 runs are complete: no-emission, no-excitation, emission-only,
exchange-0, exchange-50 and exchange-200. All six pass their short-interval
energy, cone, boundary-clearance and null-control gates. The independent
archive audit passes 105 checks, reconstructing energies within 1.78e-15.
This is a partial audit, not a completed spatial campaign.

Starting with empty fields, finite internal matter excitation produces the X
channel and the coupling produces Y. At time 4, the Y propagation-energy term
is 4.82148e-5 for chi=0, 3.20431e-5 for chi=50, and 1.12235e-5 for chi=200.
Thus increasing the response parameter does not simply increase the generated
companion channel. These terms exclude the separately recorded mixing energy;
they are not conversion fractions or the energies of local eigenmodes.

In exchange-200 the initial internal rest-energy contribution is 0.006 and
ends at 0.00469523. Matter motion also exchanges energy with the other fields.
It would be incorrect to attribute the entire field energy to internal fuel.
No-emission and no-excitation produce zero X/Y; emission-only produces zero Y.

## Remaining decision

Finish sign symmetry, time, space and rotation tests before interpreting the
spatial effect quantitatively. Then determine how much energy the evolving
source actually puts into each local propagation mode, and calculate the
resulting scalar/vector force and full light trajectories from the same
Hamiltonian. A coefficient measured at a prescribed homogeneous U cannot be
inserted into a galaxy: U, emission, composition, recoil and energy depletion
must be solved together. It also depends on wave number K, so propagation and
light-quality consequences need testing rather than assuming achromaticity.

The scalar X/Y channels are not yet derived electromagnetic photons or spin-2
gravitons. Frozen positive mode energies do not prove nonlinear stability,
long-lived swirl, attachment/detachment, a physical source lifetime, flat stellar
rotation or cluster lensing. No observed system has been fitted here.

## Reproduction and attribution

Protocol: mode-response-protocol.md (declared before execution). Source:
mode_response.py at 4e47056; immutable results: mode-response-v2/results.json.
Independent audit: run `python -B audit_mode_response.py` from this directory.
The archive records source/protocol SHA-256 hashes and the source commit.

Hamiltonian dynamics, coupled-wave eigenmodes, finite differences and the
Hellmann-Feynman matrix identity are established mathematics. The particular
mass/mixing functions and matter coupling are our candidate assumptions. Their
checked derivatives and short-time simulations are demonstrated results within
that candidate; originality of the physical mechanism is not established.
Older gravity formulas are not pass/fail targets. The full twelve-item research
goal remains active.
