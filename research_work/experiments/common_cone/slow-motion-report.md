# CC-2S: the direct swirl force is strongly speed-suppressed

Completed20 September2026. [Protocol](slow-motion-protocol.md), [calculation](slow_motion.py), [first results](slow-motion-v1/summary.json).

The declared linear-current model produces a much smaller effect at an illustrative200 km/s than in the earlier0.2c fixtures. With eta/g=1 and kappa=0.5, the directional addition to a prograde star's inward acceleration is only3.08e-8 of the scalar acceleration at twice the source-ring radius. The corresponding photon correction has magnitude4.47e-5 of the scalar bend, with a sign determined by the source orientation and light path. These are calculated weak-field ratios, not observations or full nonlinear predictions.

## What was calculated

A uniformly populated ordinary-matter ring of radius0.7 and rest mass6 sources the linear, unscreened fields. We use g=.001, source speeds0.2c,0.02c,0.002c and200 km/s, both rotation senses, and eta/g=1,10,100,1000. The stellar probe moves along+y at(x,y)=(r,0), with speed equal to the positive source-speed magnitude; source reversal leaves this probe velocity unchanged. Light travels along+x above the ring at y=b. The three radii/impact parameters are2,4,8 ring radii.

There are96 joint star/light fixtures, each repeated with4096 instead of2048 ring samples. All96 refinement comparisons pass; maximum relative change is6.68e-15. Four independent infinite-path photon quadratures pass, with maximum relative disagreement4.04e-16. No case exceeds the declared0.01 amplitude flag for the sampled |U| or |beta|. That flag is not a proof that a dominating force correction is perturbatively accurate.

The current-current interaction, velocity suppression and magnetic-type transverse force are established mechanics, not claimed inventions. They follow here from the already declared CC Hamiltonian. Writing rho_E for source energy density:

```text
U(x) = -g^2 integral rho_E(x')/[4 pi |x-x'|] d^3x'
beta(x) = -(kappa eta)^2 integral rho_E(x') v_source(x')/[4 pi |x-x'|] d^3x'
a_transverse = -grad_perp U - v_probe cross curl(beta)
```

The last equation follows from physical velocity, including its beta contribution, rather than from canonical momentum force alone. For the stellar fixture the velocity is tangential and the requested acceleration radial. The scalar/current ratios cancel the source-energy normalization. The scalar and directional light deflections were integrated over the entire straight unperturbed path, not evaluated at one crossing point.

## Stellar and light response stay linked

At200 km/s with eta/g=1, source sense+1:

| r/R or b/R | Stellar inward directional/scalar ratio | Photon directional/scalar bend ratio | Linear extrapolated eta/g for a unit stellar ratio |
|---|---:|---:|---:|
| 2 | 3.0806e-8 | -4.4689e-5 | 5,697 |
| 4 | 1.4245e-8 | -2.1184e-5 | 8,379 |
| 8 | 6.9951e-9 | -1.0465e-5 | 11,956 |

In this orientation the vector effect strengthens the star's inward acceleration but reduces the light bend on the chosen side of the ring. Reversing the source reverses both corrections. This is a directional prediction, not a uniform extra attractive mass.

The stellar fraction scales as (eta/g)^2 v_source v_star, whereas the photon fraction scales as (eta/g)^2 v_source. Thus the shared coupling cannot freely tune the two effects independently. If the linear formula is extrapolated until the stellar correction reaches unity, the same formula gives an oppositely signed photon correction of magnitude roughly1,451-1,496 times the scalar light bend for these paths. Those coupling ratios lie beyond the declared scanned maximum1000, and nonlinear terms may matter there. These figures identify a scaling tension; they are not valid strong-coupling predictions or a general impossibility theorem.

The earlier0.2c motion was useful for exposing a transient directional effect, but it cannot be presented as a model of galactic stellar speeds. The current ring is prescribed and stationary; this experiment does not establish its source budget or dynamical stability. Its exterior also remains the decreasing multipole field identified in CC-2W, not a flat outer-force profile.

## Consequence for the next candidate

Simply increasing the direct vector coupling is not an established joint solution. A promising alternative to examine is a source-generated swirl that changes the **scalar or spatial response** through a derived interaction: slow matter would then respond through its energy, rather than only through a velocity cross product. Such an interaction needs a Hamiltonian, positive and stable relevant states, a physical source budget, and the shared-cone/local tests. We have not supplied or validated that mechanism here.

This is a research direction, not a claim of novelty: modified susceptibilities, scalar/vector gravity and nonlinear gravitational field laws have precedents. Any candidate matching known MOND-like or polarization mathematics must be identified as such. No dark matter or expanding geometry was inserted, no new astronomical fit was made, and the twelve-item goal remains incomplete.

Protocol8802608 and numerical implementation827d41d precede the saved calculation. Sources and outputs are hashed in slow-motion-v1. The nonlinear CC-2 3D campaign remains separate and its equations were not changed in response to these results.
