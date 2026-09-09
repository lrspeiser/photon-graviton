# Where captured energy goes

Two explicit, energy-conserving transport examples now separate a halo's shape from its energy supply. Neither supplies a microscopic conversion, capture or binding law. These examples guide T07–T08; their prerequisite field and source-history tasks remain incomplete.

## Findings in plain language

- **Internal light can have a useful shape.** In the point-source approximation, if companions are captured quickly after creation and photon conversion is weak over a galaxy, each additional radial interval receives nearly the same deposited energy. Under the stipulated ordinary mass response, this produces a nearly flat outer circular-speed contribution.
- **The available energy is still insufficient under the archived assumptions.** Fast capture approaches the previous optimistic local-conversion budget; it cannot improve on that budget. Across the 26 matched galaxies its median shortfall remains about 889 million at 10 billion years of present luminosity and the archived conversion rate.
- **External supply has a different shape problem.** Uniform incoming companions absorbed with a constant probability per distance produce approximately uniform deposits when absorption is weak, or an outer layer when absorption is strong. Neither example gives a flat outer deposited-speed contribution. Transport inward, a capture rate varying with environment, or a changed gravity response are open extensions.

No branch is closed by these findings. They identify which assumptions a next candidate must address.

## Internal-source derivation

Assume a point source emits steady photon luminosity L, photons and companions move radially outward, and neither secondary loss nor deposit decay occurs. Let alpha be conversion per distance and beta capture per distance. P(r) and C(r) are outward luminosities; D(r) is cumulative power deposited inside r. At the center P=L and C=D=0.

```
dP/dr = -alpha P
dC/dr = alpha P - beta C
dD/dr = beta C
P + C + D = L
```

For beta unequal to alpha:

```
P/L = exp(-alpha r)
C/L = alpha [exp(-alpha r) - exp(-beta r)] / (beta-alpha)
D/L = 1 - P/L - C/L
```

For equal rates, C/L = alpha r exp(-alpha r). Zero conversion and zero capture are handled separately. The code evaluates short-path deposition with a positive integral to avoid numerical cancellation:

```
D(r)/L = integral from 0 to r:
         alpha exp(-alpha s) [1-exp(-beta (r-s))] ds
```

At radius r the deposition power per volume is q(r)=beta C(r)/(4 pi r²). If deposits stay where formed for an effective common accumulation time T, their density is rhoD=T q/c². This assumes the established steady illumination profile has operated for T; initial light-travel transients are omitted. A different companion speed or source history needs a retarded-time calculation.

In the regime **1/beta much smaller than r, and r much smaller than 1/alpha**:

```
C ≈ (alpha/beta) L
q ≈ alpha L/(4 pi r²)
MD(<r) ≈ T alpha L r / c²
vD² = G MD(<r)/r ≈ G T alpha L/c²
```

This gives an approximately flat deposited contribution. It does not prove total observed speeds are reproduced. In the opposite short-path limit, alpha r and beta r both small, D/L≈alpha beta r²/2, and vD grows approximately as sqrt(r).

In a normalized example with alpha R=10^-5, a capture length of 0.001 R gives vD(R/2)/vD(R)=0.99950; a capture length of 10 R gives 0.71295. The comparison is about shape; normalizing curves to their value at R hides the unresolved amplitude difference, which is reported separately.

The recovered local conversion coefficient is alpha=2.5216780343170007 × 10^-7 per kpc. For each of the 26 archived objects, the beta→infinity limit exactly reproduces its archived local-conversion shortfall within relative tolerance 10^-12. Four illustrative finite capture lengths—0.01, 1, 10 and 100 kpc—are recorded per galaxy. They are not measured or fitted values.

This point source is only appropriate as an illustrative exterior approximation. Actual stars occupy disks and bulges. The central density singularity is not an accepted physical prediction; finite source size, deposit dynamics and the full matter equations must replace it before galaxy fitting.

## External-source derivation

Assume a sphere of radius R receives the same companion specific intensity I from every external direction. Inside it, companions travel straight and are absorbed with constant beta; no scattering, emission, focusing or secondary energy loss is included.

For r≤R, mu the direction cosine and s the backward path to entry:

```
s(r,mu) = r mu + sqrt(R²-r²+r² mu²)
q(r) = 2 pi beta I integral[-1,1] exp[-beta s(r,mu)] dmu
incoming power = 4 pi² R² I
```

An independent incoming-ray calculation gives the total captured fraction:

```
f = integral[0,R] (2 b/R²) [1-exp(-2 beta sqrt(R²-b²))] db
```

Here b is the ray's impact parameter. This agrees with the volume integral of q for all four tested capture depths, with maximum relative numerical difference 8.93 × 10^-6. Escaping power is the incoming power times 1-f; no energy is created.

For beta R much smaller than one, q≈4 pi beta I is uniform, so MD(<r) grows as r³ and vD grows as r. At beta R=0.001, the half-radius speed ratio is 0.49995. At beta R=10, 99.5% of incoming energy is absorbed, but it is concentrated near the surface; the half-radius ratio is only 0.03973.

Efficient capture therefore does not guarantee the desired mass profile. This is a result for the specified uniform external bath and constant absorption, not for all external-supply models. Discrete sources, variable capture, scattering, deposit migration and feedback can change it.

## What is still missing

The plotted velocities include only the stipulated deposited cold-mass contribution. Photons, companions, source fuel, any supporting field, and pressure/stress contributions require a complete gravity calculation. Permanent stationary deposits need a support or binding mechanism: an energy ledger alone does not establish one. The central source's mass changes when fuel is emitted; any additional gravity claim must compare the complete before-and-after distribution without counting the emitted energy twice. Lensing remains uncomputed.

The finite-fuel and escape ledgers from the previous step remain applicable. The two spatial examples do not select a cosmic age or invent extra fuel.

## Concrete next work

1. Preserve the fast local-capture case as a candidate morphology benchmark. Replace the point source with measured disk and bulge emissivity before treating its shape as a galaxy prediction.
2. Compare candidate external-transport mechanisms only after specifying how momentum, energy and captured states are exchanged. Calculate whether they deliver energy to the required radii instead of only increasing total capture.
3. Derive a capture rate from a field interaction or scattering model. Keep all four trial capture lengths labeled illustrative until then.
4. Couple source histories and capture transport to the same finite energy supply. Do not adjust a per-galaxy normalization and call that a photon-energy derivation.
5. Derive the gravitational response and stability of the stored state, then calculate dynamics and lensing from that same state.

A user question is pending about which external-deposition extension to investigate. No changed law has been adopted.

## Evidence

- `radial-capture-results.json`: eight checks, four dimensionless morphology examples, and four capture-length scenarios for each of 26 galaxies.
- `external-capture-results.json`: four independent energy-integral comparisons and spatial diagnostics.
- `radial_capture.py`, `run_radial_capture.py`, `check_external_capture.py`: implementations and reproducible runners using the current local project layout.
- `radial-capture.png`, `external-capture.png`: visually inspected figures.

These equations are derived from the stated transport assumptions. Deriving those assumptions from first principles remains an outstanding task.
