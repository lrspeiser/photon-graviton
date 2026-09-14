# Self-illumination pilot: rotating emitters shift the momentum debt but do not remove it

13 September 2026. Protocol: [protocol.md](protocol.md), drafted before execution, with execution notes appended. This is a diagnostic of the archived RB-1 reaction under companion fields derived from declared rotating emitters. Nothing is fitted, and no co-rotation factor is assigned.

## Outcome in plain language

Stars that emit companions isotropically in their own frames create a field that partly moves with the rotating matter.

- **Receiver drag falls.** In a disk with a flat rotation curve, receivers feel κ = 0.19–0.51, against 1.742 for the external isotropic bath. Inside a rigidly rotating ring the drag is exactly zero.
- **Emitters pay instead.** Each unit of emitted energy removes (E/c^2) times the emitter's own r v. Angular momentum is conserved overall, so in the most favourable case, where every emitted companion is absorbed inside the galaxy, emitters and receivers together lose exactly what the new particles carry away. That is (absorbed energy/c^2) times the receivers' r v.
- **The floor still fails.** At the reference inventory and ideal efficiency, this floor is a median **19 times the baryons' own angular momentum** (10th to 90th percentile 3.8 to 72). It exceeds 1 in 146 of 149 galaxies; the Milky Way fiducials give 2.9 and 3.2. Self-illumination lowers the external-bath debt (median 33.7) by about 44% but does not relieve it.
- **Efficiency does not improve.** Differential rotation leaves enough relative Doppler shift that the threshold-cut spectrum still loses half of its slow products. The (v_esc/c)^3 penalty is untouched.

**For the tested RB-1 reaction, self-illumination does not supply the required momentum budget.** Other interactions are not tested here.

## Physics and provenance

Transport is steady, optically thin and along straight rays. The lab emissivity of a steady emitting flow is j_N(E) = D j'_N(E/D), with D = 1/(γ_e(1 − n·β_e)) and n the propagation direction from emitter to receiver. A uniformly moving source distribution instead has the boosted-bath form of the RB-1 revision. Isotropic emission in the emitter frame removes momentum at β_e times the emitted lab power. The reaction is the RB-1 heavy-receiver threshold process, applied beam by beam with exact kinematics and both spectrum controls. The transfer relations are established physics; applying them to companions is our hypothesis.

## Controls, all passing

- **Emission normalization.** Integrating one moving element over the sphere returns its number rate, power and momentum to 10^-15.
- **Stationary emitting shell.** Reproduces the isotropic drag: 1.74224 for the extended spectrum, and 1.73827 for the threshold-cut spectrum with its edge correction.
- **Steady co-moving flow.** A fixed shell with emitters flowing at velocity w obeys κ = (4/3 + χ/3) − (w/β)(1 + χ/3). This gives exactly 1/3 at w = β and 1.0378 at w = β/2, matched to 10^-4. A steady structure is therefore not a co-moving bath: co-rotating matter cannot give zero drag unless the whole configuration is rigid.
- **Rigid ring.** Receivers co-rotating inside a rigidly rotating ring feel zero torque to 10^-6.

## Results

**Rigid ring** (β = 0.01 and 200 km/s):
- *Co-rotating receivers inside the ring (0.5 R and 0.95 R):* κ = 0, receiver energy change 0, and bound efficiency factor 1 for both spectra. Co-rotating pairs see no Doppler shift, so the threshold edge is not smeared.
- *Exterior receiver (1.5 R):* κ = 0.556, with threshold-cut efficiency 0.95.

**Thin exponential disk** (flat v_c = 200 km/s, softening h = 0.1 R_d, emitters to 12 R_d). Angular momentum is per unit absorbed energy, in kpc·c with R_d = 1 kpc.

| R/R_d | κ_self | Threshold-cut efficiency | Receiver L_z change | Products' L_z | Emitter/receiver specific L_z |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.435 | 0.500 | −1.45×10^-4 | 3.34×10^-4 | 4.0 |
| 1 | 0.255 | 0.500 | −1.70×10^-4 | 6.67×10^-4 | 2.0 |
| 2 | 0.189 | 0.500 | −2.52×10^-4 | 1.33×10^-3 | 1.0 |
| 4 | 0.222 | 0.500 | −5.92×10^-4 | 2.67×10^-3 | 0.5 |
| 8 | 0.507 | 0.500 | −2.71×10^-3 | 5.34×10^-3 | 0.25 |

The extended-spectrum efficiency is 1.000 at every radius. Softening sensitivity is mild: κ = 0.227, 0.255 and 0.303 at 1 R_d for h = 0.05, 0.1 and 0.2, and 0.217, 0.222 and 0.233 at 4 R_d. Receiver energy changes are about ±10^-7 per absorbed energy, far below the external-bath −κβ^2 ≈ −7.7×10^-7. Receivers still lose angular momentum everywhere, because incoming companions bring less than the products carry off.

## Conservation floor across the galaxy sample

The receivers, sites and inventories are those of the RB-1 revision. All values are at ideal efficiency, with every companion absorbed internally.

| Quantity | Value |
|---|---|
| Combined emitter-plus-receiver loss / baryons' angular momentum | median 18.96 (10th to 90th percentile 3.76 to 72.5; min 0.29) |
| Galaxies above 1 | 146 of 149 |
| External-bath debt, for comparison | median 33.68 |
| Floor / external | median 0.565 |
| Milky Way fiducials, floor / external | 2.91 / 5.11 and 3.23 / 5.67 |

At the threshold-limited supply, both numbers are multiplied by the supply multiplier of 4–8×10^9.

## Decision

The protocol's rule requires the combined loss at ideal efficiency to fall below the baryons' own angular momentum. It fails in 146 of 149 galaxies, so **self-illumination does not relieve the RB-1 debt.** The lower receiver torque does not count on its own, and threshold efficiency is unchanged.

## What this means for the next candidate

For RB-1 the obstruction is momentum conservation, not the illumination geometry. When products are born moving with ordinary matter, the baryons must supply the angular momentum of the retained reservoir or of the escaping flux. A viable channel would need two things:

1. retained products born without drawing their motion from the baryons, for example a non-rotating receiver population with an adequate momentum budget, or a collective field mechanism;
2. a retained fraction far above (v_esc/c)^3.

These are requirements for the next candidate, not a proof that none exists.

## Reproduction

```sh
python research_work/results/self-illumination/pilot.py    # about 4 min on 8 workers; --canonical writes pilot-results.json here
```

`selfillum.py` holds the steady-flow beam transport, the beam-resolved RB-1 channels, bound production, emitter losses and the emission check. It reuses `capture-to-orbit/incident.py` and `kinematics.py`. The v1.5 PDF predates this result.
