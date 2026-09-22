# JR-5: gas-directed energy–companion conversion

**21 September 2026, America/Los_Angeles. Executed toy simulations followed by a real-data screen.**

## Theory first

Matter and radiation may exchange energy with a companion that propagates and populates an extended gravitational state. This candidate lets gas disrupt or enhance conversion differently along different directions. It is an implementation within the broad companion/memory portfolio, not a claim that a local scalar or any existing theory is the only possibility.

**Main result:** the toys can reduce the companion pull in a disk while increasing its projected lensing contribution. The order and clumping of gas also matter, even at fixed total column. However, the tested shared smooth-gas rules do not solve the real sample. No new real-lens or cluster success is claimed.

Baseline: `a5805452b97ef0856d364de90e394572a8fd2146`. [Protocol](PROTOCOL.md) was committed before the toy sweep and fitting at `50b36091098f302fbfeed73c52622194f1db46cc`. Original observations, uncertainties, source properties and R10 parameters remain unchanged. The scientific calculations ran locally, not on GitHub Actions: Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, Numba 0.65.1. Source files are in [code](code); compact evidence is in [results-summary.json](results-summary.json). The originating conversation also contains the full report, all array results, input snapshot, logs and figures.

## 1. Explicit energy exchange and directional transport

Let L be photon-like traveling energy, C traveling companion energy, and B bound companion energy per unit radial ray length:

    L_t + v L_r = -k_plus L + k_minus C
    C_t + v C_r =  k_plus L - k_minus C - kc C + kr B
    B_t         =                         kc C - kr B.

All rates are nonnegative; the trial L and C speeds are equal. Summing cancels internal transfers, leaving only boundary injection and escape. This is a kinetic energy ledger for postulated equal-energy channels, NOT a derivation of recoil, angular momentum, microscopic cross sections, gravitational work or a complete matter/geometry Hamiltonian.

At steady capture balance B=(kc/kr)C. The traveling companion fraction f satisfies

    df/dr = a(r,n)(1-f) - b(r,n)f,
    a = alpha*rho_g/(10^6 Msun/kpc^2),
    b = beta*(rho_g/rho_ref)*rho_g/(10^6 Msun/kpc^2),
    rho_ref = 10^7 Msun/kpc^3 = 0.01 Msun/pc^3.

Linear forward and density-quadratic reverse rates are explicit hypotheses, not established photon/gravity physics. A comparison uses both rates linear. The exact constant-cell update f_next=f_eq+(f-f_eq)exp[-(a+b)dr] keeps f in [0,1] and preserves total traveling energy. Outside gas, the mixture keeps its path history. Gas geometry determines which sequence of rates a ray encounters.

### What remains empirical

The primary source sends an identical f=1/2 mixture into every ray; this is a fixed emitter assumption. The gas-free capture/residence kernel reproduces the previously fitted R10 envelope. Holding that kernel and incoming flux fixed implies

    rho_chi,new(r,n) = 2*f(r,n)*rho_chi,R10(r).

This is conditional on the inherited kernel, NOT a derivation of the old envelope, source luminosity, physical lifetime or adequate energy supply. Equal-speed traveling total energy is unchanged by repartition, while stored energy can change through capture/residence. Merely relabeling unchanged total energy would not create extra gravity. The link to full gravitational stress, photometric source inference and matter backreaction remains open.

## 2. Toy results: same field, different directions

64 combinations compare equal gas mass in a sphere, thin disk, thick disk and ring, with alpha,beta each in {0,0.03,0.3,3}. Gas mass is 10^9 Msun and radial scale 2 kpc. Thin height is 0.2 kpc, thick height 1 kpc. The companion baseline has rc=1, rt=20, q=2; its amplitude cancels in ratios. Gas is prescribed, not hydrodynamically evolved.

The computed three-dimensional axisymmetric companion density sources ONE potential. Even Legendre multipoles give both the forces and the equal-potential weak-field lensing. No separately fitted directional or photon multiplier is used.

For the thin disk with alpha=beta=0.3, at radius/impact 2 kpc:

| Companion-only observable | Change versus gas-free companion |
|---|---:|
| Radial pull in gas-disk plane | **-12.1%** |
| Radial pull along perpendicular axis | **+0.5%** |
| Face-on light deflection | **+3.3%** |
| Edge-on deflection, in-plane impact | **+2.5%** |
| Edge-on deflection, perpendicular impact | **+14.8%** |

These are companion contributions, not total galaxy forces or measured lens angles. Local orbital gradients and path-integrated light deflections respond differently to one spatial distribution.

The ring with alpha=3,beta=0 enhances companion occupation but reduces the inner equatorial pull **8.3%**, increases the polar pull **8.7%**, and increases edge-on/perpendicular deflection **32.3%**. Additional exterior ring material can contribute an outward radial pull at an interior point. A spherical enclosed-mass shortcut misses this directional term. These are toy geometry consequences, not a claim that any specific observed galaxy has this configuration.

## 3. Same gas column, different order

Post-screen attribution controls keep the exact rate law and do not refit or reselect the data model. Take a 2.5-kpc diffuse layer at 0.002 Msun/pc^3 and a 0.1-kpc dense layer at 0.05 Msun/pc^3. Reversing their order preserves both total column (10 Msun/pc^2) and the density-squared path integral.

Starting with purely photon-like energy:

| Order | Outgoing companion energy | Outgoing photon-like energy |
|---|---:|---:|
| Diffuse then dense | **16.67%** | **83.33%** |
| Dense then diffuse | **72.31%** | **27.69%** |

The energy sum is unchanged. Local transfer matrices do not commute: the last strong layer relaxes the mixture toward its own equilibrium. For the primary 50/50 input the outputs are 16.674% and 72.314%, a factor 4.34 apart.

A clumping control at the same column changes f_out from 0.5 in a homogeneous one-kpc layer to 0.16667 in a dense 0.2-kpc clump plus 0.8 kpc of empty path. A diffuse five-kpc layer gives 0.82423. Unlike the order reversal, that comparison changes the density-squared integral. Neither total mass nor average column fully specifies this proposed interaction. Actual clumping/layer order in SPARC was not measured here.

## 4. Source-fed finite-volume evolution

Two initially empty, three-channel rays include capture, release and escape. Input is one model energy unit per time unit; no conversion to physical galaxy age is claimed. At time 400, the near-plane ray has 144.9573 bound, 15.3797 traveling and 239.6630 escaped energy; the polar ray has 111.3870, 15.2594 and 273.3537. Each sums to 400 supplied units. Maximum absolute ledger discrepancy is 6.6e-11 and no energy becomes negative.

The two bound populations reach about 91% and 87% of their steady values by time 400, not complete finite-age equilibration. Turning off a source from a stationary preparation causes bound energy to decline: near-plane 159.33 -> 96.80 -> 53.49 at elapsed times 0,100,200; polar 128.58 -> 85.28 -> 52.28. Permanent storage is not inserted.

The dynamic toy uses a bounded trial residence kernel. The observational bridge uses the inherited R10 kernel. They share kinetics and energy structure but are not one measured microscopic capture spectrum. Force/lens numbers use stationary readouts, not the still-filling time-400 state.

## 5. Real-data screen: partial improvements, no overall solution

All 149 archived SPARC galaxies and 3,152 positive-radius rows are retained. Universal coefficients are fitted on 89 galaxies, the null and four gas variants selected using 29-role mean equal-galaxy fractional-squared error, and the 31 comparison galaxies scored after selection. Sixteen optimizer starts are retained. All systems were historically examined: this is development, not fresh blind confirmation.

The input package supplies gas masses and gas-only rotation contributions, not resolved 3-D gas cubes. Exponential radial gas scales were reconstructed from the gas-only force contributions at fixed gas mass, without using total observed-rotation residuals. Height is provisionally 0.1 gas scales. The reconstruction has a median **32.2% normalized RMS error in the gas-only squared-speed contribution**; **122/149** exceed 20%, one scale hits its bound. The actual baryonic gravity still uses the published gas component. Rings, holes, clumps and real directional structures are absent from this reconstruction; that mismatch limits the screen.

| Model | All-galaxy mean velocity RMS | Median fractional RMS | Below 20% |
|---|---:|---:|---:|
| Original R10 | **17.71 km/s** | **15.05%** | **95/149** |
| Selected gas-return-only | 18.82 km/s | 15.63% | 86/149 |
| Two-way density-dependent, post-selection comparison | 19.03 km/s | 18.33% | 88/149 |

Return-only is alpha=0,beta=0.0110049186852. It wins the declared selection objective (validation mean fractional-square 0.05847 -> 0.04195 and RMSE 17.68 -> 16.74 km/s), but rescues only two previous >20% outliers while spoiling eleven previous <20% galaxies overall. This pure-return limiting control is not a claim that fundamental reversible physics has no forward channel.

The two-way fit has alpha=1.44265,beta=30, with beta on its upper bound. It reduces some large outliers, rescues **18**, but spoils **25** previously <20% curves. Its full-sample comparison was performed after selection and does not retroactively change the selection.

Examples under those SAME two-way coefficients:

| Galaxy | Original RMS | Two-way RMS |
|---|---:|---:|
| UGC07559 | 51.8% | **3.6%** |
| KK98-251 | 69.7% | **7.5%** |
| IC2574 | 57.6% | 31.2% |
| UGC07125 | 55.8% | 46.5% |
| UGC07399 | 32.6% | **46.7%, worse** |
| NGC2403 | 15.4% | **22.0%, worse** |

For the 31 comparison galaxies, mean velocity RMS is 15.09,15.32,14.72 km/s for original, selected and two-way respectively. Below-20% counts are 22,20,23. An improvement on this one mean is not an overall or statistical validation.

Angular attribution: the selected model's monopole-only approximation gives essentially unchanged physical RMSE and 87 rather than 86 curves below 20%. Directional feasibility in toys has not become a successful law in the rotation averages.

Gas-height sensitivity at FROZEN coefficients: half-height gives 79/149 below 20%, double-height 94/149, versus primary 86. No height was selected from these results. Even the 27 source reconstructions within 20% do not show a general rescue: mean velocity RMS is 18.47 baseline, 20.92 selected, 19.34 two-way. Poor gas reconstruction cannot be asserted to explain the entire failure.

## 6. What was not tested

The supplied six SLACS lenses have no usable gas maps. Setting their gas to zero and obtaining unchanged results would not validate this mechanism. No new real-lens or cluster score is claimed. Attempts to retrieve the THINGS map page yielded no usable maps in this run. Real directional validation needs resolved gas geometry, uncertain deprojection/height/clumping and distributed source locations. No hidden data edits, synthetic gas maps labeled observed, or new per-object directions were used.

Source energy supply, capture/residence derivation, gas recoil/momentum, microscopic coupling, photometric self-consistency, hydrodynamics, star formation and nested planetary behavior remain open. The equal-potential readout is an explicit trial coupling, not a derived relativistic completion.

## 7. Numerical reliability

Cell transfer vs independent ODE: fraction error 2.43e-11. Layer transfer vs independent matrix exponentials: <=2.3e-15; energy error <=2.3e-16. Toy resolution doubling changes representative disk equatorial/polar forces by 0.105%/0.474%, lens deflection by <0.083%.

Independent direct-volume/azimuth-integrated forces agree with high-resolution multipoles to 0.0345% (disk) and 0.00081% (ring); direct self-refinement is 0.119%/0.00091%, not a rigorous universal error bound. Independent face-on cylinder-mass lensing agrees within 0.028%.

At the selected data coefficients, doubling radial/angular resolution and multipole order changes velocities by at most **0.259%** across all 149 galaxies; none exceeds 1%. The conclusions/counts remain unchanged. Readback reconstructs **2,235** saved error/chi-square quantities, zero discrepancies. These are numerical checks, not validation of a physical theory.

An initial metadata-output attempt failed on a NumPy boolean before any fitting. The log is retained; only JSON serialization was fixed.

## Reproduce

Extract the original JR-1 input package from the full JR-5 conversation archive to obtain `Photon-Graviton-JR1`. With a fresh output directory:

    python code/toys.py --output new-results
    python code/data_screen.py --jr1 Photon-Graviton-JR1 --output new-results
    python code/audit.py --jr1 Photon-Graviton-JR1 --results new-results
    python code/layers.py --output new-results/layer-order.json
    python code/supplement.py --jr1 Photon-Graviton-JR1 --results new-results

The pipeline preserves prior files. Source hashes are recorded in the compact summary. The full report, raw input snapshot, all results and plots are in the originating conversation package; they are not implied to be permanently stored in Actions.

## Decision

Retain gas-directed conversion as a concrete mechanism worth developing: one field can reduce inner orbital pull and increase some projected lensing, and it can retain information about gas clumping and layer order. Do not promote the tested smooth-gas/central-emitter/R10-kernel approximation as a solution to all observations. The next physical construction is transport through measured structured gas with distributed emission and a derived capture/return kernel, not another global weakening factor or per-object lensing multiplier.

## Primary sources and project provenance

- Project JR-1 and JR-4 packages, unchanged inputs and source hashes; Alternative Gravity Branch Portfolio.
- Lelli, McGaugh & Schombert (2016), SPARC: https://arxiv.org/abs/1606.09251 ; official data portal https://astroweb.cwru.edu/SPARC/ . No fitted dark-matter halo products from the portal were used.
- Walter et al. (2008), THINGS: https://arxiv.org/abs/0810.2125 . Context for resolved HI data, not a new map dataset analyzed here.
