# Phase Junction Network

**Status:** exploratory candidate theory; internally checked at the linear electromagnetic and gravitational levels, but not a completed fundamental theory or an empirical fit.  
**Created:** 2026-09-20  
**Purpose:** preserve the rewritten Phase Junction proposal, the corrected gauge-field construction, the candidate frame-gravity branch, and reproducible zero-data checks in one self-contained research module.

## Core idea

A local matter–geometry phase difference is not itself light. For a smooth scalar phase,

\[
A_\mu=\frac{\hbar}{q}\partial_\mu\phi
\quad\Rightarrow\quad
F_{\mu\nu}=0,
\]

so the original scalar-gradient proposal is pure gauge. The corrected model moves the physical comparison to **links between neighboring junctions**:

- each link carries a compact phase comparison \(a_\ell\) and conjugate matter–geometry imbalance \(E_\ell\);
- electromagnetic field strength is the non-closing phase accumulated around a loop;
- photons are the quantized transverse normal modes of this constrained link network;
- gravity requires a richer frame-valued comparison, represented at low energy by a symmetric tensor with four first-class constraints and two transverse-traceless modes.

The resulting working statement is:

> Reality is modeled as a network of microscopic matter–geometry junctions. Electromagnetism is scalar phase holonomy on the network. Gravity is frame-and-clock holonomy on the same network. Matter appears as charged or massive network defects and endpoints, while the network constraints enforce charge and stress-energy conservation.

## What is currently established inside the candidate model

The electromagnetic construction supplies, under the stated lattice/rotor assumptions:

\[
H_A=\frac{U_A}{2}\sum_\ell E_\ell^2-K_A\sum_p\cos B_p,
\]

\[
c_\gamma=\frac{\ell\sqrt{U_AK_A}}{\hbar},
\qquad
\alpha_{\rm bare}=\frac{1}{4\pi}\sqrt{\frac{U_A}{K_A}},
\]

with two transverse modes, a long-distance Coulomb kernel, and \(E_n=\hbar\omega(n+1/2)\) after canonical quantization.

The candidate gravitational branch supplies, at the linear structural level:

- four first-class constraints acting on a symmetric spatial frame deformation;
- exactly two propagating tensor modes;
- positive transverse-traceless quadratic energy;
- linear long-wavelength dispersion;
- a synthetic lattice Green function approaching \(1/r\);
- a route to universal coupling through one shared frame;
- a common low-energy photon/gravity causal cone if both sectors use that frame.

## What is not established

This folder does **not** yet provide:

- a finite-dimensional microscopic junction Hamiltonian that derives both effective sectors;
- a derivation of the numerical impedance ratio \(Z_g/Z_A\), \(G\), or \(\alpha\) from fewer microscopic inputs;
- a derivation of fermions, chirality, generations, charge assignments, or particle masses;
- nonlinear strong-field solutions or a proof of quantum consistency;
- a solution to the vacuum-volume/cosmological-constant term;
- a demonstrated microscopic evasion of the Weinberg–Witten assumptions;
- a frozen, overconstrained empirical prediction that distinguishes the model from ordinary low-energy QED plus general relativity.

The frame-gravity branch deliberately lands in the same low-energy tensor universality class as linearized Einstein gravity. The proposed novelty is the microscopic origin of the gauge connection, frame constraints, coupling ratios, defects, and any calculable departures above the junction scale.

## Files

| File | Purpose |
|---|---|
| [`electromagnetic_derivation.md`](electromagnetic_derivation.md) | Reconstructs the original scalar proposal as a link/loop gauge network, derives its spectrum, quantization, Coulomb limit, and QED route, and separates the scalar sine-Gordon mode from the photon. |
| [`frame_gravity_derivation.md`](frame_gravity_derivation.md) | Defines the candidate frame-valued gravitational sector, its four constraints, two tensor modes, static limit, universal coupling, and parameter relations. |
| [`check_frame_gravity.py`](check_frame_gravity.py) | Reproducible no-data numerical checks for constraint counting, the tensor spectrum, the excluded negative scalar branch, and the lattice \(1/r\) Green function. |
| [`checks.json`](checks.json) | Frozen output of the current verification run. |
| [`validation_protocol.md`](validation_protocol.md) | Gates for the next microscopic derivation and the later transition to real-data tests. |
| [`manifest.json`](manifest.json) | Machine-readable scope, commands, assumptions, and result summary. |

## Reproduce the current checks

From the repository root, using the repository requirements:

```sh
pip install -r requirements.txt
python phase_junction_network/check_frame_gravity.py \
  --samples 500 \
  --seed 42 \
  --lattice-size 96 \
  --fit-r-min 4 \
  --fit-r-max 20
```

The current frozen run reports:

- spectrum failures: `0 / 500`;
- maximum relative transverse-tensor eigenvalue error: `1.1250857213591084e-15`;
- maximum constraint/gauge residual: `1.0870458088261133e-14`;
- unconstrained transverse spectrum: one negative scalar and two positive tensor modes;
- normalized RMS residual of the fitted lattice \(A/r+B\) kernel over radii 4–20 on a \(96^3\) grid: `0.002409597456672864`.

These results verify the implemented linear algebra and synthetic lattice behavior. They do not establish that nature uses the proposed microscopic interpretation.

## Immediate next task

Construct the smallest finite junction Hilbert space and local move set that simultaneously generates:

1. the electromagnetic Gauss constraint and plaquette ring exchange;
2. the three gravitational momentum constraints and one curvature/energy constraint;
3. one shared frame for all matter excitations;
4. a calculable ratio \(Z_g/Z_A\);
5. a protected low-energy matter-defect spectrum.

The model should be rejected or revised before observational fitting if that finite construction produces extra gapless scalar/vector modes, ghosts, unstable gradients, nonuniversal matter coupling, or independently tunable photon and gravity cones.
