# Validation repair and integration audit — 2026-09-20

This is the authoritative review of the earlier Stage 6B, 3E, 3I and 5A claim
boundaries. Numerical histories are retained; an audit PASS is **not** a pass of
the physical hypothesis being tested. Read the individual physical-status fields.

## Status at this checkpoint

| Track | Executed work | Physical conclusion |
|---|---|---|
| Validation baseline | Guarded quick-fit repair; independent full four-pair rerun; fixed-Hamiltonian cutoff comparison | Full finite result reproduced. One-pair and multipair comparisons must hold the Hamiltonian fixed. |
| Chain-to-3D relation | Explicit finite-link embedding and exact leakage on 2^3, 3^3, 4^3 lattices; separate full-real-space Gaussian control | Natural sheet embedding is not invariant. Stage 6B chains are reduced models, not established exact 3D sectors. |
| Shared matter loop | Same six-species domain-wall slab, all mirrors/doublers retained, both photon and frame vertices, several 3D volumes/momenta, temporal derivative and contact terms | Gauge Ward identity passes. Existing frame source leaves nonzero uniform response. No physical common-cone conclusion is licensed. |
| Local nonlinear gravity | Completed-move energy-current algebra, holonomy BCH term, cubic kinetic bootstrap, local constraint classification, bracket refinement and unprojected linear Dirac evolution | Several necessary terms/obstructions are derived. Full finite nonlinear gravitational closure remains open. |

## 1. Fixed-Hamiltonian baseline

The original `microscopic/fermionic_multipair_vacuum_results.json` is retained
byte-for-byte. Its SHA-256 is
`9fc66061e31fe53c4b8a016c5bda9eda8a5f00e0d7efa044314a443a1b5d9400`.
The fitting repair samples polarization at both quick-mode masses and rejects
single-point, repeated-abscissa, nonfinite, nonpositive and mismatched fit inputs.
Quick-mode skipped checks are `null` and explicitly listed as not tested.

`check_fixed_cutoff.py` builds the complete 6,336-state Hamiltonian once and uses

    H_N = Q_N H_full Q_N.

Rebuilding the pair operator before forming its diagonal completion is a different
model. In general,

    Q(B†B + BB†)Q != (QBQ)†(QBQ) + (QBQ)(QBQ)†.

The difference is a positive boundary diagonal, not additional hopping.
The controlled photon-gap results, in the existing energy normalization, are:

| Maximum pairs | Photon gap | Error in the small loop correction |
|---:|---:|---:|
| 1 | 0.0178699485193 | -52.7275% |
| 2 | 0.0178816292423 | -0.100742% |
| 3 | 0.0178816515913 | -0.00004966% |
| 4 | 0.0178816516023 | reference |

These percentages are **not** errors in the total gap. The former 1.3965%
comparison remains a valid comparison of rebuilt cutoff models, not a fixed-H
truncation error. The full model itself is unchanged and reproduced.

## 2. The exact 3D embedding test

For propagation along x define the explicit isometry

    V|a,b> = |E_x=0, E_y(x,y,z)=a_x, E_z(x,y,z)=b_x>,

where a and b are spin-1 chains with zero total flux. These configurations satisfy
all 3D Gauss constraints. The script computes **every** local plaquette move from
every embedded state, hence the complete H V. It does not pretend that the small
one-step support is the whole many-body Hilbert space.

The leakage operator is `(I-VV†) H V`. With plaquette hopping K=1 its norm is
6.92820323, 12.72792206 and 19.59591794 on L=2,3,4. Already the zero-flux column
has squared leakage `6 L^3`; leakage per cell does not vanish. The compressed
Hamiltonian has no chain hopping at this order. This rules out this natural
undressed embedding, **not every possible dressed embedding**.

Separately, the script evolves full three-component fields with the local 3D
Gaussian operator `curl† curl`. Its plane-wave embedding has roundoff-sized
leakage and controlled lattice errors. Equal-norm directions (3,0,0) and (2,2,1)
are compared as the 3D grid is refined. This is a conditional Gaussian-regulator
control; it supplies no missing derivation of a finite-spin interacting phase.
The earlier chain-versus-harmonic-frame speed crossing is therefore a reduced
candidate calibration, not a demonstrated common physical light cone.

## 3. Shared matter response, without a sector-specific adjustment

`check_shared_matter_ir.py` uses the existing `slab_hamiltonian` and the charge
multiplet `(-11,-5,-1,-1,9,9)`. Both responses use the same filled negative-energy
bands, full slab spectrum, energy normalization, spatial grid and boundary
conditions. Quick mode has slab width 4; full mode has width 6 and L=12,16,24.
No mirror, doubler, charge species or occupied band is discarded.

The photon vertex is the Peierls derivative of **both** kinetic and Wilson
hopping. Its second derivative supplies the contact term. The frame vertex is
exactly the existing traceless e_yy/e_zz derivative. A second control uses the
volume-preserving path `e=exp(s P)`, `P=diag(0,1,-1)/sqrt(2)`, including its derived
second-derivative contact. That contact is not an independently tuned coefficient.

The finite spectral sum is the second derivative of the free-fermion determinant
with external sources. It retains the static response, computes the temporal
coefficient analytically from energy denominators, and separately fits spatial
q^2 and q^4 terms. Three-mode versus four-mode fits and a held-out fourth momentum
expose finite-momentum sensitivity. An independent finite difference of the filled
band energy verifies all three uniform-source curvatures.

In the full run the operator gauge Ward residual is below 3e-14 and the
longitudinal static response is below 2e-12. The frame's uniform linear-path
response approaches -10.5672 in the inherited units. The determinant-one path
still gives about -21.6396. These are raw source curvatures, **not measured graviton
masses**; omitted microscopic sectors and geometry vertices must be derived before
physical pole interpretation. Nothing has been subtracted to make them vanish.

A physical speed comparison is explicitly blocked: the interacting photon/frame
kinetic residues have not been matched to this slab, no frame stress Ward identity
has been established, and the finite-spin 3D phase is still unproven. These finite
volumes also do not establish an asymptotic infrared limit; derivative fits must
not be mistaken for a converged continuum speed. The output keeps all raw data so
that further volume, slab-width and infrared-window tests can continue.

## 4. Nonlinear terms and local obstructions

The actual unit-completed moves on a three-state transport block give

    [h01,h12] = x^2 J02 - x^3 (J01+J12),
    Jij = |i><j| - |j><i|.

Thus a longer-range composite energy current is forced by the same rule. A basis
containing only nearest-neighbor currents is not closed. The energy density
contains the complete move including its diagonal partner. Exact unprojected
unitary evolution conserves the total energy and norm in this block. This is an
operator-inventory result, not an identification of these densities with the
complete gravitational constraint generators.

Holonomy multiplication independently gives `[A,B]` at quadratic order. Omitting
it leaves an O(epsilon^2) error; retaining it leaves O(epsilon^3). But a Hermitian
**unitary** completed plaquette alone begins at curvature squared, O(epsilon^4).
It does not generate the needed first-order frame/area coupling by itself.

An affine spatial-density bootstrap starts with the already selected quadratic
kinetic form `T2=tr(pi^2)-(tr pi)^2/2` and solves for all four parity-even ultralocal
h pi pi coefficients. The unique result in that declared basis is

    T3 = 2 tr(h pi^2) - tr(pi) tr(h pi)
         - (tr h) tr(pi^2)/2 + (tr h) (tr pi)^2/4.

The calculation solves a rank-four coefficient system and verifies independent
held-out identities. It assumes that g is a spatial metric and pi a weight-one
contravariant density. **Matching that transformation and interaction to the
finite junction moves remains unestablished.** No Einstein-Hilbert interaction is
inserted as a claimed microscopic answer.

At nonzero momentum, appending local determinant and trace-momentum conditions to
the four old linear constraints produces six independent conditions with bracket
rank four. Therefore there are four second-class and two remaining first-class
conditions, leaving two physical configuration modes. The two added conditions
act as gauge fixings of old directions; one cannot keep all four old constraints
first class and subtract another physical pair. At homogeneous momentum the
separate global pair leaves five shape directions, as before.

Linear Dirac-bracket evolution now preserves those gauge-fixed constraints and
energy without repeated projection. This is explicitly a **linear** evolution
check. For the candidate nonlinear Lie generator and quadratic frame self-energy,
the local central-difference brackets do not close exactly. The vector residual
scales approximately as L^-1.99; a smooth spectral control closes near roundoff.
The scalar self-energy bracket shows the same product-rule issue. The energy of
the eliminated auxiliary connection is already in the Fierz-Pauli gradient term
and is not counted a second time.

## Reproduction and claim gates

Run each executable with an explicit output path. `--quick` is a regression,
not a substitute for the full result. CI runs quick and full jobs separately and
retains their outputs. No empirical parameters or observations were fitted.

The next construction must provide a microscopic frame/area vertex, its complete
matter and regulator stress couplings, and the nonlinear local constraint
algebra. Only then can the shared infrared pole test close. The current tests
are designed to expose those failures, not turn missing requirements into green
physical claims.

## Prior-art boundary

Self-coupling/bootstrap methodology is established prior art, not a claimed
invention here: S. Deser, *Self-Interaction and Gauge Invariance*,
https://arxiv.org/abs/gr-qc/0411023 . The relevant distinction is a microscopic
junction derivation, which remains to be shown. Emergent common limiting speeds
also require a dynamical/renormalization argument; see, for comparison rather
than proof of this model, https://arxiv.org/abs/1510.07650 .
