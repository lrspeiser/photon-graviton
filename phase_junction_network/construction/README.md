# Area-reference interaction and complete candidate stress coupling

**Status:** finite operator candidate constructed; classical quadratic frame gate passes with local placement matched. Full common-cone and nonlinear quantum-gravity gates remain open. The candidate Gaussian coordinate-stress gate fails. An executable diagnostic passing is not a successful physical-theory claim.

This follows `../validation/CHECKPOINT.md` and the preregistered `PROTOCOL.md`. Historical JSON is unchanged. The inherited `x=0.13554178509861228` and, only for the classical frame oscillator, `U=0.006789913753821769` are not adjusted to match a new speed.

## 1. First-order area/curvature move

Let B be the Hermitian frame bivector and W the oriented spin holonomy at the same plaquette base. The new candidate is

    M_B = B(W-I)/(2i),
    H[M] = -x(M+M†) + x²(M†M+MM†).

Two equal coherent path amplitudes give the block `(W-I)/2`; an area-controlled quadrature supplies `B/i`. An explicit finite two-path circuit block is checked. This makes the difference a microscopic amplitude, not a subtraction of measured vacuum energy. It does **not** derive a protected low-energy ancilla sector. Equal path weights and the relative phase are additional architectural assumptions, not consequences of the completed-move rule alone.

For classical frame coefficients or commuting coordinate registers, the normalized auxiliary-spin trace is exactly

    H_area = -x tr_normalized[B Im(W)]
             + x² tr_normalized[B²(I-Re(W))].

Thus W=exp(iF) gives the first-order term `-x tr(BF)`. The complete quadratic-curvature partner is retained. At W=I the move and energy vanish for every B. A no-reference control instead has a nonzero frame-dependent flat energy.

The factorization

    H[M] = (Q†Q+QQ†-2I)/4, Q=I-2xM,

bounds each finite move below by -I/2; it is not a healthy-continuum theorem. Tests cover Hermiticity, conjugation covariance, orientation reversal of the traced density, first curvature derivative, and a 108-dimensional finite locked-coordinate realization with unprojected evolution. That representative is not a fully quantized 3D geometry.

For the inverse frame e, every sector uses

    nu = 1/det(e),
    Gamma^i = sum_a e_a^i gamma^a,
    g^ij = sum_a e_a^i e_a^j,
    B^ij = nu[Gamma^i,Gamma^j]/(2i).

Site frames and link spin transports are established techniques, not claimed inventions.

## 2. Classical quadratic kernel derived from the new move

Expand the actual B(W-I) operator at flat geometry, retaining six symmetric metric amplitudes in each real Fourier quadrature and all eighteen real connection amplitudes. The expansion is

    Im(W) = sum_r A_r + (i/2)sum_{r<s}[A_r,A_s] + O(A³),
    I-Re(W) = (sum_r A_r)²/2 + O(A³).

The four oriented insertions A_r retain the connection commutator and x² completion. No Fierz-Pauli kernel or transverse-traceless projector is inserted. Exact finite-amplitude evaluation verifies a fourth-order even Taylor remainder.

The connection Hessian has rank eighteen and is eliminated algebraically. A Fourier block of a translation-invariant quadratic expansion is exact at that order; it is **not** the rejected many-body chain embedding. Eight trigonometric nodes integrate the quadratic volume average exactly. Tests cover L=8,12,16,24,32,48,64,96,128 and directions (1,0,0), (1,1,0), (1,2,1).

### Placement and constraints

The old site-centered generator fails in mixed directions. The new coframe calculation has the backward-difference vector symbol

    d_i* = khat_i exp(-ik_i/2).

Its null directions are independently checked against the derived kernel. Tensor placement `(i+j)/2` aligns that generator with the real-khat form. The inherited kinetic metric must use the same placement. Only diagonal trace cross terms change:

    trace_at_vertex(pi) = sum_i pi_ii(x+hat_i).

This is a local integer one-link shift, not a fitted sector factor. It remains an explicitly geometrically matched kinetic candidate with inherited DeWitt coefficients, not a derivation of every finite quantum kinetic move.

With that placement, eight real linear constraints are first class, their Hamiltonian evolution closes, and unprojected linear evolution preserves them. Actual constraint reduction leaves four real configuration modes: **two polarizations per signed nonzero momentum**, with positive springs. The unreduced Hessian has two negative conformal, six gauge-zero and four positive directions. No physical ghost-free nonlinear claim follows from this linear result.

Both polarizations and all three directions extrapolate toward

    c_g²/(U x) = 1/8

in the declared normalization. The largest full-run fitted departure from 1/8 is approximately 1.1e-6. This is not a universal observed constant. Momentum-dependent splits and fit-window changes are retained. The result is a conditional classical tensor dispersion, not quantum photon/gravity speed matching.

## 3. Geometry and regulator coupling

Covariant shifts T_i carry U(1) phases and Spin(3) transport. Let D_i=T_i-I and K_i=(T_i-T_i†)/(2i). In canonically normalized site fermions:

    H_geom = H_wall + nu^(-1/2)[
        1/2 sum_i {nu Gamma^i,K_i}
        + beta/2 sum_ij D_i†(nu g^ij)D_j
    ]nu^(-1/2).

The full Wilson metric, mixed components, volume factors, spin connection, wall masses and fifth-direction hopping are included. Flat geometry recovers the earlier slab spectrum. Wilson strength and localization/slab data are inherited parameters, not predictions. Spin connections are independent first-order backgrounds; their full matter-sourced compatibility equation is not solved.

Frame stress, second-derivative contacts, charge current, spin-link current and symmetric lapse/energy density are derivatives of the complete candidate. The symmetric lapse assignment includes all terms but is not a constructed nonlinear Hamiltonian constraint. Covariance tests reconstruct the changed sources rather than defining the answer by conjugation.

### Fock completion and quartic terms

Decompose H_geom into on-site O_a and unordered-pair K_ab blocks. The candidate elementary Fock moves are

    M_aa = -c_a† O_a c_a/2,
    M_ab = -c_a† K_ab c_b, a<b.

The same completed rule applies to every move. This grouping is another explicit architectural choice, not a uniqueness or microscopic low-band matching result.

Completing a single-particle matrix and then second-quantizing misses interactions:

    (c†Ac)(c†Bc) = c†ABc + sum A_ij B_kl c_i† c_k† c_l c_j.

The full definition retains the generated four-fermion operators and their stress derivatives. A complete eight-mode/256-state Fock test verifies the one-particle restriction, nonzero quartic remainder, spin-source derivatives and energy conservation without repeated projection.

The normal-ordered quadratic part is

    H_quad = x H_geom + x² C,
    C_aa = O_a²/2 + sum_{b != a} K_ab K_ba.

Only this Gaussian part enters the large-volume one-loop calculation. Quartic corrections and frame/connection loops have **not** been calculated there.

## 4. New matter-source measurements

Charge and local-spin covariance hold near roundoff. The spin Ward bubble and contact cancel. Varying the frame while omitting required spin-link variations fails the negative control. A coordinate-unitary orbit is retained only as a reference, not as proof that physical frame deformations are gauge.

Both photon and frame responses use all six inherited species `(-11,-5,-1,-1,9,9)`, all mirrors, doublers and slab states, and the same normalization. Full mode uses width6 and L=12,16,24; quick uses width2 and L=8,12. Filling intervals are checked rather than assuming zero chemical potential after a completion shifts the spectrum.

The raw uniform shear response approaches -8.83821 in this new candidate's units. At L=24:

| Contribution | Value |
|---|---:|
| Matter bubble | -2.43670 |
| Kinetic contact | -1.50062 |
| Wilson/metric contact | -35.69334 |
| Completed-move contact | +30.79245 |
| **Total** | **-8.83821** |

These are source curvatures, not graviton masses. The Hamiltonian and normalization differ from the old uncompleted slab; a smaller raw number is not evidence of resolution. No contribution is discarded. The area interaction vanishes on the uniform flat connection.

### Leading coordinate-stress obstruction

An independent geometric coordinate source is compared with the half-density generator `D_xi={xi,K_x}/2`. Canonical measure and Wilson variations are included, and the analytic source is checked against a finite graph derivative.

For external q=2pi/L and internal p=(pi/L,pi/L,pi/L), the uncompleted geometry target has vanishing low-momentum discretization error. The direct-block completed Gaussian source instead has

    ||delta H-i[D_xi,H]||/|q| -> x² ||{O,beta}+3I||.

For the width2, unit-charge diagnostic this limit is about **0.1881**. It is derived from the completion source, not fitted. Tests reach L=256 for this small matrix-symbol calculation; this is **not** a 256³ interacting simulation.

The direct-block Gaussian candidate therefore fails the leading coordinate-stress requirement. Refining the wavelength does not remove the normalized leading defect. Internal spin/gauge covariance alone is insufficient. This does not rule out interacting corrections or a different microscopic transport/projector completion; no compensating term was inserted.

## 5. What can and cannot be matched

The physical pole gate remains false: full coordinate stress, nonlinear scalar/lapse generators and a stationary joint quantum background are not established; the finite-spin 3D photon pole is still unavailable. Temporal coefficients are Euclidean derivatives, not real-time speeds. The code retains raw constants and separately reports spatial q²/q⁴ terms and fit-window effects.

Established at declared scope: area/reference operator, derived classical quadratic frame kernel, local placement-matched linear constraints, complete candidate source definitions, exact small Fock completion, shared Gaussian response and a precise leading stress defect.

Not established: full nonlinear constraint preservation, healthy nonlinear quantum spectrum, common quantum cone, a unique microscopic model, observed matter, vacuum cancellation or empirical agreement.

The next construction must derive the regulator's transport/projector and interaction stress from complete local moves, not pick an independent coefficient to cancel the measured response. Preserve the successful area/linear tests and pass the now-explicit coordinate-stress gate.

## Reproduction

```sh
python phase_junction_network/construction/check_area_stress.py --output /tmp/area_stress_full.json
python phase_junction_network/construction/check_area_kernel.py --output /tmp/area_kernel_full.json
python phase_junction_network/construction/check_geometry_response.py --output /tmp/geometry_response_full.json
```

Each accepts `--quick`; CI keeps modes separate. Frozen outputs include source/environment hashes. Full mode means the declared finite configuration, not a complete theory.

## Prior-art boundary

- R. C. Brower et al., *Lattice Dirac Fermions on a Simplicial Riemannian Manifold*, arXiv:1610.08587; Phys. Rev. D 95, 114510 (2017). Site vierbeins and link spin transports are established methodology. https://arxiv.org/abs/1610.08587
- S. Deser, *Self-Interaction and Gauge Invariance*, Gen. Rel. Grav. 1, 9–18 (1970), archived as gr-qc/0411023. First-order self-coupling is prior art, not a microscopic derivation for this candidate. https://arxiv.org/abs/gr-qc/0411023

No originality is claimed for those component techniques. The new work here is the declared candidate and its explicit retained successes and failures.
