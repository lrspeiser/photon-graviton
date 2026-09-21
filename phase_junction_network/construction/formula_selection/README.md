# Formula-selection result: a conditional nonlinear law, not a unique gravity theory

**Question:** Can the next transport/stress correction select a distinctive gravity formula, rather than add another consistency-only calculation?

**Result:** The single-loop area move has an explicit nonlinear curvature response and a definite higher-curvature signature. But the transported-projector correction alone does not remove the actual geometric-regulator stress defect. Moreover, three finite local-support loop rules have identical weak-field gravitational kernels and different nonlinear responses. The present assumptions and successful tests do not select the original rule uniquely. This is not a proof that future nonlinear constraints cannot select it.

Parent: `b6f55a632e707ede4b872cd0092c85f485b459c5`.
Protocol committed before the executable: `e72c53296622ea25494f88d344c6b8a1a84cd87d`.
Code: `check_formula_selection.py`. Full numerical evidence: `results_full.json`.
All historical results and inherited sector coefficients are unchanged.

## 1. Derive the missing projector stress, rather than subtract the old error

Write the normal-ordered one-body completion of the existing site-block rule as

    Delta_P(A) = sum_a P_a A P_a,
    C_P(H) = Delta_P(H^2) - Delta_P(H)^2/2.

Here P_a selects a complete local site block. This exactly reproduces

    C_aa = O_a^2/2 + sum_(b != a) K_ab K_ba.

The complete derivative is

    delta C = D_H C[delta H] + D_P C[delta P],

where, writing Q=Delta_P(H),

    D_H C[dH] = Delta_P(H dH + dH H) - {Q,Delta_P(dH)}/2,
    d_P Delta_P(A) = sum_a (dP_a A P_a + P_a A dP_a),
    D_P C[dP] = d_P Delta_P(H^2) - {Q,d_P Delta_P(H)}/2.

For simultaneous transport dH=i[D,H], dP_a=i[D,P_a], the sum equals i[D,C]. The unit coefficient of the projector term follows from differentiation; it is not a value fitted to G or a speed.

Tests reconstruct both H and the projectors and compare against independent finite differences. The full run gives a maximum derivative covariance error 3.30e-16 and finite-difference error 7.87e-12. Holding the projectors fixed is a negative control, with error at least 0.332 in the sampled matrices. The inherited one-body completion is recovered to 1.22e-17.

The complete 256-state/eight-mode Fock control keeps the generated quartic operators and their derivatives. Its full derivative covariance error is 1.96e-16. The quartic interaction and quartic stress have nonzero norms, 0.78883136 and 0.42799448 respectively. They are not deleted. This is still a simultaneous-transport identity, not a construction of interacting spacetime constraints or the interacting vacuum.

## 2. Apply it to the actual geometric source: a leading defect survives

The original independently differentiated geometric regulator must not be replaced by a unitary orbit merely to force a Ward identity.

For the existing flat slab write O=H_wall+3 beta, and take p=q/2 along the tested propagation direction. Exact Fourier integration of the coordinate-generator variation gives

    (dH_orbit)_onsite = -i sin(q) beta/2,
    D_H C[dH_orbit] = -i sin(q) {O,beta}/4.

The flat C is site-independent and commutes with the spatial generator. The independently derived moving-projector source is therefore

    D_P C[dP] = +i sin(q) {O,beta}/4.

The old completion source starts at -i q({O,beta}+3I). After the projector term is added, the remaining leading completed source defect is

    -i x^2 q [3{O,beta}/4 + 3I] + O(q^2).

Thus the predicted norm divided by |q| changes from

    x^2 ||{O,beta}+3I|| = 0.188124933199

to

    x^2 ||3{O,beta}/4+3I|| = 0.154872381530.

This is partial cancellation, not restoration of the identity. A separate 108-dimensional, 3^3 graph differentiates the projectors explicitly and verifies the new source against the Fourier expression to 1.18e-16.

| Wavelength parameter L | Original error/|q| | After projector transport |
|---:|---:|---:|
| 16 | 0.18274324 | 0.15058010 |
| 64 | 0.18764630 | 0.15446651 |
| 256 | 0.18809306 | 0.15484543 |
| 1024 | 0.18812293 | 0.15487069 |
| Analytic q -> 0 limit | 0.18812493 | 0.15487238 |

These large-L entries evaluate finite Bloch matrices, not 1024^3 many-body simulations. The tests concern the minimal charge-one, width-two stress symbol; the operator differentiation identity is general. They are not a new six-species interacting loop calculation.

The reason exact simultaneous covariance is insufficient is explicit: the geometric Wilson-regulator variation is not the coordinate commutator on the full microscopic momentum range. C combines products of those microscopic moves, so low external momentum does not erase this mismatch. A full regulator-transport construction remains missing. Quartic vacuum corrections are also not computed here, so this is not a no-go theorem for an interacting completion.

## 3. The conditional single-loop formula

Keep the declared move

    M = B(W-I)/(2i),
    H_M = -x(M+M^dagger) + x^2(M^dagger M+M M^dagger).

For the normalized spin trace, its exact energy is

    E_area = -x tr_n[B Im W] + x^2 tr_n[B^2(I-Re W)].

For the aligned two-dimensional internal block B=b sigma_3, W=exp(i theta sigma_3), this becomes

    E_1(b,theta) = -x b sin(theta) + x^2 b^2 [1-cos(theta)].

The exact curvature response at fixed b is

    dE_1/dtheta = -x b cos(theta) + x^2 b^2 sin(theta),

and the area response at fixed theta is

    dE_1/db = -x sin(theta) + 2 x^2 b [1-cos(theta)].

The expansion is

    E_1 = -x b theta + (x^2 b^2/2) theta^2
          + (x b/6) theta^3 - (x^2 b^2/24) theta^4 + ... .

These are dimensionless lattice energy/source formulas in the existing normalization. theta is the holonomy eigenphase, not a radius; b is a selected area amplitude, not a source mass. This is NOT a derivation of acceleration versus distance, Newton's G, an observed strong-field correction, or a full field equation. Restoring units, source dynamics, field redefinitions and a physical observable dictionary are separate unresolved steps. No claim of literature novelty is made for trigonometric holonomy functions.

## 4. Why the weak-field results do not uniquely select that formula

Consider the modestly enlarged class of real coherent loop amplitudes

    P(W) = sum_(n>=1) a_n (W^n-I),
    M_P = B P(W)/(2i),
    sum_n n a_n = 1.

The normalization holds the original linear-curvature interaction fixed. It is not a speed fit or an independently adjusted gravity coefficient. Repeated loops have the same spatial support but different path length. Their microscopic availability and relative amplitudes have not been derived; they are alternative architectures, not consequences of the original single-loop assumption.

With the SAME completed rule and x, the exact traced energy is

    E_P = -x tr_n[B Im P(W)] + (x^2/2) tr_n[B^2 P(W)P(W)^dagger].

Since P is a polynomial in a unitary W, P commutes with P^dagger, making this trace identity valid even for noncommuting B and W in the finite coordinate/background setting.

Define moments m_j=sum_n n^j a_n. With m_1=1,

    Im P(exp(iF)) = F - m_3 F^3/6 + ...,
    P P^dagger = F^2 + (m_2^2/4-m_3/3)F^4 + ... .

All members therefore have the SAME expansion through total field order two. They inherit the same classical weak-field frame Hessian, its placement-matched two polarizations and the same conditional c_g^2/(U x)=1/8 limit. The latter ratio cannot distinguish the paths.

An independent finite-difference calculation of the complete frame/connection energy checks this statement for random directions in the 30-dimensional field block, not just for the aligned scalar example. Its maximum shared-Hessian directional error is 3.53e-15. The analytic expansion supplies the general quadratic identity; the random tests are regressions, not a proof by sampling.

The finite tests also verify Hermiticity, conjugation covariance, orientation reversal, vanishing flat energy and the same completed-move lower bound for all three alternatives. They do NOT establish full nonlinear gravitational closure, positive interacting phases or common physical poles for any of them.

## 5. A precise nonlinear signature to derive from microscopic path weights

For fixed nonzero b define the lattice-response diagnostics

    R3 = -E'''(0)/E'(0) = m_3,
    R4 = -E''''(0)/E''(0) = 4m_3 - 3m_2^2.

Primes differentiate theta, whose normalization is fixed by the same W in all models. The fixed x and b cancel. The full run obtains:

| Local coherent rule | a_1 | a_2 | R3 | R4 |
|---|---:|---:|---:|---:|
| Single loop | 1 | 0 | 1 | 1 |
| Twice-around loop, same linear normalization | 0 | 1/2 | 4 | 4 |
| Mixture | 1/2 | 1/4 | 2.5 | 3.25 |

For nonnegative a_n, the normalized mathematical weights w_n=n a_n satisfy sum w_n=1, giving the additional identity

    R4-R3 = 3 [m_3-m_2^2] = 3 Var_w(n) >= 0.

The w_n are derivative weights, not asserted Born probabilities. For signed amplitudes the algebraic moment identity survives but the variance/positivity interpretation need not.

This makes the missing microscopic information precise. The first nonlinear curvature response depends on the third path moment, which is invisible to the existing linear/quadratic gravitational tests. Specifying only locality, the common completion coefficient, a flat reference, and the normalized weak-field interaction does not determine that moment in this enlarged class.

Restricting the model to a single, degree-one loop plus reference DOES fix the conditional formula. What remains unproved is why the elementary junction dynamics selects that restriction or predicts the effective repeated-loop amplitudes after eliminating auxiliary states. The alternatives are not claimed physically inequivalent under every possible field redefinition; their lattice responses differ at fixed microscopic W, and an on-shell observable comparison remains downstream.

## Decision

Stop expanding this branch solely to obtain additional consistency passes. The transported-projector term is worth retaining, but it neither repairs the entire stress defect nor selects a unique nonlinear gravity formula.

The next formula-producing calculation would have to derive the path amplitudes a_n and complete regulator transport from a specified elementary junction Hamiltonian, with reference-channel normalization and auxiliary elimination included. Its first held-out target is m_3 (and the paired R4), not another measurement of the already shared weak-field speed ratio. If those amplitudes must simply be chosen, label the resulting law a candidate rather than a unique prediction.

All theory-level closure flags remain false. Historical full multipair output remains byte-identical. Quick/full tests are separate and independently reproduce the frozen analytic/numerical signatures; a green execution result is not a green unique-theory claim.

## Reproduction

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python phase_junction_network/construction/formula_selection/check_formula_selection.py --output /tmp/formula_selection_full.json

Add `--quick` for the smaller regression. The recorded source hashes identify the exact scripts used. Results include software versions. Numerical changes from a later microscopic construction must be versioned, not silently reconciled by editing these old results.
