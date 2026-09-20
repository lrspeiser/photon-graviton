# Candidate change after the CC-2 control campaign

20 September2026. This is an analytic design note, not an executed protocol or a claim that the requested solution exists. Finish and audit CC-2 first; preserve that candidate as the control. The twelve-item objective is unchanged.

## Evidence that motivates a change

- CC-2W's compact linear exterior has declining circular speeds and a lapse-only static deflection coefficient. Changing coupling amplitudes does not fix either structural issue.
- CC-2S's direct current contribution to stellar acceleration scales with both source and probe speed. Its photon contribution has only the source-speed factor. A strong effect in a0.2c toy is not evidence of a useful galactic effect.
- CC-1/CC-2 demonstrate how to derive a common local propagation cone and reciprocal source accounting. Preserve those properties rather than inserting independent stellar and optical corrections.

## Spatial response that modifies every sector together

A family to derive and test, using U as a dimensionless scalar response and s as a spatial-response parameter, is

```text
ds^2 = -exp(2U) dt^2 + exp(-2sU) |dx-beta dt|^2
C = exp((1+s)U), beta = C b, |b|<kappa<1
H_particle = exp(U) sqrt(m^2+exp(2sU) p^2) + beta dot p
H_field density = exp((1+3s)U) Pi^2/2
                  + exp((1-s)U) |grad F|^2/2
                  - beta dot sum_a Pi_a grad F_a + V(F)
```

The product of the field kinetic and gradient coefficients is C^2. Consequently the frozen principal cone is centered on beta with radius C, matching the photon Hamiltonian. Positive V and |b|<1 give a positive local field kinetic/gradient block. Differentiating all coefficients is essential for energy accounting.

For weak stationary U, slow matter measures grad U while light measures (1+s)grad U. Thus s=1 supplies the standard equal temporal/spatial response. This is established metric optics, not a new optical invention or proof that s=1 is dynamically required. It must be derived or declared as a constitutive assumption and tested against local observations. The field potential/constraint structure, matter regularization and nonlinear stability still require work. A conformally modified optical metric alone does not flatten an outer force law.

## A nonlinear amplitude route, with explicit reasons it may fail

One candidate hypothesis is U=g phi-lambda |A|^2/2 with lambda>0. This lets a generated vector amplitude influence the scalar response felt by slow matter. In a homogeneous, resting-matter approximation with phi held at zero and an effective positive mode stiffness Omega^2, the potential is

```text
V_local(A) = Omega^2 |A|^2/2 + rho exp(-lambda |A|^2/2)
t = lambda rho / Omega^2
V_local''(0) = Omega^2 (1-t)
|A_*|^2 = 2 log(t)/lambda                  if t>1
U_* = -log(t)
V_local(A_*) = (Omega^2/lambda) [1+log(t)]
```

These are algebraic stationarity results for a restricted local potential, not a derived galactic phase. At t>1 a small directional seed could grow without its final amplitude being proportional to source speed. Exact zero fields remain zero in an exactly symmetric deterministic calculation unless the matter current or another modeled perturbation provides a seed. Energy must come from the matter Hamiltonian; a damped or radiating state needs its lost energy recorded.

The simple formula also exposes difficulties before we invest in it. A weak response |U_*| of order1e-6 requires t to lie about1e-6 above threshold in this local approximation. It does not explain how different galaxies and clusters obtain that condition with universal constants. Nonzero |A| is not necessarily circulation: a uniform vector can have zero curl. Outside compact matter the linear tail remains a separate obstacle. The local stationary amplitude does not prove full-mode stability, long-lived vortices, or acceptable clocks and rulers.

This broad mechanism has substantial prior art: [Damour and Esposito-Farese's spontaneous scalarization](https://arxiv.org/abs/gr-qc/9602056) and [Ramazanoglu's spontaneous vectorization](https://arxiv.org/abs/1706.01056). The exponential matter coupling and threshold-growth idea must not be presented as an original theory invented here. We are not adopting those papers' neutron-star theories or inheriting their results; any examination of this branch is an attributed adaptation within our candidate, and may be rejected. No historical uniqueness is claimed.

## Decision requirements for the next protocol

1. Complete the CC-2 audit, including failures and source-radius sensitivity.
2. Declare a precise spatial-response extension and its full Hamiltonian derivatives. Test the common cone, particle convexity and the discrete energy ledger before evolution.
3. If amplitude feedback is retained, compare lambda=0 controls with both sides of the declared growth threshold. Require nonzero measured circulation, not merely vector amplitude. Include realistic dimensionless source speeds and a source-energy/angular-momentum budget.
4. Reject the branch as a joint solution if it merely needs object-by-object threshold tuning, still produces the wrong outer force, or cannot satisfy local light/matter tests. Preserve it as an unsuccessful hypothesis instead of renaming a known theory as our own.

No numerical implementation or observational success of this proposed extension is claimed in this note. It is not permission to skip the original12 requirements or the pending coupled3D validation.

Completion update,20 September2026: CC-2 is now finished and independently audited; all31 evolutions and11 comparisons pass. Its source-size sensitivity and OB-1 boundary failures remain. Proceed from that archived control, not from an assumption that the full physical solution is validated.
