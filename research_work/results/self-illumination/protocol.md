# Self-illumination pilot: does a field from rotating emitters relieve the receiver momentum debt?

Drafted before execution, 13 September 2026, after the RB-1 consistency revision. This is a diagnostic of the RB-1 reaction under a derived, not assigned, illumination. It is not a new capture branch, and nothing is fitted.

## Question

Can the companion field produced by a declared rotating source geometry reduce the total receiver torque enough to build the reference inventory, with every conserved quantity tracked in one account? The account covers the emitters, the receiving matter, bound products, escaping products and unabsorbed companions.

## Declared physics

- **Emitters.** Each emits companions isotropically in its own rest frame, with a flat number spectrum per unit rest-frame energy. Both spectrum controls apply: S1 cut at threshold, S2 extending below it.
- **Transport.** Steady, optically thin, straight rays in flat space. A steady emitting flow has lab emissivity j_N(E) = D j'_N(E/D), with D = 1/(gamma_e(1 - n.beta_e)) and n the propagation direction from emitter to receiver. This is established relativistic transfer; retardation does not enter a steady structure.
- **Receivers and reaction.** The same heavy-receiver RB-1 reaction, with the same s-wave law and exact kinematics. Each emitter beam is evaluated as a direction-resolved contribution. No co-rotation factor is assigned; the field's frame follows from the emitters' motion.
- **Emitter loss.** Isotropic emission in the emitter frame removes momentum at (emitted lab power / c^2) v_e. Unabsorbed companions carry their own momentum away.

## Geometries

1. **Controlled rigid ring** of radius R at speed v, with co-rotating receivers at 0.5 R, 0.95 R and 1.5 R. Run at beta = 0.01 (for resolution) and at 200 km/s.
2. **Declared thin disk.** Exponential emitter surface density with scale R_d and a flat rotation curve v_c = 200 km/s. The disk is a uniform slab of half-thickness h = 0.05, 0.1 or 0.2 R_d; the thickness softens the near field. Receivers sit in the midplane at 0.5–8 R_d, co-rotating at v_c. Emitters extend to 12 R_d.

## Reported quantities, per receiver

- **Drag.** The effective coefficient kappa_self = -F_t/(beta P_abs), set beside the isotropic external value of 1.742.
- **Retained efficiency.** The bound-production factor relative to the isotropic S2 uniform ball.
- **Energy.** The receiver energy change per absorbed energy.
- **Angular momentum:**
  - emitter loss per unit emitted energy;
  - receiver change, product and incident angular momentum per unit absorbed energy.

## Transfer to the galaxy sample and decision rule

The most favourable accounting absorbs every emitted companion inside the galaxy. Angular momentum conservation then makes the combined emitter-plus-receiver loss equal to the angular momentum the products carry away. Transfer that floor to each galaxy's reference inventory at ideal efficiency, using the same sites and receivers as the RB-1 revision.

Self-illumination relieves the debt only if the combined loss at ideal efficiency falls below the baryons' own angular momentum, without lowering retained efficiency. A lower receiver torque alone does not count, and reduced drag does not fix threshold efficiency.

## Execution notes, recorded after running

1. **Controls added.** A stationary emitting shell must reproduce the isotropic value. A fixed shell through which emitters flow with velocity w (a steady flow) has the analytic first-order drag kappa = (4/3 + chi/3) - (w/beta)(1 + chi/3), which is exactly 1/3 when the emitters co-move with the receiver. A uniformly moving source distribution is not a steady structure, and is equivalent instead to the boosted isotropic bath of the RB-1 revision, whose drag vanishes when co-moving. Galaxies are steady structures, so the steady-flow controls apply.
2. **Transfer computed directly.** The galaxy-sample transfer evaluates the decision rule from the conservation floor above, not from receiver torque.
