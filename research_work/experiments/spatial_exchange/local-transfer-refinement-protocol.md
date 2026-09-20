# SE-LR4: combined local transfer refinement and rotated source

Declared 20 September 2026 before execution. Keep the exact SE-LR3 equations,
epsilon2,lambda1,c_Q.5,Omega1,omega_A.2,L8,T1,dt.01 and compact radius1.2
excitation of amplitude.03. Do not tune couplings to a measured outcome.

Run n40,n48,n64 and n48 with a common proper rotation: angle.573 radians
about normalized (1,2,3). Rotate initial Q/P directions and the circular
measurement plane together. The centered radial initial envelope is unchanged.
Use the archived n32 SE-LR3 combined-space run for the previous refinement.
Do not repeat or overwrite it. Save initial/final fields and every-step ledgers.

All earlier individual gates remain: finite states, energy drift<1e-5,
relative total-angular-vector drift<1%, edge energy / initial energy<1e-5,
and128/256 circulation quadrature disagreement<1e-5 at r1,1.5,2. Final
circulation r1.5 and A angular momentum projected on the initial spin axis
must each change<5% from n48 to n64 and less than from n32 to n48. n40 is an
additional trend diagnostic, not an alternate passing pair. Use relative
denominator max(abs(finer),1e-10). Common-rotation changes must each be<2%
using the unrotated n48 value for the denominator. No change to earlier gates.

Additionally report A-sector and total canonical angular momentum projected
on the source axis outside radii1.2,1.5,2 throughout each run. Report fractions
of initial total angular momentum; signed densities are allowed. Hard radial
masks have their own grid error, so no convergence or signal threshold is
declared for this exploratory diagnostic. It measures angular momentum stored
outside the initial excitation, not a directly evaluated surface flux. Q/P
can also propagate outward and generate A locally there: do not attribute all
outer A angular momentum to an A-wave crossing the sphere.

The source remains a finite prescribed internal excitation, not a derived
ordinary-matter population. Numerical transfer does not establish long-lived
swirl, gravity, lensing or observed orbital motion. No dark matter, expansion
or comparison against old gravity predictions enters the acceptance gates.
Established local Hamiltonian calculus and discretization methods are credited.
