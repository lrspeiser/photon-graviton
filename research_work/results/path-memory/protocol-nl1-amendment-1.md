# NL-1 amendment 1: the family is linear in the roots of its amplitudes
Declared 19 September 2026 before any calculation. The protocol treats the root spectrum as a family
that is not convex in the amplitudes a_j and therefore claims no global certificate. Writing
s_j = sqrt(a_j) >= 0, the extra force is sum_j s_j sqrt(g_j(x)): linear in s with the nonnegative basis
forces sqrt(g_j) (established algebra). The galaxy loss is then stage 2's convex loss on that basis, the
cluster and lens blocks are linear in s, and the joint objective is convex, so stage 2's certificate
(two solvers, optimality, directional convexity) applies and is adopted. The declared multi-start is
kept and must agree with the certified point (gate G2 becomes: every start within 1e-6 of the certified
optimum). Nothing else changes: the family, its members, the objectives, the balanced joint combination
and the decision rule are as declared.
