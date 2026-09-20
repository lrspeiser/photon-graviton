# SV-1 amendment 1: match initial physical velocities
Declared after the original 252 runs and 36 numerical refinements completed, and after inspection of the pair-heading probes, before examining packet-comparison rankings. The original evidence and source remain unchanged.

The Hamiltonian gives qdot=B(q)p, so identical canonical momenta do not give identical physical emission speeds as eta changes. Retention differences in the original scan may therefore include changed initial speeds. Repeat the identical 252-case grid and 36 refinement cases with matched source-relative physical velocities.

Let p0 be the carrier momenta in the eta=0 packet, P0=-sum p0. The target source-relative velocity is vtarget=p0-P0/M. Set B=I-eta*w/(N-1), and solve (B+11^T/M)p=vtarget independently for each Cartesian component. Set P=-sum p and S=-sum q cross p. This preserves total canonical momentum and angular momentum at initialization and matches all source-relative physical velocities. Absolute source recoil speed and preparation energy can differ and must be reported. This is a diagnostic control, not a new fitted force law.

Store everything separately in evidence-matched-v1 and matched-evidence-sha256.json. Verify matched relative velocities, zero total canonical momentum and compensated angular momentum before integration, at absolute tolerance 1e-11. Keep the original numerical gates, durations, metrics, fixed settings and all failed cases. No observation score, cosmological source law, massless limit or causal propagation is supplied by either scan.

The pair probes already show why the distinction matters: opposite momenta can have a repulsive canonical pair force but an inward instantaneous physical separation acceleration in this particular kinetic model. Do not call the interaction an exclusive same-heading attraction or a demonstrated follower rule merely from the sign of p_i dot p_j.
