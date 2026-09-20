# TF-1 attribution register
Written with the protocol, before any calculation. Labels follow
[formula-provenance.md](../../../research_plan/formula-provenance.md). Records checked on
19 September 2026 through arXiv abstract pages or the arXiv API; primary PDFs were not re-read
here, and no result is ascribed to unread pages. No external code was copied.

| Component | Status | Source or origin |
|---|---|---|
| Lorentz-force structure a = v x B and its gravitational analogue (gravitoelectromagnetism) | Established mathematics and physics | B. Mashhoon, *Gravitoelectromagnetism: a brief review*, arXiv:gr-qc/0311030 (2003, revised 2008); abstract checked. |
| Gravitomagnetic effects are of order 1e-6 of Newtonian ones for disk galaxy rotation curves and cannot replace unseen mass in general relativity | Established literature result, used only as context | L. Ciotti, arXiv:2207.09736 (2022) and arXiv:2411.03987 (2024); K. Glampedakis and D. I. Jones, arXiv:2303.16679 (2023); A. N. Lasenby, M. P. Hobson and W. E. V. Barker, arXiv:2303.06115 (2023); abstracts checked. TF-1's vector channel is a candidate with a free coefficient, so these magnitude results do not decide it; TF-1 tests its geometry and its light-to-matter ratio instead. |
| Eikonal and Born approximations for ray deflection by a weak index gradient | Established mathematics | Standard geometrical optics; CWC-1 used the same weak-gradient comparison (ray-notes.md, 13da431). |
| Velocity-moment averages of Gaussian distributions through the Laplace identity 1/S = integral exp(-sS) ds | Established mathematics | Standard; the one-dimensional quadrature for f(beta) is written here from it. |
| First velocity moment of the collisionless Boltzmann equation (the Jeans equation) with constant anisotropy | Established mathematics | As used by CR-2 and CL-2 in this repository (results/slacs-resolved-fit/model.py). |
| One-homogeneous velocity terms in a Lagrangian (Finsler-type, Randers-type) | Established mathematics | G. Randers, Phys. Rev. 59, 195 (1941), DOI 10.1103/PhysRev.59.195, is the classical reference for the linear form; the publisher page was not fetched here, so the citation is by record only. The family L = v^2/2 + psi(x)|v|^alpha and its Euler-Lagrange consequences in G6 are derived here; originality unverified. |
| Euler's theorem for homogeneous functions, applied to the unit-speed ray argument | Established mathematics; the conclusion is derived here, originality unverified | Protocol, "A derived limitation". |
| The steering law dn/dt = g_perp P_perp(n) grad phi | Proposed here; originality unverified | Protocol. Direction-only, speed-independent forces exist in other fields (control theory, active matter); no literature search on them was made, and no novelty is claimed. |
| CWC-1's Hamiltonian, lattice, 1D fixture, receiver diagnostics and archived generated 2D field | Existing project result, reused unchanged by import and by hash | research_work/experiments/coupled_conversion (7abb52b to 9aff6b2). |
| CR-2's lens measurement interface and CL-2's LensSystem | Existing project result, reused by import | results/supported-reservoir/cr2.py, results/path-memory/cl2_sources.py. |
| The universal footprint family and CL-2's acceptance thresholds | Existing project declaration | results/path-memory/protocol-cl2.md and amendments. |
| Static Euclidean geometry PF1, D = ln(1+z)/alpha0 | Existing project premise | As RPG-1 and CL-1; the owner's no-expansion constraint. |
| The observed cosmological rate H0/c, about 2.3e-4 per Mpc at H0 = 70 km/s/Mpc | Established arithmetic on an adopted constant | Used only for context in the report; no fixture-to-galaxy mapping. |

Nothing in this stage is claimed as historically novel. Passing a gate establishes a property of the
declared candidate in the declared fixture, not priority and not empirical gravity.
