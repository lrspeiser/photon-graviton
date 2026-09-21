# Coherent junction: solve the parent, then eliminate it

Research parent: `0e8fc4087f8a676ed872daa970acc38756c78af9`.

The change here is a construction, not a revised acceptance tolerance. Two explicit
local auxiliary Hamiltonians replace arbitrary effective path weights and the
site-block re-completion. Their elimination is calculated exactly. The former
models and every archived result remain unchanged.

**Obtained:** all oriented return weights of a specified four-stage junction;
a non-polynomial transfer law; its absolute strength and induced frequency
residue; and a local coherent matter parent that removes the *constant
simultaneous-low-momentum* source error introduced by the old block completion.

**Not obtained:** a unique microscopic theory, the complete coordinate-stress
identity at arbitrary internal momentum, an interacting shared photon/gravity
cone, or an acceleration-distance law. These are two declared local mechanisms,
not yet one derived all-sector microscopic Hamiltonian. The new matter parent
is NOT equivalent to the old quartic Fock Hamiltonian.

## 1. Fix the loop amplitudes by an explicit finite Hamiltonian

Let S move a junction excitation through four directed stages around one square
plaquette. Three internal transports can be made identity by a gauge choice;
the fourth is W. Thus S is unitary and each diagonal block of S^4 is the
holonomy transported to that stage. In the chosen gauge it is W.

Summing the original completed **one-excitation** elementary link moves and a
unit auxiliary penalty gives exactly

    Q = d I - x(S+S†),    d=1+2x².

The factor 2x² is the sum of the two endpoint completions at each stage, not a
counterterm. The existing x=0.13554178509861228 is retained. Choose the decaying
factorization, determined by positivity and 0<r<1:

    alpha = [d+sqrt(d²-4x²)]/2,
    r = x/alpha,
    Q = alpha(I-rS)†(I-rS).

Now specify a **new, local chiral auxiliary Hamiltonian**

    K = [[0, A], [A†, 0]],    A=sqrt(alpha)(I-rS).

Both chiralities have K²=Q. K couples only adjacent stages plus the two chiral
states at each stage. Its auxiliary energies have magnitude at least
sqrt(alpha)(1-r)=0.8750197602. Positive and negative Dirac branches are retained;
a filled-negative-band fermion interpretation, rather than discarding the
negative branch, is available. The chiral dilation is an explicit architectural
choice. Q and K are not asserted to be the same Hamiltonian.

The lower-left zero-energy Green block at the entrance is exactly

    G_00(W) = alpha^(-1/2) (I-rho W)^(-1),    rho=r^4.

This follows by summing the powers of rS; only multiples of four return to the
entrance. It also follows by direct inversion of the finite Hermitian K. There
is no adjustable winding weight.

### The reference channel and the absolute amplitude

Use a second identical junction with identity holonomy as the existing type of
reference path. Two low ports couple to opposite auxiliary chiralities with
amplitudes x/sqrt(2), and the reference port has the opposite sign. This is an
explicit two-channel Hamiltonian, not a subtraction of an observed energy.

Its exact zero-energy Schur block is

    H_eff,RL(0) = -(x²/2)[G_00(W)-G_00(I)]
                = -kappa P(W),

    P(W) = (1-rho)(W-I)/(I-rho W),
    kappa = x² rho/[2 sqrt(alpha)(1-rho)²].

All diagonal Schur blocks vanish at zero frequency because K is chiral, not
because they were omitted. The original low-port block and both auxiliary
Hamiltonians are retained in the executable.

P has unit first curvature derivative ONLY as a reported shape convention;
**kappa is not divided away in the physical matrix.** In this construction

    rho   = 0.0003133951506985591,
    kappa = 0.000002854012469780666.

This is NOT the old area coupling x, nor the old photon/frame matching root.
Preserving the latter by rescaling kappa would be an extra intervention; none is
made here. Port layout and equal splitter amplitudes are explicit choices, not
proved uniquely selected properties of nature.

### The nonlinear function and the solved path weights

Expanding the resolvent gives

    P(W) = sum_(n>=1) a_n(W^n-I),
    a_n = (1-rho)² rho^(n-1).

The old normalization sum n a_n=1 follows, but every higher moment is now fixed
by x and the four-stage mechanism. For an eigenphase theta the oriented transfer
quadrature is

    Im P(exp(i theta)) = (1-rho)² sin(theta)
                         / [1+rho²-2rho cos(theta)].

In particular,

    m3 = -[Im P]'''(0)/[Im P]'(0)
       = (1+4rho+rho²)/(1-rho)²
       = 1.0018815500567193.

The paired transfer-power diagnostic is

    -(|P|²)''''(0)/(|P|²)''(0)
       = (1+10rho+rho²)/(1-rho)²
       = 1.0037631001134386.

These numbers are derived from a specified local Hamiltonian, not selected from
a family of polynomial weights to pass the linear-gravity tests. They describe
a microscopic transfer response, NOT an observed gravitational coefficient.
The old trigonometric area energy must not be inferred by simply completing
this effective transfer again: completion and elimination do not commute.
Doing that would define an additional model, which this work does not adopt.

### Do not drop the frequency residue

For the low ports coupled by V to heavy block K_H, the exact inverse propagator
is

    D_low(z) = z I + V(K_H-z I)^(-1)V†.

Hence its derivative at zero is

    Z_low = I + V K_H^(-2) V†,

not identity. At flat connection it is 1.018370697661399 times identity. The
actual small-phase low-band slope is consequently kappa/Z_low=
0.000002802528073848414. The code solves the full Hermitian matrix and verifies
the energy-dependent Schur equation at its actual low-band eigenvalues. The
transfer coefficient is not falsely reported as a propagation speed.

## 2. A local matter construction without the site-pinching artifact

The former Gaussian completion used the local site projectors to discard
interference terms. Instead, specify one physical and one auxiliary fermion
species per existing mode, with one common hopping rule:

    H_joint[h] = [[x h, x h], [x h, I+x h]].

Here h is the existing full geometric kinetic+Wilson+wall Hamiltonian. Every
spatial block has the same range as h, and all species inherit the same x and
unit auxiliary gap. There is no independently fitted frame/photon coefficient.
This is a NEW microscopic bilinear ansatz. It replaces, and is not secretly
equivalent to, the prior completion performed on each Fock bilinear. The old
quartic model remains committed separately; none of its terms are silently
claimed to be included in this free parent.

For every eigenvalue lambda of h, the two energies are exactly

    E_-(lambda) = x lambda + [1-sqrt(1+4x²lambda²)]/2,
    E_+(lambda) = x lambda + [1+sqrt(1+4x²lambda²)]/2.

The auxiliary branch is above 1/2 for all finite real lambda. E_- has the same
sign as lambda and its only zero is lambda=0. Therefore this enlargement does
not introduce another massless matter branch or extra zero at high momentum.
It does not, of course, solve the old chiral mirror question.

Elimination fixes the entire low-band operator:

    h_low = f(h) = x h - x²h² + x⁴h⁴ - 2x⁶h⁶ + ... .

The expansion is a small-|xh| expansion, not a globally convergent power series
on every microscopic band. The exact square-root function is used for spectra.
The inverse propagator before diagonalizing is

    D_low(z) = z-xh-x²h(z-I-xh)^(-1)h.

Where the chosen heavy pivot is singular, retain the full block rather than
using its Schur expression there. The full Hermitian model and the exact band
functions have no such singularity. This distinction matters for ultraviolet
loop integrals.

All source terms follow by differentiating this same operator. For a source J
connecting h0 and h1, delta f is its matrix divided difference, not an
independently chosen stress. In the full parent,

    J_joint = x [[J,J],[J,J]].

If R0=J-iD(h0-h1), the EXACT new residual is

    R_joint = x [[R0,R0],[R0,R0]],
    ||R_joint|| = 2x ||R0||.

The block-pinching contribution that produced the old 0.188... constant is
absent by construction. In the same simultaneous p~q->0 calculation used to
identify that error, R0=O(q^4), so ||R_joint||/|q|=O(q^3):

| L | Old pinched result /q | Coherent-parent result /q |
|---:|---:|---:|
| 16 | 0.1827432434 | 0.0087275675 |
| 64 | 0.1876462998 | 0.0001437802 |
| 256 | 0.1880930648 | 0.0000022540 |
| 1024 | 0.1881229265 | 0.00000003523 |

No threshold or old source formula was changed. An independently reconstructed
3^3 geometric graph verifies the new source. Its doubled spectrum agrees with
the analytic bands, and a complete small Fock calculation agrees with the
second-quantized bilinear parent, not the old quartic model.

### What this fixes, and what it does not

This removes the spurious **constant simultaneous-low-momentum error** of the
old site-pinched completion. It does not cure every flaw of the inherited
geometric regulator. With an internal momentum held high while external q goes
to zero, the underlying regulator still has a coordinate-response defect; the
full run records a nonzero example. Loop integrals sample those modes. Therefore
there is no claim of an exact all-momentum stress Ward identity, complete
interacting gravitational consistency, or a matched quantum light cone.

Similarly, f(h)'s derivative is not automatically the physical transition vertex
between the momentum-dependent low-band embeddings. The full doubled source and
its resolvent/residue must be used for a pole or loop calculation. The static
operator comparison is labeled accordingly.

## 3. Deliverable, not another changed test

The solved outputs are an explicit ring transfer with all a_n fixed, its
absolute coupling and residue, and an explicit local matter parent with an
analytic spectrum and inherited source rule. Historical formula-selection
results are preserved, not rewritten. No claim is made that the chosen chiral
ring and doubled matter parent have already been derived from one unique
all-sector microscopic principle.

This work does not supply g(r) or Newton's G. Converting the ring amplitude to
an area interaction in a dynamical frame theory must use its actual low-energy
Hamiltonian and residues, then vary with physical sources. It must not insert
the normalized P into the old area formula and quietly recover the old
coefficients. That would undo the point of solving the parent.

Reproduce the compact analytic/matrix check with:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python phase_junction_network/construction/coherent_junction/solve_coherent_junction.py --output /tmp/coherent_full.json

Add `--quick` for fewer sampled matrices. Results include the actual source
hashes and software versions. This study was not preregistered. The formulas
are analytic; numerical checks guard implementation, not determine parameters.

Eliminating auxiliary states and using Schur/Schrieffer-Wolff methods are prior
art, not claimed inventions: S. Bravyi, D. DiVincenzo and D. Loss,
*Schrieffer-Wolff transformation for quantum many-body systems*,
arXiv:1105.0675 (2011), doi:10.1016/j.aop.2011.06.004. A general review of why
lattice coordinate symmetry is a distinct issue is B. Bahr and B. Dittrich,
*Breaking and restoring of diffeomorphism symmetry in discrete gravity*,
arXiv:0909.5688 (2009). Neither reference establishes this proposed model.
