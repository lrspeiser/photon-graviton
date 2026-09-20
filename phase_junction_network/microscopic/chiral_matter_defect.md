# Finite chiral matter defects with protected overlap gaps

**Status:** finite prototype pass for issue #4  
**Date:** 2026-09-20  
**Calculation:** [`check_chiral_matter_defect.py`](check_chiral_matter_defect.py)  
**Frozen output:** [`chiral_matter_results.json`](chiral_matter_results.json)

## 1. What is established here

This construction supplies the first explicit finite charged-matter sector for the Phase Junction network. It is not a continuum Dirac field appended to the infrared theory. Its microscopic objects are:

1. spin-1 quantum-link lanes with three states per lane;
2. finite fermionic strand orbitals at each site and at each layer of a finite synthetic direction;
3. an exact local binding projector that turns an odd number of strands into one charged endpoint defect;
4. a gauge-covariant composite hopping move that changes the electric flux on exactly the same number of lanes as the defect charge;
5. a finite Wilson/domain-wall slab whose boundary supports chiral low-energy modes;
6. one shared frame operator for every matter species;
7. a common charge-to-localization rule that produces exponentially separated residual gaps rather than arbitrary tiny on-site energies.

The selected left-handed charge spectrum is

\[
\boxed{q=(-11,-5,-1,-1,9,9)}.
\]

It obeys both four-dimensional Abelian anomaly conditions,

\[
\sum_s q_s=0,
\qquad
\sum_s q_s^3=0,
\]

contains no pair \(q,-q\), and therefore admits no gauge-neutral left-left bilinear mass. Every charge is odd, so the bound endpoint has odd fermion parity.

At slab width \(L_s=12\), one common microscopic localization rule produces four distinct positive one-particle gaps spanning a factor

\[
\boxed{3.2452902\times 10^7}.
\]

The result is a finite, anomaly-free chiral endpoint prototype with a protected hierarchy. It is **not** a derivation of the Standard Model charge spectrum, observed generations, observed Yukawa couplings, or a completed interacting QED continuum.

## 2. Finite local Hilbert space

### 2.1 Electromagnetic link bundle

Retain the spin-1 quantum link already selected by the electromagnetic kinematic audit. Each oriented lane \(a\) on a geometric link \(\ell\) has

\[
E_{\ell a}|m\rangle=m|m\rangle,
\qquad m\in\{-1,0,1\},
\]

with raising operator \(U_{\ell a}\) satisfying

\[
[E_{\ell a},U_{\ell a}]=U_{\ell a}.
\]

A geometric link contains eleven parallel finite lanes. Its full link-bundle Hilbert space therefore has dimension

\[
3^{11}=177147.
\]

A charge-\(q\) defect uses only \(|q|\) of those lanes. Define

\[
\mathcal U^{(q)}_{xy}
=
\prod_{a=1}^{|q|}
\begin{cases}
U_{xy,a}, & q>0,\\
U_{xy,a}^{\dagger}, & q<0.
\end{cases}
\]

Then

\[
\left[
\sum_{a=1}^{|q|}E_{xy,a},
\mathcal U^{(q)}_{xy}
\right]
=q\,\mathcal U^{(q)}_{xy}.
\]

Charge is therefore not an independent label attached after the fact: it is the signed number of elementary flux lanes terminated and transported by the defect.

### 2.2 Microscopic fermionic strands

For a species of charge \(q\), introduce \(|q|\) microscopic fermionic strands

\[
f_{q,a,A}(x,s),
\qquad
 a=1,\ldots,|q|,
\qquad
 A=1,\ldots,4,
\qquad
 s=0,\ldots,L_s-1.
\]

The index \(A\) is the four-orbital Wilson/domain-wall index. Each strand orbital is a two-state fermionic Fock space. Before binding, the local dimension for one species at one spatial site and one synthetic layer is

\[
2^{4|q|}.
\]

For the largest selected charge, \(|q|=11\), this is \(2^{44}\). It is large but finite.

Add an exact local binding penalty

\[
H_{\rm bind}
=
\Lambda
\sum_{x,s,q}
\sum_{A=1}^{4}
\sum_{a<b}
\left(n_{q,a,A}(x,s)-n_{q,b,A}(x,s)\right)^2.
\]

The zero-penalty band has all strands of a species at the same occupation for each orbital. Because the composite move acts only between the all-empty and all-filled configurations, this band is an exact invariant subspace at any positive \(\Lambda\); taking \(\Lambda\) large only separates unbound strand excitations energetically. The composite operator is

\[
D_{q,A}^{\dagger}(x,s)
=
\prod_{a=1}^{|q|}f_{q,a,A}^{\dagger}(x,s).
\]

The script enumerates the one-orbital strand basis through \(|q|=11\), verifies exactly two zero-penalty states per orbital, and verifies that composite creation and annihilation never leave that band. The bound composite subspace therefore has four fermionic orbitals and dimension \(2^4=16\), independent of \(|q|\).

Because every selected \(q\) is odd, exchanging two composites crosses \(q^2\) microscopic fermion operators and gives

\[
D_q^{\dagger}(x)D_q^{\dagger}(y)
=(-1)^{q^2}
D_q^{\dagger}(y)D_q^{\dagger}(x)
=-D_q^{\dagger}(y)D_q^{\dagger}(x).
\]

Thus the endpoint exchange statistic is tied to the same odd endpoint multiplicity that defines its charge.

## 3. Exact Gauss-law compatibility

Use the site generator

\[
G_x
=
\sum_{i,a}
\left(E_{x,i,a}-E_{x-\hat i,i,a}\right)
-
\sum_{q,s,A}q\,N_{q,A}(x,s),
\]

where \(N_{q,A}=D_{q,A}^{\dagger}D_{q,A}\) in the bound subspace.

The finite spatial hopping move is

\[
H_{\rm hop}
=-t_m
\sum_{\langle xy\rangle,q,s}
D_q^{\dagger}(x,s)
\,T_{xy}[e]\,
\mathcal U_{xy}^{(q)}
D_q(y,s)
+\text{h.c.}
\]

The matter move \(y\rightarrow x\) changes \(N_x\) by \(+1\), \(N_y\) by \(-1\), and the link bundle flux by \(q\). Therefore

\[
[G_x,H_{\rm hop}]=[G_y,H_{\rm hop}]=0
\]

exactly. The regression enumerates every nonzero spin-1 bundle transition for each distinct selected charge. The largest case, \(|q|=11\), contains \(2^{11}=2048\) nonzero bundle transitions, all with zero Gauss residual.

## 4. Why these charges were selected

The script performs an exhaustive finite search over spectra satisfying:

- two through six left-handed species;
- nonzero odd integer charges;
- \(|q|\le 11\);
- primitive normalization, \(\gcd(|q_s|)=1\);
- no vectorlike pair \(q,-q\);
- \(\sum q=0\);
- \(\sum q^3=0\).

There is no solution with fewer than six species in this search class. At six species there is no solution through \(|q|\le9\). At \(|q|\le11\), the unique solution up to permutation and overall charge reversal is

\[
(-11,-5,-1,-1,9,9).
\]

The anomaly checks are exact integers:

\[
-11-5-1-1+9+9=0,
\]

\[
(-11)^3+(-5)^3+(-1)^3+(-1)^3+9^3+9^3=0.
\]

No two selected charges sum to zero. A left-left Lorentz bilinear would require \(q_i+q_j=0\), so every such mass is forbidden by the exact gauge symmetry.

This is the first nontrivial relationship among charge, chirality, and statistics:

> anomaly cancellation plus the microscopic odd-endpoint rule fixes a finite chiral charge spectrum rather than permitting independently assigned species labels.

## 5. Finite domain-wall Hamiltonian

### 5.1 Gamma algebra

Use five Hermitian \(4\times4\) matrices satisfying

\[
\{\Gamma_A,\Gamma_B\}=2\delta_{AB}.
\]

The finite synthetic direction is indexed by \(s=0,\ldots,L_s-1\). For spatial momentum \(\mathbf k\), common frame \(e_a{}^i\), and species bulk parameter \(m_5(q)\), the exact one-composite block is

\[
H_q(\mathbf k)
=
\sum_s D_{q,s}^{\dagger}
\left[
\sum_{a,i}\Gamma_a e_a{}^i\sin k_i
+
\Gamma_5
\left(
 m_5(q)+1+
 \sum_i(1-\cos k_i)
\right)
\right]
D_{q,s}
\]

\[
+
\sum_{s=0}^{L_s-2}
D_{q,s}^{\dagger}
\left(-\frac{\Gamma_5+i\Gamma_4}{2}\right)
D_{q,s+1}
+\text{h.c.}
\]

Spatial hopping in position space carries the flux-bundle operator \(\mathcal U^{(q)}\). Synthetic-direction hopping is onsite in physical space and therefore needs no electromagnetic link.

The wall-chirality operator is

\[
\chi_{\rm wall}=i\Gamma_4\Gamma_5.
\]

The four lowest states at \(\mathbf k=0\) separate into two states localized on each wall. Numerically,

\[
\langle\chi_{\rm wall}\rangle_{\rm left}=-1,
\qquad
\langle\chi_{\rm wall}\rangle_{\rm right}=+1
\]

to roundoff for every species.

This implements the domain-wall mechanism: one physical wall contains the selected left-handed spectrum, while the remote wall contains the opposite-chirality regulator partners. The finite regulator is therefore explicit about its mirror sector rather than hiding species doubling.

### 5.2 Universal frame coupling

Every species uses the identical derivative

\[
\frac{\partial H_q}{\partial e_a{}^i}
=
\sin k_i\,
I_{L_s}\otimes\Gamma_a.
\]

Charge enters only through \(\mathcal U^{(q)}\); the frame operator contains no species-dependent coefficient. The numerical audit perturbs all nine frame components and finds a maximum species-to-species derivative difference of zero at floating-point precision.

This is the required universal bare frame coupling. Whether all interacting species flow to one Lorentz cone remains part of the continuum gate in issue #8.

## 6. Chirality and species doubling

For \(-2<m_5<0\) in the stated convention, the open Wilson slab has one boundary cone at the physical Brillouin corner \(\mathbf k=0\). The Wilson term gaps the other seven corners.

For every selected species, the regression verifies:

- exactly one low-energy physical corner;
- the smallest nonzero-corner gap is larger than one in the chosen lattice units;
- the two walls have opposite chirality;
- positive and negative single-particle energies pair to better than \(1.6\times10^{-15}\).

The negative-energy branch supplies the antiparticle sector after second quantization. A left-handed particle of charge \(q\) corresponds to a right-handed antiparticle of charge \(-q\).

The regulator still contains the remote mirror wall. That is a declared finite regulator degree of freedom, not an unreported doubler. Its coupling to the physical wall is the source of the exponentially small residual gap below. A strictly local 3+1D construction that removes or symmetrically gaps the remote mirror while preserving the full interacting gauge theory is not claimed here; it remains a continuum-completion requirement under issue #8.

The result is also checked throughout a finite neighborhood rather than at one tuned matrix point. For every distinct charge, the regression shifts the local Wilson mass by \(-0.05,-0.025,0,0.025,0.05\). All twenty finite matrices remain inside the topological interval, retain exactly one light physical corner, preserve opposite wall chirality and localization, and keep every nonzero corner gapped above one lattice unit. This is a finite free-regulator stability test; radiative and interacting stability remain issue #8.

## 7. Protected mass hierarchy

### 7.1 No same-wall bilinear mass

There are two independent protections.

First, gauge invariance forbids every left-left bilinear because the selected spectrum contains no \(q_i+q_j=0\) pair.

Second, even the repeated-charge pairs \((-1,-1)\) and \((9,9)\) cannot form a same-chirality kinetic mass. In the projected Weyl space, a mass matrix would have to anticommute with all three Pauli kinetic generators. The script constructs the complete Hermitian matrix basis for two same-chirality flavors and finds nullity zero for

\[
\{M,\sigma_i\otimes I_{\rm flavor}\}=0,
\qquad i=1,2,3.
\]

A local perturbation confined to one wall may shift or mix equal-charge cones, but it cannot create a bilinear spectral gap. A gap requires coupling to opposite chirality.

### 7.2 Common microscopic localization rule

The only opposite-chirality partner is on the remote wall. Define one common endpoint-tension rule

\[
r_q=r_0+\eta(|q|-1),
\qquad
m_5(q)=-1+r_q,
\]

with

\[
r_0=0.12,
\qquad
\eta=0.04.
\]

This can be implemented as the same local synthetic-direction retention cost per extra bound endpoint strand. It uses two shared microscopic parameters for all six species, rather than one independently tuned mass per species.

For an isolated wall, \(r_q\) is the zero-momentum decay ratio into the synthetic bulk. At finite width the opposite walls overlap, giving

\[
\Delta_q(L_s)
\simeq
\frac{(1-r_q^2)r_q^{L_s}}
{1-r_q^{2L_s}}.
\]

The small number is produced by spatial separation and topology. No tiny onsite defect energy is inserted.

### 7.3 Frozen numerical spectrum

At \(L_s=12\):

| Charge | Multiplicity | \(r_q\) | \(m_5(q)\) | Numerical \(\Delta_q\) | Fitted velocity |
|---:|---:|---:|---:|---:|---:|
| \(-1\) | 2 | 0.12 | -0.88 | \(8.7876\times10^{-12}\) | 1.0000000 |
| \(-5\) | 1 | 0.28 | -0.72 | \(2.1401\times10^{-7}\) | 1.0000000 |
| \(9\) | 2 | 0.44 | -0.56 | \(4.2460\times10^{-5}\) | 1.0000000 |
| \(-11\) | 1 | 0.52 | -0.48 | \(2.8518\times10^{-4}\) | 1.0000009 |

Equal charges are exactly degenerate because they share the same microscopic rule. The four distinct gaps span

\[
\frac{\Delta_{|q|=11}}{\Delta_{|q|=1}}
=3.2452902\times10^7.
\]

For widths \(L_s=6,8,10,12,14\), each fitted slope of \(\log\Delta_q\) agrees with \(\log r_q\) within the frozen tolerance. The direct finite-matrix gaps agree with the overlap formula to at worst \(2\times10^{-5}\) relative error.

This is the second required nontrivial relationship:

\[
\boxed{
q
\longrightarrow
|q|\text{ endpoint strands}
\longrightarrow
r_q
\longrightarrow
\Delta_q\propto r_q^{L_s}.
}
\]

The relation predicts exact degeneracy for repeated charges and an ordered hierarchy from one common rule. It does not claim that this toy ordering is the observed fermion spectrum.

## 8. Stability and negative controls

The frozen script fails if any of the following occurs:

- the five gamma matrices cease to satisfy the Clifford algebra;
- the anomaly-free charge search returns a different first spectrum;
- either anomaly sum is nonzero;
- a vectorlike pair appears;
- any selected endpoint has even strand number or bosonic exchange phase;
- a flux-bundle hopping transition violates either endpoint Gauss law;
- more than one physical Brillouin corner becomes light;
- wall localization or wall chirality fails;
- the low-energy velocity becomes nonpositive or non-linear at the frozen tolerance;
- particle/antiparticle spectral pairing fails;
- the finite-width gap ceases to follow exponential overlap scaling;
- a gauge-neutral or kinetic same-wall bilinear mass appears;
- the one-cone phase fails anywhere in the frozen Wilson-mass neighborhood;
- the frame derivative becomes species dependent;
- the charge-linked gaps fail to span the frozen nontrivial hierarchy.

Every output is written to `chiral_matter_results.json`; the script exits nonzero on failure.

## 9. Acceptance-criteria map for issue #4

| Issue requirement | Construction or result |
|---|---|
| Explicit local finite Hilbert space | Spin-1 link lanes, finite strand Fock spaces, exact invariant binding band, finite \(L_s\) |
| Gauge-covariant hopping | \(D_x^\dagger T[e]\mathcal U^{(q)}_{xy}D_y\), exact Gauss commutator |
| Universal frame coupling | One identical \(e_a{}^i\Gamma_a\) operator for every species |
| Fermionic exchange statistics | Odd-strand composite gives \((-1)^{q^2}=-1\) |
| Chirality and doubling treatment | One Weyl cone per wall; seven physical corners gapped; finite mass-window stability; remote mirror retained explicitly |
| Anomaly cancellation | Exact \(\sum q=\sum q^3=0\) for the searched charge set |
| Particles and antiparticles | Exact \(E\leftrightarrow-E\) spectral pairing |
| Charge constrained microscopically | Signed number of spin-1 flux lanes; finite exhaustive anomaly search |
| Calculated dispersion and gap | Direct diagonalization of the finite slab; \(E^2=\Delta_q^2+v^2k^2+\cdots\) |
| Protected small masses | Same-wall mass forbidden; opposite-wall overlap exponentially small |
| Nontrivial charge/chirality/mass relation | odd anomaly-free charges and \(r_q=r_0+\eta(|q|-1)\), hence \(\Delta_q\propto r_q^{L_s}\) |

The finite prototype therefore passes the explicit acceptance criteria of issue #4.

## 10. Claim boundary and handoff

### Established

- one explicit finite charged endpoint construction;
- exact finite Gauss-covariant hopping through the spin-1 link sector;
- a derived fermionic exchange sign for every selected defect;
- a finite anomaly-free same-chirality charge spectrum in the declared search class;
- a Wilson/domain-wall regulator with one physical Weyl cone per wall, no light physical-corner doublers, and a finite verified topological neighborhood;
- exact universal bare frame coupling;
- particles and antiparticles;
- exponentially protected finite-slab gaps;
- a shared charge-to-gap hierarchy with fewer microscopic parameters than distinct gaps.

### Not established

- the observed Standard Model charges, color, weak isospin, or three generations;
- observed fermion masses or mixings;
- a derivation of \(r_0\) and \(\eta\) from the same move amplitudes that determine \(U_A,K_A,U_g,K_g\);
- removal or symmetric gapping of the remote mirror in the full interacting theory;
- a 3+1D deconfined electromagnetic phase with these dynamical defects;
- interacting Ward identities, radiative stability, or Lorentz recovery;
- bound-state clocks and chemistry.

Those are not reasons to reinterpret this result as failure. They are the next dependency boundaries:

- issue #3 owns common microscopic coefficients, including deriving \(r_0\) and \(\eta\);
- issue #6 owns the deconfined finite QED phase with dynamical defects;
- issue #8 owns mirror completion, interacting anomalies, unitarity, Ward identities, radiative stability, and Lorentz recovery;
- issue #9 owns any frozen empirical use of the resulting mass and charge relations.

## 11. Reproduce

From the repository root:

```sh
python phase_junction_network/microscopic/check_chiral_matter_defect.py \
  --slab-width 12 \
  --output phase_junction_network/microscopic/chiral_matter_results.json
```

## Reference used for the regulator comparison

- D. B. Kaplan, “A Method for Simulating Chiral Fermions on the Lattice,” *Physics Letters B* **288** (1992) 342–347, arXiv:hep-lat/9206013.
