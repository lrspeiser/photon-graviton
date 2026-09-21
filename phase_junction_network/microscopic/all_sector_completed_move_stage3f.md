# Stage 3F: One Completed Move Across Four Finite Sectors

**Status:** **PASS at finite algebra and one-particle scope.** The same unit completed-move coefficient works in photon, dressed-frame, minimally charged endpoint, and neutral-companion representatives.  
**Date:** 2026-09-20  
**Issues:** #3, #4, #6, and #7  
**Executable:** [`check_all_sector_completed_move.py`](check_all_sector_completed_move.py)  
**Frozen output:** [`all_sector_completed_move_results.json`](all_sector_completed_move_results.json)

## 1. Shared local rule

Every elementary move \(M\) is assigned

\[
\boxed{
H_M=-x(M+M^\dagger)+x^2(M^\dagger M+MM^\dagger)
}
\]

with

\[
x=t/\Delta=0.13554178509861228.
\]

The coefficient of the diagonal partner is one in every tested sector.

The move has the exact factorization

\[
H_M=\frac14\left(Q_M^\dagger Q_M+Q_MQ_M^\dagger-2I\right),
\qquad
Q_M=I-2xM.
\]

Positive-square factorization is standard methodology. The project-specific question is whether one Phase Junction transport algebra can use this same rule in every sector without sector-specific counterterms.

## 2. Photon sector

The Stage-3E finite photon root remains:

\[
K_A/\Delta=0.005606197571911289,
\]

\[
U/\Delta=0.006789913753821769,
\]

\[
v/K_A=-0.12924513079019587.
\]

The finite photon branch has power

\[
0.9773450471
\]

and fit residual

\[
2.91\times10^{-4}.
\]

Scanning the completion coefficient around one changes the sign of the leading electric self-energy; only the unit coefficient cancels it while preserving the ring term.

## 3. Dressed-frame sector

Use one finite frame qudit and one connection qudit with

\[
R=Z_CZ_h^\dagger,
\qquad
\overline X=X_hX_C.
\]

The completed elementary frame moves commute with the connection lock.

| Prime | Full dimension | Physical ground band | Lock gap | Maximum lock commutator | Leakage |
|---:|---:|---:|---:|---:|---:|
| 3 | 18 | 3 | 0.826993 | \(1.11\times10^{-16}\) | \(6.75\times10^{-16}\) |
| 5 | 50 | 5 | 0.821373 | \(2.22\times10^{-16}\) | \(1.17\times10^{-15}\) |

The effective frame hopping scales as

\[
K_g^{(p=3)}\propto x^{1.99533},
\qquad
K_g^{(p=5)}\propto x^{1.99626}.
\]

The complete finite dressed-shift harmonic inventory reconstructs the low-band Hamiltonian to below \(3\times10^{-15}\). The \(p=5\) model generates a small second harmonic, which is retained rather than dropped.

## 4. Minimal charged endpoint

The finite endpoint model contains:

- matter position \(L/R\);
- spin-2 electric flux;
- an intermediate charge state;
- four elementary partial moves.

The full Hilbert dimension is 30; the exact Gauss sector has dimension 6.

At the shared \(x\), the low band has

\[
K_m=0.0338029748033\,\Delta
\]

and the next band lies

\[
29.1434K_m
\]

higher.

The effective operator contains only the expected identity and \(\sigma_x\) terms to numerical precision.

### Exact constraints and response

Both endpoint Gauss generators commute exactly with the Hamiltonian.

For left endpoint number \(N_L\) and its current \(J\),

\[
\boxed{i[H,N_L]+J=0}
\]

with zero matrix residual.

In the energy basis,

\[
i(E_m-E_n)\langle m|N_L|n\rangle
+\langle m|J|n\rangle=0
\]

with maximum residual

\[
4.65\times10^{-16}.
\]

The diamagnetic and paramagnetic pure-gauge curvatures are

\[
+0.0636484940413
\]

and

\[
-0.0636484940413,
\]

leaving residual

\[
5.55\times10^{-17}.
\]

The charged hopping scales as \(x^{1.99533}\), while its residual common self-energy scales as \(x^{3.99129}\).

## 5. Neutral companion

A three-state companion conversion

\[
|A\rangle\leftrightarrow|C\rangle\leftrightarrow|B\rangle
\]

uses the same completed move.

At the shared \(x\):

\[
K_\chi=0.0174426357263\,\Delta,
\]

with next-band separation

\[
60.3840K_\chi.
\]

The conversion scales as \(x^{1.99719}\), while its residual self-energy scales as \(x^{3.99531}\).

Maximum finite conversion probability is

\[
0.998343,
\]

and forward/reverse reciprocity closes to

\[
2.78\times10^{-16}.
\]

## 6. Decision

The same coefficient one is compatible with all four finite representatives. No photon-only, frame-only, matter-only, or companion-only completion coefficient was introduced.

Stage 3F therefore passes its declared finite algebra and one-particle gate.

## 7. Claim boundary

Stage 3F does not establish:

- one simultaneous many-body Hamiltonian containing all sectors;
- a nonzero-momentum Ward identity in the spatial photon phase;
- transverse photon polarization by dynamical matter;
- virtual matter loops;
- radiative stability in an interacting continuum;
- mirror-wall completion;
- nonlinear gravity;
- physical \(\alpha\), \(G\), or formal novelty.

## 8. Next gate

Place the minimal charged endpoint in the actual finite spatial photon Hamiltonian and require:

1. exact local continuity;
2. a nonzero-momentum Ward identity;
3. exact longitudinal pure-gauge invariance;
4. nonzero transverse matter polarization;
5. finite photon pole dressing with retained transverse residue;
6. no new sector-specific counterterm.

That next gate is implemented by [`spatial_gauge_matter_stage3g.md`](spatial_gauge_matter_stage3g.md).
