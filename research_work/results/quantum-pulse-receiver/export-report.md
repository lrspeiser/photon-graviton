# Exporting the receiver energy into a moving companion

The existing receiver spectrum permits an energy/momentum-conserving export into a forward massless mode. Unconditional export leaves already processed photons' density matrix unchanged. However, full export leaves the receiver in a sharp-energy ground state, removing the coherence needed by the next conversion. This supports a finite batch followed by export at the level of asymptotic state bookkeeping; it does not yet provide continuous local void transport.

## Explicit proposed map

Use the existing receiver spectrum P_R=E_R-U in c=1 units, now restricted to E_R>=U. The full-export proposal is

    |E_R, P_R=E_R-U>_R |vacuum>_c
      -> |U, P=0>_R |E_c=E_R-U, P_c=E_R-U>_c.

At E_c=0 the outgoing state is the vacuum, not a zero-energy particle. This is a new project completion of the stipulated receiver model using known relativistic and quantum mathematics; no fundamental novelty or graviton identity is claimed. It is not a derived decay amplitude, selection rule, cross section or lifetime. The original receiver spectrum has variable invariant mass and is not the dispersion of one free fixed-mass particle.

Basis by basis,

    E_initial = E_R = U + E_c = E_final
    P_initial = E_R-U = P_c = P_final
    E_c >= 0,  E_c^2-P_c^2=0.

Distinct initial receiver energies map to distinct outgoing companion energies, so the map preserves inner products on its input subspace. A reverse channel can complete it to a unitary on an enlarged space. This is quantum state bookkeeping, not proof of a local physical interaction. Lossless travel at c is possible under a free massless propagation postulate afterward; that propagation and any coupling to the clock field remain to be specified consistently.

## Why releasing energy does not retroactively spoil the light

For any photon-receiver state rho, an isometry V acting only on the receiver and an initially empty output sector satisfies

    Tr_(R,c)[(I tensor V) rho (I tensor V_dagger)] = Tr_R[rho].

This is the standard partial-trace identity for an unconditional environment operation; see [Preskill's quantum-channel notes](https://www.preskill.caltech.edu/ph229/notes/chap3.pdf). It preserves the photon profiles and correlations that existed before export. The distinguishing receiver information has moved to the companions rather than been erased. Postselecting companion measurements would be a different operation and is not used here.

The explicit two-photon calculation uses the preceding sequence model: initial receiver mean energy 55, U=50, and two incoming photons of mean energy 10, each reduced to 5. Exported energy is 15: **10 supplied by photons plus 5 already present in the prepared receiver**. The receiver retains its baseline energy 50. Crediting all 15 to converted starlight would double-count the preparation contribution.

export.py verifies zero basis energy/momentum discrepancy and zero change in the photon density matrix after relabeling the receiver into the companion sector. The mean ledger agrees within floating-point precision. This verifies the stated isometry, not emission dynamics or arrival times.

## What happens to the next photon

If full export occurs between uses, the next photon encounters |U>, not the original coherent receiver. For the same fractional conversion q=1/S,

    |E>_gamma |U>_R -> |qE>_gamma |U+(1-q)E>_R.

Different incoming energies now leave orthogonal receiver states. Tracing the receiver removes off-diagonal photon energy coherences. A later full export preserves that reduced state; it cannot restore the lost marginal coherence by acting only on the receiver.

In the finite nine-level example, the sum of absolute off-diagonal elements is 7.2603 with the prepared coherent receiver and zero after a ground-state reset. This number is a basis-dependent diagnostic, not a measured visibility or an astronomical exclusion. The latter state has no localized Fourier pulse profile over the finite recurrence period. It does not simulate photons physically arriving uniformly through cosmic time. The result rejects automatic reuse after full discharge for this particular coherent-pulse mechanism, not every explanation of supernova timing.

## Consequence for the complete model

Energy export itself is not the obstruction: an explicit compatible output state exists. Continuous operation needs either partial export that retains a working coherent state, replenishment carrying its own energy and temporal resource, or a different local interaction. We have not derived any of those. A finite batch-and-export map also leaves energy temporarily in the receiver; its dwell time is unknown, so it does not establish the requested absence of long-term void accumulation.

Next work must connect local conversion, receiver maintenance and export in one process before fitting an astronomical rate. The pulse map still does not derive causal event transport, the distance dependence of S, capture at wells, or a lensing/rotation field. No new observational fit or final holdout was used. All six objectives remain open.
