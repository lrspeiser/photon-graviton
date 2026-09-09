# Ordinary gravitons: small-energy emission during photon scattering

## Why this calculation is needed

We have tested hypothetical scalar companion waves. An ordinary graviton is a different kind of wave, so those results cannot simply be renamed. Here we ask whether ordinary gravity supplies the small energy transfers our redshift idea needs. If we skip this, neither the graviton identity nor its production strength is supported.

The useful distinction is between **many very small packets** and **enough total energy**. Ordinary gravity allows very low-energy radiation during scattering. Its packet count grows as we include ever smaller packets, but its carried energy does not grow without limit. This calculation makes that distinction quantitative.

This is an optional comparison under ordinary weak-field gravity. It does not adopt expansion, a Big Bang, dark matter, a new target population, or a special void clock. No complete astronomical theory or measured graviton signal is claimed.

## Source basis and process selection

[Cachazo and Strominger, section 1, equations 2–4](https://arxiv.org/html/1404.4091) state the leading soft-graviton factor and its connection to momentum conservation. We restore the gravitational coupling below. We use only the leading term, not the paper's subleading conjecture.

[Ahmadiniaz et al.](https://arxiv.org/html/1908.03425) construct amplitudes with a massive scalar target, photons and one graviton, including contributions from photon lines. Their one-photon photoproduction channel removes the incoming photon; it does not by itself redshift a surviving photon. The two-photon channel permits the process of interest, gamma + T -> gamma' + g + T'. We have not implemented their complete finite-energy amplitude.

The angular bounds, numerical diagnostics and throughput comparison below are our conditional calculations from the leading soft factor. They are not results or observational claims attributed to those papers.

## 1. A physical coupling and its domain

Take Einstein gravity, Maxwell light and an ordinary target with a specified elastic scattering interaction. In units hbar=c=1, write the gravitational action coefficient as 2/kappa^2 with kappa^2=32 pi G and expand the metric as g=eta+kappa h. The first-order graviton coupling has magnitude kappa h_mu_nu T^mu_nu/2. It couples to the full conserved stress tensor, including the photon and target; selecting only the photon contribution would generally break the required conservation identity.

For a hard elastic event gamma(k)+T(p) -> gamma(k')+T(p'), adding a graviton q=omega(1,n) with omega much smaller than the hard scales gives

    M_with_g = S_g M_hard + O(omega^0)
    S_g = (kappa/2) sum_a eta_a p_a^mu p_a^nu epsilon*_mu_nu / (p_a.q)

Here eta=-1 for incoming and +1 for outgoing legs, signature is +---, and each physical graviton polarization has unit tensor norm. Hard momenta conserve energy and momentum at the expansion point. The finite-omega event requires correspondingly shifted final momenta; using unchanged hard momenta is only the leading asymptotic calculation.

This is a derivation of a universal low-energy emission factor from a chosen ordinary interaction framework. It is not a first-principles selection of that framework as the law of the fictional universe.

## 2. Conservation and the two tensor polarizations

Set the incident photon energy E=1 temporarily and place the target at rest with mass M. If the photon scatters by theta, exact elastic kinematics are

    E' = E / [1 + (E/M)(1-cos theta)]
    p' = p+k-k'.

Define qhat=(1,n) and J^mu_nu=sum eta p^mu p^nu/(p.qhat). Contracting with qhat gives sum eta p=0. This identity checks that unphysical polarization choices cannot create radiation.

For a unit graviton direction n, the spatial transverse projector is P=I-nn^T. The physical polarization sum for a symmetric spatial J is

    K(n) = ||P J P||_F^2 - [Tr(P J P)]^2/2 >= 0.

The code checks this expression against explicit plus/cross tensors, verifies the conservation identity in 100 random finite-mass events, and checks the recoiling target's mass shell. These tests validate the implemented soft factor, not the complete interaction probability.

For a real emission the full ledger is E+M=E_gamma'+E_g+E_target'. A changed target kinetic energy is part of the ledger. Graviton energy cannot automatically be identified with all photon energy loss, and target recoil is not a graviton deposit.

## 3. Packet count, carried energy and a bounded angular factor

Combining the leading amplitude with one-graviton phase space gives, conditional on a hard event,

    dP_g = [G/(2 pi^2)] K(n) dOmega domega/omega.

Let I(theta)=integral K(n)/E^2 dOmega and C(theta)=G E^2 I(theta)/(2 pi^2). Then dP_g=C domega/omega. C is dimensionless and proportional to (E/E_Planck)^2, with the non-reduced Planck energy E_Planck=sqrt(hbar c^5/G).

In the strict heavy-target limit, the target TT contributions vanish and E'=E. Each photon leg's TT norm is E(1+n.n_photon)/sqrt(2), at most sqrt(2)E. The triangle inequality therefore gives

    K <= 8 E^2; I <= 32 pi; C <= (16/pi)(E/E_Planck)^2.

This is a deliberately loose bound on the **leading soft term for this heavy elastic event**, not an upper bound on every possible graviton process. Finite-target corrections, finite-energy amplitudes, coherent macroscopic sources and different incident states are outside that bound.

At exactly zero photon deflection the incoming and outgoing hard states coincide and the leading factor cancels. This does not prove that every subleading forward process vanishes. For backscatter, K/E^2=2 cos^2(theta_g), giving the independent analytic check I=8 pi/3.

Numerical angular integrals at 256 polar nodes and 512 azimuth nodes:

| Photon scattering angle | I in the heavy-target limit |
|---|---:|
| 0 degrees | 0 |
| 10 degrees | 1.9267 |
| 30 degrees | 9.6640 |
| 90 degrees | 21.6095 |
| 180 degrees | 8.3776 |

Doubling both angular resolutions from the coarser grid changes these integrals by less than 0.005 times max(1,I). At 90 degrees, finite masses M/E=100,1000,10000 approach the heavy result with differences 0.215,0.0216,0.00216. No hard-event angular distribution is assumed; these are conditional factors, not a total target cross section.

## 4. Why an infrared packet count cannot solve the energy budget

For fractional graviton energies f=omega/E between f_min and delta much smaller than one, the leading expressions are

    P_resolved = C log(delta/f_min)
    <E_g>/E = C(delta-f_min)
    <(E_g/E)^2> = C(delta^2-f_min^2)/2.

These are real-emission contributions per hard event, including its usually zero graviton output. They are not moments conditional on a graviton having been observed. Taking f_min toward zero increases the first expression but leaves finite energy moments. An inclusive prediction needs virtual corrections and unresolved/multiple-soft emission; a diverging lowest-order packet count must not be treated as a physical probability larger than one.

The code integrates these three moments independently in log frequency for three lower cutoffs. Reducing the ceiling delta to demand smaller packets also reduces the energy carried by that band. The soft theorem alone does not prohibit emission above the selected ceiling or supply a mechanism that filters it out.

## 5. Quantitative implication for our redshift benchmark

For an illustrative 1 eV photon and soft ceiling delta=0.001, ordinary constants give

    C <= 3.42e-56
    mean soft-graviton energy / initial photon energy <= 3.42e-59

within the leading heavy-target comparison. These are expectation values per hard scattering, not the size of every emitted packet. Delta is a chosen diagnostic boundary, not an observed cutoff; finite-soft corrections have not been bounded by a full amplitude calculation.

To supply our reused empirical alpha=0.0002488993/Mpc through this band alone, a necessary leading-order throughput would be at least alpha/(C delta), about 7.28e54 hard events per Mpc. No target density or cross section has been selected to create that rate. Extrapolating independent transparent propagation to such a rate would be invalid; this is a diagnostic of how weak the ordinary soft channel is, not a viable repeated-scattering construction.

For a fixed hard angle in the heavy limit and a fixed fractional band, C scales as E^2. The full propagation rate also contains the energy dependence and angles of the hard scattering cross section. Consequently, achromatic redshift is not implied by the soft theorem, nor have we proved its impossibility for all interactions.

## Consequences and next work

Ordinary gravitons now have a specified small-energy production route, and its energy cannot be increased merely by counting more nearly energy-free packets. This ordinary soft-scattering branch does not supply the desired conversion strength under transparent independent scattering. It also inherits a hard scattering event that can redirect light.

Do not adopt it as the final law or generalize this result to all companion fields. A different coupling, a derived coherent background or a joint material response would require its own strength, energy source, spectrum and angular calculation. Finite-energy ordinary emission remains a separate calculation, not something that can be inferred by extending the soft formula to omega=E. Capture, storage, lensing, event timing and the original 32 observational areas remain unresolved by this pass.

Reproduce with `python -X utf8 research_work/results/soft-graviton/check_soft_graviton.py`, or run the full repository diagnostic runner. Results are saved in `soft-graviton-results.json`; fresh runs default to an ignored generated directory. Mathematical checks do not establish astronomical agreement.
