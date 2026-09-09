# A companion/time field that reacts to photon energy transfer

## Purpose and scope

The agreed receiver is the companion/graviton sector. The question here is whether an explicit interaction can transfer photon energy into that sector while conserving momentum and allowing field disturbances to propagate. This optional classical scalar-field candidate extends the earlier matched-wave receiving field beyond its homogeneous calculation. It is not an identification with ordinary spin-2 gravitons or a finished gravitational theory.

Unlike the prescribed void profile in the cumulative-time example, the field now evolves under the same Hamiltonian as the photon packet. Its initial motion and energy are recorded. There is no added receiver term chosen afterward to close the energy ledger.

## 1. A closed effective Hamiltonian

Work in one periodic spatial dimension of length L and reference c=1 units. Expand a positive field as

    n(x)=n0 + sum_j [qc_j cos(k_j x)+qs_j sin(k_j x)],
    k_j=2*pi*j/L.

The homogeneous mode has inertia M0=K L and each sine/cosine mode has inertia M=K L/2. A Gaussian spatial smoothing, w_j=exp[-(sigma_x*k_j)^2/2], defines the field n_bar(X) sampled by a finite photon packet. This form factor is a declared coarse-grained interaction, not a derived quantum cutoff.

The candidate Hamiltonian is

    H_field = P0^2/(2M0)
              + sum_j [(Pc_j^2+Ps_j^2)/(2M)
                        + M*v^2*k_j^2*(qc_j^2+qs_j^2)/2],
    H_photon = P_packet / n_bar(X),
    H_total = H_field + H_photon.

All terms are nonnegative for n_bar>0 and positive packet momentum. A fixed-number photon packet is represented by a ray degree of freedom; this is not a single-photon quantum emission calculation. For a field not being driven, the nonzero modes satisfy omega_j=v*k_j and can carry disturbances independently of the photon ray.

Hamilton's equations give

    X_dot=1/n_bar,
    P_packet_dot=P_packet*(partial_X n_bar)/n_bar^2,
    n0_dot=P0/M0;   P0_dot=P_packet/n_bar^2,
    qc_dot=Pc/M;    qs_dot=Ps/M,
    Pc_dot=-M*v^2*k^2*qc + (P_packet/n_bar^2)*w*cos(kX),
    Ps_dot=-M*v^2*k^2*qs + (P_packet/n_bar^2)*w*sin(kX).

The photon drives both the mean field and its spatial modes. In return their values change photon propagation. Energy conservation follows because H_total has no explicit time dependence. The instantaneous photon loss is

    dH_photon/dt=-(P_packet/n_bar^2)*(partial_t n_bar),
    dH_field/dt=-dH_photon/dt.

Thus the companion sector receives exactly the lost photon energy through the interaction itself. A reversible Hamiltonian does not guarantee loss rather than gain in every possible field state; all sampled trajectories here show monotone photon-energy loss over the tested interval. Long-term behavior must not be inferred from this short finite example.

## 2. Momentum also balances

The Fourier truncation and smoothing preserve continuous translations. The conserved momentum is

    P_total=P_packet + sum_j k_j(qc_j*Ps_j-qs_j*Pc_j).

Differentiating the field contribution gives minus the photon force. No fixed lattice or external boundary force absorbs the missing momentum. In addition, the nonzero-mode energy obeys E_modes>=v*|P_field|, which is checked numerically. This is momentum in the declared effective rest-frame model, not a proof of relativistic covariance or a gravitational stress-energy completion.

## 3. What the calculations show

The main scan uses L=8, v=0.5, smoothing width 0.2, initial n0=1, and initial homogeneous rolling rate 0.05 for duration 3, all in dimensionless units. Eighteen runs vary field inertia, initial packet energy and 16/32/64 modes. Initially all nonzero modes are empty.

For K=1 and 64 modes:

| Initial photon-packet energy | Fraction of photon energy transferred | Fraction of field energy gain in nonzero-mode quadratic energy |
| --- | ---: | ---: |
| 0.001 | 13.179% | 0.784% |
| 0.01 | 14.382% | 7.103% |
| 0.1 | 24.932% | 36.767% |

The initial rolling field has energy 0.01 in these cases. That pre-existing energy is not credited as photon-produced companion energy. In the middle case the photon loses 0.00143815; the homogeneous kinetic energy increases by 0.00133600 and the nonzero-mode energy increases by 0.000102152. Their sum equals the loss.

Radiation loading affects the field. Increasing the packet energy changes its fractional loss; increasing inertia suppresses this response. Consequently an intensity-independent distance-only stretching coefficient is not automatically derived by this candidate. Packet energy here is a classical loading parameter; this scan does not distinguish a change in photon number from a change in individual photon frequency and is not a measured color-dependence test.

With no initial rolling motion, the three corresponding fractional losses are about 0.175%, 1.721% and 15.053%. The photon itself excites a receiving field and spatial disturbances; approximately 68% of that field energy is in the nonzero-mode quadratic contribution. These controls establish a possible self-driven transfer in the toy dynamics, not the required astronomical rate or a universal exponential redshift law.

## 4. Energy flow must include the background cross term

The mode-energy column is **not** the fraction able to travel to another location. Fourier modes are a global decomposition, and energy flow against a moving background includes interference terms. Let n=1+u*t+delta_n, with initial rolling rate u. The energy density above the background and its flux are

    delta_e=K*u*delta_n_t
             + K/2[(delta_n_t)^2+v^2*(delta_n_x)^2],
    F=-K*v^2*(u+delta_n_t)*delta_n_x.

Omitting the terms linear in u would omit real energy transport. Therefore the follow-up diagnostic measures full energy flow through the boundaries of a fixed region around the initial photon position. It checks

    final excess field energy - initial excess field energy
      = work transferred inside the region - net boundary energy outflow.

For K=1, initial packet energy 0.01 and 64 modes, work deposited inside [-0.1,0.9] is 0.000264338, net boundary outflow is 0.0000994257, and final excess field energy is 0.000164912. The residual is below 4e-12. Energy through the left boundary is positive and outward, while the photon moves right. This demonstrates energy propagation away from the photon route in this one-dimensional candidate; it does not demonstrate capture by a gravity well. The net outward flux can be smaller than one boundary's outward contribution because the other boundary can admit energy produced elsewhere.

The local check uses the actual smoothed source and the full field energy flux. It does not count only nonzero Fourier energy as a portable reservoir. Periodic boundaries, background evolution and interaction support remain part of the model; no outgoing-wave fraction at infinity is measured.

## 5. Verification, limitations and next work

Across the main scan and zero-rolling controls, maximum relative total-energy drift is below 4e-15 and momentum drift below 2e-15. Changing from 32 to 64 modes changes the final photon energy by less than 6e-16 of its initial energy in the tested smooth cases. Translating the photon within the periodic box leaves energy and momentum results unchanged to about 1e-19. The largest local transport-balance residual is below 6e-9 of initial packet energy. The no-photon control creates no wave energy and leaves the rolling field energy unchanged. These are finite-domain numerical checks, not universal accuracy or stability bounds.

The useful advance is a closed candidate interaction with photon feedback, a specified energy receiver, spatial field response, energy flux and momentum conservation. It goes beyond attaching an energy ledger to a prescribed redshift profile. But it does not yet derive the whole-signal void law or the observed rate. The background, positive-n domain, scalar identity, finite packet smoothing and one-dimensional geometry are explicit assumptions.

Next test how a localized field state can yield the required nearly common frequency/duration stretch for separate weak signals in the presence of ambient radiation, while preserving source/detector clock standards. Its void dependence, matter coupling, possible reverse transfer, companion capture, supported storage and joint rotation/lensing response must be derived. The companion receiver remains fixed; what is being tested is whether this candidate implements all the required behavior.

```sh
python research_work/results/companion-backreaction/check.py
```

Requires NumPy and SciPy. The offline diagnostic suite includes this calculation. Green checks establish mathematical consistency only; the full research goal remains open.
