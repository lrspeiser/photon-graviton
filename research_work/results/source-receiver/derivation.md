# Use the converter's spectrum as the receiver's input

## Why this calculation is necessary

We have separately tested a converter that produces companion waves and a receiver that absorbs them. Their favorable properties cannot be chosen independently. Here we use the converter's predicted energy distribution to drive the receiver and account for energy remaining in transit, received and escaped. Without this link, a capture model could rely on a spectrum the source never produces.

This is a conditional spectral-transport comparison, not a new successful redshift theory. The original point-oscillator converter still has its demonstrated color dependence and broad photon-angle problem. Nothing in this capture calculation repairs those production failures. It develops a reusable source-to-receiver test for the larger program.

## 1. Specify the regime in which separate responses can be linked

Take the [point-oscillator conversion amplitude](../oscillator-response/derivation.md) and the [scalar receiving response](../broadband-capture/derivation.md). Treat emission and reception as spatially separated, incoherent, weak processes with an effectively isotropic companion bath. This is a kinetic approximation, not multiplication of coherent collective amplitudes. Field-mediated correlations, near-field interactions and multiple-scattering changes to transport are omitted and would have to be small for this approximation to apply.

The source photon energy E is fixed for each comparison; no stellar spectrum, physical target abundance or absolute conversion normalization is inferred. A distribution of stellar photon energies would require integrating the differential source against its number flux. We keep frequency-independent escape and no travel redshifting as explicit diagnostic assumptions, not established companion properties.

The one-zone escape time below is prescribed. Especially at strong capture, it is not a derived spherical escape law, geometry or transport solution. A spatially resolved model must replace it before galaxy predictions. No independent dark-matter component, expansion or cosmic age is assumed.

## 2. Source energy spectrum follows the conversion rate

For one incident photon energy E, the earlier differential conversion rate is proportional to

    d sigma/d omega proportional to E(E-omega)^3 omega / D_em(omega)
    D_em=(Omega_em^2-omega^2)^2+gamma_em^2 omega^2.

Multiplying by companion energy omega gives its spectral energy supply. Dropping an overall normalization constant at fixed E,

    S(omega) = A omega^2 (1-omega/E)^3 / D_em(omega), 0<omega<E.

S is energy supplied per volume, time and companion-frequency interval. It is not companion number density. A contains photon number flux and the conversion amplitude; it has not been calibrated to the stellar energy budget. For variable source photon energies, the omitted E-dependent factor must be restored before averaging: the source includes an integral of photon number flux times E(E-omega)^3.

This is the emitted companion channel. If extra damping channels are present in the converter, their direct receiving energy and altered branching fractions must also be included in the full photon ledger. We do not assume all photon loss becomes this S.

The receiver has cross section sigma_abs proportional to 1/D_rec, where D_rec has its own resonance and total damping. Its continuum receives energy but has not been shown to store it permanently. Equal total response widths below need not mean equal branching fractions or equal microscopic couplings at emitter and receiver.

## 3. A frequency-resolved kinetic energy ledger

Let Gamma_abs(omega)=n_receiver c sigma_abs(omega). For a homogeneous volume with escape time t_esc,

    du_omega/dt = S_omega - [Gamma_abs(omega)+1/t_esc] u_omega
    dU_received,omega/dt = Gamma_abs(omega) u_omega
    dU_escaped,omega/dt = u_omega/t_esc.

Starting empty under steady supply, u_omega=S_omega[1-exp(-lambda_omega t)]/lambda_omega, with lambda_omega=Gamma_abs+1/t_esc. Thus

    u_omega + U_received,omega + U_escaped,omega = S_omega t.

The code verifies this time-dependent identity at four companion energies and compares the reservoir evolution with its analytic solution. The prescribed source remains an external input; a closed cosmic calculation still needs finite stellar fuel, source evolution and boundary conditions.

At steady state, u_omega=S_omega t_esc/[1+Gamma_abs t_esc]. Absorption therefore reshapes the bath when it is strong. Using the unattenuated production spectrum at all capture strengths would be inconsistent.

Define a=Gamma_abs(Omega_rec)t_esc and L_rec(omega)=sigma_abs(omega)/sigma_abs(Omega_rec). Then

    received fraction = integral S a L_rec/(1+a L_rec) / integral S
    escaped fraction  = integral S /(1+a L_rec) / integral S.

Their sum is one. The calculations record power received by the continuum, not permanent gravitating deposits. Heating, return emission, storage saturation and moving-receiver corrections are still absent from this energy-only steady state.

## 4. A matched-response identity

The low-speed absorption slowing coefficient derived previously is

    C = integral u sigma_abs [4+dlnsigma_abs/dlnomega]
        / [3 integral u sigma_abs].

In the weak-capture limit the bath spectrum is proportional to S. If emitter and receiver have equal resonance frequency and total width, their response denominators match. Write sigma for their shared frequency shape and q=(1-omega/E)^3. The weight is W=omega^2 q sigma^2. Integration by parts gives

    integral W (omega sigma'/sigma)
      = -3/2 integral W + 3/2 integral W omega/(E-omega),

with vanishing endpoints. Hence

    C = 5/6 + (1/2) weighted_average_W[omega/(E-omega)] > 5/6

for finite E. As E/Omega_em tends to infinity, C tends to 5/6. This does not require a thermal spectrum. The code checks the identity for three damping ratios, two finite photon energies and the asymptotic mathematical reference. The infinite-energy reference is not a physical infinite-energy source or an absolute rate extrapolation.

Thus matching the two resonances does not supply the earlier cancellation. Making the resonance narrow alone does not send this matched coefficient to zero. The single-energy cancellation and the integrated source-receiver spectrum are different limits; one cannot substitute for the other.

## 5. Capture strength changes efficiency and the spectrum together

For E=100 Omega_em, equal resonance frequencies and damping ratios gamma/Omega=0.1:

| a: resonant capture rate times escape time | Received power fraction | Absorption slowing C |
|---|---:|---:|
| 0.01 | 0.00500 | 0.83714 |
| 1 | 0.29452 | 0.75321 |
| 100 | 0.89674 | 0.43080 |
| 10000 | 0.98169 | 0.35417 |

Strong capture receives more of the supplied energy and changes the spectral weighting, but slowing remains positive in these examples. The a=0 entries in the saved data are formal weak-capture coefficient limits; at exactly zero capture, the actual absorbed power and force are zero.

Reradiated scalar scattering shares the receiver's spectral shape in the constant-damping comparison and contributes additional force. The table normalizes the absorption contribution to receiving power only. Full force per permanently stored energy requires the scattering and storage branching ratios. Strong scattering can also invalidate the prescribed escape-time approximation.

## 6. Detuning is a real change with multiple consequences

An exploratory scan also changes the receiver resonance. At a=1 and the same source example, receiver/emitter frequency ratios 0.5, 0.8, 1, 1.2 and 2 give different received fractions and slowing coefficients. Some coefficients are negative, so the source-weighted response can accelerate an already moving receiver rather than slow it in this driven approximation. Energy is supplied by the companion bath; this is not free acceleration or a thermal equilibrium prediction.

This scan holds each receiver's resonant capture rate times escape time fixed, not necessarily its physical number density. A real oscillator's resonant cross section changes with frequency and branching fractions. Therefore the table is a parameter comparison, not evidence that one can retune a single material at fixed abundance and obtain all those efficiencies. No detuning is selected as the final law, and no parameter is fitted to astronomical data.

## What is now joined, and what is not

The emitted companion energy spectrum, receiver frequency response, depletion by absorption and escape, and the resulting local low-speed drift coefficient are now evaluated together under explicit kinetic assumptions. This is stronger than assigning an independent favorable capture spectrum.

The production vertex's color/angle failures remain. Absolute source strength, a real stellar spectrum, causal spatial transport, receiver identity and capacity, all emitted and reverse channels, orbital support, self-gravity and lensing are not established. A candidate that resolves production can be tested with this frequency-resolved framework; a successful capture coefficient alone cannot promote the present converter to a complete theory.

Run `python -X utf8 research_work/results/source-receiver/check_source_receiver.py`. Saved numerical evidence is in `source-receiver-results.json`. All 20 tasks and 32 observational requirements remain in scope; these checks are mathematical diagnostics, not astronomical validation.
