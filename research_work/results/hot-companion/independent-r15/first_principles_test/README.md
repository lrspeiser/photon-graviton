# Motion-opened radiation from ordinary matter's internal oscillations

## Status and scope

Independent exploratory calculation performed on 24 September 2026. This is a test of a possible microscopic source mechanism, not a derivation or validation of the complete streaming-gravity law. The astrophysical regression suite was not run and the user's repository, blog, and adopted law were not changed.

The round-14 report was checked against GitHub commit `242eb1dcdf1f06b91ada4e428bfade333c5fcd0a`, including `coherent_force_v13.py`, `singer_v14.py`, and the relevant BLOG.md section. The existing power-balance implementation sets a target phase using an arcsine of a clipped supply-minus-loss expression. This is an energy-motivated closure, not a dynamically integrated reservoir.

## Three calculations

### 1. A genuine gain-reservoir emitter changes from supplying to absorbing

Let a be a real, phase-locked oscillation envelope and n its gain-reservoir energy. Let s be an incident wave amplitude. The local equations used were

    da/dt = (g*n - gamma_i - gamma_e)*a + sqrt(2*gamma_e)*s
    dn/dt = P - gamma_R*n - 2*g*n*a*a
    s_out = s - sqrt(2*gamma_e)*a.

They imply

    d(a*a+n)/dt = P - gamma_R*n - 2*gamma_i*a*a + s*s - s_out*s_out.

Positive outgoing-minus-incoming power means supplying the wave; negative means absorbing. No feed/absorb switch or target-phase arcsine appears.

With P=4, g=1, gamma_R=1, gamma_i=0.1, gamma_e=1, input amplitudes 0, 1, 3, 4, 10 gave net outgoing powers 2.6364, 2.7242, 0.2446, -2.0489, -29.7909. The analytic crossover amplitude is 3.122499. This result does not establish attraction, the astrophysical screening factor, or a physical relation between background acceleration and incident amplitude.

### 2. Why a fixed-supply version is not enough

A pumped quiet mode and three radiating modes, with internal losses, produced approximately fourfold additional radiation per doubling of detuning. It did so by redirecting power that otherwise went to internal losses. Its chosen cold radiative output was only 0.00029996 of a pump of 4 (one part in about 13,335).

A negative control removed every nonradiative loss. With fixed pump 0.04, the stationary total radiative output stayed 0.04 across every tested detuning, although power moved between radiation channels. This agrees with the balance equation: constant pump plus no accumulation or internal loss fixes the total luminosity.

Thus this version cannot simultaneously identify its total pump with the old cold escaping power and obtain arbitrarily larger steady hot output.

### 3. Finite internal-energy reservoir: the better source candidate

Two oscillations can couple to an outgoing wave through their sum. Their antisymmetric combination D cancels in that channel, while the symmetric combination B radiates. Opposite small frequency shifts mix D and B linearly. Three radiating channels model three directions.

    dD/dt   = -gamma_0*D - i*sum(delta_j*B_j)
    dB_j/dt = -gamma*B_j - i*delta_j*D

Here D is internal energy of the ordinary source, not additional unseen matter. The field amplitudes are normalized so that

    E = |D|^2 + sum(|B_j|^2)
    dE/dt = -2*gamma_0*|D|^2 - 2*gamma*sum(|B_j|^2).

The energy leaving the modes is explicitly accumulated as radiation. There is no pump or unrecorded internal loss in this refinement. The carrier frequency is removed in this envelope description; this does not calculate the translational work of the material particles.

For small detuning and slowly depleted D, B_j approximately equals -i*delta_j*D/gamma. If delta_j = chi*w_j and the three w components have variance sigma^2, then

    P_motion / P_cold approximately = 3*chi^2*sigma^2/(gamma_0*gamma)
                                  = 3*sigma^2/u_eff^2,
    u_eff^2 = gamma_0*gamma/chi^2.

The exponent arises from linear amplitude mixing followed by quadratic energy; it is not a sigma-squared term inserted into a power source.

For exponentially correlated detuning with direction-randomization rate nu,

    P_motion(nu)/P_motion(0) approximately = gamma/(gamma+nu).

This limit assumes motion changes detuning without directly destroying internal phase coherence. It does not justify treating every collisional medium identically.

## New numerical results

`finite_reservoir_test.py` used gamma=1, gamma_0=1e-6, initial stored energy 1, 8 seeds and 64 independent local realizations per seed (512 per condition). The default run lasts 80 response-time units, with the first 10 discarded for rate measurements and dt=0.005. The equations have no gravitational force term. Detuning is prescribed; particle trajectories are not evolved.

q is the one-component rms detuning divided by gamma. It is proportional to relative speed only under the stated linear detuning assumption. It is not the old code's re-timing speed parameter.

| q | Speed relative to q=0.001 | Additional radiated power | Normalized additional power | Ratio to previous doubled speed |
|---|---:|---:|---:|---:|
| 0.001 | 1 | 5.94556e-6 | 1.00000 | — |
| 0.002 | 2 | 2.37503e-5 | 3.99463 | 3.99463 |
| 0.004 | 4 | 9.44931e-5 | 15.89306 | 3.97861 |
| 0.008 | 8 | 3.70031e-4 | 62.23649 | 3.91595 |
| 0.016 | 16 | 1.36396e-3 | 229.40857 | 3.68608 |
| 0.032 | 32 | 4.10051e-3 | 689.67701 | 3.00633 |

Normalization is to measured output at q=0.001. Raw JSON and CSV are authoritative; the last two normalized values are rounded narrative values. The larger-speed runs deplete their reservoirs more during the same observation window. The small-perturbation law is not an unlimited luminosity law.

At q=0.004, direction-randomization rates nu/gamma = 1, 10, 20 left 50.70%, 9.23%, and 4.84% of the freely moving bright-channel power. The weak-coupling predictions are 50%, 9.09%, and 4.76%. These comparisons have finite-ensemble uncertainty; the fixed free velocity sample is reused in paired speed comparisons.

After stopping detuning at time 40, the extra bright output decayed below 1e-55 in the tested runs. Some internal energy was permanently spent; this is ordinary fuel depletion, not persistent phase locking. A cold state at the remaining source energy is not identical in total mass/energy to the original source.

Checks included a matrix-exponential comparison (maximum state error 3.22e-15), timestep refinement (0.01 to 0.0025 for a collision case, about 0.14% change in mean bright power with independently sampled paths), and integrated modal-energy accounting (worst residual across finite-reservoir runs and refinement about 1.7e-12 of initial energy).

A separate 100,000-pair kinematic test used detuning proportional to radial relative velocity n dot (v_1-v_2). A uniform velocity boost changed it by at most 1.78e-15 and rigid rotation gave at most 2.18e-15. Independent Gaussian velocities with unit one-component dispersion gave radial relative variance 1.99953, against 2. This verifies the algebra of the chosen coupling; it does not derive that coupling or show that differential rotation, shear, or a multi-stream distribution is treated correctly by real matter.

## Meaning for the streaming-gravity theory

This is a possible derivation of the heat-source exponent and collision suppression, plus a separate dynamic example of feeding weak waves and absorbing strong ones. It is not one finished microscopic interaction delivering the whole theory.

The finite-reservoir version changes the interpretation of the law's power: ell would be a cold leakage rate per unit mass, while motion opens additional conversion of internal matter energy. Keeping ell as an absolute fixed total supply is incompatible with calling the larger stationary luminosity merely a release of that same already-escaping supply.

If internal energy is identified with M*c^2, then gamma_0=ell/(2*c^2). The coefficient u_eff still depends on an independently unknown detuning coupling and response rate; its equality to the separately fitted transport speed 169 km/s has not been established. Nor have attraction, its radial and mass dependence, the spatial source integral, light bending, strong-field screening, or collision transport memory been recovered.

The natural next integration test is to use emitted waves from these evolving internal modes to drive test emitters, derive forces from a common local interaction, and measure force, phase, energy and momentum separately. Do not convert sqrt(radiated power) into a force in the implementation.

## Reproduce

Install the versions in `requirements.txt` in an appropriate Python environment, then run:

    python test_dark_bright.py --out results --part static
    python test_dark_bright.py --out results --part collisions
    python test_dark_bright.py --out results --part controls
    python test_dark_bright.py --out results --part no-sink
    python finite_reservoir_test.py

`results/no_internal_sink.json` is the original no-sink control output; the `--part no-sink` command writes the same experiment in `results/no-sink.json` with the general script's wrapper. JSON files preserve the measured numbers and seeds. `manifest.json` records provenance. No original private repository code is redistributed.

## Background, not evidence for gravity

The cancellation and mode-coupling construction uses general wave-system mathematics. Relevant primary examples are Koshelev et al., Physical Review Letters 121, 193903 (2018), arXiv:1809.00330, and Fan, Suh and Joannopoulos, Journal of the Optical Society of America A 20, 569–572 (2003). Neither supports this candidate as a theory of gravity; they support the availability of the mathematical ingredients.
