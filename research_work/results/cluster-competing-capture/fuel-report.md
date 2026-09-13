# Fuel accounting without fixing the universe age

## Result

For light funded by hydrogen-to-helium fusion, each unit of deposited energy E/c^2 requires at least 140.48 units of burned hydrogen mass in the limiting case where all nuclear release becomes target deposits. At 10% net delivery this becomes 1,404.79 units. The fuel can be anywhere in the contributing source population; it is not restricted to stars inside the receiving cluster.

This does not establish an actual energy shortage. We have not inventoried the contributing hydrogen or reconstructed its history. It does establish an age-independent cost that hypothetical old emission must pay. A longer duration cannot increase the energy released by burning the same fuel from the same initial to final nuclear state.

## Measured input and known derivation

The NIST isotopic tables give neutral atomic masses [hydrogen-1](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=H) 1.00782503223 u and [helium-4](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=He) 4.00260325413 u. Using neutral endpoints consistently includes the associated electron/positron accounting in the net reaction. The established mass-defect relation gives

    epsilon_HHe=(4 m_H-m_He)/(4 m_H)=0.00711851608.

This is approximately 0.712% of the mass of the hydrogen actually processed, not that percentage of every star's entire mass. The result is total available reaction energy; treating all of it as photons is an optimistic upper bound because other channels, including neutrinos, carry energy.

Let F be the fraction of this available energy ultimately retained in the selected target through the photon/companion channel. It includes radiative fraction, conversion, geometric interception, competing capture and retention, so 0<=F<=1. Then known energy conservation requires

    E_deposit <= F epsilon_HHe M_H_burned c^2
    M_H_burned >= (E_deposit/c^2)/(F epsilon_HHe).

This is known nuclear physics and bookkeeping applied to the hypothetical transport, not a unique photon-to-gravity law. The E/c^2 quantity is energy-equivalent mass, not automatically an observed lensing mass. A proposed enhanced gravitational response needs its own consistent field equations and cannot be inferred from this inequality.

| Net target fraction F | Minimum burned hydrogen / deposited energy-equivalent mass |
|---|---:|
| 1 | 140.48 |
| 0.5 | 280.96 |
| 0.1 | 1,404.79 |
| 0.01 | 14,047.87 |

F=1 is a limiting bound, not a measured or generally achievable efficiency. The test explicitly accounts for helium products and energy outside the target; it does not make undelivered energy disappear. No capture efficiency is fitted here.

## Connection to the previous history test

The earlier constructed histories had identical current bolometric radiation despite different stored energy. Under this hydrogen-burning funding assumption, their minimum fuel costs differ too. The reference's two emitted energy units require 280.96 mass units (where the mass unit is one energy unit/c^2). The most extreme alternative's 9,959.23 emitted units require at least 1,399,060.11 mass units. Matching today's radiation therefore does not eliminate the histories' vastly different material requirements.

Fuel products may be retained in stars, ejected, or processed into heavier nuclei. Counting helium alone without tracking those channels is insufficient. No primordial helium abundance, Big Bang origin, fixed initial composition, universe age or universe size is imposed here. Initial composition and imported/exported matter remain explicit unknowns.

This bound applies specifically to the H-to-He funding channel. Later nuclear burning, gravitational contraction, accretion or another energy source require their own budgets; they are not excluded by the table. Recycling helium back into hydrogen requires energy and cannot create an additional free release while retaining the first release in permanent deposits. The system must not count both the source's unreduced initial mass and its exported energy as new total energy.

## Next data requirement

The external-source question should be addressed by an energy inventory over the contributing population, with present stellar/remnant/gas estimates, possible compositions and physically allowed histories, coupled to the competing-receiver transport. Receiver cross-section and source abundance alone do not supply that inventory. Actual motion/lensing targets then constrain the same stored field. No actual cluster fuel sufficiency, mass fit, or full first-principles companion model is established. All six objectives remain open and observational holdouts remain unopened.

Run `python research_work/results/cluster-competing-capture/fuel.py`. Results retain the laboratory inputs, source links, four delivery fractions, energy balances and costs of the previous synthetic histories.
