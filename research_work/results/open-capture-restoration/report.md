# Open capture can restore a clock state, but does not yet give a universal redshift law

We implemented a capture/restoration rule with an explicit available-energy condition and let the radiation leave the domain. The receiver's clock state can return close to its initial value while the well gains deposited energy. This resolves one operational question: restoration can be modeled without simply resetting the endpoint by hand. It does not resolve the model's dependence on source position, emission time and initial field seed.

## Rule and energy derivation

Postulate stored state energy M=K(n-1), a traveling receiving sector T, and reference radiation energy R=G+H. Let

    r = Gamma*(n-1)*max(0,1-R/(K*n))
    n_t = T-r
    h = n_t/n.

The max factor is an analyst-chosen energy-availability law, not a fit or a unique physical derivation. It prevents this restoration channel from spending more than the specified state sector can supply to radiation. Use the same field state for radiation propagation and local clocks as before.

The open continuity laws are

    G_t+(vG)_x = -hG,   H_t+(vH)_x = -hH
    T_t+(vT)_x = (R/n-K)*T-Gamma*T
    D_t = Gamma*T+(K-R/n)*r,   v=1/n.

Adding M_t=K(T-r) cancels the local transfer terms exactly. The total reference energy changes only through the boundary radiation/receiving flux. The restoration deposit term is nonnegative by construction. These conservation and positivity statements are mathematical consequences of the postulates, not an action or momentum derivation. Stationary field memory is included as energy, not left free.

## Demonstration setup

We retain a [0,10] domain, initial radiation G=0.2,H=0.02 and a 0.01-amplitude receiving seed near x=2. Incoming radiation switches off smoothly near time 2, allowing the region to clear. A prescribed well rate rises near x=8 and remains high to the receiving boundary. It is an imposed test well, not gravity derived from the deposited energy. Compare K=0.05 and 0.5 with restoration on/off, plus the empty-seed control. All values are nondimensional and none is calibrated to real redshift.

## Finest-grid results

| Case | Early boundary probe, x=0 to 10, t_e=1 | Probe beginning near seed, x=2 to 9, t_e=0.1 | Late boundary probe, x=0 to 10, t_e=5 |
|---|---:|---:|---:|
| Low state-energy cost, no restoration | -0.003255 | +0.107063 | -0.003255 |
| Low cost, restoration | -0.003255 | +0.107063 | -0.002841 |
| High cost, no restoration | -0.000097 | +0.027786 | -0.000097 |
| High cost, restoration | -0.000198 | +0.027552 | -0.000087 |

Entries are measured z with endpoint clocks included. The low-cost rows use 1280 cells; high-cost rows use 640. Negative values are blueshifts. The zero-seed case gives zero for all probes. Identical EM/GW arrivals are imposed by the shared ray law.

The positive interior-probe result is not evidence for a general redshift mechanism. Its emission occurs while the probe overlaps the seed, unlike the previously tested later interior emission. Position, timing, memory cost and source history differ from older runs; changes cannot be attributed to capture alone. Within the paired low-cost runs, restoration has effectively no effect on either early probe because sufficient radiation remains to close the restoration channel then. The late probe changes, but remains blueshifted.

At t=24 the low-cost receiver's n is 1.00326574 without restoration and 1.00000002545 with restoration. Total deposited energy rises from 0.0173424 to 0.0188044. Stored memory drops from 0.00442780 to 0.00293724. The remainder is accounted for through outgoing radiation/receiving energy; the deposit increment is not a second copy of the full released memory. Restoration operates after the receiver's radiation has largely left. Do not substitute the final t=24 clock state for its value at the earlier signal arrival.

## Remaining mismatch with the intended process

Some energy still remains in stationary memory away from the well after the traveling sector has left. Calling it 'memory' does not satisfy the user's requirement that converted energy travel onward instead of permanently collecting in voids. This provisional closure makes that issue visible. It needs a physical memory-transport or relaxation channel with its energy and force accounted for, not a relabeling of the residual storage.

An energy-funded reset therefore improves the specification, but does not yet explain the widespread distance-redshift trend. The next postulate to examine is how field state and its stored energy move or relax together, and how the initial seed is supplied across real paths. A revision must preserve the observable clock factors and the energy accounting rather than optimizing selected probe paths. The well's gravity, lensing and momentum response are still not calculated by this one-dimensional test.

## Verification and limits

Seventeen integrations cover the prescribed controls and refinements. The two low-cost cases initially failed the 0.003 absolute-z refinement gate at 640 cells and were continued to 1280. Their largest new change is 0.000956, passing the original gate; high-cost and empty-seed cases already pass at 640. Failed earlier comparisons remain archived. This is an initial absolute-error check, not a claim of high relative precision for tiny z values.

Reference energy including memory and boundary fluxes closes to 2.05e-14 absolute. Energy sectors stay nonnegative within the declared tolerance; the deposit-rate formula is analytically nonnegative for nonnegative T and n>=1. Independent proper-clock intervals agree with the phase/event relation within 6.23e-6 relative. This does not establish a full conserved stress-energy description or observational agreement.

Reproduce with `python research_work/results/open-capture-restoration/run.py`, followed by the same command with `--refine-low`. Both commands are required for the reported finest low-cost cases. No observational parameters, alpha, unseen data or source-energy-supply conclusion changed. This advances the operational energy/restoration analysis, while all nine research goals remain incomplete.
