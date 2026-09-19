# RUT-1 stage 8, amendment 1: a replacement for the ring-limit gate, declared after the declared one failed

**This amendment is written after a declared gate failed, and it says so.** Stage 7's amendment came before
its campaign. This one does not: parts L and P of stage 8 have been run as declared and their archive is
pushed (aaf0bfe), with gate L3 recorded as **failed** and the numerical verification status as failed. That
record is not changed by anything here. `protocol-rut8.md` is not edited.

Part M — the three new simulations — does not depend on L3 and proceeds as declared.

## What L3 required, and how it fell

Three narrow cold annuli against the exact-history ring theory of stage 5, evaluated at each annulus's own
mean radius and total writing rate:

| required | result |
|---|---|
| the relative difference in m = 2 growth rate falls at every step | 23.35%, 4.67%, 0.63%: **met** |
| the last difference is under 2% | 0.63%: **met** |
| the negative control — the ring at stage 5's radius, R = 1 — fails the 2% | it is off by 72%: **rejected, as required** |
| the fitted order of convergence in the radial width lies between 1.5 and 2.5 | **2.507: not met** |

So L3 failed as declared, by 0.007 in its fourth criterion, with the two formulations agreeing *better* than
the gate required.

## Why the criterion was wrong, and whose mistake it was

Mine. I prototyped the ring limit on a coarser orbit quadrature (64 × 32) than the one I then declared for
the stage (192 × 64). There the differences were 23.35%, 4.70% and 1.03%, the three-point order was 2.16, and
the orders between successive pairs of a four-annulus sequence were 2.2, 2.2 and 2.5 — which I wrote into the
protocol, and then declared 2.5 as a ceiling anyway. At the declared quadrature the narrowest annulus agrees
with the ring better (0.63% against 1.03%), so the fitted order rose past the ceiling.

Two things are wrong with the criterion itself, beyond how it was sized:

* **A ceiling on a convergence order tests nothing.** Converging faster than expected is not a defect of
  either calculation. What a limit check has to establish is that the difference falls, that it ends small,
  and that it falls at least as fast as a consistent method should.
* **At the 0.5% level the difference is no longer only the finite-width effect.** The narrowest annulus is
  0.0066 wide on a radial grid of spacing 0.00046, fourteen nodes across its width, so the population
  solver's own discretization error is of the order of the difference being fitted. An order fitted through
  that point is not an asymptotic order.

This is the third time in two stages that I have sized a check on something other than what it was then run
on — stage 7's two negative controls on the warm prototype, its mode window on the cold ring, and now this.
The habit that would have prevented all three is the same: run the declared configuration itself before
declaring a number about it.

## L3b, the replacement — declared now, before it is run

**Four** annuli: the declared three, and a fourth at (dL, dE) = (0.005, 0.0000625) on 4,800 radial nodes,
solved like the others for a 9.99% mean support, at the declared discretization. The fourth has **not been
run at the declared quadrature**; in the 64 × 32 prototype it differed from the ring by 0.18%.

| gate | requirement | negative control the gate must reject |
|---|---|---|
| L3b the ring limit | the relative difference in m = 2 growth rate against the exact-history ring, at the annulus's own mean radius and total writing rate, falls at every one of the three steps; the last is under 0.5%; and the order fitted through all four is at least 1.5 | the ring evaluated at R = 1 must differ from the narrowest annulus by more than 2%, as before |

L3b is stricter than L3 where L3 was meaningful — a fourth, narrower annulus, and a final difference four
times smaller — and drops only the ceiling.

**How the two are reported.** L3 stays in the archive as declared and as failed. L3b is reported beside it.
The numerical verification status names both: L3 failed as declared, superseded by amendment 1; L3b, its
outcome. The status is "passed" only if every other gate passes and L3b passes, and whenever it is quoted the
L3 failure is quoted with it. If L3b fails, the status is failed and nothing further is declared in its
place without the owner.

## Declared thresholds for L3b (read by the driver)

```json
{
  "L3b": {"ring_limit": [[0.04, 0.004, 600], [0.02, 0.001, 1200], [0.01, 0.00025, 2400], [0.005, 0.0000625, 4800]],
          "last": 0.005, "order_min": 1.5, "control_min": 0.02}
}
```
