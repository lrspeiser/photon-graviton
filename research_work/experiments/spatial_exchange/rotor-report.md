# SE-Rot: finite rotor-to-curl interaction checks

All 96 local and 48 grid interaction fixtures pass.

| Check | Maximum error |
|---|---:|
| Local mass derivative | 1.00348004e-10 |
| Local rotation/torque identity | 3.47621661e-18 |
| Joint grid/position/internal derivative | 9.14527334e-11 |
| Discrete curl adjoint identity | 8.22387426e-18 |

The positive internal energy uses K=P-epsilon Bbar cross Q, Bbar the sampled
curl of A. Its derivative supplies a curl-shaped field source and reciprocal
position/internal reactions. Canonical source spin is Q cross P. This adds
a finite, declared angular reservoir rather than creating rotation for free.

This is an interaction derivation and implementation check, not a full field
evolution. Existing SE-1/SE-2 simulations have not been changed. The next step
is to embed this mass in the common Hamiltonian, include internal angular
momentum in the full ledger, validate all derivatives and propagation, then
evolve from empty fields with spin-reversal and zero-spin controls.

Its derivative coupling can alter field response near the regularized source;
old cone/stability tests do not automatically validate it. More circulation
would also not by itself establish an inward force or stronger cluster lens.

## Attribution

| Ingredient | Status |
|---|---|
| Shifted-momentum square and Hamiltonian derivatives | Established mathematics; not claimed new |
| Reciprocal spin/rotor angular exchange | Prior art, including Banerjee and Fehske |
| Spin/field torque accounting | Prior art, including Turcati et al. |
| This regularized curl coupling inside the candidate internal mass | Specific hypothesis; originality unestablished |
| Successful derivative/torque checks | Demonstrated for the declared fixtures only |

Primary abstracts reviewed: [Banerjee and Fehske](https://arxiv.org/abs/2604.23768v2),
[Turcati et al.](https://arxiv.org/abs/1907.12196v2). These related non-gravity
models are attribution context, not gravity acceptance targets. No claim is
made that those papers use our specific equations.

Raw results and source/protocol hashes: rotor-checks-v1/results.json.
All twelve original research requirements remain active.
