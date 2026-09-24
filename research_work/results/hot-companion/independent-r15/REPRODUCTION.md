# Independent calculation, round 15: reproduction record

`first_principles_test/` is the independent calculation "Motion-opened radiation from ordinary matter's internal
oscillations" (24 September 2026), kept exactly as received: its README, scripts, requirements, manifest, figure and
result files. Every file's SHA-256 matches the manifest.

## Reproduction in this repository's environment

Python 3.11, numpy and scipy as installed here, numba 0.67.0 (the calculation lists 0.65.1).

    cd first_principles_test
    python finite_reservoir_test.py
    python test_dark_bright.py --out results --part static
    python test_dark_bright.py --out results --part collisions
    python test_dark_bright.py --out results --part controls
    python test_dark_bright.py --out results --part no-sink

The outputs are in `reproduced/`.

| Output | Against the received file |
|---|---|
| `finite_reservoir.json` | byte-identical |
| `static.json` | byte-identical |
| `collisions.json` | byte-identical |
| `controls.json` | byte-identical |
| `no-sink.json` | every number identical; the wrapper differs, as the calculation's README says |

The finite-reservoir run takes about 17 seconds, and the four other parts under 30 seconds each.

## What round 15 did with it

Three parts are in `../README.md` §25:

1. `code/wave_dark_v15.py` tests the calculation's postulated coupling against emitters that talk only through the wave.
2. `code/receivers_v15.py` tests which test bodies a passing wave pulls, including the calculation's gain-reservoir
   oscillator.
3. `code/reservoir_force_v15.py` is the integration test the calculation's README asks for. The pieces' internal modes
   radiate real waves. Test bodies feel only the local wave. Force, timing, energy and momentum are each measured on
   their own.
