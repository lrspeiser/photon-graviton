# First execution: output serialization failure

Source commit 1fcc6e2 executed all fixtures but failed before writing results.json:
`TypeError: Object of type int64 is not JSON serializable`.
The category counts inherited NumPy integer types. No result file was emitted,
and no numerical pass is claimed for this invocation. Preserve this directory.
The next execution uses mode-response-v2 and explicit Python int conversion
for three counts. Equations, fixtures, finite-difference step and gates are unchanged.
