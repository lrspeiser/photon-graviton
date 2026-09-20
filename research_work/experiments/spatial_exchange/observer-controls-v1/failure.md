# First execution: JSON output failure

Source a88eb7f completed the calculations but failed before writing results.json:
`TypeError: Object of type bool is not JSON serializable` (NumPy boolean).
The uniform control's pass flag was not converted to a native Python bool.
No numerical pass is claimed for this unarchived execution. The second execution
uses observer-controls-v2 and changes only that cast and the output directory.
Equations, initial fixtures, tolerances and gates remain unchanged.
