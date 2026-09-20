#!/usr/bin/env python3
"""Run the Stage 6B detuned finite-electromagnetic dynamics campaign."""
from __future__ import annotations

import finite_em_dynamics_impl as implementation

_original_spectral_residue = implementation.spectral_residue


def _spectral_residue_with_compatibility_key(*args, **kwargs):
    result = _original_spectral_residue(*args, **kwargs)
    result["scalar_residue_in_photon_branch"] = result[
        "scalar_residue_in_first_photon_branch"
    ]
    return result


implementation.spectral_residue = _spectral_residue_with_compatibility_key


if __name__ == "__main__":
    raise SystemExit(implementation.main())
