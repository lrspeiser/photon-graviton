"""Opt-in corrected WarmAnnulus for new runs; historical rut6 imports stay untouched."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

from .sampling import VERSION, build_distribution, sample, verify_distribution

ARCHIVED_PATH = Path(__file__).resolve().parents[1] / 'results/path-memory/equilibrium.py'
_name = '_photon_graviton_equilibrium_for_sampling_v2'
if _name not in sys.modules:
    _spec = importlib.util.spec_from_file_location(_name, ARCHIVED_PATH)
    if _spec is None or _spec.loader is None:
        raise ImportError(f'Cannot load archived equilibrium module: {ARCHIVED_PATH}')
    _module = importlib.util.module_from_spec(_spec)
    sys.modules[_name] = _module
    try:
        _spec.loader.exec_module(_module)
    except Exception:
        sys.modules.pop(_name, None)
        raise
ArchivedWarmAnnulus = sys.modules[_name].WarmAnnulus


class WarmAnnulus(ArchivedWarmAnnulus):
    """Same equilibrium solver and density; corrected sample() for follow-on experiments.

    Inherited cutoffs, source normalization and continuum-convergence limitations
    are NOT repaired here. Do not substitute into historical rut6 replays.
    """
    sampling_version = VERSION

    def sampling_distribution(self, n_L=160, n_vr=96, reach=5.):
        return build_distribution(self, n_L, n_vr, reach)

    def verify_sampling(self, n_L=160, n_vr=96, reach=5.):
        return verify_distribution(self, n_L=n_L, n_vr=n_vr, reach=reach)

    def sample(self, n, seed=0, n_L=160, n_vr=96, reach=5.):
        return sample(self, n, seed, n_L, n_vr, reach)
