#!/usr/bin/env python3
"""Load the exact JR-1 R10 forward model and expose a single-lens nuisance interface.

No fitting or observational scoring happens here. This module only reconstructs
the already published R10 response (JR-1 run.py plus the population_followup.py
patch) and adds a uniform per-lens nuisance handle so that every one of the six
systems can be treated identically, instead of the JR-1 split where only the four
optimizer-fitted lenses carried individual offsets.

The ten universal field/source coefficients and the two shared population scales
are never re-fitted by this package. They stay at their published R10 values.
"""
from __future__ import annotations
import importlib.util, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
JRDIR = ROOT / 'research_work/results/joint-response-iteration'


def _load_jr1():
    """Reproduce JR-1 run.py under the R10 spheroid-population patch."""
    sys.path.insert(0, str(JRDIR))
    spec = importlib.util.spec_from_file_location('run', JRDIR / 'run.py')
    J = importlib.util.module_from_spec(spec)
    sys.modules['run'] = J
    spec.loader.exec_module(J)
    base_source, base_lens = J.source_parameters, J.Lens
    J.DEFAULTS['logusph'] = 0.
    J.BOUNDS['logusph'] = (-.3, .4)

    def source(Mstar, Mgas, Re, p, saturation=False, spheroid_fraction=0.):
        f = np.asarray(spheroid_fraction)
        u = 10 ** p.get('logusph', 0.)
        factor = 1 + (u - 1) * f
        return base_source(np.asarray(Mstar) * factor, Mgas, Re, p, saturation, u * f / factor)

    class R10Lens(base_lens):
        """R10 lens with one explicit (stellar offset, anisotropy) nuisance pair."""

        def local_parameters(self, p, spec):
            local = dict(p)
            local['logu'] = p['logu'] + p['logusph'] + p.get('dm_local', 0.)
            local['beta'] = p.get('beta_local', 0.)
            local['logusph'] = 0.
            return local

    J.source_parameters = source
    J.Lens = R10Lens
    return J


J = _load_jr1()
SPEC = dict(J.SPECS[9])
SPEC['name'] = 'R10_frozen_universal'
SPEC['free'] = []
RESULTS = json.loads((JRDIR / 'RESULTS.json').read_text())
# The ten universal field/source coefficients plus the two shared population
# scales. The four published per-lens dm_*/b_* entries are deliberately dropped:
# this package re-derives an equivalent pair from one observable at a time.
UNIVERSAL_KEYS = ('logA', 'logc', 'p', 'dA', 'dc', 'q', 'logt', 'logsph',
                  'core_sph', 'q_sph', 'logu', 'logusph')
R10_UNIVERSAL = {k: RESULTS['R10_parameters'][k] for k in UNIVERSAL_KEYS}
LENSES = list(J.LENSES)
def _jr1_nuisance():
    """JR-1's own convention: individual offsets for the four optimizer-fitted
    lenses, and dm=0 with the mean fitted anisotropy for the two excluded ones."""
    key = lambda n: n.replace('-', '_').replace('+', 'p')
    fitted = {n: dict(dm=RESULTS['R10_parameters']['dm_' + key(n)],
                      beta=RESULTS['R10_parameters']['b_' + key(n)]) for n in J.TRAIN_LENSES}
    mean_beta = float(np.mean([v['beta'] for v in fitted.values()]))
    out = dict(fitted)
    for n in J.TRANSFER_LENSES:
        out[n] = dict(dm=0., beta=mean_beta)
    return out


JR1_NUISANCE = _jr1_nuisance()
BETA_PRIOR_SIGMA = 0.30          # JR-1 residual penalty width on orbital anisotropy
BETA_BOUNDS = (-1.0, 0.35)       # JR-1 bounds


def params(dm=0., beta=0., companion=True):
    """Universal R10 coefficients plus one lens's ordinary nuisance pair."""
    p = J.DEFAULTS | R10_UNIVERSAL
    p['dm_local'], p['beta_local'] = float(dm), float(beta)
    return p


def make_lens(name, n=1601, order=64, deproj_order=128):
    return J.Lens(name, n, order, deproj_order)


def vrms(lens, dm, beta, companion=True):
    """Seeing-convolved annular V_rms in km/s, or None where the moment is invalid."""
    try:
        v, _, _ = lens.prediction(params(dm, beta), SPEC, baryons_only=not companion)
    except ValueError:
        return None
    return None if not np.all(np.isfinite(v)) else v


def einstein_angle(lens, dm, companion=True):
    """Predicted Einstein angle in arcsec. Independent of orbital anisotropy."""
    return float(lens.angle(params(dm, 0.), SPEC, baryons_only=not companion))
