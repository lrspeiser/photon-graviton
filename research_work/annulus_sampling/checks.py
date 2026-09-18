"""Database-free numerical checks. Run from repository root:

    python -m research_work.annulus_sampling.checks

Numerical failures exit nonzero; no archive, model score or physical gate is
rewritten. The historical H6 remains failed and is not evaluated here.
"""
from __future__ import annotations

import copy
from dataclasses import replace
import hashlib
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.special import gamma, gammainc

from .adapter import ARCHIVED_PATH, ArchivedWarmAnnulus, WarmAnnulus
from .sampling import build_distribution, radial_weights, sample, verify_distribution


class AnalyticMeasure:
    """Quadrature fixture with an independently known Sigma proportional to 1/R.

    Two L nodes at +/-1 have equal energy in Phi=-1/(2R^2), and E_circ=0.
    The radial-velocity integral has an incomplete-gamma closed form. This is
    an analytic sampling fixture, not a self-consistent physical annulus.
    """
    def __init__(self, nonuniform=False):
        self.r = np.geomspace(.4, 3., 41) if nonuniform else np.linspace(.4, 3., 41)
        self.C = np.zeros_like(self.r)
        self.sigma = np.zeros_like(self.r)
        self.L0, self.dL, self.dE, self.w, self.tau_keep = 0., 1., .7, .2, 10.

    def _grids(self, C, n_L=2, n_vr=64, reach=5., sigma=None):
        if n_L != 2:
            raise ValueError('the analytic fixture has exactly two L nodes')
        xg, wg = np.polynomial.legendre.leggauss(n_vr)
        return np.array([-1., 1.]), np.zeros(2), xg, wg

    def phi_total(self, r, C):
        return -.5/np.asarray(r)**2

    def surface_density(self, C, n_L=2, n_vr=64, reach=5., sigma=None):
        integral_vr = .5*(8*self.dE**2)**.25*gamma(.25)*gammainc(.25, reach**2/2)
        return 4*np.exp(-.5)*integral_vr/self.r


def fingerprint(a):
    h = hashlib.sha256()
    for name in ('r', 'C', 'sigma'):
        x = np.asarray(getattr(a, name))
        h.update(name.encode()); h.update(x.dtype.str.encode()); h.update(x.tobytes())
    h.update(repr((a.L0, a.dL, a.dE, a.w, a.tau_keep)).encode())
    return h.hexdigest()


class SamplerChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a = AnalyticMeasure()
        cls.kw = dict(n_L=2, n_vr=64, reach=5.)
        cls.nodes = build_distribution(cls.a, **cls.kw)

    def test_01_closed_form_density_marginal_and_mass(self):
        result = verify_distribution(self.a, self.nodes, **self.kw)
        self.assertTrue(result['numerical_verification_passed'], result)
        r = self.a.r
        expected = radial_weights(r)/(r[-1]-r[0])
        np.testing.assert_allclose(self.nodes.radial_marginal(), expected, rtol=0, atol=1e-12)
        self.assertAlmostEqual(self.nodes.moments()['radius']['mean'], (r[0]+r[-1])/2, places=12)

    def test_02_nonuniform_radial_grid(self):
        a = AnalyticMeasure(nonuniform=True)
        nodes = build_distribution(a, **self.kw)
        self.assertTrue(verify_distribution(a, nodes, **self.kw)['numerical_verification_passed'])
        np.testing.assert_allclose(nodes.radial_marginal(), radial_weights(a.r)/(a.r[-1]-a.r[0]),
                                   atol=1e-12, rtol=0)

    def test_03_deliberate_old_radius_bias_is_caught(self):
        biased = replace(self.nodes, node_mass=self.nodes.node_mass*self.nodes.node_radius)
        result = verify_distribution(self.a, biased, **self.kw)
        self.assertFalse(result['numerical_verification_passed'])
        self.assertGreater(result['radial_marginal_max_absolute_error'], 1e-3)
        moment = self.nodes.moments()['radius']
        measured = biased.moments()['radius']['mean'] - moment['mean']
        self.assertAlmostEqual(measured, moment['variance']/moment['mean'], places=12)

    def test_04_fixed_seed_draw_moments(self):
        n = 100_000
        x, v = self.nodes.draw(n, seed=1349)
        r = np.linalg.norm(x, axis=1)
        vr = np.sum(x*v, axis=1)/r
        ell = x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]
        energy = .5*np.sum(v*v, axis=1) + self.a.phi_total(r, None)
        draws = {'radius': r, 'angular_momentum': ell, 'energy': energy,
                 'radial_speed_squared': vr**2}
        for name, arr in draws.items():
            m = self.nodes.moments()[name]
            se = np.sqrt(m['variance']/n)
            self.assertLessEqual(abs(float(arr.mean())-m['mean']), 6*se+1e-12, name)
        self.assertLess(abs(float(vr.mean())), 6*np.sqrt(float(np.mean(vr**2))/n))
        self.assertLess(abs(float(np.mean(x[:, 0]/r))), 6/np.sqrt(2*n))
        self.assertLess(abs(float(np.mean(x[:, 1]/r))), 6/np.sqrt(2*n))

    def test_05_repeated_seed_and_particle_integrals(self):
        x, v = self.nodes.draw(512, seed=10)
        x2, v2 = self.nodes.draw(512, seed=10)
        np.testing.assert_array_equal(x, x2); np.testing.assert_array_equal(v, v2)
        ell = x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0]
        np.testing.assert_allclose(np.abs(ell), np.ones(len(x)), atol=2e-15, rtol=0)
        other, _ = self.nodes.draw(512, seed=11)
        self.assertFalse(np.array_equal(x, other))

    def test_06_source_is_not_mutated(self):
        before = fingerprint(self.a)
        build_distribution(self.a, **self.kw)
        sample(self.a, 32, seed=7, **self.kw)
        verify_distribution(self.a, **self.kw)
        self.assertEqual(before, fingerprint(self.a))
        self.assertFalse(self.nodes.node_mass.flags.writeable)

    def test_07_empty_draw(self):
        x, v = self.nodes.draw(0)
        self.assertEqual(x.shape, (0, 2)); self.assertEqual(v.shape, (0, 2))

    def test_08_invalid_particle_counts(self):
        for n in (-1, 2.5, True, np.nan):
            with self.subTest(n=n), self.assertRaises(ValueError):
                self.nodes.draw(n)
        self.assertEqual(self.nodes.draw(np.int64(3))[0].shape, (3, 2))

    def test_09_invalid_source_arrays(self):
        for field, value in (('C', np.zeros(3)), ('sigma', np.full(41, -1.)),
                             ('r', np.zeros(41)), ('r', np.arange(41)[::-1]+1.),
                             ('C', np.full(41, np.nan))):
            a = copy.deepcopy(self.a); setattr(a, field, value)
            with self.subTest(field=field), self.assertRaises(ValueError):
                build_distribution(a, **self.kw)

    def test_10_invalid_width_and_quadrature_parameters(self):
        for field in ('dL', 'dE', 'w', 'tau_keep'):
            for value in (0., -1., np.nan, np.inf):
                a = copy.deepcopy(self.a); setattr(a, field, value)
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    build_distribution(a, **self.kw)
        for overrides in ({'n_L': 1}, {'n_vr': 1}, {'reach': 0.}, {'reach': np.nan}):
            kw = dict(self.kw, **overrides)
            with self.subTest(overrides=overrides), self.assertRaises(ValueError):
                build_distribution(self.a, **kw)

    def test_11_nonuniform_L_grid_is_explicitly_rejected(self):
        a = copy.deepcopy(self.a)
        a._grids = lambda C, n_L, n_vr, reach: (np.array([-1., 0., 2.]), np.zeros(3),
                                               *np.polynomial.legendre.leggauss(n_vr))
        with self.assertRaisesRegex(ValueError, 'uniform'):
            build_distribution(a, n_L=3, n_vr=64)

    def test_12_no_live_phase_space(self):
        a = copy.deepcopy(self.a)
        a.phi_total = lambda r, C: np.full(len(r), 1e6)
        with self.assertRaisesRegex(ValueError, 'no positive probability'):
            build_distribution(a, **self.kw)

    def test_13_existing_WarmAnnulus_integration(self):
        # Real archived density and grids; this is a current-potential test, not a fixed-point claim.
        a = WarmAnnulus(dL=.12, dE=.03, n_r=96)
        kw = dict(n_L=48, n_vr=32, reach=5.)
        before = fingerprint(a)
        result = a.verify_sampling(**kw)
        self.assertTrue(result['numerical_verification_passed'], result)
        x, v = a.sample(100, seed=4, **kw)
        xx, vv = sample(a, 100, seed=4, **kw)
        np.testing.assert_array_equal(x, xx); np.testing.assert_array_equal(v, vv)
        self.assertEqual(before, fingerprint(a))
        self.assertIs(WarmAnnulus.surface_density, ArchivedWarmAnnulus.surface_density)
        self.assertIs(WarmAnnulus.solve, ArchivedWarmAnnulus.solve)
        self.assertIsNot(WarmAnnulus.sample, ArchivedWarmAnnulus.sample)

    def test_14_original_module_stays_unchanged(self):
        before = hashlib.sha256(ARCHIVED_PATH.read_bytes()).hexdigest()
        WarmAnnulus(n_r=48).sample(32, n_L=24, n_vr=16)
        self.assertEqual(before, hashlib.sha256(ARCHIVED_PATH.read_bytes()).hexdigest())

    def test_15_report_safety_and_status(self):
        from .audit import ROOT, new_output_path, write_new_report
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'report.json'
            write_new_report(path, {'numerical_verification': 'failed', 'physical_stability': 'not_evaluated'})
            before = path.read_bytes()
            with self.assertRaises(FileExistsError):
                write_new_report(path, {'numerical_verification': 'passed'})
            self.assertEqual(path.read_bytes(), before)
            bad = Path(tmp)/'bad.json'
            with self.assertRaises(ValueError):
                write_new_report(bad, {'value': float('nan')})
            self.assertFalse(bad.exists())
        with self.assertRaises(ValueError):
            new_output_path(ROOT/'research_work/results/new-sampler-result.json')


if __name__ == '__main__':
    unittest.main(verbosity=2)
