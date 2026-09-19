"""Offline algebra/provenance checks. These do not certify originality or gravity.

Run from any working directory: python -B research_plan/prior-art/verify.py
Optional --write-report PATH creates a NEW report; existing paths are refused.
No scientific simulation or archive is imported, executed or overwritten.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import unittest

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent


class AlgebraChecks(unittest.TestCase):
    def test_published_interpolation_identity(self):
        x = sp.symbols('x', positive=True)
        s = sp.sqrt(1+4*x)
        self.assertEqual(sp.simplify((s-1)/(s+1)-4*x/(1+s)**2), 0)

    def test_spherical_inverse(self):
        z = sp.symbols('z', positive=True)  # sqrt(y); avoids branch ambiguity
        x = z*z+z
        mu = z/(1+z)
        self.assertEqual(sp.simplify(x*mu-z*z), 0)

    def test_known_functional_derivative(self):
        t = sp.symbols('t', positive=True)
        x, mu, F = t*t+t, t/(1+t), t**4+sp.Rational(2,3)*t**3
        self.assertEqual(sp.simplify(sp.diff(F,t)-2*x*mu*sp.diff(x,t)), 0)

    def test_cascade_elimination(self):
        tf, tk, s = sp.symbols('tf tk s', positive=True)
        eliminated = tf*s*s+(1+tf/tk)*s+1/tk
        self.assertEqual(sp.simplify(eliminated-(1+s*tf)*(1+s*tk)/tk), 0)

    def test_kernel_laplace(self):
        tf, tk, s = sp.symbols('tf tk s', positive=True)
        transform = tk/(tk-tf)*(1/(s+1/tk)-1/(s+1/tf))
        H = tk/((1+s*tk)*(1+s*tf))
        self.assertEqual(sp.simplify(transform-H), 0)
        self.assertEqual(sp.simplify(H.subs(s,0)-tk), 0)

    def test_equal_time_kernel_limit(self):
        tf, tk, u = sp.symbols('tf tk u', positive=True)
        h = tk/(tk-tf)*(sp.exp(-u/tk)-sp.exp(-u/tf))
        self.assertEqual(sp.simplify(sp.limit(h,tf,tk)-(u/tk)*sp.exp(-u/tk)), 0)

    def test_static_is_not_dynamic_equivalence(self):
        tk, tf = 10.0, 3.0
        H = lambda s: tk/((1+s*tk)*(1+s*tf))
        self.assertEqual(H(0), tk)
        self.assertGreater(abs(H(1j)-tk), 1.0)

    def test_stationary_field_elimination(self):
        W = sp.Matrix([[1,sp.Rational(1,3)],[sp.Rational(1,3),2]])
        m = sp.Matrix([2,3])
        alpha, tk = sp.symbols('alpha tk', positive=True)
        h = alpha*tk*W*m
        total = -(m.T*W*h)[0]+(h.T*h)[0]/(2*alpha*tk)
        target = -alpha*tk*(m.T*W*W*m)[0]/2
        self.assertEqual(sp.simplify(total-target), 0)

    def test_reciprocal_dissipation_identity(self):
        alpha, tf, tk, gamma = sp.symbols('alpha tf tk gamma', positive=True)
        h1,h2,v1,v2,b1,b2 = sp.symbols('h1 h2 v1 v2 b1 b2', real=True)
        h,hd,b = sp.Matrix([h1,h2]),sp.Matrix([v1,v2]),sp.Matrix([b1,b2])
        hdd = (alpha*b-gamma*hd-h/tk)/tf
        matter = -(b.T*hd)[0]
        field = tf/alpha*(hd.T*hdd)[0]+(h.T*hd)[0]/(alpha*tk)
        self.assertEqual(sp.simplify(matter+field+gamma/alpha*(hd.T*hd)[0]), 0)

    def test_phase_space_jacobian(self):
        R,L = sp.symbols('R L', positive=True)
        self.assertEqual(sp.simplify(R*sp.diff(L/R,L)), 1)

    def test_gaussian_pair_force(self):
        x,y,w,k,m1,m2 = sp.symbols('x y w k m1 m2', real=True, positive=True)
        K = sp.exp(-(x-y)**2/(2*w*w))
        U = -k*m1*m2*K
        self.assertEqual(sp.simplify(-sp.diff(U,x)/m1-k*m2*sp.diff(K,x)), 0)
        self.assertEqual(sp.simplify(sp.diff(U,x)+sp.diff(U,y)), 0)

    def test_gaussian_ring_identity(self):
        # Two independently evaluated representations; quadrature nodes are not source-ring data.
        from scipy.special import i0e
        from numpy.polynomial.legendre import leggauss
        z,weights = leggauss(256)
        theta = np.pi*(z+1)
        for r,R,w in ((1.,1.,.1),(.7,1.,.2),(1.4,.9,.3)):
            direct = .5*np.sum(weights*np.exp(-(r*r+R*R-2*r*R*np.cos(theta))/(2*w*w)))
            closed = np.exp(-(r-R)**2/(2*w*w))*i0e(r*R/w**2)
            self.assertLess(abs(direct-closed), 2e-13)


class RegisterChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sources = json.loads((HERE/'sources.json').read_text(encoding='utf-8'))
        cls.claims = json.loads((HERE/'claims.json').read_text(encoding='utf-8'))

    def test_unique_sources_and_inspection(self):
        rows = self.sources['sources']
        self.assertEqual(len(rows), len({s['id'] for s in rows}))
        for s in rows:
            self.assertTrue(s['url'].startswith('https://'))
            self.assertTrue(s['inspection'])
            self.assertTrue(s['limit_of_comparison'])

    def test_claim_sources_resolve(self):
        ids = {s['id'] for s in self.sources['sources']}
        for c in self.claims['claims']:
            self.assertLessEqual(set(c['sources']), ids)
            self.assertTrue(c['permitted_wording'])
            self.assertTrue(c['prohibited_overclaim'])

    def test_exact_claim_has_equation_source(self):
        src = {s['id']:s for s in self.sources['sources']}
        for c in self.claims['claims']:
            if c['status'] in ('EXACT_FUNCTION','EXACT_ALGEBRAIC_LIMIT'):
                self.assertTrue(any('equations' in src[i]['inspection'] for i in c['sources']))

    def test_preprint_and_unretrieved_leads_separated(self):
        src = {s['id']:s for s in self.sources['sources']}
        self.assertIn('preprint',src['LV2026']['publication_status'])
        leads = {s['id'] for s in self.sources['followup_leads']}
        self.assertTrue(leads.isdisjoint(src))
        for c in self.claims['claims']:
            self.assertTrue(leads.isdisjoint(c['sources']))

    def test_markdown_reference_ids(self):
        text = (HERE/'audit.md').read_text(encoding='utf-8')
        ids = {s['id'] for s in self.sources['sources']}
        definitions = set(re.findall(r'^\[([A-Za-z][A-Za-z0-9]*)\]: https://',text,re.M))
        self.assertEqual(ids, definitions)
        self.assertNotIn('\ue200',text)  # no tool-only citation syntax in GitHub documents

    def test_package_manifest(self):
        manifest = json.loads((HERE/'file-manifest.json').read_text(encoding='utf-8'))
        for rel,expected in manifest['sha256'].items():
            path = (HERE/rel).resolve()
            self.assertTrue(path.is_relative_to(HERE.resolve()))
            self.assertTrue(path.is_file(),rel)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),expected,rel)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',type=Path)
    args = parser.parse_args()
    if args.write_report and args.write_report.exists():
        parser.error('Report already exists; refusing to overwrite it.')
    suite = unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__])
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    report = dict(audit_date='2026-09-19',baseline='3a80fec8ae7d45cc6906c20cb9c537c3ac2c9698',
        checks_run=result.testsRun,passed=result.testsRun-len(result.failures)-len(result.errors)-len(result.skipped),
        failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),
        scope='Algebraic and register-integrity checks only; no global novelty certificate and no gravity simulation verification.',
        literature_records=len(RegisterChecks.sources['sources']),unretrieved_leads=len(RegisterChecks.sources['followup_leads']),
        historical_science_suite='not_run',science_files_changed=False,
        environment=dict(python=sys.version.split()[0],numpy=np.__version__,sympy=sp.__version__))
    if args.write_report:
        args.write_report.parent.mkdir(parents=True,exist_ok=True)
        args.write_report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())
