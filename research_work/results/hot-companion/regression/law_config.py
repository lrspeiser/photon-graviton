"""The law under test: its constants and any candidate amendments, in one dictionary.

Every test in the suite reads the law from this dictionary and nothing else, so a change to the
law is a change to one JSON file in candidates/ (or to the default below).

The adopted law (round 11) is the round-3 law with one amendment, the gradual release adopted for
Cassini in round 9: the companion leaving an emitter is released over L = 30,000 AU (0.15 pc, about
720 years of travel at u), R(r) = 1 - exp(-r/L) (ADOPTED below). Its three constants were refitted in
round 11 in the project's own (static) distances, with X-COP's stellar profiles deprojected
(run-xcop-static-v11): a = 6.547e-11 m/s^2, g_d = 2.107e-10 m/s^2, u = 162.6 km/s.
'round9' loads the same law with the round-3 constants (fitted in the release's LCDM units, with
X-COP's projected stars); 'round3' loads the round-3 law, released at once.

The law (round 3; README section 9):
    h = g_N + exp(-|g_N|/g_d) sqrt(a (|g_N| + S)) (g_N + g_hot)/(|g_N| + |g_hot|),   lap Phi = -div h
    S = G int k rho_free / d^2,  g_hot = G int k rho_free (x'-x)/|x'-x|^3,  k = 3 sigma^2/u^2 (stars), 0 (gas)
    memory: the companion keeps its emitter's velocity; a fresh sphere of radius u t forms around stopped gas
    reach: u t_age (13 Gyr: 2.16 Mpc with the round-11 constants)

Keys of the law dictionary (code units: kpc, km/s, Msun; accelerations in (km/s)^2/kpc):
    a_code, a_SI        deep-regime strength a
    lam                 g_d / a  (the release scale in units of a)
    g_d_SI              release scale g_d = lam a
    u_kms               the companion's speed u
    reach_kpc           u x 13 Gyr, where the companion's pull stops
    release_length_au   the companion is released gradually over a length L around each emitter,
                        R(r) = 1 - exp(-r/L); proposed in round 7, adopted in round 9 at 30,000 AU
                        (0 = at once, the round-3 law)
    external_hold       candidate amendment (round 7): how strongly a subsystem's companion follows an
                        outside galaxy's pull (dwarfs, the Sun, wide binaries); 1 = fully (the law as is)
    gd_scale            candidate amendment (round 7): g_d multiplied by this; already folded into lam and g_d_SI
    refit               constants refitted on their home data after the amendments were applied
Candidate files (JSON) may set: name, description, base ('round3', 'round11' or a path to a results.json
with a 'constants' block), gd_scale, release_length_au, external_hold, refit (list: 'a' on SPARC with g_d
held, 'u' on X-COP), and explicit overrides a_SI / lam / u_kms. Amendments a candidate does not set take
their adopted values (ADOPTED). The candidates in candidates/ were written in rounds 7-9 and keep the
round-3 constants (base 'round3', the default), so they reproduce the tests they were made for; a new
candidate that changes the adopted law sets base 'round11'.
"""
from __future__ import annotations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
K_SI = 1e6 / 3.0856775814913673e19          # (km/s)^2/kpc -> m/s^2
KPC_PER_KMS_GYR = 1.0227121650537077         # 1 km/s for 1 Gyr, in kpc
AGE_GYR = 13.0
ADOPTED = dict(release_length_au=30000.0)    # round 9: gradual release over 0.15 pc (Cassini)


def base_constants(base='round3'):
    if base == 'round11':          # refitted in the project's own distances (run-xcop-static-v11)
        c = json.loads((RESULTS / 'run-xcop-static-v11/xcop_static_v11.json').read_text())['joint_refit_static'][-1]
        return dict(a_code=c['a_SI'] / K_SI, lam=c['lam'], u_kms=c['u_kms'])
    path = RESULTS / 'run-v3/results.json' if base == 'round3' else Path(base)
    c = json.loads(path.read_text())['constants']
    return dict(a_code=c['a_code'], lam=c['lam'], u_kms=c['u_kms'])


def make_law(a_code, lam, u_kms, name='custom', description='', release_length_au=0.0, external_hold=1.0,
             gd_scale=1.0, refit=(), base='round3'):
    return dict(name=name, description=description, base=base,
                a_code=float(a_code), a_SI=float(a_code * K_SI), lam=float(lam), g_d_SI=float(lam * a_code * K_SI),
                u_kms=float(u_kms), reach_kpc=float(u_kms * KPC_PER_KMS_GYR * AGE_GYR),
                release_length_au=float(release_length_au), external_hold=float(external_hold),
                gd_scale=float(gd_scale), refit=list(refit))


def load_law(spec=None):
    """spec: None or 'round11' (the adopted law), 'round9' (the round-3 constants with round 9's
    amendment), 'round3' (before round 9's amendment), a candidate name in candidates/, or a path to a
    JSON file."""
    if spec == 'round9':
        c = base_constants()
        return make_law(c['a_code'], c['lam'], c['u_kms'], name='round9',
                        description='the round-3 constants (run-v3) with gradual release over 0.15 pc, adopted in round 9',
                        release_length_au=ADOPTED['release_length_au'])
    if spec in (None, '', 'round11'):
        c = base_constants('round11')
        return make_law(c['a_code'], c['lam'], c['u_kms'], name='round11', base='round11',
                        description='the round-9 law with its three constants refitted in the project\'s own (static) distances: '
                                    'a and g_d on SPARC, u on X-COP converted to the static law (round 11)',
                        release_length_au=ADOPTED['release_length_au'])
    if spec == 'round3':
        c = base_constants()
        return make_law(c['a_code'], c['lam'], c['u_kms'], name='round3',
                        description='the round-3 law with its published constants (run-v3), released at once')
    p = Path(spec)
    if not p.exists():
        p = HERE / 'candidates' / (spec if spec.endswith('.json') else spec + '.json')
    cfg = json.loads(p.read_text())
    c = base_constants(cfg.get('base', 'round3'))
    a = cfg['a_SI'] / K_SI if 'a_SI' in cfg else c['a_code']
    lam = cfg.get('lam', c['lam'])
    u = cfg.get('u_kms', c['u_kms'])
    gs = float(cfg.get('gd_scale', 1.0))
    return make_law(a, lam * gs, u, name=cfg.get('name', p.stem), description=cfg.get('description', ''),
                    release_length_au=cfg.get('release_length_au', ADOPTED['release_length_au']), external_hold=cfg.get('external_hold', 1.0),
                    gd_scale=gs, refit=cfg.get('refit', []), base=cfg.get('base', 'round3'))


def with_constants(law, a_code=None, lam=None, u_kms=None):
    """A copy of the law with new constants (amendments kept)."""
    a = law['a_code'] if a_code is None else a_code
    lm = law['lam'] if lam is None else lam
    u = law['u_kms'] if u_kms is None else u_kms
    out = make_law(a, lm, u, name=law['name'], description=law['description'], release_length_au=law['release_length_au'],
                   external_hold=law['external_hold'], gd_scale=law['gd_scale'], refit=law['refit'], base=law['base'])
    for k in ('refit_log',):
        if k in law: out[k] = law[k]
    return out


def describe(law):
    parts = [f"a = {law['a_SI']:.4e} m/s^2", f"g_d = {law['g_d_SI']:.4e} m/s^2 (lambda {law['lam']:.3f})",
             f"u = {law['u_kms']:.1f} km/s", f"reach {law['reach_kpc'] / 1000:.2f} Mpc"]
    if law['release_length_au'] > 0: parts.append(f"gradual release over {law['release_length_au']:.0e} AU")
    if law['external_hold'] != 1.0: parts.append(f"external hold x{law['external_hold']:g}")
    if law['gd_scale'] != 1.0: parts.append(f"g_d x{law['gd_scale']:g}")
    if law['refit']: parts.append('refit: ' + ', '.join(law['refit']))
    return '; '.join(parts)
