"""The law under test: its constants and any candidate amendments, in one dictionary.

Every test in the suite reads the law from this dictionary and nothing else, so a change to the
law is a change to one JSON file in candidates/ (or to the default below).

The adopted law (round 12) is the round-3 law with one amendment, the gradual release adopted for
Cassini in round 9: the companion leaving an emitter is released over L = 30,000 AU (0.15 pc, about
840 years of travel at u), R(r) = 1 - exp(-r/L) (ADOPTED below). Every data set is in the project's own
(static) distances, 1 + z = e^(alpha D), with alpha = 2.3645e-4 per Mpc (x0.95 of rounds 10-11; H0-like
70.9), fitted jointly to Pantheon+ supernovae and SPARC's Hubble-flow galaxies, which are now placed at
those distances too (run-distance-scale-v12). The three constants were refitted there, with X-COP's
stellar profiles deprojected: a = 6.298e-11 m/s^2, g_d = 2.027e-10 m/s^2, u = 169.4 km/s.
'round11' loads the round-11 constants (alpha x1, SPARC at its published distances: a = 6.547e-11,
g_d = 2.107e-10, u = 162.6); 'round9' the round-3 constants (fitted in the release's LCDM units, with
X-COP's projected stars); 'round3' the round-3 law, released at once.

The law (round 3; README section 9):
    h = g_N + exp(-|g_N|/g_d) sqrt(a (|g_N| + S)) (g_N + g_hot)/(|g_N| + |g_hot|),   lap Phi = -div h
    S = G int k rho_free / d^2,  g_hot = G int k rho_free (x'-x)/|x'-x|^3,  k = 3 sigma^2/u^2 (stars), 0 (gas)
    memory: the companion keeps its emitter's velocity; a fresh sphere of radius u t forms around stopped gas
    reach: u t_age (13 Gyr: 2.25 Mpc with the round-12 constants)

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
    alpha_per_Mpc       the static distance law's scale, read by every conversion to the project's distances
                        (collisions_v10.ALPHA, set by common.apply_distances); 2.489e-4 before round 12
    sparc_distances     'published' (Lelli et al. 2016) or 'static' (the Hubble-flow galaxies at ln(1 + z)/alpha)
    heat_exponent       candidate amendment (round 14): the heat weight k = 3 (sigma/u)^p; the law has p = 2. Set
                        in code/law.py (HEAT_P) by common.apply_distances before any test or refit
    hot_geometry        candidate amendment (round 19): how the hot matter's extra glow is heard. 'two_way' (the law:
                        every shell, inside and outside the receiver), 'one_way' (inner shells only: a medium whose
                        wave travels outward) or 'one_way_vector' (the inner shells' net flux, weight 1). Set in
                        code/law.py (HOT_GEOMETRY); spherical sums only (the collision maps' 3D sums stay two-way).
                        'stream' (with stream_kappa_per_Mpc): a stream absorbing inward-travelling waves, which also
                        filters the clusters' cold glow (run_v3.cluster_M3)
    distance_variant    registered comparison (round 21): 'fixed' (the adopted law, D_A = D) or 'metric' (D_A = D/(1 + z),
                        so that D_L = (1 + z)^2 D_A); set in code/collisions_v10.py (VARIANT). The SLACS test keeps its own
                        conversion (lenses_t35's project distances) and does not see it
    eta_path            registered comparison (round 21): the path factor sqrt(1 + eta z/(1 + z)) on D (collisions_v10.ETA_PATH)
    collision_field     registered comparison (round 24): the field on the collision maps. 'law' (the adopted law's
                        bullet_v4.kappa_map_v4) or a mode of code/eft_field_v24.py (MODES: 'eft', the Casimir-EFT note's
                        local equation; 'eft_memory', 'eft_scalar', ...), installed by common.apply_distances
Candidate files (JSON) may set: name, description, base ('round3', 'round11', 'round12' or a path to a results.json
with a 'constants' block), gd_scale, release_length_au, external_hold, refit (list: 'a' on SPARC with g_d
held, 'u' on X-COP), and explicit overrides a_SI / lam / u_kms. Amendments a candidate does not set take
their adopted values (ADOPTED). The candidates in candidates/ were written in rounds 7-9 and keep the
round-3 constants (base 'round3', the default), so they reproduce the tests they were made for; a new
candidate that changes the adopted law sets base 'round12' (which also brings its distance scale and SPARC's
static distances).
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
ALPHA_ROUND10 = 2.488993286382367e-4         # per Mpc: the static distance law's scale, 1 + z = e^(alpha D), rounds 10-11


def base_constants(base='round3'):
    if base == 'round11':          # refitted in the project's own distances (run-xcop-static-v11)
        c = json.loads((RESULTS / 'run-xcop-static-v11/xcop_static_v11.json').read_text())['joint_refit_static'][-1]
        return dict(a_code=c['a_SI'] / K_SI, lam=c['lam'], u_kms=c['u_kms'])
    if base == 'round12':          # the distance scale fitted jointly, SPARC in the static law too (run-distance-scale-v12)
        d = json.loads((RESULTS / 'run-distance-scale-v12/adopted_v12.json').read_text())
        c = d['joint_refit'][-1]
        return dict(a_code=c['a_SI'] / K_SI, lam=c['lam'], u_kms=c['u_kms'], alpha_per_Mpc=d['alpha_per_Mpc'])
    path = RESULTS / 'run-v3/results.json' if base == 'round3' else Path(base)
    c = json.loads(path.read_text())['constants']
    return dict(a_code=c['a_code'], lam=c['lam'], u_kms=c['u_kms'])


def make_law(a_code, lam, u_kms, name='custom', description='', release_length_au=0.0, external_hold=1.0,
             gd_scale=1.0, refit=(), base='round3', alpha_per_Mpc=ALPHA_ROUND10, sparc_distances='published',
             heat_exponent=2.0, hot_geometry='two_way', stream_kappa_per_Mpc=0.0, distance_variant='fixed', eta_path=0.0,
             collision_field='law'):
    """alpha_per_Mpc: the static distance law's scale, used by every conversion to the project's distances
    (collisions_v10.ALPHA). sparc_distances: 'published' (Lelli et al. 2016) or 'static' (the Hubble-flow
    galaxies at D = ln(1 + z)/alpha, round 12)."""
    return dict(name=name, description=description, base=base,
                a_code=float(a_code), a_SI=float(a_code * K_SI), lam=float(lam), g_d_SI=float(lam * a_code * K_SI),
                u_kms=float(u_kms), reach_kpc=float(u_kms * KPC_PER_KMS_GYR * AGE_GYR),
                release_length_au=float(release_length_au), external_hold=float(external_hold),
                gd_scale=float(gd_scale), refit=list(refit), alpha_per_Mpc=float(alpha_per_Mpc), sparc_distances=sparc_distances,
                heat_exponent=float(heat_exponent), hot_geometry=hot_geometry,
                stream_kappa_per_Mpc=float(stream_kappa_per_Mpc), distance_variant=distance_variant, eta_path=float(eta_path),
                collision_field=collision_field)


def load_law(spec=None):
    """spec: None or 'round12' (the adopted law), 'round11' (its constants before the distance scale was fitted
    jointly), 'round9' (the round-3 constants with round 9's amendment), 'round3' (before round 9's amendment), a
    candidate name in candidates/, or a path to a JSON file."""
    if spec == 'round9':
        c = base_constants()
        return make_law(c['a_code'], c['lam'], c['u_kms'], name='round9',
                        description='the round-3 constants (run-v3) with gradual release over 0.15 pc, adopted in round 9',
                        release_length_au=ADOPTED['release_length_au'])
    if spec == 'round11':
        c = base_constants('round11')
        return make_law(c['a_code'], c['lam'], c['u_kms'], name='round11', base='round11',
                        description='the round-9 law with its three constants refitted in the project\'s own (static) distances: '
                                    'a and g_d on SPARC, u on X-COP converted to the static law (round 11)',
                        release_length_au=ADOPTED['release_length_au'])
    if spec in (None, '', 'round12'):
        c = base_constants('round12')
        return make_law(c['a_code'], c['lam'], c['u_kms'], name='round12', base='round12',
                        description='the round-11 law with the distance scale fitted jointly (supernovae and SPARC\'s Hubble-flow '
                                    'galaxies; alpha x0.95, H0-like 70.9), SPARC in the static law too, and a, g_d and u refitted there (round 12)',
                        release_length_au=ADOPTED['release_length_au'], alpha_per_Mpc=c['alpha_per_Mpc'], sparc_distances='static')
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
                    gd_scale=gs, refit=cfg.get('refit', []), base=cfg.get('base', 'round3'),
                    alpha_per_Mpc=cfg.get('alpha_per_Mpc', c.get('alpha_per_Mpc', ALPHA_ROUND10)),
                    sparc_distances=cfg.get('sparc_distances', 'static' if cfg.get('base') == 'round12' else 'published'),
                    heat_exponent=cfg.get('heat_exponent', 2.0), hot_geometry=cfg.get('hot_geometry', 'two_way'),
                    stream_kappa_per_Mpc=cfg.get('stream_kappa_per_Mpc', 0.0),
                    distance_variant=cfg.get('distance_variant', 'fixed'), eta_path=cfg.get('eta_path', 0.0),
                    collision_field=cfg.get('collision_field', 'law'))


def with_constants(law, a_code=None, lam=None, u_kms=None):
    """A copy of the law with new constants (amendments kept)."""
    a = law['a_code'] if a_code is None else a_code
    lm = law['lam'] if lam is None else lam
    u = law['u_kms'] if u_kms is None else u_kms
    out = make_law(a, lm, u, name=law['name'], description=law['description'], release_length_au=law['release_length_au'],
                   external_hold=law['external_hold'], gd_scale=law['gd_scale'], refit=law['refit'], base=law['base'],
                   alpha_per_Mpc=law.get('alpha_per_Mpc', ALPHA_ROUND10), sparc_distances=law.get('sparc_distances', 'published'),
                   heat_exponent=law.get('heat_exponent', 2.0), hot_geometry=law.get('hot_geometry', 'two_way'),
                   stream_kappa_per_Mpc=law.get('stream_kappa_per_Mpc', 0.0),
                   distance_variant=law.get('distance_variant', 'fixed'), eta_path=law.get('eta_path', 0.0),
                   collision_field=law.get('collision_field', 'law'))
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
    if law.get('alpha_per_Mpc', ALPHA_ROUND10) != ALPHA_ROUND10:
        parts.append(f"distance scale alpha x{law['alpha_per_Mpc'] / ALPHA_ROUND10:.4f} (H0-like {law['alpha_per_Mpc'] * 299792.458:.1f})")
    if law.get('sparc_distances', 'published') != 'published': parts.append(f"SPARC distances {law['sparc_distances']}")
    if law.get('heat_exponent', 2.0) != 2.0: parts.append(f"heat weight k = 3 (sigma/u)^{law['heat_exponent']:g}")
    if law.get('hot_geometry', 'two_way') == 'stream':
        parts.append(f"heard through a stream absorbing inward waves at {law['stream_kappa_per_Mpc']:g}/Mpc")
    elif law.get('hot_geometry', 'two_way') != 'two_way': parts.append(f"hot matter heard {law['hot_geometry'].replace('_', ' ')}")
    if law.get('distance_variant', 'fixed') != 'fixed': parts.append(f"distance geometry {law['distance_variant']} (D_A = D/(1 + z))")
    if law.get('eta_path', 0.0): parts.append(f"path factor sqrt(1 + {law['eta_path']:g} z/(1 + z))")
    if law.get('collision_field', 'law') != 'law': parts.append(f"collision maps: {law['collision_field']} (code/eft_field_v24.py)")
    return '; '.join(parts)
