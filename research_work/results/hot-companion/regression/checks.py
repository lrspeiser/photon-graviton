"""One check = one comparison of the law with one measurement, scored the same way everywhere.

status  pass   within 2 standard errors of the measurement (or inside the required range)
        close  within 3
        fail   further out
        info   a number we track but do not grade (disputed data, or a quantity with no measurement)
        error  the test itself crashed
score   lower is better: the distance from the measurement in standard errors (z), or for a limit
        the value over the limit. Scores let the runner see a test drift before its status changes.
"""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
import math

RANK = {'pass': 0, 'close': 1, 'fail': 2}
PASS_Z, CLOSE_Z = 2.0, 3.0


@dataclass
class Check:
    id: str
    group: str
    title: str
    value: float | None
    unit: str = ''
    target: str = ''
    status: str = 'info'
    score: float | None = None
    note: str = ''
    refs: str = ''
    detail: dict = field(default_factory=dict)

    def as_dict(self):
        d = asdict(self)
        for k in ('value', 'score'):
            if d[k] is not None:
                d[k] = float(d[k])
        return d


def _grade(z):
    z = abs(z)
    return ('pass' if z <= PASS_Z else 'close' if z <= CLOSE_Z else 'fail'), z


def z_check(value, obs, err):
    """A measurement obs +- err."""
    return _grade((value - obs) / err)


def range_check(value, lo, hi, err):
    """Inside [lo, hi] passes; outside, graded by the distance to the nearest edge in units of err."""
    if lo <= value <= hi:
        return 'pass', 0.0
    return _grade((lo - value) / err if value < lo else (value - hi) / err)


def floor_check(value, floor, err):
    """A lower bound (e.g. Clowe et al.'s kappas): at or above the floor passes; below, graded in err."""
    return ('pass', 0.0) if value >= floor else _grade((floor - value) / err)


def at_most(value, limit, close_limit=None):
    """An upper limit (e.g. 'the extra pull must be below 1e-10 of Newton')."""
    s = abs(value) / limit if limit else math.inf
    if abs(value) <= limit:
        return 'pass', s
    if close_limit is not None and abs(value) <= close_limit:
        return 'close', s
    return 'fail', s


def rms_z(values, obs, errs):
    """Several bins at once: root-mean-square of the z scores, graded like a single z."""
    zs = [(v - o) / e for v, o, e in zip(values, obs, errs)]
    return _grade(math.sqrt(sum(z * z for z in zs) / len(zs)))


def make(id, group, title, value, crit=None, **kw):
    """Build a Check; crit is the (status, score) pair from one of the helpers above (None = info)."""
    st, sc = crit if crit is not None else ('info', None)
    return Check(id=id, group=group, title=title, value=None if value is None else float(value), status=st, score=sc, **kw)


def compare(new, old, tol_abs=0.1, tol_rel=0.05):
    """How a check moved relative to the baseline: regressed / improved / worse / better / same /
    new / known (a fail already in the baseline) / changed (the number moved by > 0.1% while its grade and
    score did not, e.g. a tracked number, or a value moving inside a passing range) / error."""
    if new['status'] == 'error':
        return 'error'
    if old is None:
        return 'new'
    rn, ro = RANK.get(new['status']), RANK.get(old['status'])
    if rn is not None and ro is not None:
        if rn > ro: return 'regressed'
        if rn < ro: return 'improved'
    s1, s0 = new.get('score'), old.get('score')
    if s1 is not None and s0 is not None:
        d = s1 - s0
        if abs(d) > max(tol_abs, tol_rel * abs(s0)):
            return 'worse' if d > 0 else 'better'
    if new['status'] == 'fail' and old['status'] == 'fail':
        return 'known'
    v1, v0 = new.get('value'), old.get('value')
    if v1 is not None and v0 is not None and abs(v1 - v0) > 1e-3 * max(abs(v0), 1e-300):
        return 'changed'
    return 'same'
