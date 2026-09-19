"""Shared static-space experiment utilities; no geometry or hidden-matter model."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PM = ROOT / 'research_work/results/path-memory'
sys.path.insert(0, str(PM))

INPUTS = set()

def read_json(path):
    path = Path(path)
    INPUTS.add(path)
    return json.loads(path.read_text(encoding='utf-8'))

def track(path):
    path = Path(path)
    INPUTS.add(path)
    return path

def plain(x):
    if isinstance(x, (bool, np.bool_)): return bool(x)
    if isinstance(x, dict): return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)): return [plain(v) for v in x]
    if isinstance(x, (float, np.floating)):
        if not np.isfinite(x): raise ValueError('Nonfinite result cannot be silently saved')
        return float(x)
    if isinstance(x, (int, np.integer)): return int(x)
    return x

def save(path, data):
    Path(path).write_text(json.dumps(plain(data), indent=2, allow_nan=False)+'\n',
                          encoding='utf-8', newline='\n')

def provenance():
    paths = INPUTS | set(HERE.glob('*.py')) | set(HERE.glob('*.md'))
    paths |= {PM / 'steady_field.py', PM / 'cluster_lensing.py'}
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}

def inv_cdf(r, density, u):
    c = cumulative_trapezoid(np.maximum(density, 0), r, initial=0)
    if c[-1] <= 0: raise ValueError('Empty emission profile')
    return np.interp(u, c/c[-1], r)

def directions(u, v):
    z = 2*u-1
    s = np.sqrt(np.maximum(0, 1-z*z))
    return np.column_stack((s*np.cos(2*np.pi*v), s*np.sin(2*np.pi*v), z))

def fit_shape(model, y, err):
    norm = max(float(np.sum(model*y/err**2)/np.sum((model/err)**2)), 0)
    pred = norm*model
    return dict(chi2=np.sum(((pred-y)/err)**2), nuisance=norm, prediction=pred)
