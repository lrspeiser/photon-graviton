"""TF-1 evidence bookkeeping, after CWC-1's common.py (existing project code, reused)."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CWC = ROOT/'research_work/experiments/coupled_conversion'
PM = ROOT/'research_work/results/path-memory'


def plain(x):
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (tuple, list, np.ndarray)):
        return [plain(v) for v in x]
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        if not np.isfinite(x):
            raise ValueError('Nonfinite evidence')
        return float(x)
    return x


def save(path, data):
    Path(path).write_text(json.dumps(plain(data), indent=1, allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hashes(paths):
    return {Path(p).resolve().relative_to(ROOT).as_posix(): sha(p) for p in paths}


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))
