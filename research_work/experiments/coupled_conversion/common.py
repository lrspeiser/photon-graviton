"""Evidence IO for the declared coupled-conversion campaign."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def plain(x):
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list,np.ndarray)):return [plain(v) for v in x]
    if isinstance(x,(np.integer,int)):return int(x)
    if isinstance(x,(np.floating,float)):
        if not np.isfinite(x):raise ValueError('Nonfinite evidence')
        return float(x)
    return x

def save(path,data):
    Path(path).write_text(json.dumps(plain(data),indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def hashes(paths):
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))
