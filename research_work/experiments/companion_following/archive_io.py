"""Lossless byte chunking for large immutable NumPy archives; no data reduction."""
from __future__ import annotations
import hashlib
import io
import json
from pathlib import Path
import numpy as np

def open_npz(path):
    path = Path(path)
    if path.exists():
        return np.load(path, allow_pickle=False)
    index = json.loads(path.with_name(path.name+".parts.json").read_text())
    payload=[]
    for part in index["parts"]:
        raw=path.with_name(part["name"]).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=part["sha256"]:
            raise ValueError("Archive part checksum mismatch: "+part["name"])
        payload.append(raw)
    raw=b"".join(payload)
    if len(raw)!=index["bytes"] or hashlib.sha256(raw).hexdigest()!=index["sha256"]:
        raise ValueError("Reassembled archive checksum mismatch")
    return np.load(io.BytesIO(raw), allow_pickle=False)

def package(root):
    for path in sorted(Path(root).rglob("*.npz")):
        if path.stat().st_size <= 48*1024**2: continue
        indexpath=path.with_name(path.name+".parts.json")
        if indexpath.exists(): raise FileExistsError(indexpath)
        index=dict(original=path.name,bytes=path.stat().st_size,parts=[])
        digest=hashlib.sha256()
        with path.open("rb") as source:
            i=0
            while raw:=source.read(32*1024**2):
                digest.update(raw)
                part=path.with_name(path.name+f".part{i:03d}")
                with part.open("xb") as dest: dest.write(raw)
                index["parts"].append(dict(name=part.name,bytes=len(raw),
                     sha256=hashlib.sha256(raw).hexdigest()))
                i+=1
        index["sha256"]=digest.hexdigest()
        with indexpath.open("x",encoding="utf8") as target:
            json.dump(index,target,indent=2);target.write("\n")
        print(f"Losslessly packaged {path}: {len(index['parts'])} parts",flush=True)

if __name__=="__main__":
    package(Path(__file__).resolve().parent/"evidence"/"screen-v1")
