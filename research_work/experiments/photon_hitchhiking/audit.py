"""Read archived evidence and pinned Git objects; never rerun simulations."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding="utf8"))
def main():
    checks=[];hash_count=0;source_count=0;trajectory_arrays=0
    def check(name,ok,detail=None):checks.append(dict(name=name,passed=bool(ok),detail=detail))
    for directory,names in [
      (HERE.parent/"companion_following",["evidence-sha256.json"]),
      (HERE.parent/"matter_streams",["evidence-sha256.json"]),
      (HERE,["evidence-sha256.json","curvature-evidence-sha256.json","numerical-evidence-sha256.json"])]:
        for filename in names:
            manifest=read(directory/filename)
            bad=[]
            for rel,expected in manifest.items():
                p=directory/rel
                actual=hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
                if actual!=expected:bad.append(rel)
                hash_count+=1
            check(str((directory/filename).relative_to(REPO))+" exact bytes",not bad,bad)
    for directory in [HERE.parent/"matter_streams",HERE]:
        for manifest in directory.glob("evidence*/manifest.json"):
            m=read(manifest);entries=m.get("sources",m.get("source_sha256",{}))
            bad=[]
            for rel,expected in entries.items():
                path=rel if rel.startswith("research_work/") else (directory/rel).relative_to(REPO).as_posix()
                try:
                    data=subprocess.check_output(["git","show",m["git_head"]+":"+path],cwd=REPO,stderr=subprocess.PIPE)
                    actual=hashlib.sha256(data).hexdigest()
                except subprocess.CalledProcessError:actual=None
                if actual!=expected:bad.append(path)
                source_count+=1
            check(str(manifest.relative_to(REPO))+" source commit",not bad,bad)
        for archive in directory.glob("evidence*/*/trajectories.npz"):
            with np.load(archive) as z:
                state=z["state"]
                check(str(archive.relative_to(REPO))+" finite state",np.isfinite(state).all())
                if directory==HERE and state.shape[-1]==9:
                    # Archived canonical ray state, not JSON summary.
                    root=archive.parts[-3]
                    curvature="curvature" in root or "curvature" in archive.parts[-2]
                    bound=2 if curvature else 1
                    p=state[:,:,3];moment=state[:,:,4:6]
                    err=float(np.max(np.linalg.norm(moment,axis=-1)-bound*p))
                    check(str(archive.relative_to(REPO))+" occupation/memory",p.min()>=-1e-8 and p.max()<=1+1e-8 and err<1e-8,
                          dict(min_p=float(p.min()),max_p=float(p.max()),moment_excess=err))
                trajectory_arrays+=1
    gf=read(HERE.parent/"companion_following"/"audit.json")
    check("GF1 completed audit",gf["all_evidence_checks_passed"])
    num=read(HERE/"evidence-numerical-v1"/"summary.json")
    check("HH follow-up controls",num["controls_passed"])
    for root in ["evidence-v1","evidence-curvature-v1","evidence-numerical-v1"]:
        check(root+" completion",read(HERE/root/"complete.json")["complete"])
    result=dict(all_passed=all(c["passed"] for c in checks),hashes_checked=hash_count,
        pinned_sources_checked=source_count,trajectory_arrays_checked=trajectory_arrays,checks=checks,
        note="Integrity and implementation checks do not override scientific or timestep failures.")
    (HERE/"integrity-audit.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in result.items() if k!="checks"}))
    for c in checks:
        if not c["passed"]:print(c)
    if not result["all_passed"]:raise SystemExit(1)
if __name__=="__main__":main()
