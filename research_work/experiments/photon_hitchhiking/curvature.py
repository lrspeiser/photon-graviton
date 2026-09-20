"""HH-2 curvature-carrying attachment, compared with preserved HH-1."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
import hitch as h
from campaign import dump
HERE=Path(__file__).resolve().parent

def curvature(x,pitch):
    r=np.linalg.norm(x,axis=1)
    phi=np.arctan2(x[:,1],x[:,0]);beta=np.arctan(pitch)
    psi=phi+beta
    amplitude=np.divide(np.sin(beta),r,out=np.zeros_like(r),where=r!=0)
    return np.stack([-np.sin(psi)*amplitude,np.cos(psi)*amplitude],axis=-1)

def rhs(t,y,p):
    out=h.mrhs(t,y,p)
    psi,rate=h.local(t,y,p)
    capture=rate*(1-y[:,3])
    C=curvature(y[:,:2],p["pitch"])
    out[:,4:6]=capture[:,None]*C-p["off"][:,None]*y[:,4:6]
    return out

def controls():
    x=np.array([[1.,2.],[-.8,1.4],[3.,-2.],[.7,-.4]])
    pitch=np.array([.5,1.,2.,4.])
    def tangent(x):
        psi=np.arctan2(x[:,1],x[:,0])+np.arctan(pitch)
        return np.stack([np.cos(psi),np.sin(psi)],axis=-1)
    t=tangent(x);eps=1e-6
    derivative=(tangent(x+eps*t)-tangent(x-eps*t))/(2*eps)
    error=float(np.max(np.abs(derivative-curvature(x,pitch))))
    if error>1e-8:raise AssertionError(error)
    return dict(curvature_finite_difference_error=error,passed=True)

def main():
    out=HERE/"evidence-curvature-v1";out.mkdir(exist_ok=False)
    source=[HERE/"curvature.py",HERE/"hitch.py",HERE/"protocol-curvature.md",h.STREAM/"campaign.py"]
    dump(out/"manifest.json",dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=h.REPO,text=True).strip(),
        sources={p.relative_to(h.REPO).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in source}))
    dump(out/"controls.json",controls())
    primary=h.mean_trace(h.configs(),out/"mean-screen",.02,fun=rhs,moment_bound=2.)
    fine=h.mean_trace(h.REFINE,out/"mean-refinement",.01,fun=rhs,moment_bound=2.)
    m0=h.stochastic(h.REPS,out/"events-dt001",.01,190621,carried="curvature")
    m1=h.stochastic(h.REPS,out/"events-dt0005",.005,190622,carried="curvature")
    coarse={(r["id"],r["impact"]):r for r in primary["results"]}
    refs=[]
    for r in fine["results"]:
        q=coarse[(r["id"],r["impact"])]
        both=r["arrived"] and q["arrived"]
        da=abs(r["angle"]-q["angle"]) if both else None
        dy=abs(r["observer_y"]-q["observer_y"]) if both else None
        refs.append(dict(id=r["id"],impact=r["impact"],angle_difference=da,position_difference=dy,
            passed=r["arrived"]==q["arrived"] and (da is None or da<.01) and (dy is None or dy<.01)))
    ranks=[]
    for r in h.configs():
        pairs=[a for a in primary["pairs"] if a["id"]==r["id"] and a["impact"] in (1,2,4)]
        if all(p["arrived"] for p in pairs):
            ranks.append(dict(id=r["id"],inward=float(np.mean([p["inward"] for p in pairs])),
                             sideways=float(np.mean([abs(p["sideways"]) for p in pairs]))))
    ranks.sort(key=lambda r:(-r["inward"],r["id"]))
    stats=[]
    for a,b in zip(m0["ensembles"],m1["ensembles"]):
        se=np.hypot(a["se_angle"],b["se_angle"])
        diff=abs(a["mean_angle"]-b["mean_angle"])
        stats.append(dict(id=b["id"],impact=b["impact"],difference=diff,combined_se=float(se),
                          difference_in_standard_errors=diff/se if se>0 else None))
    audit=dict(controls_passed=True,occupation_primary=primary["occupation_pass"],
        occupation_refined=fine["occupation_pass"],eligible_cases=len(ranks),
        inward_cases=sum(r["inward"]>0 for r in ranks),ranking=ranks,
        refined_passing=sum(r["passed"] for r in refs),refined_count=len(refs),
        refinements=refs,stochastic_refinement=stats)
    dump(out/"audit.json",audit);dump(out/"complete.json",dict(complete=True))
    dump(HERE/"curvature-evidence-sha256.json",{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
        for p in out.rglob("*") if p.is_file()})
    print(f"HH-2: {audit['inward_cases']}/{len(ranks)} eligible cases bend inward; {audit['refined_passing']}/60 comparisons pass.",flush=True)
if __name__=="__main__":main()
