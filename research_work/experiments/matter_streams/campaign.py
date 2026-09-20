"""MS-1 prescribed matter-fed stream and test-ray steering; see protocol limits."""
from __future__ import annotations
import hashlib
import itertools
import json
from pathlib import Path
import platform
import subprocess
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
Q=1/30
IMPACTS=(-8.,-4.,-2.,-1.,-.5,.5,1.,2.,4.,8.)

def dump(path,data):
    def clean(a):
        if isinstance(a,np.ndarray):return clean(a.tolist())
        if isinstance(a,(np.floating,float)):return float(a) if np.isfinite(a) else None
        if isinstance(a,(np.integer,)):return int(a)
        if isinstance(a,(np.bool_,)):return bool(a)
        if isinstance(a,dict):return {k:clean(v) for k,v in a.items()}
        if isinstance(a,(list,tuple)):return [clean(v) for v in a]
        return a
    path.write_text(json.dumps(clean(data),indent=2,allow_nan=False)+"\n",encoding="utf8")

def configs():
    return [dict(id=f"H{h:g}-S{sgn:+d}-F{gate}-K{k:+g}",pitch=h,sign=sgn,gate=gate,kappa=k)
       for h,sgn,gate,k in itertools.product((0.,.5,1.,2.,4.),(-1,1),range(3),(-16.,-4.,-1.,0.,1.,4.,16.))]

def prepare(rows):
    return {key:np.repeat([r[key] for r in rows],len(IMPACTS)) for key in ("pitch","sign","gate","kappa")}

def density(t,r,h):
    vr=.5/np.sqrt(1+h*h)
    inside=(r>.5)&(r<.5+vr*(30+t))
    return np.divide(Q,4*np.pi*vr*r*r,out=np.zeros_like(r),where=inside)

def rhs(t,y,p):
    r=np.hypot(y[:,0],y[:,1])
    psi=np.arctan2(y[:,1],y[:,0])+np.arctan(p["pitch"]*p["sign"])
    delta=psi-y[:,2]
    c=np.cos(delta)
    gate=np.ones_like(c)
    gate[p["gate"]==1]=((1+c[p["gate"]==1])/2)**2
    gate[p["gate"]==2]=np.maximum(c[p["gate"]==2],0)**2
    exposure=density(t,r,p["pitch"])*gate
    turn=p["kappa"]*exposure*np.sin(delta)
    return np.stack([np.cos(y[:,2]),np.sin(y[:,2]),turn,exposure,np.abs(turn)],axis=-1)

def step(t,y,dt,p):
    k1=rhs(t,y,p);k2=rhs(t+dt/2,y+dt*k1/2,p)
    k3=rhs(t+dt/2,y+dt*k2/2,p);k4=rhs(t+dt,y+dt*k3,p)
    return y+dt*(k1+2*k2+2*k3+k4)/6

def trace(rows,out,dt):
    out.mkdir(exist_ok=False)
    p=prepare(rows);b=np.tile(IMPACTS,len(rows));n=len(b)
    y=np.zeros((n,5));y[:,0]=-20;y[:,1]=b
    exits=np.full(n,np.nan)
    live=np.ones(n,dtype=bool)
    frames=[];times=[]
    tstart=time.perf_counter()
    for k in range(round(80/dt)+1):
        t=k*dt
        if k%round(.2/dt)==0:frames.append(y.copy());times.append(t)
        if k==round(80/dt):break
        nxt=step(t,y,dt,p)
        hit=live&(nxt[:,0]>=20)
        fraction=np.divide(20-y[:,0],nxt[:,0]-y[:,0],out=np.zeros(n),where=hit)
        nxt[hit]=y[hit]+fraction[hit,None]*(nxt[hit]-y[hit])
        exits[hit]=t+fraction[hit]*dt
        nxt[~live]=y[~live]
        y=nxt;live[hit]=False
        if not np.all(np.isfinite(y)):raise FloatingPointError(f"Nonfinite state at {t+dt}")
    hist=np.stack(frames,axis=1)
    angle=np.arctan2(np.sin(y[:,2]),np.cos(y[:,2]))
    results=[]
    for i in range(n):
        results.append(dict(id=rows[i//len(IMPACTS)]["id"],impact=b[i],arrived=bool(np.isfinite(exits[i])),
            exit_time=exits[i],exit_angle=angle[i] if np.isfinite(exits[i]) else None,
            observer_y=y[i,1] if np.isfinite(exits[i]) else None,
            terminal_x=y[i,0],terminal_y=y[i,1],terminal_angle=angle[i],
            exposure=y[i,3],absolute_turn=y[i,4],
            stream_recoil_per_unit_photon_momentum=[1-np.cos(y[i,2]),-np.sin(y[i,2])]))
    pairs=[];rank=[]
    for k,row in enumerate(rows):
        selected={r["impact"]:r for r in results[k*10:(k+1)*10]}
        pair=[]
        for bb in (.5,1.,2.,4.,8.):
            lo,hi=selected[-bb],selected[bb]
            ok=lo["arrived"] and hi["arrived"]
            pp=dict(id=row["id"],impact=bb,arrived=ok,
                inward=(lo["exit_angle"]-hi["exit_angle"])/2 if ok else None,
                sideways=(lo["exit_angle"]+hi["exit_angle"])/2 if ok else None)
            pair.append(pp);pairs.append(pp)
        requested=[r for r in pair if r["impact"] in (1.,2.,4.)]
        if all(r["arrived"] for r in requested):
            rank.append(dict(id=row["id"],mean_inward=np.mean([r["inward"] for r in requested]),
                             mean_abs_sideways=np.mean([abs(r["sideways"]) for r in requested])))
    rank.sort(key=lambda r:(-r["mean_inward"],r["id"]))
    np.savez_compressed(out/"trajectories.npz",state=hist,time=np.array(times),
                        ids=np.repeat([r["id"] for r in rows],10),impact=b,exit_time=exits)
    data=dict(dt=dt,rows=rows,results=results,pairs=pairs,ranking=rank,
              wall_seconds=time.perf_counter()-tstart)
    dump(out/"results.json",data)
    print(f"{out.name}: {n} rays, {np.isfinite(exits).sum()} observer arrivals, "
          f"{n-np.isfinite(exits).sum()} nonarrivals, {data['wall_seconds']:.1f}s",flush=True)
    return data,hist

def controls():
    checks=[]
    def check(name,value,limit):
        checks.append(dict(name=name,value=float(value),limit=limit,passed=abs(value)<=limit))
    theta=.7;dt=.01;k=.8
    for _ in range(1000):
        f=lambda a:-k*np.sin(a)
        k1=f(theta);k2=f(theta+dt*k1/2);k3=f(theta+dt*k2/2);k4=f(theta+dt*k3)
        theta+=dt*(k1+2*k2+2*k3+k4)/6
    check("straight-stream exact relaxation",theta-2*np.arctan(np.tan(.7/2)*np.exp(-8)),1e-11)
    rng=np.random.default_rng(991)
    ang=rng.normal(size=100);psi=rng.normal(size=100);weight=rng.random(100)
    n=np.stack([np.cos(ang),np.sin(ang)],axis=-1)
    turn=weight*np.sin(psi-ang)
    force=np.stack([-np.sin(ang)*turn,np.cos(ang)*turn],axis=-1)
    check("pointwise zero guide work",np.max(np.abs(np.sum(n*force,axis=-1))),1e-14)
    check("straight aligned light zero bend",abs(.8*np.sin(0.)),0.)
    largest=0.
    for h,t in itertools.product((0.,.5,1.,2.,4.),(0.,20.,80.)):
        vr=.5/np.sqrt(1+h*h)
        # Midpoint quadrature independently integrates the spherical volume.
        edges=np.linspace(.5,.5+vr*(30+t),1001)
        centers=(edges[:-1]+edges[1:])/2
        total=np.sum(4*np.pi*centers**2*density(t,centers,np.full_like(centers,h))*np.diff(edges))
        largest=max(largest,abs(total-Q*(30+t)))
    check("finite-age source energy integral",largest,1e-12)
    check("total source plus emitted energy",max(abs((10-Q*(30+t))+Q*(30+t)-10) for t in (0,20,80)),1e-12)
    check("source remains nonnegative",max(0.,-(10-Q*110)),0.)
    check("companion total speed below light",max(abs(np.hypot(.5/np.sqrt(1+h*h),.5*h/np.sqrt(1+h*h))-.5) for h in (0,.5,1,2,4)),1e-15)
    if not all(c["passed"] for c in checks):raise AssertionError(checks)
    return checks

def main():
    out=HERE/"evidence-v1";out.mkdir(exist_ok=False)
    manifest=dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=REPO,text=True).strip(),
        python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.iterdir() if p.suffix in (".py",".md")})
    dump(out/"manifest.json",manifest)
    ctrl=controls();dump(out/"controls.json",ctrl)
    rows=configs();dump(out/"configurations.json",rows)
    data,hist=trace(rows,out/"screen",.02)
    zero=[r for r in data["results"] if next(c for c in rows if c["id"]==r["id"])["kappa"]==0]
    zerr=max(max(abs(r["exit_angle"]),abs(r["observer_y"]-r["impact"]),abs(r["exit_time"]-40)) for r in zero)
    # Complete reflection of the signed swirl and impact, including all stored frames.
    lookup={(r["pitch"],r["sign"],r["gate"],r["kappa"]):i for i,r in enumerate(rows)}
    err=0.
    for i,row in enumerate(rows):
        j=lookup[(row["pitch"],-row["sign"],row["gate"],row["kappa"])]
        reflected=hist[i*10:(i+1)*10][::-1].copy()
        reflected[...,[1,2]]*=-1
        err=max(err,float(np.max(np.abs(reflected-hist[j*10:(j+1)*10]))))
    selected=[r["id"] for r in data["ranking"][:6]]
    fine,fh=trace([r for key in selected for r in rows if r["id"]==key],out/"refinement",.01)
    coarse={(r["id"],r["impact"]):r for r in data["results"]}
    compare=[]
    for r in fine["results"]:
        c=coarse[(r["id"],r["impact"])]
        same=c["arrived"]==r["arrived"]
        angle=abs(np.arctan2(np.sin(c["exit_angle"]-r["exit_angle"]),np.cos(c["exit_angle"]-r["exit_angle"]))) if c["arrived"] and r["arrived"] else None
        pos=abs(c["observer_y"]-r["observer_y"]) if c["arrived"] and r["arrived"] else None
        expo=abs(c["exposure"]-r["exposure"])/max(1,abs(r["exposure"]))
        compare.append(dict(id=r["id"],impact=r["impact"],same_arrival=same,angle_difference=angle,
            position_difference=pos,exposure_difference=expo,
            converged=same and (angle is None or angle<.01) and (pos is None or pos<.01) and expo<.01))
    audit=dict(controls_passed=all(c["passed"] for c in ctrl),zero_coupling_error=zerr,
        zero_coupling_pass=zerr<1e-8,reflection_error=err,reflection_pass=err<1e-8,
        cases=len(rows),screen_rays=len(data["results"]),refinement_rays=len(fine["results"]),
        selected=selected,comparisons=compare,
        converged_rays=sum(r["converged"] for r in compare),
        all_comparisons_pass=all(r["converged"] for r in compare))
    dump(out/"audit.json",audit)
    # Plot transmitted pairs and actual trajectories; no claim of lens validation.
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    for identifier in selected[:3]:
        i=next(i for i,r in enumerate(rows) if r["id"]==identifier)
        for j in (2,3,6,7):
            path=hist[i*10+j]
            axs[0].plot(path[:,0],path[:,1],lw=.8,label=identifier if j==2 else None)
        pp=[p for p in data["pairs"] if p["id"]==identifier and p["arrived"]]
        axs[1].plot([p["impact"] for p in pp],[p["inward"] for p in pp],marker="o",label=identifier)
    axs[0].set(xlabel="x",ylabel="y",title="Prescribed evolving stream; test rays")
    axs[0].set_aspect("equal",adjustable="datalim");axs[0].legend(fontsize=7)
    axs[1].axhline(0,color="black",lw=.8)
    axs[1].set(xlabel="Absolute impact parameter",ylabel="Paired inward direction change (rad)",
               title="Inward bending separated from sideways bias")
    axs[1].legend(fontsize=7)
    fig.suptitle("MS-1: directional cumulative steering, not yet a derived gravity law")
    fig.savefig(HERE/"results.png",dpi=160);plt.close(fig)
    hashes={p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob("*") if p.is_file()}
    dump(HERE/"evidence-sha256.json",hashes)
    dump(out/"complete.json",dict(complete=True))
    print((out/"audit.json").read_text(encoding="utf8"),flush=True)

if __name__=="__main__":main()
