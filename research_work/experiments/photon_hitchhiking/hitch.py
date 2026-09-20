"""HH-1 mean occupation and discrete attachment/release trajectories."""
from __future__ import annotations
import hashlib,itertools,json,platform,subprocess,sys,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
STREAM=HERE.parent/"matter_streams"
sys.path.insert(0,str(STREAM))
from campaign import density,dump,IMPACTS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def row(h,on,off,k):
    return dict(id=f"H{h:g}-A{on:g}-D{off:g}-K{k:g}",pitch=h,on=on,off=off,kappa=k)
def configs():return [row(h,a,b,k) for h,a,b,k in itertools.product((.5,2.,4.),(1.,10.,100.),(.1,1.,10.),(1.,4.,16.))]
REPS=[row(2,10,1,4),row(4,100,.1,16),row(.5,1,10,1)]
REFINE=REPS+[row(.5,100,.1,16),row(4,1,10,1),row(2,10,.1,4)]
def pars(rows,repeat):return {k:np.repeat([r[k] for r in rows],repeat) for k in ("pitch","on","off","kappa")}

def local(t,y,p):
    psi=np.arctan2(y[:,1],y[:,0])+np.arctan(p["pitch"])
    gate=((1+np.cos(psi-y[:,2]))/2)**2
    rate=p["on"]*density(t,np.hypot(y[:,0],y[:,1]),p["pitch"])*gate
    return psi,rate

def mrhs(t,y,p):
    psi,rate=local(t,y,p)
    capture=rate*(1-y[:,3])
    direction=np.stack([np.cos(psi),np.sin(psi)],axis=-1)
    turn=p["kappa"]*(np.cos(y[:,2])*y[:,5]-np.sin(y[:,2])*y[:,4])
    z=np.zeros_like(y)
    z[:,0]=np.cos(y[:,2]);z[:,1]=np.sin(y[:,2]);z[:,2]=turn
    z[:,3]=capture-p["off"]*y[:,3]
    z[:,4:6]=capture[:,None]*direction-p["off"][:,None]*y[:,4:6]
    z[:,6]=np.abs(turn);z[:,7]=capture;z[:,8]=p["off"]*y[:,3]
    return z

def rk(t,y,dt,p):
    k1=mrhs(t,y,p);k2=mrhs(t+dt/2,y+dt*k1/2,p)
    k3=mrhs(t+dt/2,y+dt*k2/2,p);k4=mrhs(t+dt,y+dt*k3,p)
    return y+dt*(k1+2*k2+2*k3+k4)/6

def observer(y,nxt,t,dt,live,exits):
    hit=live&(nxt[:,0]>=20)
    frac=np.divide(20-y[:,0],nxt[:,0]-y[:,0],out=np.zeros(len(y)),where=hit)
    nxt[hit]=y[hit]+frac[hit,None]*(nxt[hit]-y[hit])
    exits[hit]=t+frac[hit]*dt
    nxt[~live]=y[~live];live[hit]=False
    return nxt

def mean_trace(rows,out,dt):
    out.mkdir(exist_ok=False)
    p=pars(rows,10);imp=np.tile(IMPACTS,len(rows));n=len(imp)
    y=np.zeros((n,9));y[:,0]=-20;y[:,1]=imp
    live=np.ones(n,dtype=bool);exits=np.full(n,np.nan)
    hist=[];times=[];max_p=0.;min_p=0.;moment_violation=0.
    for step in range(round(80/dt)+1):
        t=step*dt
        max_p=max(max_p,float(y[:,3].max()));min_p=min(min_p,float(y[:,3].min()))
        moment_violation=max(moment_violation,float(np.max(np.linalg.norm(y[:,4:6],axis=1)-y[:,3])))
        if step%round(.2/dt)==0:hist.append(y.copy());times.append(t)
        if step==round(80/dt):break
        y=observer(y,rk(t,y,dt,p),t,dt,live,exits)
        if not np.isfinite(y).all():raise FloatingPointError(f"nonfinite at {t}")
    results=[]
    for i in range(n):
        angle=np.arctan2(np.sin(y[i,2]),np.cos(y[i,2]))
        results.append(dict(id=rows[i//10]["id"],impact=imp[i],arrived=bool(np.isfinite(exits[i])),
            exit_time=exits[i],angle=float(angle),observer_y=y[i,1] if np.isfinite(exits[i]) else None,
            terminal_occupation=y[i,3],absolute_turn=y[i,6],expected_captures=y[i,7],
            expected_releases=y[i,8],stream_recoil_per_photon_momentum=[1-np.cos(angle),-np.sin(angle)]))
    pairs=[]
    for i,r in enumerate(rows):
        by={q["impact"]:q for q in results[i*10:(i+1)*10]}
        for b in (.5,1,2,4,8):
            a,c=by[-b],by[b];ok=a["arrived"] and c["arrived"]
            pairs.append(dict(id=r["id"],impact=b,arrived=ok,inward=(a["angle"]-c["angle"])/2 if ok else None,
                              sideways=(a["angle"]+c["angle"])/2 if ok else None))
    data=dict(rows=rows,results=results,pairs=pairs,dt=dt,min_p=min_p,max_p=max_p,
              moment_violation=moment_violation,occupation_pass=min_p>=-1e-8 and max_p<=1+1e-8 and moment_violation<1e-8)
    np.savez_compressed(out/"trajectories.npz",state=np.stack(hist,axis=1),time=np.array(times),
        ids=np.repeat([r["id"] for r in rows],10),impact=imp,exit_time=exits)
    dump(out/"results.json",data)
    print(f"{out.name}: {n} mean rays, {np.isfinite(exits).sum()} arrivals, occupation pass {data['occupation_pass']}",flush=True)
    return data

def fly(y,dt,bound,heading,kappa):
    def f(q):
        turn=np.where(bound,kappa*np.sin(heading-q[:,2]),0.)
        return np.stack([np.cos(q[:,2]),np.sin(q[:,2]),turn],axis=-1)
    k1=f(y);k2=f(y+dt*k1/2);k3=f(y+dt*k2/2);k4=f(y+dt*k3)
    return y+dt*(k1+2*k2+2*k3+k4)/6

def stochastic(rows,out,dt,seed):
    out.mkdir(exist_ok=False)
    count=128;bs=(-2.,-1.,1.,2.)
    p=pars(rows,4*count)
    imp=np.tile(np.repeat(bs,count),len(rows));n=len(imp)
    y=np.zeros((n,3));y[:,0]=-20;y[:,1]=imp
    bound=np.zeros(n,dtype=bool);heading=np.zeros(n)
    captures=np.zeros(n,dtype=int);releases=np.zeros(n,dtype=int)
    live=np.ones(n,dtype=bool);exits=np.full(n,np.nan)
    rng=np.random.default_rng(seed)
    hist=[];attachments=[];times=[]
    for step in range(round(80/dt)+1):
        t=step*dt
        if step%round(.2/dt)==0:
            hist.append(y.copy());attachments.append(bound.copy());times.append(t)
        if step==round(80/dt):break
        y=observer(y,fly(y,dt/2,bound,heading,p["kappa"]),t,dt/2,live,exits)
        psi,rate=local(t+dt/2,y,p)
        old=bound.copy()
        attach=live&~old&(rng.random(n)<-np.expm1(-rate*dt))
        release=live&old&(rng.random(n)<-np.expm1(-p["off"]*dt))
        bound[attach]=True;heading[attach]=psi[attach];captures[attach]+=1
        bound[release]=False;releases[release]+=1
        y=observer(y,fly(y,dt/2,bound,heading,p["kappa"]),t+dt/2,dt/2,live,exits)
        if not np.isfinite(y).all():raise FloatingPointError(f"nonfinite MC at {t}")
    angle=np.arctan2(np.sin(y[:,2]),np.cos(y[:,2]))
    ensembles=[]
    for i,r in enumerate(rows):
        for j,b in enumerate(bs):
            sl=slice((i*4+j)*count,(i*4+j+1)*count)
            hit=np.isfinite(exits[sl])
            vals=angle[sl][hit]
            avg=float(vals.mean()) if len(vals) else None
            sd=float(vals.std(ddof=1)) if len(vals)>1 else None
            ensembles.append(dict(id=r["id"],impact=b,n=count,arrivals=int(hit.sum()),
                mean_angle=avg,sd_angle=sd,se_angle=sd/np.sqrt(len(vals)) if sd is not None else None,
                angular_blur_to_mean=sd/abs(avg) if sd is not None and abs(avg)>1e-15 else None,
                mean_captures=float(captures[sl].mean()),max_captures=int(captures[sl].max()),
                mean_delay=float(np.mean(exits[sl][hit]-40)) if hit.any() else None))
    np.savez_compressed(out/"trajectories.npz",state=np.stack(hist,axis=1),bound=np.stack(attachments,axis=1),
        time=np.array(times),ids=np.repeat([r["id"] for r in rows],4*count),impact=imp,
        exit_time=exits,heading=heading,captures=captures,releases=releases)
    data=dict(dt=dt,seed=seed,per_ensemble=count,rows=rows,ensembles=ensembles,
              total_histories=n,arrivals=int(np.isfinite(exits).sum()))
    dump(out/"results.json",data)
    print(f"{out.name}: {n} event histories, {data['arrivals']} arrivals",flush=True)
    return data

def controls():
    checks=[]
    def add(name,value,limit):checks.append(dict(name=name,value=float(value),limit=limit,passed=abs(value)<limit))
    a,b=.7,.2;p=0.;m=np.zeros(2);direction=np.array([.6,.8]);dt=.02
    z=np.zeros(3)
    def rhs(q):return np.r_[a*(1-q[0])-b*q[0],a*(1-q[0])*direction-b*q[1:]]
    for _ in range(250):
        k1=rhs(z);k2=rhs(z+dt*k1/2);k3=rhs(z+dt*k2/2);k4=rhs(z+dt*k3)
        z+=dt*(k1+2*k2+2*k3+k4)/6
    exact=a/(a+b)*(1-np.exp(-(a+b)*5))
    add("uniform occupation exact solution",z[0]-exact,1e-9)
    add("uniform carried direction equals p*t",np.max(np.abs(z[1:]-exact*direction)),1e-9)
    # Actual ray equation controls at a location inside the source field.
    y=np.zeros((1,9));y[0,:2]=[1.,1.]
    param={k:np.array([v]) for k,v in dict(pitch=2.,on=0.,off=1.,kappa=4.).items()}
    dd=mrhs(0,y,param)
    add("zero capture and empty state stays free",np.max(np.abs(dd[0,2:])),1e-14)
    y[0,3]=.5;y[0,4:6]=[.3,.4];param["on"][:]=10;param["kappa"][:]=0
    add("zero guide coupling leaves heading unchanged",mrhs(0,y,param)[0,2],1e-14)
    if not all(c["passed"] for c in checks):raise AssertionError(checks)
    return checks

def main():
    out=HERE/"evidence-v1";out.mkdir(exist_ok=False)
    files=list(HERE.glob("*.py"))+list(HERE.glob("*.md"))+[STREAM/"campaign.py",STREAM/"protocol.md"]
    dump(out/"manifest.json",dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=REPO,text=True).strip(),
        python=platform.python_version(),numpy=np.__version__,
        sources={p.relative_to(REPO).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}))
    dump(out/"controls.json",controls())
    primary=mean_trace(configs(),out/"mean-screen",.02)
    fine=mean_trace(REFINE,out/"mean-refinement",.01)
    mc0=stochastic(REPS,out/"events-dt001",.01,190619)
    mc1=stochastic(REPS,out/"events-dt0005",.005,190620)
    coarse={(r["id"],r["impact"]):r for r in primary["results"]}
    refinements=[]
    for r in fine["results"]:
        c=coarse[(r["id"],r["impact"])]
        both=r["arrived"] and c["arrived"]
        da=abs(r["angle"]-c["angle"]) if both else None
        dy=abs(r["observer_y"]-c["observer_y"]) if both else None
        refinements.append(dict(id=r["id"],impact=r["impact"],angle_difference=da,position_difference=dy,
            passed=r["arrived"]==c["arrived"] and (da is None or da<.01) and (dy is None or dy<.01)))
    stats=[]
    for c,r in zip(mc0["ensembles"],mc1["ensembles"]):
        se=np.hypot(c["se_angle"],r["se_angle"])
        difference=abs(c["mean_angle"]-r["mean_angle"])
        stats.append(dict(id=r["id"],impact=r["impact"],absolute_mean_difference=difference,
            combined_standard_error=float(se),difference_in_standard_errors=difference/se if se>0 else None))
    inward_ids=[]
    for rr in configs():
        pairs=[r for r in primary["pairs"] if r["id"]==rr["id"] and r["impact"] in (1,2,4)]
        if all(r["arrived"] for r in pairs) and np.mean([r["inward"] for r in pairs])>0:inward_ids.append(rr["id"])
    audit=dict(mean_cases=81,mean_primary_rays=810,mean_refinement_rays=60,event_histories=3072,
        occupation_primary=primary["occupation_pass"],occupation_refined=fine["occupation_pass"],
        refined_passing=sum(r["passed"] for r in refinements),refined_count=len(refinements),
        refinements=refinements,stochastic_refinement=stats,positive_mean_inward_cases=inward_ids)
    dump(out/"audit.json",audit);dump(out/"complete.json",dict(complete=True))
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout="constrained")
    for rr in REPS:
        pp=[r for r in primary["pairs"] if r["id"]==rr["id"] and r["arrived"]]
        axs[0].plot([r["impact"] for r in pp],[r["inward"] for r in pp],marker="o",label=rr["id"])
        ens=[r for r in mc1["ensembles"] if r["id"]==rr["id"]]
        axs[1].errorbar([r["impact"] for r in ens],[r["mean_angle"] for r in ens],
                        yerr=[r["sd_angle"] for r in ens],marker="o",capsize=3,label=rr["id"])
    axs[0].axhline(0,color="black",lw=.8)
    axs[0].set(xlabel="Absolute impact",ylabel="Paired inward bend (rad)",title="Mean-state attachment: conditional ray test")
    axs[1].set(xlabel="Signed impact",ylabel="Exit angle (rad)",title="Discrete captures: mean +/- one standard deviation")
    for ax in axs:ax.legend(fontsize=7)
    fig.suptitle("HH-1: temporary coupling can turn light; scatter is a separate prediction")
    fig.savefig(HERE/"results.png",dpi=160);plt.close(fig)
    dump(HERE/"evidence-sha256.json",{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
        for p in out.rglob("*") if p.is_file()})
    print(f"HH-1 complete: {len(inward_ids)}/81 mean cases have positive paired inward bending; "
          f"{audit['refined_passing']}/60 timestep comparisons pass.",flush=True)
if __name__=="__main__":main()
