"""SV-1 amendment 1: unchanged Hamiltonian, matched initial physical velocities."""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1";os.environ["OMP_NUM_THREADS"]="1"
import json,hashlib,platform,subprocess,time
from datetime import datetime,timezone
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import mechanics as M
from run_campaign import save,digest
HERE=M.HERE
def matched_initial(spec):
    y,_=M.initial(spec["packet"]);q,p=M.unpack(y)
    target=p[:-1]-p[-1]/M.M
    _,w=M.pair(q[:-1],spec["ell"])
    matrix=np.eye(M.N)-spec["eta"]/(M.N-1)*w+np.ones((M.N,M.N))/M.M
    mom=np.linalg.solve(matrix,target)
    spin=-np.cross(q[:-1],mom).sum(axis=0)
    return M.pack(q,np.vstack((mom,-mom.sum(axis=0)))),spin
def matching_controls():
    rows=[]
    for packet in M.PACKETS:
        for ell in (.25,1.,3.):
            for eta in (-.8,-.4,0.,.2,.4,.6,.8):
                spec=dict(packet=packet,ell=ell,eta=eta,mu=0.)
                y,spin=matched_initial(spec);original,_=M.initial(packet)
                _,p0=M.unpack(original);v,_=M.unpack(M.rhs(0,y,spec))
                mom,ang=M.invariants(y,spin)
                vals=dict(relative_velocities=np.max(abs(v[:-1]-v[-1]-(p0[:-1]-p0[-1]/M.M))),
                    total_momentum=np.linalg.norm(mom),total_angular_momentum=np.linalg.norm(ang))
                for name,val in vals.items():
                    rows.append(dict(name=f"{packet} ell={ell} eta={eta} {name}",value=float(val),limit=1e-11,passed=bool(val<=1e-11)))
    return rows
def integrate(spec,fine=False):
    y,spin=matched_initial(spec);e0=M.energy(y,spec);comp=M.components(y,spec)
    scale=max(1,float(np.sum(abs(comp))));mom0,ang0=M.invariants(y,spin)
    t0=time.monotonic()
    sol=solve_ivp(lambda t,y:M.rhs(t,y,spec),(0,20),y,method="DOP853",
        rtol=1e-10 if fine else 1e-8,atol=1e-12 if fine else 1e-10,
        max_step=.05 if fine else .1,t_eval=np.linspace(0,20,101))
    series=[];energies=[];md=[];ad=[];metrics=[]
    for state in sol.y.T:
        met=M.metric(state,spec);mom,ang=M.invariants(state,spin)
        metrics.append(met);energies.append(M.energy(state,spec));md.append(np.linalg.norm(mom-mom0));ad.append(np.linalg.norm(ang-ang0))
    series={k:[m[k] for m in metrics] for k in metrics[0]}
    ed=max(abs(np.array(energies)-e0))/scale
    min_eig=min(series["minimum_kinetic_eigenvalue"])
    gates=dict(integration=bool(sol.success and len(sol.t)==101),energy=bool(ed<=1e-6),
        momentum=bool(max(md)<=1e-7),angular=bool(max(ad)<=1e-6),kinetic=bool(min_eig>=.2-1e-12))
    return dict(spec=spec,fine=fine,success=sol.success,message=sol.message,nfev=sol.nfev,
        seconds=time.monotonic()-t0,energy_initial=e0,energy_scale=scale,energy_drift=ed,
        momentum_drift=max(md),angular_drift=max(ad),minimum_kinetic_eigenvalue=min_eig,
        preparation_energy=e0+M.N*spec["mu"]/M.EPS+np.dot(spin,spin)/(2*M.INERTIA),
        compensating_source_spin=spin,spin_energy=np.dot(spin,spin)/(2*M.INERTIA),
        gates=gates,passed=all(gates.values()),times=sol.t,energy_series=energies,
        metrics=series,late={k:float(np.mean(v[-20:])) for k,v in series.items()}),sol.y.T
def main():
    out=HERE/"evidence-matched-v1";out.mkdir(exist_ok=False)
    save(out/"manifest.json",dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=M.ROOT,text=True).strip(),
        utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        hashes={p.relative_to(M.ROOT).as_posix():digest(p) for p in list(HERE.glob("*.py"))+list(HERE.glob("*.md"))}))
    controls=matching_controls();save(out/"controls.json",controls)
    if not all(r["passed"] for r in controls):raise AssertionError("Preserved control failure")
    print(f"{len(controls)} initial-velocity matching controls pass",flush=True)
    rows=[];trajectories={};fine_rows=[];refinements=[];fine_traj={}
    for packet in M.PACKETS:
        for mu in (0.,1.):
            for ell in (.25,1.,3.):
                for eta in (-.8,-.4,0.,.2,.4,.6,.8):
                    spec=dict(id=f"{packet}-mu{mu:g}-ell{ell:g}-eta{eta:g}",packet=packet,mu=mu,ell=ell,eta=eta)
                    row,path=integrate(spec);rows.append(row);trajectories[spec["id"]]=path
                print(f"{len(rows)}/252 cases; failures so far {sum(not r['passed'] for r in rows)}",flush=True)
                save(out/"runs.json",rows)
    np.savez_compressed(out/"trajectories.npz",**trajectories)
    for old in rows:
        s=old["spec"]
        if s["ell"]!=1 or s["eta"] not in (-.8,0.,.8):continue
        row,path=integrate(s,True);fine_rows.append(row);fine_traj[s["id"]]=path
        coarse=trajectories[s["id"]]
        n=M.N+1
        error=float(np.max(np.linalg.norm((coarse[:,:3*n]-path[:,:3*n]).reshape(-1,3),axis=1))/max(1,np.max(np.linalg.norm(path[:,:3*n].reshape(-1,3),axis=1))))
        refinements.append(dict(id=s["id"],relative_position_error=error,limit=.001,
            passed=bool(error<=.001 and row["passed"]),fine_numerical_passed=row["passed"]))
    save(out/"refined-runs.json",fine_rows);save(out/"refinements.json",refinements)
    np.savez_compressed(out/"refined-trajectories.npz",**fine_traj)
    comparisons=[]
    for r in rows:
        s=r["spec"]
        if s["eta"]==0:continue
        b=next(x for x in rows if x["spec"]["packet"]==s["packet"] and x["spec"]["mu"]==s["mu"] and x["spec"]["ell"]==s["ell"] and x["spec"]["eta"]==0)
        comparisons.append(dict(id=s["id"],packet=s["packet"],eta=s["eta"],mu=s["mu"],ell=s["ell"],
            both_numerically_passed=bool(r["passed"] and b["passed"]),
            late_radius_ratio=r["late"]["rms_radius"]/b["late"]["rms_radius"],
            late_width_ratio=r["late"]["packet_width"]/b["late"]["packet_width"],
            retained_difference=r["late"]["retained"]-b["late"]["retained"],
            circulation_difference=r["late"]["circulation"]-b["late"]["circulation"],
            heading_difference=r["late"]["heading"]-b["late"]["heading"],
            neighbor_weight=r["late"]["neighbor_weight"]))
    save(out/"comparisons.json",comparisons)
    save(out/"complete.json",dict(runs=len(rows),passed_runs=sum(r["passed"] for r in rows),
        refined_runs=len(fine_rows),passed_refinements=sum(r["passed"] for r in refinements),
        maximum_scaled_energy_drift=max(r["energy_drift"] for r in rows),
        maximum_momentum_drift=max(r["momentum_drift"] for r in rows),
        maximum_angular_drift=max(r["angular_drift"] for r in rows),
        physical_solution=False,observational_fit=False,causal_field_derived=False,
        production_derived=False,massless_limit_derived=False))
    save(HERE/"matched-evidence-sha256.json",{p.relative_to(HERE).as_posix():digest(p) for p in sorted(out.iterdir())})
    print(f"Finished {len(rows)} cases and {len(fine_rows)} refinements; {sum(r['passed'] for r in rows)} primary numerical passes, {sum(r['passed'] for r in refinements)} refinement passes",flush=True)
if __name__=="__main__":main()
