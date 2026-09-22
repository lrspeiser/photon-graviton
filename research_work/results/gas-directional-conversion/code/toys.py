"""Execute toy geometry sweep before looking at observational residuals."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse import lil_matrix,csc_matrix,bmat
from scipy.sparse.linalg import expm_multiply,spsolve
from scipy.integrate import solve_ivp
from core import AxisModel,transfer,gas_density,baseline

def energy_evolution(n=64,alpha=.3,beta=.3,mu=.05,duration=400.):
    """Initially empty, source-fed three-channel finite-volume system.
    Total channel energy plus boundary escape equals integrated input.
    Fixed speed v=1 model kpc/time; no years or physical source sufficiency claim.
    """
    rf=np.linspace(0.,16.,n+1);r=(rf[:-1]+rf[1:])/2;dx=rf[1]-rf[0];adv=1/dx
    rho=gas_density(r,np.array([mu]),1e9,2.,.1)[:,0]
    ap=alpha*rho/1e6;am=beta*(rho/1e7)*rho/1e6
    _,density,_=baseline(r,1.,1.,20.,2.)
    # A bounded toy capture/residence kernel, not an astrophysical lifetime.
    residence=1+20*(4*np.pi*r*r*density)/(4*np.pi*r*r*density).max()
    capture=.2;release=capture/residence
    N=3*n+2;gen=lil_matrix((N,N));source=np.zeros(N);source[0]=source[n]=.5
    escape=3*n
    for i in range(n):
        L=i;T=n+i;B=2*n+i
        for a,b,k in [(L,T,ap[i]),(T,L,am[i]),(T,B,capture),(B,T,release[i])]:
            gen[b,a]+=k;gen[a,a]-=k
        for a,offset in [(L,0),(T,n)]:
            dest=a+1 if i<n-1 else escape
            gen[dest,a]+=adv;gen[a,a]-=adv
    gen[:,-1]=source[:,None];gen=csc_matrix(gen)
    state=np.zeros(N);state[-1]=1
    hist=expm_multiply(gen,state,start=0,stop=duration,num=33,endpoint=True)
    times=np.linspace(0,duration,len(hist));total=hist[:,:-1].sum(1)
    Q=gen[:3*n,:3*n];steady=spsolve(-Q,source[:3*n])
    fc=steady[n:2*n]/(steady[:n]+steady[n:2*n]);expect=[];prev=.5
    for aa,bb in zip(ap,am):
        prev=(prev+aa*dx)/(1+(aa+bb)*dx);expect.append(prev)
    shut=expm_multiply(Q,steady,start=0,stop=200.,num=3)
    return dict(cells=n,alpha=alpha,beta=beta,ray_mu=mu,final_time=duration,
       maximum_energy_ledger_absolute_error=float(np.max(abs(total-times))),
       minimum_energy=float(hist[:,:-1].min()),
       steady_fraction_recurrence_error=float(np.max(abs(fc-expect))),
       final_bound_fraction_of_steady=float(hist[-1,2*n:3*n].sum()/steady[2*n:3*n].sum()),
       final_stored_energy=float(hist[-1,2*n:3*n].sum()),
       final_traveling_energy=float(hist[-1,:2*n].sum()),final_escaped_energy=float(hist[-1,escape]),
       bound_after_source_off_0_100_200=[float(s[2*n:3*n].sum()) for s in shut],
       elapsed_time_units='model units; no inferred astrophysical lifetime',
       microscopic_limit='Energy bookkeeping for phenomenological transitions, not a recoil/quantum/gravity-work Hamiltonian')

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    target=a.output/'toy-results.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    geometries=[('sphere','sphere',1.),('thin_disk','disk',.1),('thick_disk','disk',.5),('ring','ring',.1)]
    opacities=[0.,.03,.3,3.]
    for name,kind,height in geometries:
      base=AxisModel(1.,1.,20.,2.,nr=384,nmu=64,lmax=20).solve(0.,2.,0.,0.)
      bends={s:base.bend(2.,s) for s in ['face','edge','edge_vertical']}
      for ep in opacities:
       for dp in opacities:
        model=AxisModel(1.,1.,20.,2.,nr=384,nmu=64,lmax=20).solve(1e9,2.,ep,dp,height=height,kind=kind)
        res=model.diagnostic();res.update(gas_geometry=name,alpha=ep,beta=dp,
          bend_ratio={s:model.bend(2.,s)/bends[s] for s in bends})
        rows.append(res)
      print('TOY',name,'complete',flush=True)
    # Independent local characteristic ODE versus constant-cell exact update.
    r=np.linspace(0,8,101);rho=np.full((100,1),4e7)
    f=transfer(r,rho,.3,.3)[:,0]
    sol=solve_ivp(lambda x,y:[12*(1-y[0])-48*y[0]],(0,8),[.5],t_eval=r[1:],rtol=1e-10,atol=1e-12)
    energy=[energy_evolution(64,.3,.3,m) for m in [.05,.8]]
    refinements=[]
    for kind in ['disk','ring','sphere']:
      coarse=AxisModel(1,1,20,2,nr=384,nmu=64,lmax=20).solve(1e9,2,.3,.3,kind=kind)
      fine=AxisModel(1,1,20,2,nr=768,nmu=128,lmax=40).solve(1e9,2,.3,.3,kind=kind)
      pts=np.array([.5,1,2,4,8])
      refinements.append(dict(kind=kind,equatorial_force_relative=float(np.max(abs(coarse.radial(pts)/fine.radial(pts)-1))),
        polar_force_relative=float(np.max(abs(coarse.radial(pts,1)/fine.radial(pts,1)-1))),
        lens_relative={s:coarse.bend(2,s)/fine.bend(2,s)-1 for s in bends}))
    result=dict(experiment='JR-5 toy stage',baseline='a5805452b97ef0856d364de90e394572a8fd2146',
        geometry_combinations=len(rows),gas_mass_Msun=1e9,gas_radial_scale_kpc=2.,
        rows=rows,local_ode_error=float(np.max(abs(f-sol.y[0]))),energy_evolutions=energy,refinements=refinements)
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('ENERGY',json.dumps(energy),flush=True);print('REFINEMENT',json.dumps(refinements),flush=True)
if __name__=='__main__':main()
