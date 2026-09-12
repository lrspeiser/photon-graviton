"""Spherical outer-law diagnostic. All quantities dimensionless; no data fit."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

HERE=Path(__file__).resolve().parent

def grid(n,outer):
    return np.r_[0.,np.geomspace(1e-6,outer,n-1)]

def gravity(r,rho):
    mass=4*np.pi*cumulative_trapezoid(rho*r*r,r,initial=0)
    shells=4*np.pi*cumulative_trapezoid((rho*r)[::-1],r[::-1],initial=0)[::-1]*-1
    depth=shells.copy()
    depth[1:]+=mass[1:]/r[1:]
    return mass,depth

def run(p,C,outer,n,rtol,atol):
    r=grid(n,outer);wb=1/np.sqrt(1+r*r)
    def rhs(s,rho):
        _,wd=gravity(r,rho)
        W=wb+wd
        return C*(W/(1+W))**p
    sol=solve_ivp(rhs,(0,1),np.zeros(n),rtol=rtol,atol=atol)
    assert sol.success,sol.message
    rho=sol.y[:,-1]
    assert np.min(sol.y)>=-1e-12 and np.max(sol.y-C*sol.t[None,:])<1e-8
    mass,wd=gravity(r,rho)
    vcs=[float(np.interp(x,r,mass)/x) for x in (1,3,10)]
    return dict(p=p,C=C,outer=outer,nodes=n,evaluations=sol.nfev,
                mass=float(mass[-1]),central_depth=float(wd[0]),vc2=vcs,
                edge_density=float(rho[-1]),edge_capture_fraction=float(rho[-1]/C),
                diagnostics=[float(mass[-1]),float(wd[0])]+vcs)

results=[];checks=[]
for p in (4,6):
    for C in (.01,1.,100.):
        for outer in (30.,100.,300.,1000.):
            coarse=run(p,C,outer,1000,1e-7,1e-10)
            fine=run(p,C,outer,2000,2e-9,1e-12)
            change=float(np.max(np.abs(np.array(coarse['diagnostics'])/fine['diagnostics']-1)))
            checks.append(dict(p=p,C=C,outer=outer,max_relative_change=change,passes=change<.005))
            results.append(dict(coarse=coarse,refined=fine))
            print(p,C,outer,fine['mass'],change,flush=True)

analytic=[]
for outer in (30.,1000.):
    r=grid(2000,outer);rho=3/(4*np.pi)*(1+r*r)**-2.5
    m,w=gravity(r,rho)
    exactm=r**3/(1+r*r)**1.5
    exactw=1/np.sqrt(1+r*r)-(1+outer*outer)**-1.5
    error=max(float(abs(m[-1]/exactm[-1]-1)),float(np.max(abs(w/exactw-1))))
    analytic.append(dict(outer=outer,max_relative_error=error,passes=error<.005))

boundary=[]
for p in (4,6):
    for C in (.01,1.,100.):
        pair=[next(x['refined'] for x in results if x['refined']['p']==p and x['refined']['C']==C and x['refined']['outer']==R) for R in (300.,1000.)]
        ratios=np.array(pair[1]['diagnostics'])/pair[0]['diagnostics']
        boundary.append(dict(p=p,C=C,ratios_1000_over_300=ratios.tolist(),
                             finite_range_convergence=bool(np.max(abs(ratios-1))<.01)))
out=dict(results=results,refinement_checks=checks,analytic_checks=analytic,boundary_checks=boundary,
         status='Synthetic outer-law sensitivity; no observations fitted or validated.')
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
assert all(x['passes'] for x in checks+analytic),'Retained failed numerical gate; refine before interpreting'
