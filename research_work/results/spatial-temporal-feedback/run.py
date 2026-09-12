from pathlib import Path
import json,sys
import numpy as np
from scipy.integrate import solve_ivp,quad
HERE=Path(__file__).resolve().parent
results=[]
first_order='--first-order' in sys.argv
for cells in [80,160,320,640]:
    dx=10/cells;x=(np.arange(cells)+.5)*dx
    for case in ['zero_seed','travel','capture']:
        seed=np.zeros(cells) if case=='zero_seed' else .01*np.exp(-((x-2)/.4)**2)
        gamma=np.exp(-((x-8)/.5)**2) if case=='capture' else np.zeros(cells)
        initial=np.concatenate([np.ones(cells),np.full(cells,.2),np.full(cells,.02),seed,np.zeros(cells),[0.,0.]])
        def rhs(t,y):
            n,G,H,T,dep=y[:-2].reshape(5,cells);v=1/n;h=T/n
            def transport(u,inflow):
                if first_order: right=v*u
                else:
                    slope=np.zeros(cells);left=u[1:-1]-u[:-2];rightdiff=u[2:]-u[1:-1]
                    slope[1:-1]=np.where(left*rightdiff>0,np.sign(left)*np.minimum(abs(left),abs(rightdiff)),0.)
                    facev=np.r_[2/(n[:-1]+n[1:]),v[-1]]
                    right=facev*(u+.5*slope)
                flux=np.r_[inflow,right];return -np.diff(flux)/dx,flux[-1]
            dG,oG=transport(G,.2);dH,oH=transport(H,.02);dT,oT=transport(T,0.)
            return np.r_[T,dG-h*G,dH-h*H,dT+h*(G+H)-gamma*T,gamma*T,.22,oG+oH+oT]
        sol=solve_ivp(rhs,[0,24],initial,method='RK23',rtol=2e-7,atol=1e-10,max_step=.35*dx,dense_output=True)
        assert sol.success
        sample=sol.sol(np.linspace(0,24,241));fields=sample[:-2].reshape(5,cells,-1)
        ledger=dx*fields[1:].sum(axis=(0,1))+sample[-1]-sample[-2]
        err=float(np.max(abs(ledger-ledger[0])));assert err<1e-8 and fields[1:].min()>-1e-10
        def state(s,t):
            z=sol.sol(t)[:-2].reshape(5,cells)
            return np.interp(s,x,z[0]),np.interp(s,x,z[3])
        def ray(a,b,te):
            def flow(s,y):
                n,nt=state(s,y[0]);return [n,nt]
            answer=solve_ivp(flow,[a,b],[te,0.],rtol=1e-11,atol=1e-13,max_step=dx/4)
            assert answer.success and answer.y[0,-1]<24
            return answer.y[:,-1]
        probes=[]
        for a,b in [(0.,10.),(2.,9.)]:
            te=1.;to,logJ=ray(a,b,te);ne,_=state(a,te);no,_=state(b,to);S=np.exp(logJ)*ne/no
            eps=.001;early=ray(a,b,te-eps)[0];late=ray(a,b,te+eps)[0]
            din=quad(lambda t:1/state(a,t)[0],te-eps,te+eps,epsabs=1e-12)[0]
            dout=quad(lambda t:1/state(b,t)[0],early,late,epsabs=1e-12)[0]
            rel=float(abs(dout/din/S-1));assert rel<2e-4
            probes.append(dict(start=a,end=b,arrival=float(to),coordinate_stretch=float(np.exp(logJ)),endpoint_clock_factor=float(ne/no),measured_redshift=float(S-1),proper_interval_stretch=float(dout/din),clock_check_relative_error=rel,shared_messenger_delay=0.))
        final=fields[:,:,-1]
        results.append(dict(cells=cells,case=case,max_reference_energy_error=err,min_energy=float(fields[1:].min()),max_n=float(fields[0].max()),final_deposit=float(dx*sum(final[4])),final_receiving_energy=float(dx*sum(final[3])),initial_receiving_energy=float(dx*sum(seed)),total_inflow=float(sample[-2,-1]),total_outflow=float(sample[-1,-1]),probes=probes))
comparisons=[]
for case in ['zero_seed','travel','capture']:
    for low,high in [(80,160),(160,320),(320,640)]:
        coarse=next(r for r in results if r['case']==case and r['cells']==low);fine=next(r for r in results if r['case']==case and r['cells']==high)
        differences=[abs(a['measured_redshift']-b['measured_redshift']) for a,b in zip(coarse['probes'],fine['probes'])]
        comparisons.append(dict(case=case,coarse_cells=low,fine_cells=high,absolute_redshift_refinement_differences=differences,passes_initial_gate=max(differences)<.003))
result=dict(runs=results,refinement=comparisons,transport='first-order upwind' if first_order else 'minmod MUSCL',scope='Source-driven prescribed-response spatial kinematics, reference-energy ledger only')
(HERE/('first-order-results.json' if first_order else 'results.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(result,indent=2))
