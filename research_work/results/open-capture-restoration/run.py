from pathlib import Path
import json,sys
import numpy as np
from scipy.integrate import solve_ivp,quad
HERE=Path(__file__).resolve().parent
refine_only='--refine-low' in sys.argv
runs=json.loads((HERE/'results.json').read_text())['runs'] if refine_only else []
runs=[r for r in runs if r['cells']!=1280]
cases=[('zero_seed',.05,True,False),('low_cost_no_reset',.05,False,True),('low_cost_reset',.05,True,True),('high_cost_no_reset',.5,False,True),('high_cost_reset',.5,True,True)]
for cells in ([1280] if refine_only else [160,320,640]):
    dx=10/cells;x=(np.arange(cells)+.5)*dx;gamma=1/(1+np.exp(-(x-8)/.25))
    for label,K,restore,seeded in cases:
        if refine_only and not label.startswith('low_cost'):continue
        seed=.01*np.exp(-((x-2)/.4)**2) if seeded else np.zeros(cells)
        initial=np.r_[np.ones(cells),np.full(cells,.2),np.full(cells,.02),seed,np.zeros(cells),0.,0.]
        def response(n,R,T):
            r=gamma*(n-1)*np.maximum(0,1-R/(K*n)) if restore else np.zeros(cells)
            return T-r,r
        def rhs(t,y):
            n,G,H,T,D=y[:-2].reshape(5,cells);R=G+H;ndot,r=response(n,R,T);h=ndot/n
            facev=np.r_[2/(n[:-1]+n[1:]),1/n[-1]]
            def adv(u,inflow):
                slope=np.zeros(cells);a=u[1:-1]-u[:-2];b=u[2:]-u[1:-1]
                slope[1:-1]=np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.)
                flux=np.r_[inflow,facev*(u+.5*slope)];return -np.diff(flux)/dx,flux[-1]
            source=(1-np.tanh((t-2)/.1))/2
            dG,oG=adv(G,.2*source);dH,oH=adv(H,.02*source);dT,oT=adv(T,0.)
            return np.r_[ndot,dG-h*G,dH-h*H,dT+(R/n-K)*T-gamma*T,gamma*T+(K-R/n)*r,.22*source,oG+oH+oT]
        sol=solve_ivp(rhs,[0,24],initial,method='RK23',rtol=2e-7,atol=1e-10,max_step=.35*dx,dense_output=True);assert sol.success
        sample=sol.sol(np.linspace(0,24,241));f=sample[:-2].reshape(5,cells,-1);memory=K*(f[0]-1)
        balance=dx*(f[1:].sum(axis=(0,1))+memory.sum(axis=0))+sample[-1]-sample[-2]
        error=float(max(abs(balance-balance[0])));assert error<1e-8 and min(f[1:].min(),memory.min())>-1e-10
        def state(s,t):
            n,G,H,T,D=sol.sol(t)[:-2].reshape(5,cells);nt,_=response(n,G+H,T)
            return np.interp(s,x,n),np.interp(s,x,nt)
        def ray(a,b,te):
            def rhsray(s,z):
                n,nt=state(s,z[0]);return [n,nt]
            rr=solve_ivp(rhsray,[a,b],[te,0.],rtol=1e-11,atol=1e-13,max_step=dx/4);assert rr.success and rr.y[0,-1]<24;return rr.y[:,-1]
        probes=[]
        for a,b,te in [(0.,10.,1.),(2.,9.,.1),(0.,10.,5.)]:
            to,logJ=ray(a,b,te);ne,_=state(a,te);no,_=state(b,to);S=np.exp(logJ)*ne/no
            eps=.001;early=ray(a,b,te-eps)[0];late=ray(a,b,te+eps)[0]
            din=quad(lambda t:1/state(a,t)[0],te-eps,te+eps,epsabs=1e-12)[0];dout=quad(lambda t:1/state(b,t)[0],early,late,epsabs=1e-12)[0]
            cerr=abs(dout/din/S-1);assert cerr<2e-4
            probes.append(dict(start=a,end=b,emission=te,arrival=float(to),coordinate_stretch=float(np.exp(logJ)),endpoint_factor=float(ne/no),measured_z=float(S-1),clock_error=float(cerr),shared_messenger_delay=0.))
        final=f[:,:,-1]
        runs.append(dict(cells=cells,case=label,K=K,restoration=restore,energy_error=error,final_memory=float(dx*memory[:,-1].sum()),final_deposit=float(dx*final[4].sum()),final_receiver=float(dx*final[3].sum()),final_radiation=float(dx*(final[1]+final[2]).sum()),energy_out=float(sample[-1,-1]),energy_in=float(sample[-2,-1]),final_receiver_clock_n=float(final[0,-1]),probes=probes))
ref=[]
for label,*_ in cases:
    for a,b in [(160,320),(320,640)]:
        lo=next(r for r in runs if r['cells']==a and r['case']==label);hi=next(r for r in runs if r['cells']==b and r['case']==label)
        diff=[abs(p['measured_z']-q['measured_z']) for p,q in zip(lo['probes'],hi['probes'])]
        ref.append(dict(case=label,coarse=a,fine=b,absolute_z_difference=diff,passes=max(diff)<.003))
    if any(r['case']==label and r['cells']==1280 for r in runs):
        lo=next(r for r in runs if r['cells']==640 and r['case']==label);hi=next(r for r in runs if r['cells']==1280 and r['case']==label)
        diff=[abs(p['measured_z']-q['measured_z']) for p,q in zip(lo['probes'],hi['probes'])]
        ref.append(dict(case=label,coarse=640,fine=1280,absolute_z_difference=diff,passes=max(diff)<.003))
result=dict(runs=runs,refinement=ref,scope='Open reference-energy/memory closure; no full momentum or proper-energy completion')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(result,indent=2))
