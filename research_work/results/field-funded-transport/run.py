from pathlib import Path
import json,sys
import numpy as np
from scipy.integrate import solve_ivp,quad

HERE=Path(__file__).resolve().parent
refine='--refine' in sys.argv
runs=json.loads((HERE/'results.json').read_text())['runs'] if refine else []
runs=[r for r in runs if r['cells']!=640]
K=2.;sigma=.5
for cells in ([640] if refine else [80,160,320]):
    dx=10/cells;x=(np.arange(cells)+.5)*dx
    for label,capture,seeded in [('no_capture',False,True),('capture',True,True),('zero_seed',True,False)]:
        if refine and not seeded:continue
        gamma=1/(1+np.exp(-(x-8)/.25)) if capture else np.zeros(cells)
        seed=.1*np.exp(-((x-2)/.4)**2) if seeded else np.zeros(cells)
        initial=np.r_[np.full(cells,.22),seed,np.zeros(2*cells),0.,0.,0.,0.]
        def fields_rhs(t,y):
            R,T,M,D=y[:-4].reshape(4,cells);n=1+M/K;b=R/(K*n);d=1-b;Q=T-M
            assert n.min()>0 and d.min()>0 and (n*d-sigma).min()>0
            facev=np.r_[2/(n[:-1]+n[1:]),1/n[-1]]
            def adv(u,inflow,factor=1.):
                slope=np.zeros(cells);a=u[1:-1]-u[:-2];b0=u[2:]-u[1:-1]
                slope[1:-1]=np.where(a*b0>0,np.sign(a)*np.minimum(abs(a),abs(b0)),0.)
                flux=np.r_[inflow,factor*facev*(u+.5*slope)]
                return -np.diff(flux)/dx,flux[-1]
            incoming=.22*(1-np.tanh((t-2)/.1))/2
            aR,oR=adv(R,incoming);aT,oT=adv(T,0.);aM,oM=adv(M,0.,sigma)
            Mt=(aM+Q)/d
            return np.r_[aR-b*Mt,aT-Q-gamma*T,Mt,gamma*T,incoming,oR+oT+oM,dx*np.sum(b*Mt),oR]
        sol=solve_ivp(fields_rhs,[0,24],initial,method='RK23',rtol=2e-7,atol=1e-10,max_step=.3*dx,dense_output=True)
        assert sol.success
        f=sol.y[:-4].reshape(4,cells,-1);nall=1+f[2]/K;dall=1-f[0]/(K*nall)
        minimum=float(f.min());min_d=float(dall.min());min_gap=float((nall*dall-sigma).min())
        balance=dx*f.sum(axis=(0,1))+sol.y[-3]-sol.y[-4]
        radiation_balance=dx*f[0].sum(axis=0)+sol.y[-1]+sol.y[-2]-sol.y[-4]
        radiation_error=float(max(abs(radiation_balance-radiation_balance[0])))
        energy_error=float(max(abs(balance-balance[0])))
        assert minimum>-1e-9 and energy_error<1e-8 and radiation_error<1e-8
        def state(s,t):
            assert 0<=t<=24
            y=sol.sol(t);n=1+y[:-4].reshape(4,cells)[2]/K
            nt=fields_rhs(t,y)[:-4].reshape(4,cells)[2]/K
            return np.interp(s,x,n),np.interp(s,x,nt)
        def ray(a,b,te):
            def dr(s,y):
                n,nt=state(s,y[0]);return [n,nt]
            rr=solve_ivp(dr,[a,b],[te,0.],rtol=1e-10,atol=1e-12,max_step=dx/3)
            assert rr.success;return rr.y[:,-1]
        probes=[]
        for a,b,te in [(0.,10.,1.),(2.,9.,.1),(0.,10.,5.)]:
            to,logJ=ray(a,b,te);ne=state(a,te)[0];no=state(b,to)[0];S=np.exp(logJ)*ne/no
            eps=.001;early=ray(a,b,te-eps)[0];late=ray(a,b,te+eps)[0]
            din=quad(lambda t:1/state(a,t)[0],te-eps,te+eps,epsabs=1e-12)[0]
            dout=quad(lambda t:1/state(b,t)[0],early,late,epsabs=1e-12)[0]
            error=float(abs(dout/din/S-1));assert error<2e-4
            probes.append(dict(start=a,end=b,emission=te,arrival=float(to),coordinate_stretch=float(np.exp(logJ)),endpoint_factor=float(ne/no),measured_z=float(S-1),clock_error=error,shared_messenger_delay=0.))
        result=dict(cells=cells,case=label,K=K,sigma=sigma,minimum_energy=minimum,min_d=min_d,min_speed_gap=min_gap,
            energy_error=energy_error,final_energies=dict(zip(['R','T','M','D'],map(float,dx*f[:,:,-1].sum(axis=1)))),
            energy_in=float(sol.y[-4,-1]),energy_out=float(sol.y[-3,-1]),
            net_radiation_to_M=float(sol.y[-2,-1]),radiation_out=float(sol.y[-1,-1]),radiation_ledger_error=radiation_error,probes=probes)
        runs.append(result);print(json.dumps(result),flush=True)
refinement=[]
for label in ['no_capture','capture','zero_seed']:
    pairs=[(80,160),(160,320)]+([(320,640)] if any(r['cells']==640 and r['case']==label for r in runs) else [])
    for a,b in pairs:
        lo=next(r for r in runs if r['cells']==a and r['case']==label)
        hi=next(r for r in runs if r['cells']==b and r['case']==label)
        diff=[abs(p['measured_z']-q['measured_z']) for p,q in zip(lo['probes'],hi['probes'])]
        refinement.append(dict(case=label,coarse=a,fine=b,absolute_z_difference=diff,passes=max(diff)<.003))
result=dict(runs=runs,refinement=refinement,empty_T_boundary_derivative=.1,
    empty_T_boundary_explanation='At T=0 and T_x=0, T_t=-Q=M/tau, independent of M_x. For prior example M=0.1,tau=1 it is +0.1.',
    scope='Finite-domain synthetic reference-energy and clock test; no global domain proof or observational fit.')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
