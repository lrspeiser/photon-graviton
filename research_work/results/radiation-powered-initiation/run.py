from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp,quad

HERE=Path(__file__).resolve().parent
K=2.;sigma=.5
cases=[('no_initiation',0.,10.),('initiation_no_capture',.01,0.),('initiation_capture',.01,10.),('stronger_initiation_capture',.05,10.)]
runs=[]
for cells in [80,160,320]:
    dx=10/cells;x=(np.arange(cells)+.5)*dx
    for label,epsilon,strength in cases:
        gamma=strength*(np.exp(-x*x)+np.exp(-(x-10)**2))
        initial=np.r_[np.full(cells,.22),np.zeros(3*cells),np.zeros(5)]
        def rhs(t,y):
            R,T,M,D=y[:-5].reshape(4,cells);n=1+M/K;b=R/(K*n);d=1-b;Q=T-M;I=epsilon*R
            assert n.min()>0 and d.min()>0 and (n*d-sigma).min()>0
            facev=np.r_[2/(n[:-1]+n[1:]),1/n[-1]]
            def adv(u,inflow,factor=1.):
                slope=np.zeros(cells);a=u[1:-1]-u[:-2];bb=u[2:]-u[1:-1]
                slope[1:-1]=np.where(a*bb>0,np.sign(a)*np.minimum(abs(a),abs(bb)),0.)
                flux=np.r_[inflow,factor*facev*(u+.5*slope)]
                return -np.diff(flux)/dx,flux[-1]
            aR,oR=adv(R,.22);aT,oT=adv(T,0.);aM,oM=adv(M,0.,sigma);Mt=(aM+Q)/d
            return np.r_[aR-b*Mt-I,aT-Q-gamma*T+I,Mt,gamma*T,.22,oR+oT+oM,dx*np.sum(b*Mt),oR,dx*np.sum(I)]
        sol=solve_ivp(rhs,[0,40],initial,method='RK23',rtol=2e-7,atol=1e-10,max_step=.3*dx,dense_output=True)
        assert sol.success
        f=sol.y[:-5].reshape(4,cells,-1);n=1+f[2]/K;d=1-f[0]/(K*n)
        total=dx*f.sum(axis=(0,1))+sol.y[-4]-sol.y[-5]
        rad=dx*f[0].sum(axis=0)+sol.y[-2]+sol.y[-3]+sol.y[-1]-sol.y[-5]
        energy_error=float(max(abs(total-total[0])));rad_error=float(max(abs(rad-rad[0])))
        assert f.min()>-1e-9 and energy_error<1e-8 and rad_error<1e-8
        def state(s,t):
            assert 0<=t<=40
            y=sol.sol(t);n=1+y[:-5].reshape(4,cells)[2]/K;nt=rhs(t,y)[:-5].reshape(4,cells)[2]/K
            return np.interp(s,x,n),np.interp(s,x,nt)
        def ray(te):
            def dr(s,y):
                n,nt=state(s,y[0]);return [n,nt]
            rr=solve_ivp(dr,[0,10],[te,0.],rtol=1e-10,atol=1e-12,max_step=dx/3)
            assert rr.success;return rr.y[:,-1]
        probes=[]
        for te in [.1,5.,20.]:
            to,logJ=ray(te);ne=state(0,te)[0];no=state(10,to)[0];S=np.exp(logJ)*ne/no
            h=.001;early=ray(te-h)[0];late=ray(te+h)[0]
            din=quad(lambda t:1/state(0,t)[0],te-h,te+h,epsabs=1e-12)[0]
            dout=quad(lambda t:1/state(10,t)[0],early,late,epsabs=1e-12)[0]
            err=float(abs(dout/din/S-1));assert err<2e-4
            probes.append(dict(emission=te,arrival=float(to),coordinate_stretch=float(np.exp(logJ)),endpoint_clock_factor=float(ne/no),measured_z=float(S-1),photon_survival=float(np.exp(-epsilon*(to-te))),clock_error=err,shared_messenger_delay=0.))
        final=dx*f[:,:,-1].sum(axis=1)
        result=dict(cells=cells,case=label,epsilon=epsilon,capture_strength=strength,minimum_energy=float(f.min()),min_d=float(d.min()),min_speed_gap=float((n*d-sigma).min()),
            energy_error=energy_error,radiation_error=rad_error,final_energies=dict(zip(['R','T','M','D'],map(float,final))),
            photon_energy_in=float(sol.y[-5,-1]),total_energy_out=float(sol.y[-4,-1]),temporal_work=float(sol.y[-3,-1]),photon_energy_out=float(sol.y[-2,-1]),whole_photon_conversion=float(sol.y[-1,-1]),
            net_radiation_transfer=float(sol.y[-3,-1]+sol.y[-1,-1]),final_clock_range=[float(n[:,-1].min()),float(n[:,-1].max())],probes=probes)
        runs.append(result);print(json.dumps(result),flush=True)
refinement=[]
for label,*_ in cases:
    for a,b in [(80,160),(160,320)]:
        lo=next(r for r in runs if r['cells']==a and r['case']==label);hi=next(r for r in runs if r['cells']==b and r['case']==label)
        diff=[abs(p['measured_z']-q['measured_z']) for p,q in zip(lo['probes'],hi['probes'])]
        refinement.append(dict(case=label,coarse=a,fine=b,absolute_z_difference=diff,passes=max(diff)<.003))
(HERE/'results.json').write_text(json.dumps(dict(runs=runs,refinement=refinement,scope='Synthetic seeded-by-photon-removal branch; surviving-photon redshift and extinction distinguished; no observational fit.'),indent=2)+'\n',encoding='utf8',newline='\n')
