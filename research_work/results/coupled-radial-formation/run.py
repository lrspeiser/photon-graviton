from pathlib import Path
import sys,json
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'depleted-capture-feedback'))
from run import Model

def simulate(Cdot,n,nq,dt,interval,epsilon):
    radiation=Model(30.,n,nq);rgrid=radiation.r
    x=np.array([]);v=np.array([]);mass=np.array([])
    absorbed_sum=transmitted_sum=incident_sum=rest_sum=bulk_sum=0.
    kinetic_insertion=potential_insertion=0.
    max_residual=max_beta=max_speed=max_opacity_rate=0.
    old_k=None;outputs=[];diagnostic_flags={}
    def mechanics():
        if not len(x):return np.array([]),0.,0.
        r=abs(x);order=np.argsort(r);mr=mass[order];rr=r[order]
        enclosed=np.cumsum(mr)-.5*mr
        adep=-np.sign(x[order])*enclosed*rr/(rr*rr+epsilon*epsilon)**1.5
        acceleration=np.empty(len(x));acceleration[order]=adep
        acceleration-=x/(1+x*x)**1.5
        K=float(mass@(.5*v*v))
        U=-float(np.sum(mr*enclosed/np.sqrt(rr*rr+epsilon*epsilon)))-float(mass@(1/np.sqrt(1+x*x)))
        return acceleration,K,U
    def depth():
        if not len(x):return np.zeros(n)
        order=np.argsort(abs(x));rr=abs(x[order]);mm=mass[order]
        cumulative=np.r_[0.,np.cumsum(mm)]
        outer=np.r_[np.cumsum((mm/np.sqrt(rr*rr+epsilon*epsilon))[::-1])[::-1],0.]
        idx=np.searchsorted(rr,rgrid,side='right')
        return cumulative[idx]/np.sqrt(rgrid*rgrid+epsilon*epsilon)+outer[idx]
    stride=round(interval/dt);steps=round(3/dt)
    accel=np.array([])
    for step in range(steps):
        if step%stride==0:
            _,Kbefore,Ubefore=mechanics()
            W=radiation.wb+depth();k=.1*(W/(1+W))**6
            if old_k is not None:
                max_opacity_rate=max(max_opacity_rate,float(max(abs(np.log(k/old_k))))/interval*(60/1000))
            old_k=k
            shell,total,out=radiation.absorption(k)
            energy_equiv=interval*(Cdot/.1)*shell
            q=shell/radiation.V
            cumulative=np.cumsum(shell)-shell+shell*(radiation.r**3-radiation.a**3)/(radiation.b**3-radiation.a**3)
            beta=k*cumulative/(4*np.pi*rgrid*rgrid*q)
            assert np.all(beta>=0) and np.all(beta<1)
            rest=energy_equiv*np.sqrt(1-beta*beta)
            new_v=-1000*beta
            max_beta=max(max_beta,float(max(beta)))
            if max_beta>.03: diagnostic_flags.setdefault('injected_beta_above_0.03',step*dt)
            if max_opacity_rate>.1: diagnostic_flags.setdefault('opacity_change_per_crossing_above_0.1',step*dt)
            x=np.r_[x,rgrid];v=np.r_[v,new_v];mass=np.r_[mass,rest]
            accel,Kafter,Uafter=mechanics()
            kinetic_insertion+=Kafter-Kbefore;potential_insertion+=Uafter-Ubefore
            absorbed_sum+=float(energy_equiv.sum());rest_sum+=float(rest.sum())
            bulk_sum+=float(np.sum(energy_equiv-rest))
            transmitted_sum+=interval*(Cdot/.1)*out
            incident_sum+=interval*(Cdot/.1)*np.pi*30**2
        v+=.5*dt*accel;x+=dt*v
        accel,K,U=mechanics();v+=.5*dt*accel
        K=float(mass@(.5*v*v));H=K+U
        residual=abs(H-kinetic_insertion-potential_insertion)/max(K+abs(U),1e-30)
        max_residual=max(max_residual,residual);max_speed=max(max_speed,float(max(abs(v)))/1000)
        if max_speed>.03: diagnostic_flags.setdefault('speed_above_0.03c',(step+1)*dt)
        if step+1 in [round(t/dt) for t in (1,2,3)]:
            fractions=[float(mass[abs(x)<r].sum()/mass.sum()) for r in (.1,1,3,10)]
            outputs.append(dict(T=(step+1)*dt,mass=float(mass.sum()),enclosed_fractions=fractions,
                                outside_capture_fraction=float(mass[abs(x)>30].sum()/mass.sum()),
                                max_speed_over_c_to_date=max_speed,max_opacity_log_change_per_crossing_to_date=max_opacity_rate,max_mechanical_residual_to_date=max_residual,
                                K=K,U=U,H=H,kinetic_insertion=kinetic_insertion,gravitational_insertion_work=potential_insertion,
                                central_deposit_depth=float(np.sum(mass/np.sqrt(x*x+epsilon*epsilon)))))
    ledgers=dict(radiation_relative_error=abs(absorbed_sum+transmitted_sum-incident_sum)/incident_sum,
                 local_conversion_relative_error=abs(rest_sum+bulk_sum-absorbed_sum)/absorbed_sum,
                 retained_mass_relative_error=abs(mass.sum()-rest_sum)/rest_sum,
                 max_mechanical_residual=max_residual)
    assert max(ledgers[k] for k in ('radiation_relative_error','local_conversion_relative_error','retained_mass_relative_error'))<1e-10
    result=dict(Cdot=Cdot,n=n,impact_nodes=nq,dt=dt,injection_interval=interval,epsilon=epsilon,
                cohort_count=len(x),outputs=outputs,ledgers=ledgers,diagnostic_flags=diagnostic_flags,max_injected_beta=max_beta,max_speed_over_c=max_speed,
                max_opacity_log_change_per_crossing=max_opacity_rate,absorbed_energy_equivalent=absorbed_sum,
                rest_mass=rest_sum,bulk_energy_equivalent=bulk_sum,numerical_energy_gate=max_residual<.005)
    return result

Cdot=float(sys.argv[1]);assert Cdot in (.1,1.)
cases=[]
core='--core' in sys.argv
configs=((128,8,.001,.01,.0125),(128,8,.001,.01,.00625)) if core else ((64,4,.002,.02,.05),(128,8,.001,.01,.05),(128,8,.001,.01,.025))
for config in configs:
    result=simulate(Cdot,*config);cases.append(result)
    print(Cdot,config,result['outputs'][-1]['mass'],result['numerical_energy_gate'],result['diagnostic_flags'],flush=True)
    (HERE/f"{'core' if core else 'results'}-C{Cdot:g}.json").write_text(json.dumps(cases,indent=2)+'\n',encoding='utf8',newline='\n')
