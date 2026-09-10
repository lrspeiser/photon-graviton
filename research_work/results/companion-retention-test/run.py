"""Autonomous companion-wave energy exchange in the proposed equal-speed field."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import solve_ivp,simpson,cumulative_trapezoid
HERE=Path(__file__).resolve().parent

def simulate(rolling,amplitude,cells,length=16.):
    L=length;dx=L/cells;x=np.arange(cells)*dx;K=1.;c0=1.
    k=2*np.pi*np.fft.fftfreq(cells,d=dx);k[cells//2]=0
    def D(f):return np.fft.ifft(1j*k*np.fft.fft(f)).real
    def acceleration(n):return c0*c0*D(D(np.log(n)))/n
    delta=(x-L/4+L/2)%L-L/2
    wave=amplitude*np.exp(-.5*(delta/.5)**2)*np.cos(2*np.pi*delta)
    wave-=wave.mean()
    y0=np.r_[1+wave,rolling-c0*D(wave)+rolling*wave/2]
    def rhs(t,y):
        n,p=y[:cells],y[cells:]
        if min(n)<=0:raise ValueError('Nonpositive field')
        return np.r_[p/K,K*acceleration(n)]
    times=np.linspace(0,6,301)
    sol=solve_ivp(rhs,[0,6],y0,t_eval=times,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.25*dx)
    assert sol.success
    homogeneous=[];traveling=[];total=[];means=[];centroid=[];work=[]
    for y in sol.y.T:
        n,p=y[:cells],y[cells:];pm=p.mean()
        density=(p-pm)**2/(2*K)+K*c0*c0*D(np.log(n))**2/2
        eh=L*pm*pm/(2*K);ew=dx*density.sum()
        homogeneous.append(eh);traveling.append(ew);total.append(eh+ew);means.append(n.mean())
        z=np.sum(density*np.exp(2j*np.pi*x/L));centroid.append(np.angle(z))
        work.append(L*pm*acceleration(n).mean())
    homogeneous,traveling,total,means=map(np.asarray,(homogeneous,traveling,total,means))
    displacement=(np.unwrap(centroid)-centroid[0])*L/(2*np.pi)
    geometric=cumulative_trapezoid(c0/means,times,initial=0)
    relative_total=float(max(abs(total-total[0]))/total[0])
    # Report conservation relative to wave energy too, so a large background cannot hide error.
    relative_wave=float(max(abs(total-total[0]))/traveling[0])
    closure=float(abs((traveling[0]-traveling[-1])-(homogeneous[-1]-homogeneous[0]))/traveling[0])
    transfer=float(simpson(work,x=times))
    assert relative_wave<1e-6 and closure<1e-6 and sol.y[:cells].min()>0
    result={'rolling_rate':rolling,'initial_wave_amplitude':amplitude,'cells':cells,'box_length':L,
       'initial_wave_energy':float(traveling[0]),'final_wave_energy':float(traveling[-1]),
       'wave_energy_retained_fraction':float(traveling[-1]/traveling[0]),
       'background_energy_gain':float(homogeneous[-1]-homogeneous[0]),'initial_background_energy':float(homogeneous[0]),
       'wave_loss_background_gain_residual_over_initial_wave':closure,'energy_error_over_initial_total':relative_total,
       'energy_error_over_initial_wave':relative_wave,'integrated_background_work':transfer,
       'work_integral_error_over_initial_wave':float(abs(transfer-(homogeneous[-1]-homogeneous[0]))/traveling[0]),
       'mean_field_final':float(means[-1]),'homogeneous_adiabatic_retention_prediction':float(means[0]/means[-1]),
       'wave_energy_centroid_displacement':float(displacement[-1]),'integral_shared_speed_displacement':float(geometric[-1]),
       'minimum_field':float(sol.y[:cells].min()),
       'history':{'time':times[::10].tolist(),'wave_energy_fraction':(traveling[::10]/traveling[0]).tolist(),'mean_field':means[::10].tolist(),'centroid_displacement':displacement[::10].tolist()}}
    assert result['work_integral_error_over_initial_wave']<1e-5
    return result

def main():
    runs=[];refinements=[]
    for rolling,amp in [(0.,.001),(.05,.001),(.1,.001),(.05,.003)]:
        pair=[]
        for cells in [256,512]:
            r=simulate(rolling,amp,cells);runs.append(r);pair.append(r)
            print(rolling,amp,cells,r['wave_energy_retained_fraction'],r['background_energy_gain'],flush=True)
        changes={key:abs(pair[1][key]-pair[0][key]) for key in ['wave_energy_retained_fraction','wave_energy_centroid_displacement','mean_field_final']}
        assert max(changes.values())<1e-5
        refinements.append({'rolling_rate':rolling,'amplitude':amp,'absolute_changes':changes})
    larger=simulate(.05,.001,512,length=32.)
    reference=next(x for x in runs if x['rolling_rate']==.05 and x['initial_wave_amplitude']==.001 and x['cells']==256)
    box_change=abs(larger['wave_energy_retained_fraction']-reference['wave_energy_retained_fraction'])
    assert box_change<1e-4
    out={'scope':'Autonomous nonlinear field, no photons or external driving during evolution. Wave/background energy defined in one reference Hamiltonian. No capture or astronomical inference.',
       'equations':'H=integral [Pi^2/(2K)+K*c0^2/2*(partial_x ln n)^2] dx; n_t=Pi/K; Pi_t=K*c0^2/n*partial_xx ln n.',
       'runs':runs,'refinements':refinements,'double_box_control':larger,'double_box_retention_change':box_change,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()
