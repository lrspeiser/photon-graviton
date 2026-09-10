"""Finite signal packets and their shared receiver; no astronomical calibration."""
from pathlib import Path
import json,hashlib,importlib.util
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
prior=HERE.parent/'weak-signal-timing/check.py'
spec=importlib.util.spec_from_file_location('weak_probe',prior);weak=importlib.util.module_from_spec(spec);spec.loader.exec_module(weak)

def simulate(energy,modes=32,color=1.,rolling=.05,reference_field=None,return_field=False,driver_energy=.01):
    L=8.;K=1.;v=.5;mass0=K*L;mass=mass0/2
    k=2*np.pi*np.arange(1,modes+1)/L;window=np.exp(-.5*(.2*k)**2)
    end=4.5;source=.4;detector=2.4;off=2+4*modes
    frequencies=np.array([1.,color,color,color]);weights=np.array([driver_energy,energy/3/color,energy/3/color,energy/3/color])
    def optical(y,x):
        x=np.atleast_1d(x);co=np.cos(x[:,None]*k);si=np.sin(x[:,None]*k)
        qc=y[2:2+modes];qs=y[2+modes:2+2*modes];pc=y[2+2*modes:2+3*modes];ps=y[2+3*modes:off]
        n=y[0]+(co*window)@qc+(si*window)@qs
        nx=(-si*window*k)@qc+(co*window*k)@qs
        nt=y[1]/mass0+((co*window)@pc+(si*window)@ps)/mass
        return n,nx,nt,co,si
    def rhs(t,y):
        x=y[off:off+4];p=y[off+4:off+8]
        n,nx,_,co,si=optical(y,x)
        if min(n)<=0:raise ValueError('Nonpositive optical factor')
        loading=weights*p/n**2; spring=mass*v*v*k*k
        clock_n=optical(y,[source,detector])[0]
        return np.r_[y[1]/mass0,loading.sum(),y[2+2*modes:off]/mass,
           -spring*y[2:2+modes]+window*(loading@co),
           -spring*y[2+modes:2+2*modes]+window*(loading@si),1/n,p*nx/n**2,1/clock_n]
    initial=np.zeros(off+10);initial[0]=1;initial[1]=mass0*rolling
    initial[off:off+4]=[.4,.35,.15,-.05];initial[off+4:off+8]=frequencies
    sol=solve_ivp(rhs,[0,end],initial,dense_output=True,method='DOP853',rtol=2e-11,atol=1e-13,max_step=.02)
    assert sol.success
    def field(t,x):
        if not 0<=t<=end:raise ValueError('No extrapolation')
        a=optical(sol.sol(t),[x]);return tuple(float(z[0]) for z in a[:3])
    def ledger(y):
        qc=y[2:2+modes];qs=y[2+modes:2+2*modes];pc=y[2+2*modes:2+3*modes];ps=y[2+3*modes:off]
        zero=y[1]**2/(2*mass0)
        wave=np.sum((pc*pc+ps*ps)/(2*mass)+mass*v*v*k*k*(qc*qc+qs*qs)/2)
        photons=weights*y[off+4:off+8]/optical(y,y[off:off+4])[0]
        momentum=float(weights@y[off+4:off+8]+np.sum(k*(qc*ps-qs*pc)))
        return zero,wave,photons,momentum
    led=[ledger(sol.sol(t)) for t in np.linspace(0,end,181)]
    lower_bounds=[]
    for t in np.linspace(0,end,181):
        y=sol.sol(t)
        lower_bounds.append(y[0]-np.sum(window*np.hypot(y[2:2+modes],y[2+modes:2+2*modes])))
    assert min(lower_bounds)>0
    total=np.array([z+w+p.sum() for z,w,p,_ in led]);mom=np.array([j for _,_,_,j in led])
    ee=float(max(abs(total-total[0]))/total[0]);me=float(max(abs(mom-mom[0]))/max(mom[0],1e-12))
    assert ee<1e-8 and me<1e-8
    rows=[]
    for i in range(1,4):
        te=brentq(lambda t:sol.sol(t)[off+i]-source,0,end,xtol=1e-13)
        to=brentq(lambda t:sol.sol(t)[off+i]-detector,te,end,xtol=1e-13)
        ye,yo=sol.sol(te),sol.sol(to);ne=field(te,source)[0];no=field(to,detector)[0]
        freq_e=ye[off+4+i]/ne;freq_o=yo[off+4+i]/no
        check=weak.probe(field,source,detector-source,te,freq_e,end)
        assert abs(check['arrival_time']-to)<2e-8
        assert abs(check['frequency_stretch']-freq_e/freq_o)<2e-8
        rows.append({'packet':i,'source_crossing_time':float(te),'detector_crossing_time':float(to),'reference_frequency_stretch':float(freq_e/freq_o),
          'source_q1_clock':float(ye[-2]),'detector_q1_clock':float(yo[-1]),
          'q1_frequency_stretch':float(freq_e/freq_o*ne/no),
          'independent_frozen_field_jacobian':check['local_event_stretch'],
          'source_frequency':float(freq_e),'detector_frequency':float(freq_o),
          'initial_packet_energy':float(weights[i]*frequencies[i])})
        if reference_field is not None:
            matched=weak.probe(reference_field,source,detector-source,te,freq_e,end)
            rows[-1]['weak_background_stretch_at_same_source_time']=matched['frequency_stretch']
            rows[-1]['feedback_stretch_increment']=float(freq_e/freq_o-matched['frequency_stretch'])
    intervals=[]
    for a,b in zip(rows[:-1],rows[1:]):
        sr=(b['detector_crossing_time']-a['detector_crossing_time'])/(b['source_crossing_time']-a['source_crossing_time'])
        sq=(b['detector_q1_clock']-a['detector_q1_clock'])/(b['source_q1_clock']-a['source_q1_clock'])
        intervals.append({'reference_duration_stretch':sr,'q1_duration_stretch':sq,'relative_to_first_carrier_stretch':sr/a['reference_frequency_stretch']-1})
    z,w,p,_=led[-1];z0,w0,p0,_=led[0]
    result={'signal_total_initial_energy':energy,'driver_initial_energy':driver_energy,'modes':modes,'signal_carrier_multiplier':color,'initial_rolling_rate':rolling,
      'sampled_global_optical_lower_bound':float(min(lower_bounds)),
      'packets':rows,'intervals':intervals,'maximum_relative_energy_error':ee,'maximum_relative_momentum_error':me,
      'initial_field_energy':float(z0+w0),'final_field_energy':float(z+w),'field_energy_gain':float(z+w-z0-w0),
      'driver_energy_lost':float(p0[0]-p[0]),'signal_energy_lost':float(p0[1:].sum()-p[1:].sum()),
      'photon_field_exchange_residual':float((p0-p).sum()-(z+w-z0-w0))}
    return (result,field) if return_field else result

def main():
    runs=[];reference_field=None
    for energy in [0.,1e-5,.001,.01,.1]:
        for modes in [32,64]:
            if energy==0 and modes==64:
                r,reference_field=simulate(energy,modes,return_field=True)
            else:r=simulate(energy,modes,reference_field=reference_field)
            runs.append(r)
            print('E',energy,'modes',modes,'stretch',[x['reference_frequency_stretch'] for x in r['packets']],flush=True)
    color=simulate(.01,64,3.)
    analytic=simulate(0.,64,driver_energy=0.)
    for packet in analytic['packets']:
        assert abs(packet['reference_frequency_stretch']-np.exp(.1))<1e-8
        assert abs(packet['q1_frequency_stretch']-1)<1e-8
    for interval in analytic['intervals']:
        assert abs(interval['reference_duration_stretch']-np.exp(.1))<1e-8
        assert abs(interval['q1_duration_stretch']-1)<1e-8
    reference=next(r for r in runs if r['signal_total_initial_energy']==.01 and r['modes']==64)
    color_error=max(abs(a['reference_frequency_stretch']-b['reference_frequency_stretch']) for a,b in zip(color['packets'],reference['packets']))
    assert color_error<1e-8
    changes=[]
    for a,b in zip(runs[::2],runs[1::2]):
        delta=max(abs(x['reference_frequency_stretch']-y['reference_frequency_stretch']) for x,y in zip(a['packets'],b['packets']))
        assert delta<1e-6
        changes.append({'energy':a['signal_total_initial_energy'],'max_stretch_change':delta})
    result={'scope':'Finite signal feedback in the existing periodic Hamiltonian. Reference and passive q=1/n clock controls; no self-consistent matter or cosmic fit.',
      'runs':runs,'color_control':color,'analytic_homogeneous_control':analytic,'fixed_energy_color_stretch_difference':color_error,'mode_refinement':changes,
      'known_equations':'H=field kinetic+gradient energy+sum_j N_j k_j/n(X_j). Fixed photon weights N_j distinguish packet loading from carrier frequency.',
      'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'run.py',prior]}}
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()
