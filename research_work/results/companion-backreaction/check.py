"""Closed spectral companion-field/ray dynamics with exact translation symmetry."""
from pathlib import Path
import json,hashlib,os
import numpy as np
from scipy.integrate import solve_ivp,simpson

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'companion-backreaction'


def simulate(modes,inertia,packet_energy,rolling,protocol,x0=None):
    cfg=protocol['parameters'];length=cfg['length'];speed=cfg['wave_speed']
    mass0=inertia*length;mass=mass0/2
    wave_number=2*np.pi*np.arange(1,modes+1)/length
    window=np.exp(-.5*(cfg['packet_smoothing_sigma']*wave_number)**2)
    # State: n0, P0, cosine coordinates, sine coordinates, cosine/sine momenta, ray X/P.
    def unpack(y):return y[0],y[1],y[2:2+modes],y[2+modes:2+2*modes],y[2+2*modes:2+3*modes],y[2+3*modes:2+4*modes],y[-2],y[-1]
    def rhs(t,y):
        n0,p0,qc,qs,pc,ps,x,p=unpack(y)
        co=np.cos(wave_number*x);si=np.sin(wave_number*x)
        n=n0+np.sum(window*(qc*co+qs*si))
        nx=np.sum(window*wave_number*(-qc*si+qs*co))
        if n<=0:raise ValueError('Outside positive propagation domain')
        source=p/n**2
        return np.r_[p0/mass0,source,pc/mass,ps/mass,
                     -mass*speed**2*wave_number**2*qc+source*window*co,
                     -mass*speed**2*wave_number**2*qs+source*window*si,
                     1/n,p*nx/n**2]
    initial=np.zeros(4*modes+4);initial[0]=cfg['initial_n0'];initial[1]=mass0*rolling
    initial[-2]=cfg['initial_x'] if x0 is None else x0;initial[-1]=packet_energy*initial[0]
    times=np.linspace(0,cfg['duration'],121)
    sol=solve_ivp(rhs,[0,cfg['duration']],initial,t_eval=times,method='DOP853',rtol=2e-11,atol=1e-13,max_step=.02)
    assert sol.success
    eh=[];ew=[];ep=[];mom=[];field_mom=[];min_n=[]
    region_energy=[];region_source=[];region_outflux=[];left_flux=[];region_cross=[]
    grid=np.linspace(0,length,256,endpoint=False)
    cg=np.cos(wave_number[:,None]*grid);sg=np.sin(wave_number[:,None]*grid)
    region=np.linspace(initial[-2]-.5,initial[-2]+.5,129)
    cr=np.cos(wave_number[:,None]*region);sr=np.sin(wave_number[:,None]*region)
    for y in sol.y.T:
        n0,p0,qc,qs,pc,ps,x,p=unpack(y)
        n=n0+np.sum(window*(qc*np.cos(wave_number*x)+qs*np.sin(wave_number*x)))
        eh.append(p0*p0/(2*mass0))
        ew.append(np.sum((pc*pc+ps*ps)/(2*mass)+mass*speed**2*wave_number**2*(qc*qc+qs*qs)/2))
        ep.append(p/n)
        pm=np.sum(wave_number*(qc*ps-qs*pc));field_mom.append(pm);mom.append(p+pm)
        min_n.append(np.min(n0+qc@cg+qs@sg))
        nt=p0/mass0+(pc@cr+ps@sr)/mass
        nx=(-wave_number*qc)@sr+(wave_number*qs)@cr
        excess_density=inertia*(rolling*(nt-rolling)+.5*((nt-rolling)**2+speed**2*nx**2))
        flux=-inertia*speed**2*nt*nx
        kernel=(1+2*np.sum(window[:,None]*np.cos(wave_number[:,None]*(region-x)),axis=0))/length
        region_energy.append(simpson(excess_density,x=region))
        region_source.append(simpson(p/n**2*kernel*nt,x=region))
        region_outflux.append(flux[-1]-flux[0]);left_flux.append(-flux[0])
        region_cross.append(simpson(inertia*rolling*(nt-rolling),x=region))
    eh,ew,ep,mom,field_mom=map(np.asarray,(eh,ew,ep,mom,field_mom))
    total=eh+ew+ep;scale=max(total[0],1e-15)
    energy_error=float(np.max(abs(total-total[0]))/scale)
    momentum_error=float(np.max(abs(mom-mom[0]))/max(packet_energy,1e-15))
    assert energy_error<protocol['checks']['maximum_relative_energy_drift']
    assert momentum_error<protocol['checks']['maximum_momentum_drift_over_initial_packet_momentum']
    assert min(min_n)>0
    # A nonzero-mode scalar wave obeys E_wave >= v*|P_wave|.
    assert np.min(ew-speed*abs(field_mom))>-1e-11
    lost=float(ep[0]-ep[-1]);gain=float(eh[-1]-eh[0]+ew[-1])
    region_work=float(simpson(region_source,x=times));region_out=float(simpson(region_outflux,x=times))
    region_error=abs(region_energy[-1]-region_energy[0]-region_work+region_out)
    assert region_error/max(packet_energy,1e-15)<protocol['checks']['maximum_local_energy_balance_error_over_initial_packet_energy']
    return {'modes':modes,'inertia':inertia,'initial_packet_energy':packet_energy,'initial_rolling_rate':rolling,
            'initial_x':float(initial[-2]),'initial_homogeneous_field_energy':float(eh[0]),
            'photon_energy_lost':lost,'photon_fraction_lost':lost/packet_energy if packet_energy else None,
            'homogeneous_field_energy_gain':float(eh[-1]-eh[0]),'nonzero_wave_energy_gain':float(ew[-1]),
            'nonzero_wave_share_of_received_energy':float(ew[-1]/gain) if gain>1e-15 else None,
            'final_packet_energy':float(ep[-1]),'final_field_momentum':float(field_mom[-1]),
            'photon_to_field_energy_balance_error':abs(lost-gain),'maximum_relative_total_energy_drift':energy_error,
            'maximum_relative_total_momentum_drift':momentum_error,'minimum_field_on_sampled_grid':float(min(min_n)),
            'photon_energy_increase_steps':int(np.sum(np.diff(ep)>1e-11*max(packet_energy,1e-15))),
            'local_transport':{'interval':[float(region[0]),float(region[-1])],
                               'source_work_in_region':region_work,'net_energy_out_through_region_boundaries':region_out,
                               'energy_out_through_left_boundary':float(simpson(left_flux,x=times)),
                               'final_excess_field_energy_in_region':float(region_energy[-1]),
                               'final_background_cross_energy_in_region':float(region_cross[-1]),
                               'local_energy_balance_absolute_error':float(region_error),
                               'interpretation':'Uses the full field energy flux including its cross term with the initial rolling background. Nonzero Fourier-mode energy alone is not an outgoing-energy fraction.'},
            'final_packet_x':float(sol.y[-2,-1])}


def main():
    protocol=json.loads((HERE/'protocol.json').read_text());p=protocol['parameters'];g=protocol['grid'];runs=[]
    for modes in g['modes']:
        for inertia in g['inertia']:
            for energy in g['initial_packet_energy']:
                result=simulate(modes,inertia,energy,p['initial_rolling_rate'],protocol);runs.append(result)
                print(f"modes={modes} K={inertia} E={energy}: loss={result['photon_fraction_lost']:.6g}, wave share={result['nonzero_wave_share_of_received_energy']:.6g}",flush=True)
    convergence=[]
    for inertia in g['inertia']:
        for energy in g['initial_packet_energy']:
            pair=[next(r for r in runs if r['modes']==n and r['inertia']==inertia and r['initial_packet_energy']==energy) for n in [32,64]]
            change=abs(pair[1]['final_packet_energy']-pair[0]['final_packet_energy'])/energy
            assert change<protocol['checks']['maximum_32_to_64_mode_final_energy_change_over_initial_packet_energy']
            convergence.append({'inertia':inertia,'packet_energy':energy,'relative_final_energy_change':change})
    controls=[simulate(64,1,e,0,protocol) for e in g['initial_packet_energy']]
    shifted=simulate(64,1,.01,p['initial_rolling_rate'],protocol,x0=1.7)
    reference=next(r for r in runs if r['modes']==64 and r['inertia']==1 and r['initial_packet_energy']==.01)
    translation_error=max(abs(shifted[key]-reference[key]) for key in ['photon_energy_lost','nonzero_wave_energy_gain','final_field_momentum'])
    assert translation_error<1e-10
    empty=simulate(32,1,0,p['initial_rolling_rate'],protocol)
    assert empty['nonzero_wave_energy_gain']==0 and abs(empty['homogeneous_field_energy_gain'])<1e-14
    result={'scope':protocol['scope'],'checks_pass':True,'runs':runs,'mode_convergence':convergence,
            'zero_initial_rolling_controls':controls,'translation_control':shifted,'maximum_translation_energy_momentum_difference':translation_error,
            'no_photon_control':empty,'energy_receiver':'The candidate companion field; homogeneous and propagating-mode contributions are reported separately.',
            'full_void_time_model_derived':False,'traveling_graviton_identity_established':False,'capture_or_gravity_response_established':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'check.py']}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'companion-backreaction-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'checks_pass':True,'max_convergence_change':max(r['relative_final_energy_change'] for r in convergence),'zero_rolling':controls},indent=2))


if __name__=='__main__':main()
