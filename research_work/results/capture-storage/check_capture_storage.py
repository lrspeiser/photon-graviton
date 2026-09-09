"""Resonance area, reversible three-level storage, capacity and energy accounting."""
from pathlib import Path
import json
import os
import numpy as np
from scipy.constants import hbar, electron_volt, proton_mass, c
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'capture-storage'


def generator(kc,ks,kd,nc,ns=0.,nr=0.,eb=1.,ed=.9):
    """States g,b,d; bosonic bath occupations and secular weak-coupling rates.

    Last six rows accumulate energy entering/leaving through each bath. Finite
    source/receiver inventories can be supplied from these counters; occupations
    here are prescribed reservoir conditions, not a globally closed galaxy model.
    """
    matrix=np.zeros((9,9))
    edges=[(0,1,kc*nc,kc*(nc+1),eb),
           (2,1,ks*ns,ks*(ns+1),eb-ed),
           (0,2,kd*nr,kd*(nr+1),ed)]
    for j,(low,high,up,down,gap) in enumerate(edges):
        matrix[high,low]+=up; matrix[low,low]-=up
        matrix[low,high]+=down; matrix[high,high]-=down
        matrix[3+2*j,low]=gap*up
        matrix[4+2*j,high]=gap*down
    return matrix


def evolve(matrix,initial,times,eb=1.,ed=.9):
    y0=np.r_[initial,np.zeros(6)]
    sol=solve_ivp(lambda t,y:matrix@y,(0,float(times[-1])),y0,t_eval=times,
                  method='DOP853',rtol=2e-10,atol=1e-12)
    assert sol.success
    final_exact=expm(matrix*times[-1])@y0
    error=float(np.max(abs(sol.y[:,-1]-final_exact)))
    assert error<2e-7*max(1,float(np.max(abs(final_exact))))
    internal=eb*sol.y[1]+ed*sol.y[2]
    net_in=np.sum(sol.y[[3,5,7]],axis=0)-np.sum(sol.y[[4,6,8]],axis=0)
    ledger=float(np.max(abs(internal-internal[0]-net_in)))
    assert ledger<2e-8*max(1,float(np.max(abs(net_in))))
    assert np.max(abs(sol.y[:3].sum(axis=0)-1))<1e-9
    assert np.min(sol.y[:3])>-1e-9
    return sol.y, {'energy_ledger_max_error':ledger,'ODE_vs_matrix_exponential_error':error}


def main():
    areas=[]
    # Normalize cross sections to sigma_unit, energy offsets to total width/2.
    for gc,gs in [(1.,0.),(1.,.1),(1.,1.),(1.,10.),(1e-6,0.)]:
        total=gc+gs
        inclusive=quad(lambda x:gc/total/(1+x*x)*total/2,-np.inf,np.inf)[0]
        stored=quad(lambda x:gc*gs/total**2/(1+x*x)*total/2,-np.inf,np.inf)[0]
        assert np.isclose(inclusive,np.pi*gc/2,rtol=1e-10)
        assert np.isclose(stored,np.pi*gc*gs/(2*total),rtol=1e-10,atol=1e-15)
        areas.append({'Gamma_c':gc,'Gamma_s':gs,'inclusive_area_over_sigma_unit':inclusive,
                      'shelving_area_over_sigma_unit':stored,'shelving_peak_over_sigma_unit':gc*gs/total**2})

    lifetimes=[]
    for years in [1e9,1e10]:
        seconds=years*365.25*86400
        # Retain at least 90% of a cohort over the illustrative interval.
        width=-(hbar/electron_volt)*np.log(.9)/seconds
        band=1e-6
        lifetimes.append({'retention_interval_years':years,'required_survival':.9,
                          'maximum_direct_width_eV':width,'illustrative_bandwidth_eV':band,
                          'maximum_flat_band_average_over_sigma_unit':np.pi*width/(2*band)})

    cases=[]
    for kc,ks,kd,nc in [(1.,9.,0.,.01),(1.,9.,.001,.01),(1.,9.,.1,.01),
                         (1.,1.,0.,.1),(1.,9.,0.,10.)]:
        matrix=generator(kc,ks,kd,nc)
        values,checks=evolve(matrix,[1.,0.,0.],np.linspace(0,3000,301))
        final=values[:3,-1]
        if kd==0:
            expected=np.array([0.,0.,1.])
            assert abs(final[2]-1)<1e-8
        else:
            w=kc*nc
            b=1/(1+(w+kc+ks)/w+ks/kd)
            expected=np.array([(w+kc+ks)*b/w,b,ks*b/kd])
            assert np.max(abs(final-expected))<1e-8
        cases.append({'kc':kc,'ks':ks,'kd':kd,'incident_occupation':nc,
                      'final_g_b_d':final.tolist(),'analytic_steady_g_b_d':expected.tolist(),
                      'stored_energy_per_receiver':float(.9*final[2]),
                      'incident_energy_counter':float(values[3,-1]),
                      'returned_companion_energy_counter':float(values[4,-1]),
                      'shelving_radiation_energy_counter':float(values[6,-1]),
                      'deposit_leak_energy_counter':float(values[8,-1]),'checks':checks})

    # Isolated initially caught packet: shelving competes with companion reemission.
    y,check=evolve(generator(1.,9.,0.,0.),[0.,1.,0.],np.linspace(0,20,101))
    assert np.allclose([y[2,-1],y[4,-1],y[6,-1]],[.9,.1,.09],atol=1e-9)
    # Eb=1, Ed=.9: final stored .81, returned .10, shelving radiation .09.
    isolated={'stored_fraction':float(.9*y[2,-1]),'returned_fraction':float(y[4,-1]),
              'shelving_radiation_fraction':float(y[6,-1]),'checks':check}

    reverse=[]
    for ns in [1e-4,.01,1.]:
        kc,ks=1.,9.
        upward=ks*ns; downward=ks*(ns+1); a=kc+downward
        slow=2*upward*kc/(a+upward+np.sqrt((a+upward)**2-4*upward*kc))
        matrix=generator(kc,ks,0.,0.,ns=ns)
        # Exact exponentials avoid explicit timestepping over widely separated rates.
        initial=np.r_[[0.,0.,1.],np.zeros(6)]
        final=expm(matrix*(10/slow))@initial
        internal=final[1]+.9*final[2]
        net_in=final[[3,5,7]].sum()-final[[4,6,8]].sum()
        assert abs(internal-.9-net_in)<1e-8
        eigenvalues=np.linalg.eigvals(matrix[:3,:3])
        assert np.min(abs(eigenvalues+slow))<1e-10
        assert final[2]<1e-4
        reverse.append({'shelving_bath_occupation':ns,'effective_slow_decay_rate':float(slow),
                        'dark_population_after_ten_slow_lifetimes':float(final[2]),
                        'ledger_error':float(abs(internal-.9-net_in))})

    # Incoherent broadband reservoir master equation: a Gibbs distribution is
    # stationary when all three mode occupations correspond to the same T.
    thermal=[]
    for temperature in [.05,.5,5.]:
        gaps=np.array([1.,.1,.9])
        occupations=1/np.expm1(gaps/temperature)
        matrix=generator(1.,9.,.01,float(occupations[0]),ns=float(occupations[1]),nr=float(occupations[2]))
        gibbs=np.exp(-np.array([0.,1.,.9])/temperature);gibbs/=gibbs.sum()
        error=float(np.max(abs(matrix[:3,:3]@gibbs)))
        assert error<1e-12
        thermal.append({'temperature_in_Eb_units':temperature,'Gibbs_stationarity_error':error})

    rest_eV=proton_mass*c*c/electron_volt
    results={'scope':'Conditional resonant absorber and reversible three-level receiver. No material identity, cosmic abundance, cross-section normalization, ordinary-graviton capture rate or halo is established.',
             'resonance_integrals':areas,'direct_storage_width_examples':lifetimes,
             'driven_storage_cases':cases,'single_caught_packet':isolated,
             'reverse_shelving_bath':reverse,'thermal_detailed_balance':thermal,
             'capacity_example':{'assumption':'One 1 eV storage slot per proton; illustrative, not an identified material.',
                                 'maximum_extra_rest_mass_fraction':float(1/rest_eV),
                                 'energy_per_proton_eV_for_equal_added_rest_mass':float(rest_eV)},
             'unresolved':['Actual receiver and its microscopic transition matrix elements',
                           'Companion spectrum and absorption bandwidth',
                           'Temperature, collisions and all reverse channels',
                           'Finite receiver abundance, storage capacity and spatial distribution',
                           'Supply, recoil, gravitational support and simultaneous lensing',
                           'Emission spectrum and destination of capture-release energy']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'capture-storage-results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'isolated_capture':isolated,'capacity_example':results['capacity_example'],
                       'direct_storage_width_examples':lifetimes},indent=2))


if __name__=='__main__':main()
