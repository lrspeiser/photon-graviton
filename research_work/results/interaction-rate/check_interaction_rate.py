"""Necessary interaction-rate conditions; no microscopic graviton rate claimed."""
from pathlib import Path
import json
import os
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'interaction-rate'
C=299792.458


def dot(a,b):
    return a[0]*b[0]-np.dot(a[1:],b[1:])


def vertex(k,kp,e,ep):
    return dot(k,kp)*dot(e,ep)-dot(k,ep)*dot(kp,e)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    fitted=json.loads((ROOT/'research_work/results/conversion-first/results.json').read_text(encoding='utf-8'))
    alpha=fitted['fitted']['alpha_per_mpc']
    # The scalar F^2 vertex: Ward identities and real, collinear vacuum limit.
    rng=np.random.default_rng(2026090903)
    ward_errors=[]; collinear=[]; noncollinear_q2=[]
    for _ in range(100):
        energy=rng.uniform(.5,5); remaining=rng.uniform(.1,.9)*energy
        theta=rng.uniform(.01,np.pi-.01)
        k=np.array([energy,0,0,energy])
        kp=np.array([remaining,remaining*np.sin(theta),0,remaining*np.cos(theta)])
        e=np.array([0,1,0,0]); ep=np.array([0,np.cos(theta),0,-np.sin(theta)])
        ward_errors.extend([abs(vertex(k,kp,k,ep)),abs(vertex(k,kp,e,kp))])
        q=k-kp
        q2=dot(q,q)
        assert q2<0
        noncollinear_q2.append(float(q2))
        kp0=np.array([remaining,0,0,remaining])
        for ep0 in [np.array([0,1,0,0]),np.array([0,0,1,0])]:
            collinear.append(abs(vertex(k,kp0,e,ep0)))
    assert max(ward_errors)<1e-12 and max(collinear)<1e-12
    # Compare absolute vs fractional energy loss at one calibrated reference line.
    reference_z=.1; reference_energy=1.
    tau=np.log1p(reference_z); distance=tau/alpha
    absolute_loss=reference_energy*reference_z/(1+reference_z)
    chromatic=[]
    for energy in [.5,1.,2.,4.]:
        fixed_final=energy-absolute_loss
        fractional_final=energy*np.exp(-tau)
        chromatic.append(dict(initial_energy_ev=energy,fixed_quantum_mean_z=energy/fixed_final-1,
                              fractional_loss_z=energy/fractional_final-1,
                              absolute_loss_ev=absolute_loss, path_mpc=distance))
    assert np.ptp([r['fractional_loss_z'] for r in chromatic])<1e-14
    # Energy-independent jump rates on an interior energy ladder. A birth event
    # debits the environment; a death event credits it. Boundary probabilities
    # are checked before using the infinite-ladder first/second moment formulas.
    step=.001; energies=np.arange(1,2002)*step; initial_index=999
    initial=np.zeros(len(energies)); initial[initial_index]=1
    length=100.; coupling_rate=1.
    bath=[]
    for occupation in [0.,1.,10.]:
        down=coupling_rate*(occupation+1); up=coupling_rate*occupation
        diagonal=np.full(len(energies),-(up+down))
        diagonal[0]=-up; diagonal[-1]=-down
        generator=diags([np.full(len(energies)-1,up),diagonal,np.full(len(energies)-1,down)],[-1,0,1],format='csc')
        probability=expm_multiply(generator*length,initial,traceA=float(generator.diagonal().sum()*length))
        mean=float(probability@energies)
        variance=float(probability@((energies-mean)**2))
        predicted_mean=1-step*coupling_rate*length
        predicted_variance=step**2*coupling_rate*(2*occupation+1)*length
        assert abs(probability.sum()-1)<1e-10 and probability.min()>-1e-12
        assert probability[0]+probability[-1]<1e-15
        assert abs(mean-predicted_mean)<1e-9 and abs(variance-predicted_variance)<1e-10
        bath.append(dict(occupation=occupation, mean_photon_energy_ev=mean,
                         net_energy_received_by_environment_ev=1-mean, variance_ev2=variance,
                         predicted_variance_ev2=predicted_variance,
                         added_energy_rms_as_km_s=float(C*np.sqrt(variance)/mean),
                         probability_at_boundaries=float(probability[0]+probability[-1])))
    # For scale-free fractional jumps at fixed illustrative width, the calibrated
    # alpha implies a rate and mean path between events, not a measured coupling.
    width_conditions=[]
    for shift in [1.1,2.]:
        for width in [1.,10.,100.]:
            eps=np.log1p((width/C)**2)/np.log(shift)
            event_rate=alpha/eps
            width_conditions.append(dict(shift=shift,illustrative_added_width_km_s=width,
                                         max_fraction_per_event=float(eps),
                                         min_events_per_mpc=float(event_rate),
                                         max_mean_event_spacing_pc=float(1e6/event_rate)))
    result=dict(scope='Necessary constraints and stipulated transport comparisons; no new astronomical fit or physical rate',
                source_alpha_per_mpc=alpha,
                scalar_vertex=dict(signature='+---',samples=100,max_ward_error=float(max(ward_errors)),
                                   max_collinear_vertex=float(max(collinear)),
                                   all_noncollinear_companion_four_momenta_spacelike=True),
                absolute_vs_fractional=chromatic,
                reversible_bath_ladder=bath,
                fractional_event_requirements=width_conditions,
                checks=dict(scalar_gauge_and_kinematic_identities=True,
                            fixed_fraction_achromatic=True,
                            ladder_distribution_moments_and_environment_ledger=True),
                assumptions=['Scalar F^2 vertex is one comparison, not an ordinary graviton identity.',
                             'Standard massless photon dispersion and vacuum momentum conservation for vertex result.',
                             'Ladder down/up rates are specified, not derived; environment occupancy fixed in reservoir limit.',
                             'Environment energy is explicitly debited/credited; its finite capacity must be derived for a real model.',
                             'A source of recoil, medium dispersion, time-dependent background or different action requires new calculations.'])
    (OUT/'rate-constraints.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],chromatic=chromatic,bath=bath,
                         z1_width10=[r for r in width_conditions if r['shift']==2 and r['illustrative_added_width_km_s']==10][0]),indent=2))


if __name__=='__main__':
    main()
