from pathlib import Path
import numpy as np
from scipy.linalg import expm
import json,csv,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
source=ROOT/'research_work/results/electromagnetic-audit/protocol.json'
alpha=json.loads(source.read_text())['alpha_per_mpc']
def spectrum(nu):return np.exp(-.5*((nu-1.4)/.015)**2)+.6*np.exp(-.5*((nu-2.4)/.015)**2)
def integrate(y,x):return float(np.trapezoid(y,x))
rows=[];max_unitary=0;max_analytic=0
for D in [30.660139,40.7,100,1000]:
    S=np.exp(alpha*D);fraction=1-1/S;theta=np.arcsin(np.sqrt(fraction))
    U=expm(-1j*np.array([[0,theta],[theta,0]]))
    max_unitary=max(max_unitary,float(np.max(abs(U.conj().T@U-np.eye(2)))))
    analytic=np.array([[np.cos(theta),-1j*np.sin(theta)],[-1j*np.sin(theta),np.cos(theta)]])
    max_analytic=max(max_analytic,float(np.max(abs(U-analytic))))
    trials=[]
    for count in [25001,50001]:
        nu=np.linspace(.5,3,count);N=spectrum(nu)
        photon=abs(U[0,0])**2*N;gravity=abs(U[1,0])**2*N
        shifted=S*spectrum(S*nu)
        n0=integrate(N,nu);e0=integrate(nu*N,nu)
        ng=integrate(photon,nu);eg=integrate(nu*photon,nu)
        ec=integrate(nu*gravity,nu);ns=integrate(shifted,nu);es=integrate(nu*shifted,nu)
        trial=np.array([ng/n0,eg/e0,ec/e0,(eg/ng)/(e0/n0),ns/n0,es/e0,(es/ns)/(e0/n0)])
        assert abs((eg+ec)/e0-1)<1e-12
        assert np.max(abs(trial-np.array([1/S,1/S,fraction,1,1,1/S,1/S])))<1e-11
        trials.append(trial)
    assert np.max(abs(trials[0]-trials[1]))<1e-12
    rows.append(dict(distance_mpc=D,target_redshift=S-1,chosen_conversion_fraction=fraction,
        chosen_mixing_angle=float(theta),photon_number_remaining=trials[1][0],
        photon_energy_remaining=trials[1][1],gravity_energy_fraction=trials[1][2],
        surviving_mean_frequency_ratio=trials[1][3],required_mean_frequency_ratio=trials[1][6],
        predicted_conversion_line_redshift=0.,event_interval_ratio=1.,
        required_event_interval_ratio=float(S),extra_relative_delay_common_speed_seconds=0.))
assert max_unitary<1e-12 and max_analytic<1e-12
with (HERE/'comparison.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
result=dict(status='Synthetic spectral discriminator; conversion angle matches energy by construction, not an astrophysical fit.',
    alpha_per_mpc=alpha,cases=rows,checks=dict(max_unitarity_error=max_unitary,max_expm_trigonometric_error=max_analytic,
        energy_integrals_and_grid_refinement='pass'),
    sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,HERE/'protocol.md']})
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
