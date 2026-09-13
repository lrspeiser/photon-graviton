"""Local energy-momentum ledger for massless companion absorption by a receiver."""
from pathlib import Path
import json,math

OUT=Path(__file__).resolve().parent

def absorb(e,j):
    # Receiver initially at rest, Mc^2=1; incident E=e and |p|c=j*e.
    final_rest=math.sqrt(1+2*e+(1-j*j)*e*e)
    kinetic_fraction=j*j*e/(1+e+final_rest)
    stored_fraction=1-kinetic_fraction
    speed=j*e/(1+e)
    # Independent direct sum of two opposite null four-vectors and receiver.
    Eplus=e*(1+j)/2;Eminus=e*(1-j)/2
    energy=1+Eplus+Eminus;momentum=Eplus-Eminus
    invariant_residual=abs(energy*energy-momentum*momentum-final_rest*final_rest)/(energy*energy)
    ledger_residual=abs((1+e*stored_fraction+e*kinetic_fraction)-energy)/energy
    momentum_residual=abs(speed*energy-momentum)/max(e,1e-300)
    rest_gain_residual=abs((final_rest-1)-e*stored_fraction)/max(1,e)
    assert max(invariant_residual,ledger_residual,momentum_residual,rest_gain_residual)<2e-13
    assert 0<=stored_fraction<=1 and 0<=speed<1
    return dict(incident_energy_over_receiver_rest_energy=e,net_momentum_fraction=j,
        internal_rest_energy_gain_over_incident_energy=stored_fraction,
        recoil_kinetic_energy_over_incident_energy=kinetic_fraction,receiver_speed_over_c=speed,
        invariant_scaled_residual=invariant_residual,energy_ledger_scaled_residual=ledger_residual,
        momentum_scaled_residual=momentum_residual)

rows=[absorb(e,j) for e in [1e-12,1e-6,.001,.1,1.,10.,1000.] for j in [0.,.1,1.]]
thresholds=[]
for b in [1e-4,.001,.01]:
    limit=b/(1-b)
    lower=absorb(.99*limit,1.);upper=absorb(1.01*limit,1.)
    assert lower['receiver_speed_over_c']<b<upper['receiver_speed_over_c']
    thresholds.append(dict(assumed_escape_speed_over_c=b,maximum_single_direction_energy_over_receiver_rest_energy=limit))

result=dict(scope='Special-relativistic local absorption into a receiver initially at rest, no escaping products. Companion p=E/c assumed. Kinematic allowance, not an interaction rate, lifetime or microscopic mechanism.',
    provenance='Known energy-momentum conservation and invariant-mass relation applied to hypothetical companions. Numerical ratios are synthetic, not measured receiver masses.',
    formulas=['M_final*c^2=sqrt((M*c^2+E)^2-j^2*E^2)',
    'Delta_internal=(M_final-M)*c^2; K_recoil=E-Delta_internal',
    'v/c=j*E/(M*c^2+E); j=|sum(p)|*c/E',
    'For a single direction and weak external well, v<v_escape requires E/(M*c^2)<b/(1-b), b=v_escape/c'],
    cases=rows,weak_well_retention_examples=thresholds)
(OUT/'capture-ledger-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(dict(single_direction_cases=[r for r in rows if r['net_momentum_fraction']==1.],thresholds=thresholds),indent=2))
