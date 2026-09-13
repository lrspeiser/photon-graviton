"""Conditional local capture kinematics, not cross sections or a formation fit."""
from pathlib import Path
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
source = HERE / 'circular-reservoir-support-refined.json'
support = json.loads(source.read_text())
rows = []
for row in support['rows']:
    eps = row['kinetic_to_rest_energy']
    # Recoil/input <= eps is a comparison budget, not a physical capture criterion.
    q_composite = 2*eps/(1-eps)**2
    q_separate = 2*eps/(1-eps**2)
    beta = row['max_circular_speed_kms']/299792.458
    rows.append(dict(label=row['label'],kinetic_budget_fraction=eps,
        composite_min_seed_rest_energy_over_packet=1/q_composite,
        separate_min_seed_rest_energy_over_packet=1/q_separate,
        pair_max_energy_imbalance=beta,
        equal_energy_isotropic_pair_slow_fraction=beta**2,
        constant_cross_section_collision_weighted_fraction=2*beta**2-beta**4,
        speed_ceiling_kms=row['max_circular_speed_kms']))

checks=[]
for q in (1e-12,1e-9,1e-6,.001,.1,1.,10.,1000.):
    # Set seed rest energy = 1 and c = 1; incident null packet (q, q).
    mf=math.sqrt(1+2*q)
    stored=q*2/(mf+1)
    recoil=q*q/(1+q+mf)
    # Separate new particle at rest, unchanged seed invariant mass.
    seed_final_energy=math.hypot(1,q)
    separate_recoil=q*q/(seed_final_energy+1)
    separate_stored=q-separate_recoil
    assert abs((stored+recoil)/q-1)<2e-13
    assert abs((separate_stored+separate_recoil)/q-1)<2e-13
    assert math.isclose((1+q)**2-q*q,mf*mf,rel_tol=1e-10)
    # Seed carries all momentum q; the new particle has zero momentum.
    assert separate_stored>0
    checks.append(dict(q=q,composite_stored_fraction=stored/q,
        composite_recoil_fraction=recoil/q,separate_stored_fraction=separate_stored/q,
        separate_recoil_fraction=separate_recoil/q))
for row in rows:
    eps=row['kinetic_budget_fraction']
    q=1/row['separate_min_seed_rest_energy_over_packet']
    assert math.isclose(q/(math.hypot(1,q)+1),eps,rel_tol=1e-12)
    q=1/row['composite_min_seed_rest_energy_over_packet']
    assert math.isclose(q/(1+q+math.sqrt(1+2*q)),eps,rel_tol=1e-12)
    # Equal-energy two-null-packet invariant: beta^2=(1+cos theta)/2.
    beta=row['speed_ceiling_kms']/299792.458
    cos_boundary=2*beta*beta-1
    assert abs((cos_boundary+1)/2-beta*beta)<1e-16
    for delta in (0.,beta/2,beta,2*beta):
        fraction=max(0.,(beta*beta-delta*delta)/(1-delta*delta))
        assert 0<=fraction<=beta*beta*(1+1e-12)
result=dict(scope='Local kinematics only. Optional receiver and pair branches; no new fits, interaction rate, mass spectrum, or orbital formation.',
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),rows=rows,checks=checks,
    verification='Energy and invariant identities checked across eight incident/seed ratios; analytic budget inverses and pair-angle boundary checked for all ten support profiles.')
(HERE/'capture-recoil-threshold-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
for row in rows:
    print(row['label'],f"seed/packet >= {row['separate_min_seed_rest_energy_over_packet']:.3g}",f"pair fraction {row['equal_energy_isotropic_pair_slow_fraction']:.3g}")
