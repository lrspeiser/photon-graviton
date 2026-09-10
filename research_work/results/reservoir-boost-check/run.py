"""Late-state center-of-energy test for the closed canonical recoil model."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=HERE.parent/'bounded-profile-momentum/moving-results.json'

def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    d=json.loads(SOURCE.read_text())
    for path,expected in d['source_hashes'].items():assert digest(ROOT/path)==expected
    rows=[]
    for r in d['runs']:
        mu=r['support_rest_energy_ratio'];Eg=r['outgoing_wave_energy_fraction']
        Er=r['reservoir_energy_fraction'];K=r['support_momentum_fraction']
        Es=mu+r['support_recoil_energy_fraction']
        total=Eg+Er+Es;momentum=Eg+K
        speed_support=K/Es
        # Both outgoing profiles translate at c=1 after leaving the region.
        energy_flow=Eg+Er+Es*speed_support
        centroid_speed=energy_flow/total
        required_speed=momentum/total
        deficit=energy_flow-momentum
        assert abs(deficit-Er)<1e-12 and abs(momentum-1)<1e-10
        # Independently form finite-time energy centroids after all interactions.
        times=np.array([0.,.25,1.])
        positions=np.array([22+times,22+times,r['final_support_position']+speed_support*times])
        center=np.array([Eg,Er,Es])@positions/total
        finite_speed=(center[-1]-center[0])/(times[-1]-times[0])
        assert abs(finite_speed-centroid_speed)<1e-12
        rows.append(dict(support_rest_energy_ratio=mu,tight=r['tight'],
            total_energy_including_support_rest_energy=total,total_canonical_momentum=momentum,
            outgoing_reservoir_energy=Er,energy_centroid_speed=centroid_speed,
            required_relativistic_centroid_speed=required_speed,
            centroid_speed_difference=centroid_speed-required_speed,
            relative_center_speed_mismatch=centroid_speed/required_speed-1,
            energy_flow_minus_canonical_momentum=deficit,
            added_momentum_if_reservoir_relabelled_massless=Er))
    # No-conversion control satisfies the relation with an initially stationary support.
    for mu in [10.,100.,10000.]:
        Eg,Er,K,Es=1.,0.,0.,mu
        assert (Eg+Er+K)/(Eg+Er+Es)==(Eg+K)/(Eg+Er+Es)
    result=dict(scope='Necessary symmetric-stress/boost condition for the modeled late state; not a general no-go theorem.',
        assumed_closed_system='Only outgoing EM, outgoing advected reservoir and the previously modeled support.',
        known_identity='dX_energy/dt=c^2*P_total/E_total for an isolated system with conserved symmetric stress-energy.',
        rows=rows,all_nonzero_conversion_cases_fail=True,
        omitted_background_or_preferred_frame_ruled_out=False,holdouts_opened=False,
        source_hashes={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),SOURCE]})
    assert all(r['centroid_speed_difference']>0 for r in rows)
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
