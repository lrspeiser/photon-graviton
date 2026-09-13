"""Capture area versus storage capacity in the existing exact-third profile."""
from pathlib import Path
import hashlib
import json
import math
from scipy.integrate import quad
from scipy.special import gamma, gammainc

P=Path(__file__).resolve().parent
OLD=P.parent/'isotropic-galaxy-transfer'
BASE=P.parents[2]/'temporal_candidate_audit/data'
fitfile=OLD/'third-radiation-retention-results.json'
predfile=OLD/'third-radiation-retention-predictions.json'
fit=json.loads(fitfile.read_text())
model=fit['models']['attenuated']
saved=json.loads(predfile.read_text())
for name,digest in fit['input_sha256'].items():
    assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest
meta={}
for line in (BASE/'SPARC_Lelli2016c.mrt').read_text().splitlines():
    f=line.split()
    if len(f)==19:
        try: meta[f[0]]=float(f[11])
        except ValueError: pass
C=model['C_Msun_kpc3']; k=model['k0_per_kpc']; scale=model['scale_to_disk']
out=dict(scope='Capacity inferred from the existing deposited-density shape at retention one; conditional storage interpretation, no fitted changes',
    C_Msun_kpc3=C,k0_per_kpc=k,scale_to_disk=scale,
    input_sha256=fit['input_sha256'],
    model_input_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [fitfile,predfile]},rows=[])
for row in saved:
    a=scale*meta[row['galaxy']]
    T=math.pi*k*a/2
    area=math.pi*a*a*(T**(2/3)*gamma(1/3)*gammainc(1/3,T)+math.expm1(-T))
    area_quad=math.pi*a*a*quad(lambda y:-math.expm1(-T*y**1.5)/y**2,0,1,epsabs=1e-9)[0]
    assert abs(area_quad/area-1)<1e-7
    capacity=C/k*area
    unattenuated_capacity=C*math.pi**2*a**3
    fraction=capacity/unattenuated_capacity
    assert 0 < fraction < 1
    area_per_capacity=area/capacity
    assert abs(area_per_capacity/(k/C)-1)<1e-12
    out['rows'].append(dict(galaxy=row['galaxy'],split=row['split'],a_kpc=a,
        center_chord_optical_depth=T,capture_area_kpc2=area,
        reference_capacity_mass_Msun=capacity,
        unattenuated_capacity_mass_Msun=unattenuated_capacity,
        attenuated_to_unattenuated_capacity=fraction,
        area_per_reference_capacity_kpc2_per_Msun=area_per_capacity,
        area_per_unattenuated_capacity_kpc2_per_Msun=area/unattenuated_capacity,
        cross_section_relative_check=abs(area_quad/area-1)))
out['ranges']={key:[min(r[key] for r in out['rows']),max(r[key] for r in out['rows'])]
    for key in ['a_kpc','center_chord_optical_depth','capture_area_kpc2','reference_capacity_mass_Msun','attenuated_to_unattenuated_capacity','area_per_reference_capacity_kpc2_per_Msun']}
assert len(out['rows'])==149
(P/'capture-capacity-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['ranges'],indent=2))
