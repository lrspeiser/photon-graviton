"""Direction-only frame diagnostic. No candidate outcome table is read."""
import csv
import hashlib
import json
import math
from pathlib import Path
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

HERE=Path(__file__).resolve().parent
C=299792.458
speed=369.82
b=speed/C
gamma=1/math.sqrt(1-b*b)
apex=SkyCoord(l=264.021*u.deg,b=48.253*u.deg,frame='galactic').icrs
axis=apex.cartesian.xyz.value


def factor(n):
    # n points toward source in the solar frame; photon travels along -n.
    return 1/(gamma*(1-b*np.dot(axis,n)))


def boost(energy,momentum,velocity):
    beta=np.linalg.norm(velocity)
    if beta==0:return energy,momentum.copy()
    direction=velocity/beta;g=1/math.sqrt(1-beta*beta)
    projection=np.dot(direction,momentum)
    return g*(energy+beta*projection),momentum+((g-1)*projection+g*beta*energy)*direction


features=json.loads((HERE/'elves-field-feature-audit.json').read_text())
decisions=json.loads((HERE/'elves-identity-decisions.json').read_text())
status={r['target_name']:r['decision'] for r in decisions['records']}
records=[]
for row in features['records']:
    if not row['provisional_eligible']:continue
    sky=SkyCoord(ra=row['ra_deg']*u.deg,dec=row['dec_deg']*u.deg,frame='icrs')
    n=sky.cartesian.xyz.value;f=factor(n)
    energy,momentum=boost(1.,-n,b*axis)
    assert math.isclose(f,1/energy,rel_tol=1e-14)
    eback,pback=boost(energy,momentum,-b*axis)
    assert abs(eback-1)<1e-14 and np.max(abs(pback+n))<1e-14
    records.append(dict(target_name=row['target_name'],decision=status[row['target_name']],
        ra_deg=row['ra_deg'],dec_deg=row['dec_deg'],solar_to_cmb_spectral_factor=f,
        zero_z_equivalent_cdz_kms=C*(f-1),first_order_projection_kms=speed*np.dot(axis,n)))
assert len(records)==29
assert math.isclose(factor(axis),math.sqrt((1+b)/(1-b)),rel_tol=1e-14)
assert math.isclose(factor(-axis),math.sqrt((1-b)/(1+b)),rel_tol=1e-14)
transverse=np.cross(axis,[0,0,1]);transverse/=np.linalg.norm(transverse)
assert math.isclose(factor(transverse),1/gamma,rel_tol=1e-14)
with (HERE/'elves-frame-factors.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');writer.writeheader();writer.writerows(records)
pending=[r for r in records if r['decision']=='pending_freshness_and_host_audit']
assert len(pending)==26
summary=dict(source='https://doi.org/10.1051/0004-6361/201833880',
    speed_kms=speed,apex_galactic_deg=[264.021,48.253],apex_icrs_deg=[apex.ra.deg,apex.dec.deg],
    astropy_version=astropy.__version__,rows=29,pending_rows=26,
    pending_zero_z_equivalent_cdz_range_kms=[min(r['zero_z_equivalent_cdz_kms'] for r in pending),max(r['zero_z_equivalent_cdz_kms'] for r in pending)],
    feature_sha256=hashlib.sha256((HERE/'elves-field-feature-audit.json').read_bytes()).hexdigest(),
    protocol_sha256=hashlib.sha256((HERE/'elves-frame-protocol.md').read_bytes()).hexdigest(),
    target_outcome_values_read=False,
    status='Conditional frame sensitivity; not an applied catalog correction, independent galaxy motion measurement or fitted redshift model.')
(HERE/'elves-frame-results.json').write_text(json.dumps(summary,indent=2)+'\n',newline='\n')
print(json.dumps(summary,indent=2))
