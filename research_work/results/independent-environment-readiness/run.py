"""Feature-only catalog audit; no new target redshifts read or fitted."""
from pathlib import Path
import csv, json, hashlib, math
from collections import Counter
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
raw=HERE/'features.tsv'
lines=raw.read_text().splitlines()
header=next(i for i,s in enumerate(lines) if s.startswith('Name\t'))
rows=list(csv.DictReader([lines[header]]+[s for s in lines[header+3:] if s.strip() and not s.startswith('#')], delimiter='\t'))
assert len(rows)==869
allowed={'TRGB','Cep','geom','SBF','CMD','HB','RR','PNLF','SN'}
def angle(s,ra=False):
    vals=s.split(); sign=-1 if vals[0].startswith('-') else 1
    return sign*(abs(float(vals[0]))+float(vals[1])/60+float(vals[2])/3600)*(15 if ra else 1)
def vector(ra,dec):
    ra,dec=np.deg2rad([ra,dec])
    return np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)])
features=[]
for r in rows:
    method=r['f_Dist'].strip(); d=float(r['Dist'])
    status='indicator_candidate' if method in allowed else 'pending_calibration_or_membership'
    if method in {'h',"h'"}:status='excluded_redshift_distance'
    if method=='txt':status='excluded_unmeasured_distance'
    features.append(dict(name=r['Name'].strip(), ra_deg=angle(r['RAJ2000'],True),
        dec_deg=angle(r['DEJ2000']), distance_mpc=d, method=method,
        status=status, inside_nominal_11mpc=d<=11,
        exposure='distance_and_position_viewed; no_velocity_requested; not_certified_fresh'))
with (HERE/'distance-feature-register.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(features[0]),lineterminator='\n');w.writeheader();w.writerows(features)
tracers=[r for r in features if r['status']=='indicator_candidate' and r['distance_mpc']<=11]
xyz=np.array([vector(r['ra_deg'],r['dec_deg'])*r['distance_mpc'] for r in tracers])

def occupied(n,L,radius):
    # Known ray/sphere intersection and interval-union geometry.
    center=xyz@n
    perpendicular2=np.sum(xyz*xyz,axis=1)-center**2
    good=perpendicular2<=radius**2
    half=np.sqrt(np.maximum(0,radius**2-perpendicular2[good]))
    intervals=sorted((max(0,a),min(L,b)) for a,b in zip(center[good]-half,center[good]+half) if min(L,b)>max(0,a))
    total=0; end=0
    for a,b in intervals:
        total+=max(0,b-max(a,end));end=max(end,b)
    return total

target_path=ROOT/'redshift_paper/all_164_groups.csv'
# Only features are converted; all 164 targets were already exposed historically.
targets=list(csv.DictReader(target_path.open()))
out=[]
for t in targets:
    d=float(t['catalog_distance_mpc']); L=min(d,11)
    n=vector(float(t['ra_deg']),float(t['dec_deg']))
    for radius in [.5,1,2]:
        V=L-occupied(n,L,radius)
        out.append(dict(pgc=t['pgc'],distance_mpc=d,radius_mpc=radius,
            nominal_covered_fraction=L/d,unmapped_path_mpc=d-L,
            local_no_neighbor_path_mpc=V,full_proxy_lower_mpc=V,
            full_proxy_upper_mpc=V+d-L))
with (HERE/'path-coverage.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]),lineterminator='\n');w.writeheader();w.writerows(out)

# Independent midpoint integration verifies interval-union implementation.
maxerr=0
for t in targets[::20]:
    d=float(t['catalog_distance_mpc']); L=min(d,11)
    n=vector(float(t['ra_deg']),float(t['dec_deg']))
    grid=(np.arange(20000)+.5)*L/20000
    for radius in [.5,1,2]:
        near=np.zeros(len(grid),dtype=bool)
        for p in xyz:
            near |= np.sum((grid[:,None]*n-p)**2,axis=1)<radius**2
        error=abs(float(near.mean()*L)-occupied(n,L,radius));maxerr=max(maxerr,error)
        assert error<.01
coverage=np.array([min(float(t['catalog_distance_mpc']),11)/float(t['catalog_distance_mpc']) for t in targets])
result=dict(source_rows=len(rows),unique_names=len(set(r['name'] for r in features)),
    method_counts=dict(Counter(r['method'] for r in features)),
    status_counts=dict(Counter(r['status'] for r in features)),
    eligible_indicator_tracers_within_11mpc=len(tracers),targets=len(targets),
    nominal_covered_fraction_min_median_max=[float(x) for x in [min(coverage),np.median(coverage),max(coverage)]],
    targets_with_whole_path_inside_11mpc=int(sum(coverage==1)),
    targets_more_than_half_path_outside=int(sum(coverage<.5)),
    midpoint_union_max_error_mpc=maxerr,
    scope='Catalog-proximity geometry, not physical void reconstruction. Nominal radial coverage is optimistic; angular/stellar-population selection and distance uncertainties are unresolved. No new redshifts accessed or model fit.',
    sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw,HERE/'ReadMe',HERE/'protocol.md',target_path]})
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
