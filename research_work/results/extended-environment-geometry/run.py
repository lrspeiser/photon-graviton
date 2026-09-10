from pathlib import Path
import csv, json, hashlib
from collections import Counter
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
source=ROOT/'temporal_candidate_audit/data/cf4_table2.dat'
metadata=ROOT/'temporal_candidate_audit/data/cf4_ReadMe.txt'
target_source=ROOT/'redshift_paper/all_164_groups.csv'
methods={'SNIA':(41,47,48,52),'TF':(53,59,60,64),'FP':(65,71,72,76),
         'SBF':(77,83,84,89),'SNII':(90,96,97,101),'TRGB':(102,107,108,112),
         'CEP':(113,119,120,125),'MASER':(126,131,132,136)}
preference=['MASER','TRGB','CEP','SBF']
targets=list(csv.DictReader(target_source.open()))
target_ids={int(r['pgc']) for r in targets}
measurements=[]; selected=[]; lines=source.read_text().splitlines()
assert len(lines)==55877
for line in lines:
    pgc=int(line[:7]); ra=float(line[137:145]);dec=float(line[146:154]); choices={}
    for method,(a,b,e,f) in methods.items():
        if not line[a:b].strip():continue
        mu=float(line[a:b]);sigma=float(line[e:f]) if line[e:f].strip() else 0
        d=10**((mu-25)/5)
        row=dict(pgc=pgc,method=method,ra_deg=ra,dec_deg=dec,modulus=mu,sigma_mag=sigma,
                 distance_mpc=d,minus_1sigma_mpc=d-10**((mu-sigma-25)/5),
                 plus_1sigma_mpc=10**((mu+sigma-25)/5)-d,
                 linear_sigma_mpc=np.log(10)/5*d*sigma)
        measurements.append(row)
        if sigma>0:choices[method]=row
    if pgc not in target_ids:
        for method in preference:
            if method in choices:selected.append(choices[method]);break
assert not target_ids & {r['pgc'] for r in selected}
def write_csv(name, rows):
    with (HERE/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
# The 56k source table is already tracked; retain selected features rather than another full dump.
write_csv('selected-tracers.csv',selected)
def directions(ra,dec):
    a,b=np.deg2rad(ra),np.deg2rad(dec)
    return np.column_stack([np.cos(b)*np.cos(a),np.cos(b)*np.sin(a),np.sin(b)])
n=directions([float(t['ra_deg']) for t in targets],[float(t['dec_deg']) for t in targets])
tn=directions([r['ra_deg'] for r in selected],[r['dec_deg'] for r in selected])
mu=np.array([r['modulus'] for r in selected]);se=np.array([r['sigma_mag'] for r in selected])
td=np.array([float(t['catalog_distance_mpc']) for t in targets])
tmu=np.array([float(t['sbf_modulus_mag']) for t in targets]);tse=np.array([float(t['sbf_modulus_error_mag']) for t in targets])
assert np.allclose(10**((tmu-25)/5),td,rtol=1e-12)
dist=10**((mu-25)/5);xyz=tn*dist[:,None]
radii=[1,3,5]
def vacant(direction,L,positions,radius):
    q=positions@direction
    perp=np.maximum(0,np.sum(positions*positions,axis=1)-q*q)
    valid=perp<=radius*radius
    width=np.sqrt(radius*radius-perp[valid])
    starts=np.maximum(0,q[valid]-width);ends=np.minimum(L,q[valid]+width)
    intervals=sorted((a,b) for a,b in zip(starts,ends) if b>a)
    covered=0.;last=0.
    for a,b in intervals:covered+=max(0,b-max(a,last));last=max(last,b)
    answer=L-covered
    assert -1e-10<=answer<=L+1e-10
    return max(0,answer)
def features(ds,positions):
    return np.array([[vacant(direction,L,positions,r) for direction,L in zip(n,ds)] for r in radii])
central=features(td,xyz)
rng=np.random.default_rng(2026091004);draws=[];draw_distances=[]
for i in range(64):
    ds=10**((tmu+tse*rng.normal(size=len(td))-25)/5)
    pos=tn*10**((mu+se*rng.normal(size=len(mu))-25)/5)[:,None]
    draws.append(features(ds,pos));draw_distances.append(ds)
draws=np.array(draws);draw_distances=np.array(draw_distances)
calibration=features(td*10**(.1/5),xyz*10**(.1/5))
def design(d,v):
    residual=v-d*(d@v)/(d@d)
    ratio=np.linalg.norm(residual)/np.linalg.norm(v)
    X=np.column_stack([d/np.linalg.norm(d),v/np.linalg.norm(v)])
    singular=np.linalg.svd(X,compute_uv=False)
    return float(ratio),float(singular[-1]/singular[0])
assert design(td,td)[0]<1e-14
maxerr=0
for i in range(0,len(td),20):
    grid=(np.arange(20000)+.5)*td[i]/20000
    nearest=np.full(len(grid),np.inf)
    for p in xyz:
        nearest=np.minimum(nearest,np.sum((grid[:,None]*n[i]-p)**2,axis=1))
    for j,r in enumerate(radii):
        err=abs(float(np.mean(nearest>r*r)*td[i])-central[j,i]);maxerr=max(maxerr,err)
        assert err<.02
paths=[];summaries=[]
for j,r in enumerate(radii):
    lo,hi=np.percentile(draws[:,j,:],[16,84],axis=0)
    for i,t in enumerate(targets):paths.append(dict(pgc=t['pgc'],distance_mpc=td[i],radius_mpc=r,
        catalog_no_neighbor_mpc=central[j,i],sampled_p16_mpc=lo[i],sampled_p84_mpc=hi[i],
        common_point1mag_shift_mpc=calibration[j,i]))
    identification=np.array([design(ds,v) for ds,v in zip(draw_distances,draws[:,j,:])])
    beta=(td@central[j])/(td@td)
    projected_draws=draws[:,j,:]-beta*draw_distances
    projected_lo,projected_hi=np.percentile(projected_draws,[16,84],axis=0)
    summaries.append(dict(radius_mpc=r,median_catalog_no_neighbor_fraction=float(np.median(central[j]/td)),
        median_p84_minus_p16_mpc=float(np.median(hi-lo)),
        fixed_distance_projection_beta=float(beta),
        median_projected_p84_minus_p16_mpc=float(np.median(projected_hi-projected_lo)),
        median_absolute_common_shift_mpc=float(np.median(abs(calibration[j]-central[j]))),
        residual_fraction_and_singular_ratio=design(td,central[j]),
        sampled_design_p16_p84=np.percentile(identification,[16,84],axis=0).tolist()))
write_csv('path-sensitivity.csv',paths)
shells=[]
for lower,upper in zip([0,11,25,50,100,200],[11,25,50,100,200,np.inf]):
    count=int(np.sum((dist>=lower)&(dist<upper)))
    shells.append(dict(lower_mpc=lower,upper_mpc=upper if np.isfinite(upper) else None,count=count,
        sampling_density_per_mpc3=count/(4*np.pi/3*(upper**3-lower**3)) if np.isfinite(upper) else None))
summary=dict(source_rows=len(lines),individual_method_counts=dict(Counter(r['method'] for r in measurements)),
    selected_tracers=len(selected),selected_method_counts=dict(Counter(r['method'] for r in selected)),
    selected_distance_min_median_max=np.percentile(dist,[0,50,100]).tolist(),
    selected_linear_radial_sigma_p16_p50_p84_mpc=np.percentile([r['linear_sigma_mpc'] for r in selected],[16,50,84]).tolist(),
    selected_over_25mpc_sigma_above_1mpc=int(sum(r['distance_mpc']>25 and r['linear_sigma_mpc']>1 for r in selected)),
    selected_over_25mpc=sum(r['distance_mpc']>25 for r in selected),
    shells=shells,geometry=summaries,midpoint_max_error_mpc=maxerr,
    interpretation='Feature and uncertainty diagnostics only; not a complete void map, observational fit or fresh holdout.',
    sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,metadata,target_source,HERE/'protocol.md']})
(HERE/'results.json').write_text(json.dumps(summary,indent=2)+'\n',newline='\n')
print(json.dumps(summary,indent=2))
