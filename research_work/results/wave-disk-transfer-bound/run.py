"""Necessary spherical-source mass bound; no new wave-profile fit or photon budget."""
from pathlib import Path
import json,hashlib,zipfile,io
import numpy as np
H=Path(__file__).resolve().parent;ROOT=H.parents[2];R=H.parent
meta_file=ROOT/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt'
rotation_file=ROOT/'temporal_candidate_audit/data/Rotmod_LTG.zip'
old_file=R/'joint-galaxy-audit/galaxy-rotation-predictions.json'
selection=R/'wave-lens-training/parameter-grid/selected-calibration.json'
masscal=R/'stellar-mass-degeneracy/frozen-calibration.json'
fsource=json.loads(selection.read_text())['source_to_stellar_mass_ratio']
lam=json.loads(masscal.read_text())['lambda_mass']
metadata={}
for line in meta_file.read_text().splitlines():
    v=line.split()
    if len(v)!=19:continue
    try:metadata[v[0]]={'distance_Mpc':float(v[2]),'distance_error_Mpc':float(v[3]),'distance_method':int(v[4]),
         'inclination_deg':float(v[5]),'luminosity_1e9_Lsun':float(v[7]),'luminosity_error_1e9_Lsun':float(v[8]),
         'disk_scale_kpc':float(v[11]),'quality':int(v[17])}
    except ValueError:pass
assert len(metadata)==175
old=[r for r in json.loads(old_file.read_text()) if r['model']=='baryons']
previous={(r['galaxy'],r['R_kpc']):r for r in old}
names={r['galaxy'] for r in old};assert len(names)==149 and len(old)==3150
G=4.30091727003628e-6
rows=[];verification=[]
with zipfile.ZipFile(rotation_file) as archive:
    for name in sorted(names):
        meta=metadata[name]
        raw=np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name+'_rotmod.dat'))))
        gas=raw[:,3]*abs(raw[:,3]);star=.5*raw[:,4]*abs(raw[:,4])+.7*raw[:,5]*abs(raw[:,5])
        vbar2=gas+star
        good=np.isfinite(raw).all(axis=1)&(raw[:,0]>0)&(raw[:,1]>0)&(raw[:,2]>0)&(vbar2>0)
        raw=raw[good];gas=gas[good];star=star[good];vbar2=vbar2[good]
        assert len(raw)>=5 and meta['quality']<=2 and meta['inclination_deg']>=30
        for rr,vb in zip(raw,vbar2):verification.append(abs(previous[name,float(rr[0])]['predicted_kms']-np.sqrt(vb)))
        order=np.argsort(raw[:,0]);raw=raw[order];gas=gas[order];star=star[order]
        for scale,label in [(1.,'archived_stellar_baseline'),(lam,'shared_stellar_mass_increase')]:
            # .7 times TOTAL light is a generous upper bound for the adopted .5 disk/.7 bulge split.
            upper_stellar=scale*.7*meta['luminosity_1e9_Lsun']*1e9
            assert upper_stellar>0
            vb2=gas+scale*star;r=raw[:,0];v=raw[:,1];err=raw[:,2]
            required=np.maximum(v*v-vb2,0)*r/G
            required_low=np.maximum(np.maximum(v-err,0)**2-vb2,0)*r/G
            limit=np.sqrt(np.maximum(vb2+G*fsource*upper_stellar/r,0))
            rows.append({'galaxy':name,'historical_split':previous[name,float(r[-1])]['split'],
                'scenario':label,'stellar_mass_scale':scale,**meta,
                'stellar_mass_upper_Msun':upper_stellar,'source_mass_upper_Msun':fsource*upper_stellar,
                'last_radius_kpc':float(r[-1]),'last_observed_speed_km_s':float(v[-1]),'last_speed_error_km_s':float(err[-1]),
                'last_baryonic_speed_km_s':float(np.sqrt(max(vb2[-1],0))),
                'last_maximum_speed_with_source_km_s':float(limit[-1]),
                'minimum_f_from_last_radius':float(required[-1]/upper_stellar),
                'minimum_f_from_last_radius_speed_minus_one_error':float(required_low[-1]/upper_stellar),
                'minimum_f_from_all_radii':float(np.max(required)/upper_stellar),
                'fails_last_radius_necessary_bound':bool(required[-1]>fsource*upper_stellar),
                'fails_last_radius_after_speed_minus_one_error':bool(required_low[-1]>fsource*upper_stellar),
                'fails_any_radius_necessary_bound':bool(np.max(required)>fsource*upper_stellar)})
assert len(verification)==3150 and max(verification)<1e-9
summary={'classification':'Conditional transfer of lens f to exposed SPARC data; generous spherical mass bound, not a wave solution',
         'f_from_lens_training':fsource,'galaxies':149,'radial_rows_verified':len(verification),
         'max_archived_baryonic_speed_difference_km_s':max(verification),'scenarios':{},
         'input_sha256':{str(file.relative_to(ROOT)):hashlib.sha256(file.read_bytes()).hexdigest() for file in [meta_file,rotation_file,old_file,selection,masscal]}}
for label in ['archived_stellar_baseline','shared_stellar_mass_increase']:
    ss=[r for r in rows if r['scenario']==label]
    groups={}
    for group,lo,hi in [('L_under_1e9',0,1),('L_1e9_to_1e10',1,10),('L_at_least_1e10',10,float('inf'))]:
        sub=[r for r in ss if lo<=r['luminosity_1e9_Lsun']<hi]
        groups[group]={'n':len(sub),'fail_last_radius':sum(r['fails_last_radius_necessary_bound'] for r in sub),
                       'median_minimum_f_last_radius':float(np.median([r['minimum_f_from_last_radius'] for r in sub]))}
    summary['scenarios'][label]={'fail_last_radius':sum(r['fails_last_radius_necessary_bound'] for r in ss),
        'fail_last_radius_speed_minus_one_error':sum(r['fails_last_radius_after_speed_minus_one_error'] for r in ss),
        'fail_any_radius':sum(r['fails_any_radius_necessary_bound'] for r in ss),
        'median_minimum_f_last_radius':float(np.median([r['minimum_f_from_last_radius'] for r in ss])),
        'luminosity_groups':groups,
        'cepheid_or_trgb_distance_subset':{'n':sum(r['distance_method'] in [2,3] for r in ss),
            'fail_last_radius':sum(r['distance_method'] in [2,3] and r['fails_last_radius_necessary_bound'] for r in ss)},
        'historical_roles':{split:{'n':sum(r['historical_split']==split for r in ss),
             'fail_last_radius':sum(r['historical_split']==split and r['fails_last_radius_necessary_bound'] for r in ss)} for split in ['train','validation','test']}}
for file,data in [('results.json',summary),('galaxy-bounds.json',rows)]:
    (H/file).write_text(json.dumps(data,indent=2)+'\n',newline='\n')
print(json.dumps(summary,indent=2))
