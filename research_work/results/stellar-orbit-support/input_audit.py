"""Recover omitted input provenance and flag distance/parallax disagreement.

No distances are replaced, no stars are dropped, and no held-out velocities
are scored. Overlap flags are diagnostics, not independent sigma tests.
"""
import numpy as np
import pandas as pd
from astropy.table import Table
from run import CACHE,OUT,HERE,save,digest

def clean(value):
    if pd.isna(value):return ''
    return value.decode().strip() if isinstance(value,bytes) else str(value).strip()

source=CACHE/'APOGEE_DR17_EDR3_StarHorse_v2.fits'
sh=Table.read(source)[['APOGEE_ID','DR3_source_id','StarHorse_INPUTFLAGS','dist05','dist95']].to_pandas()
sh=sh.rename(columns={'DR3_source_id':'source_id'})
sh['APOGEE_ID']=sh.APOGEE_ID.map(clean)
sh['StarHorse_INPUTFLAGS']=sh.StarHorse_INPUTFLAGS.map(clean)
columns=['source_id','APOGEE_ID','gaia_quality_candidate','dist16','dist50','dist84','parallax','parallax_error','ruwe','FE_H','ALPHA_M']
parent=pd.read_parquet(CACHE/'matched-with-gaia-covariance.parquet',columns=columns)
parent['APOGEE_ID']=parent.APOGEE_ID.map(clean)
metadata=parent.merge(sh,on=['APOGEE_ID','source_id'],how='left',validate='one_to_one')
assert len(metadata)==len(parent) and metadata.dist05.notna().all()
metadata['starhorse_used_parallax']=metadata.StarHorse_INPUTFLAGS.str.contains('PARALLAX',case=False,regex=False)
assert (metadata.dist05>0).all() and (metadata.dist95>=metadata.dist05).all()
metadata['implied_parallax_low_mas']=1/metadata.dist95-.037
metadata['implied_parallax_high_mas']=1/metadata.dist05+.003
for width in [3,5,10]:
    metadata[f'parallax_interval_disjoint_{width}']=((metadata.parallax-width*metadata.parallax_error>metadata.implied_parallax_high_mas)|
                    (metadata.parallax+width*metadata.parallax_error<metadata.implied_parallax_low_mas))
fidelity=Table.read(OUT/'gaia-astrometric-fidelity.fits').to_pandas()
metadata=metadata.merge(fidelity,on='source_id',how='left',validate='one_to_one')
roles=pd.read_parquet(CACHE/'stellar-spatial-holdouts.parquet',columns=['source_id','holdout_role','sky_pixel'])
metadata=metadata.merge(roles,on='source_id',how='left',validate='one_to_one')
assert metadata.loc[metadata.gaia_quality_candidate,'holdout_role'].notna().all()
metadata['distance_processing_route']=np.where(metadata.starhorse_used_parallax,
    'included_parallax_joint_approximation','omitted_parallax_needs_separate_likelihood_review')
metadata.loc[metadata.parallax_interval_disjoint_5,'distance_processing_route']='distance_parallax_disagreement_review'
metadata.loc[(metadata.fidelity_v1<.5)|(metadata.fidelity_v2<.5),'distance_processing_route']='astrometric_fidelity_review'
metadata.loc[metadata.fidelity_v1.isna()|metadata.fidelity_v2.isna(),'distance_processing_route']='missing_fidelity_review'
metadata['inference_review_required']=metadata.distance_processing_route.ne('included_parallax_joint_approximation')
metadata.to_parquet(OUT/'stellar-input-provenance.parquet',index=False)
training=pd.read_parquet(CACHE/'stellar-errors-conditional_minus_0.017.parquet',filters=[('holdout_role','==','training')])
d=training.merge(metadata.drop(columns=['holdout_role','sky_pixel']),on='source_id',validate='one_to_one')
selected=d.mean_R_kpc.between(.5,9)&(d.mean_z_kpc.abs()<=1.5)&(d.FE_H>=-.5)&(d.ALPHA_M<.15)
summary={'parent_sources':len(metadata),'quality_candidate_sources':int(metadata.gaia_quality_candidate.sum()),
         'metadata_file_sha256':digest(OUT/'stellar-input-provenance.parquet'),
         'raw_starhorse_sha256':digest(source),'samples':{}}
for name,part in [('all_training',d),('orbit_support_training',d[selected])]:
    rows=[]
    for used in [False,True]:
        s=part[part.starhorse_used_parallax==used]
        rows.append(dict(starhorse_used_parallax=used,stars=len(s),
            fidelity_v1_below_half=int((s.fidelity_v1<.5).sum()),
            fidelity_v2_below_half=int((s.fidelity_v2<.5).sum()),
            interval_disjoint={str(n):int(s[f'parallax_interval_disjoint_{n}'].sum()) for n in [3,5,10]}))
    summary['samples'][name]=dict(stars=len(part),groups=rows)
summary['selected_training_processing_routes']={str(k):int(v) for k,v in d[selected].distance_processing_route.value_counts().items()}
region_masks={'bulge_plane':(d.mean_R_kpc<3.5)&(d.mean_z_kpc.abs()<.2),
              'bulge_offplane':(d.mean_R_kpc<3.5)&d.mean_z_kpc.abs().between(.5,1.5),
              'disk_plane':(d.mean_R_kpc>=5)&(d.mean_z_kpc.abs()<.2),
              'disk_offplane':(d.mean_R_kpc>=5)&d.mean_z_kpc.abs().between(.5,1.5)}
summary['selected_training_regions']={name:dict(stars=int((selected&mask).sum()),
    starhorse_parallax_used=int(d.loc[selected&mask,'starhorse_used_parallax'].sum())) for name,mask in region_masks.items()}
# Record the four high-speed training representatives that prompted the audit.
launch=pd.read_parquet(OUT/'expanded-launches.parquet')
audited=launch.merge(metadata[['source_id','StarHorse_INPUTFLAGS','starhorse_used_parallax','fidelity_v1','fidelity_v2',
    'parallax_interval_disjoint_5','distance_processing_route','inference_review_required']],on='source_id',validate='one_to_one')
audited.to_parquet(OUT/'expanded-launches-audited.parquet',index=False)
summary['launches_requiring_input_review']=int(audited.inference_review_required.sum())
summary['audited_launches_sha256']=digest(OUT/'expanded-launches-audited.parquet')
speed=np.linalg.norm(launch[['vx_kms','vy_kms','vz_kms']].to_numpy(),axis=1)
extreme=launch.loc[speed>500,['source_id']].merge(d,on='source_id',validate='one_to_one')
examples=[]
for _,r in extreme.iterrows():
    v=np.array([r.mean_vR_kms,r.mean_vphi_kms,r.mean_vz_kms]);speed=float(np.linalg.norm(v))
    examples.append(dict(source_id=str(int(r.source_id)),APOGEE_ID=r.APOGEE_ID,
        StarHorse_INPUTFLAGS=r.StarHorse_INPUTFLAGS,starhorse_used_parallax=bool(r.starhorse_used_parallax),
        distance_median_kpc=float(r.dist50),distance_05_95_kpc=[float(r.dist05),float(r.dist95)],
        gaia_parallax_mas=float(r.parallax),gaia_parallax_error_mas=float(r.parallax_error),
        inverse_raw_parallax_kpc=float(1/r.parallax),
        implied_tangential_scale_ratio=float(r.dist50*r.parallax),
        fidelity_v1=float(r.fidelity_v1),fidelity_v2=float(r.fidelity_v2),
        prior_approximate_galactocentric_speed_kms=speed,
        disjoint_at_5_error_widths=bool(r.parallax_interval_disjoint_5)))
summary['fast_launch_examples']=examples
summary['interpretation']='An overlap/provenance audit, not a formal sigma significance test. Distances may share parallax inputs and have non-Gaussian systematic errors. High fidelity does not prove the crossmatch or the distance model is correct.'
summary['raw_catalogs_or_holdout_roles_changed']=False
summary['new_distance_estimates_computed']=False
save('input-audit.json',summary)
print(summary['samples'],flush=True)
