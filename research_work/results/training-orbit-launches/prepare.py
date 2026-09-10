"""Rebuild provisional orbit seeds using pre-existing association flags only."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache'
OUT=CACHE/'training-orbit-launches';OUT.mkdir(parents=True,exist_ok=True)
SOURCE=HERE.parent/'stellar-orbit-support/run.py'
spec=importlib.util.spec_from_file_location('old_orbit_support',SOURCE)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

paths=dict(moments=CACHE/'stellar-catalogs/stellar-errors-conditional_minus_0.017.parquet',
           parent=CACHE/'stellar-catalogs/matched-with-gaia-covariance.parquet',
           screen=CACHE/'stellar-association-window/training-screen.parquet',
           previous_launches=CACHE/'stellar-orbit-support/expanded-launches.parquet',
           provenance=CACHE/'stellar-orbit-support/stellar-input-provenance.parquet')


def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    hashes={name:digest(path) for name,path in paths.items()}
    d=pd.read_parquet(paths['moments'],filters=[('holdout_role','==','training')])
    assert len(d)==77927 and d.source_id.is_unique and (d.holdout_role=='training').all()
    ids=d.source_id.tolist()
    chem=pd.read_parquet(paths['parent'],columns=['source_id','FE_H','ALPHA_M'],filters=[('source_id','in',ids)])
    d=d.merge(chem,on='source_id',validate='one_to_one').sort_values('source_id')
    d=d[d.mean_R_kpc.between(.5,9)&(d.mean_z_kpc.abs()<=1.5)&(d.FE_H>=-.5)&(d.ALPHA_M<.15)].reset_index(drop=True)
    assert len(d)==27884
    def state(frame):
        y=np.column_stack([frame['mean_'+key].to_numpy(float) for key in old.LABELS])
        y[:,2]-=np.radians(old.P['bar_angle_degrees'])
        return y
    before,strata=old.anchors(state(d),d.source_id.to_numpy())
    original=pd.read_parquet(paths['previous_launches'])
    assert set(d.iloc[before].source_id)==set(original.source_id)
    screen=pd.read_parquet(paths['screen'],columns=['source_id','APOGEE_ID','window_residual_arcsec'])
    d=d.merge(screen,on='source_id',validate='one_to_one')
    prov=pd.read_parquet(paths['provenance'],columns=['source_id','starhorse_used_parallax','parallax_interval_disjoint_5'])
    d=d.merge(prov,on='source_id',validate='one_to_one')
    assert d.source_id.is_monotonic_increasing
    d['seed_eligibility']=np.where(d.window_residual_arcsec.isna(),'unassessed',
                                 np.where(d.window_residual_arcsec>.5,'flagged_position','passes_position_screen'))
    original_audit=original[['source_id']].merge(d[['source_id','APOGEE_ID','window_residual_arcsec','seed_eligibility']],on='source_id',validate='one_to_one')
    eligible=d[d.seed_eligibility=='passes_position_screen'].reset_index(drop=True)
    y=state(eligible)
    picked,strata=old.anchors(y,eligible.source_id.to_numpy())
    x,v=old.cartesian(y[picked])
    seeds=eligible.iloc[picked][['source_id','APOGEE_ID','sky_pixel','FIELD','FE_H','ALPHA_M',
                                'window_residual_arcsec','starhorse_used_parallax','parallax_interval_disjoint_5']].copy()
    for i,name in enumerate(['x_kpc','y_kpc','z_kpc','vx_kms','vy_kms','vz_kms']):seeds[name]=np.c_[x,v][:,i]
    seeds['R_stratum']=[s[0] for s in strata];seeds['height_stratum']=[s[1] for s in strata]
    seeds['holdout_role']='training';seeds['provisional_orbit_seed']=True
    assert len(seeds)==72 and seeds.source_id.is_unique and (seeds.window_residual_arcsec<=.5).all()
    path=OUT/'launches.parquet';seeds.to_parquet(path,index=False)
    # Preserve the full selected training sample and each exclusion reason in a
    # sidecar. This does not delete stars from any observational likelihood.
    side=d[['source_id','seed_eligibility','window_residual_arcsec']].copy()
    side['selected_as_new_seed']=side.source_id.isin(seeds.source_id)
    side.to_parquet(OUT/'seed-eligibility.parquet',index=False)
    records=[]
    for row in seeds.to_dict('records'):
        row={k:(v.decode('utf-8') if isinstance(v,bytes) else v) for k,v in row.items()}
        row['source_id']=str(row['source_id'])
        records.append(row)
    (HERE/'launches.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    oldids=set(original.source_id);newids=set(seeds.source_id)
    oldflag=original_audit[original_audit.seed_eligibility!='passes_position_screen'].copy()
    oldflag['source_id']=oldflag.source_id.astype(str)
    cells=[]
    for (ir,iz),part in seeds.groupby(['R_stratum','height_stratum']):
        speeds=np.linalg.norm(part[['vx_kms','vy_kms','vz_kms']].to_numpy(),axis=1)
        cells.append(dict(R_stratum=int(ir),height_stratum=int(iz),seeds=len(part),
                          positive_z=int((part.z_kpc>0).sum()),negative_z=int((part.z_kpc<0).sum()),
                          maximum_seed_speed_kms=float(speeds.max())))
    result=dict(scope='Provisional training-only orbit launch preparation; no orbit fit or new gravity score.',
        selected_training_stars=len(d),eligibility_counts={str(k):int(v) for k,v in d.seed_eligibility.value_counts().items()},
        original_seed_count=len(original),original_seed_eligibility={str(k):int(v) for k,v in original_audit.seed_eligibility.value_counts().items()},
        original_flagged_or_unassessed_seeds=oldflag.where(pd.notna(oldflag),None).to_dict('records'),
        new_seed_count=len(seeds),retained_original_seeds=len(oldids&newids),newly_selected_seeds=len(newids-oldids),
        dropped_original_seed_ids=[str(i) for i in sorted(oldids-newids)],
        velocity_clipping_applied=False,gravity_residual_selection_applied=False,
        new_seed_parallax_disagreement_count=int(seeds.parallax_interval_disjoint_5.sum()),
        new_seed_starhorse_parallax_used_count=int(seeds.starhorse_used_parallax.sum()),
        cells=cells,source_hashes=hashes,launches_sha256=digest(path),
        previous_72_seeds_exactly_reproduced_before_screen=True,
        approximate_distance_moments_unchanged=True,position_screen_not_proof_of_correct_association=True,
        holdouts_opened=False,orbits_integrated=False,stellar_likelihood_fitted=False,
        code_hashes={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),SOURCE]})
    # NaN association residuals are represented explicitly as null in audit JSON.
    for row in result['original_flagged_or_unassessed_seeds']:
        if pd.isna(row['window_residual_arcsec']):row['window_residual_arcsec']=None
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    for name,path in paths.items():assert digest(path)==hashes[name]
    print(json.dumps({k:result[k] for k in ['eligibility_counts','original_seed_eligibility','retained_original_seeds','newly_selected_seeds','new_seed_parallax_disagreement_count','cells']},indent=2))


if __name__=='__main__':main()
