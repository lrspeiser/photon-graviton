"""Ensure the provenance repair preserves identifiers, roles and original inputs."""
import json
import numpy as np
import pandas as pd
from astropy.table import Table
from run import CACHE,OUT,HERE,save,digest

protocol=json.loads((HERE.parent/'stellar-holdout-audit/protocol.json').read_text())
assert digest(CACHE/protocol['parent_catalog'])==protocol['parent_sha256']
parent=pd.read_parquet(CACHE/protocol['parent_catalog'],columns=['source_id','APOGEE_ID'])
roles=pd.read_parquet(CACHE/'stellar-spatial-holdouts.parquet',columns=['source_id','holdout_role'])
m=pd.read_parquet(OUT/'stellar-input-provenance.parquet')
assert len(m)==len(parent) and m.source_id.is_unique and m.source_id.dtype.kind in 'iu'
assert set(m.source_id)==set(parent.source_id)
check=roles.merge(m[['source_id','holdout_role']],on='source_id',validate='one_to_one',suffixes=('_old','_new'))
assert (check.holdout_role_old==check.holdout_role_new).all()
assert (m.starhorse_used_parallax==m.StarHorse_INPUTFLAGS.str.contains('PARALLAX',case=False,regex=False)).all()
assert (m.parallax_interval_disjoint_10<=m.parallax_interval_disjoint_5).all()
assert (m.parallax_interval_disjoint_5<=m.parallax_interval_disjoint_3).all()
f=Table.read(OUT/'gaia-astrometric-fidelity.fits').to_pandas()
assert set(f.source_id)==set(parent.source_id)
a=pd.read_parquet(OUT/'expanded-launches-audited.parquet')
assert len(a)==72 and a.source_id.is_unique and (a.holdout_role=='training').all()
assert set(a.source_id).issubset(set(roles.loc[roles.holdout_role=='training','source_id']))
save('input-verification.json',dict(original_parent_hash_unchanged=True,role_assignments_unchanged=True,
    all_metadata_identifiers_unique_and_exact=True,fidelity_sources=len(f),
    audited_training_launches=len(a),inference_review_required=int(a.inference_review_required.sum()),
    distance_or_velocity_values_replaced=False,heldout_kinematic_scores_evaluated=False))
print('Input provenance, exact identifiers and original holdout roles verified.')
