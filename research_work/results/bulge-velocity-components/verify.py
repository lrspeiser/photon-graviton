from run import HERE,CACHE,OUT,sha
import json
import numpy as np
import pandas as pd
r=json.loads((HERE/'results.json').read_text())
d=pd.read_parquet(OUT/'training-components.parquet')
assert len(d)==27884 and d.source_id.is_unique and d.holdout_role.eq('training').all()
roles=pd.read_parquet(CACHE/'stellar-catalogs/stellar-spatial-holdouts.parquet',columns=['source_id','holdout_role'])
assert set(d.source_id)<=set(roles.loc[roles.holdout_role.eq('training'),'source_id'])
np.testing.assert_allclose(d.vz_median_distance,d.rv_vertical+d.pm_vertical+7.25,rtol=0,atol=1e-9)
assert r['coordinate_sum_max_error_kms']<1e-9
masks={'bulge_plane':(d.mean_R_kpc<3.5)&(d.mean_z_kpc.abs()<.2),
       'bulge_offplane':(d.mean_R_kpc<3.5)&d.mean_z_kpc.abs().between(.5,1.5),
       'disk_plane':(d.mean_R_kpc>=5)&(d.mean_z_kpc.abs()<.2),
       'disk_offplane':(d.mean_R_kpc>=5)&d.mean_z_kpc.abs().between(.5,1.5)}
for name,mask in masks.items():
    q=d[mask]
    for kind,part in [('all',q),('unflagged_or_unassessed',q[~(q.window_residual_arcsec>.5)])]:
        expected=r['regions'][name][kind]
        assert len(part)==expected['stars']
        for key,value in expected['stds_kms'].items():np.testing.assert_allclose(part[key].astype(float).std(ddof=0),value,rtol=1e-12)
        x=expected['vertical_variance_decomposition_kms2']
        np.testing.assert_allclose(x['rv']+x['proper_motion_distance']+x['twice_cross_covariance'],x['total'],rtol=1e-12)
assert sum(c['flagged']['stars'] for c in r['common_bulge_cells'])==23
assert sum(c['other']['stars'] for c in r['common_bulge_cells'])==74
paths={'parent':CACHE/'stellar-catalogs/matched-with-gaia-covariance.parquet','moments':CACHE/'stellar-catalogs/stellar-errors-conditional_minus_0.017.parquet','association_screen':CACHE/'stellar-association-window/training-screen.parquet'}
for k,p in paths.items():assert sha(p)==r['source_hashes'][k]
result={'training_identifiers_only':True,'velocity_component_sum':True,'per_region_spreads_reconstructed':True,'variance_covariance_identity':True,'common_cell_counts':True,'input_hashes_unchanged':True,'final_holdout_opened':False}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
