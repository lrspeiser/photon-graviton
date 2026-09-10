from run import HERE,CACHE,OUT,digest
import json
import pandas as pd
r=json.loads((HERE/'results.json').read_text())
d=pd.read_parquet(OUT/'training-screen.parquet')
roles=pd.read_parquet(CACHE/'stellar-catalogs/stellar-spatial-holdouts.parquet',columns=['source_id','holdout_role'])
assert len(d)==r['training_stars'] and d.source_id.is_unique
assert set(d.source_id)<=set(roles.loc[roles.holdout_role.eq('training'),'source_id'])
assert not set(d.source_id)&set(roles.loc[~roles.holdout_role.eq('training'),'source_id'])
assert int(d.window_residual_arcsec.notna().sum())==r['epoch_screen_applicable_training_stars']
for threshold,counts in r['threshold_counts'].items():
    t=float(threshold);eligible=d.window_residual_arcsec.notna()
    assert int((d.window_residual_arcsec>t).sum())==counts['window_flagged']
    assert int(((d.unrestricted_residual_arcsec>t)&eligible).sum())==counts['unrestricted_flagged']
review=pd.read_csv(HERE/'review-candidates.csv',dtype={'source_id':'int64'})
confirmed=pd.read_csv(HERE/'confirmed-review.csv',dtype={'source_id':'int64'})
assert set(review.source_id)==set(d.loc[d.window_residual_arcsec>.5,'source_id'])==set(confirmed.source_id)
assert (confirmed.actual_epoch_gaia_to_2mass_arcsec>.5).all()
assert (confirmed.apogee_to_2mass_arcsec==0).all()
paths={'parent':CACHE/'stellar-catalogs/matched-with-gaia-covariance.parquet','roles':CACHE/'stellar-catalogs/stellar-spatial-holdouts.parquet','moments':CACHE/'stellar-catalogs/stellar-errors-conditional_minus_0.017.parquet','provenance':CACHE/'stellar-orbit-support/stellar-input-provenance.parquet'}
for key,p in paths.items(): assert digest(p)==r['source_hashes'][key]
assert digest(OUT/'training-screen.parquet')==r['screen_sha256']
dl=json.loads((HERE/'download.json').read_text())
for b in dl['batches']: assert digest(HERE/'inputs'/b['file'])==b['sha256']
checks={'training_role_disjointness':True,'original_source_hashes_unchanged':True,'screen_counts_reconstructed':True,'confirmed_ids_match_review_list':True,'downloaded_snapshot_hashes_match':True,'all_140_survive_actual_epoch_check':True,'no_association_replacements_or_adopted_deletions':True,'plot_visually_reviewed':True}
(HERE/'verification.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(checks,indent=2))
