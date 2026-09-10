"""Check new paths and compare coverage on identical non-launch training stars."""
from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('expansion_driver', HERE/'run.py')
run = importlib.util.module_from_spec(spec); spec.loader.exec_module(run)
cov, CACHE, OUT, ROOT = run.cov, run.CACHE, run.OUT, run.ROOT

def main():
    integration = json.loads((HERE/'integration.json').read_text())
    assert integration['completed']
    selection = json.loads((HERE/'selection.json').read_text())
    assert selection['protocol_sha256'] == run.digest(HERE/'protocol.md')
    for path, sha in selection['input_hashes'].items(): assert run.digest(ROOT/path) == sha
    assert integration['driver_sha256'] == run.digest(HERE.parent/'full-bar-orbits/run.py')
    assert integration['launches_sha256'] == run.digest(OUT/'launches.parquet')
    assert integration['code_sha256'] == run.digest(HERE/'run.py')
    for path, sha in integration['input_field_hashes'].items(): assert run.digest(ROOT/path) == sha
    seeds = pd.read_parquet(OUT/'launches.parquet')
    full = run.load('expansion_field_audit', HERE.parent/'full-bar-completion/run.py')
    field_paths = [full.CACHE/name for name in ['axisymmetric-reference-80-L128-N512.npz','fine.npz','finer.npz']]
    axis, lower, higher = [full.load_field(p) for p in field_paths]
    checks, accepted, accepted_ids = [], [], []
    for row in integration['rows']:
        j = row['seed_index']; path = OUT/f'orbit-{j}.npz'
        item = dict(seed_index=j, numerical_pass=row['numerical_pass'], field_pass=False)
        if row['error'] is None:
            with np.load(path) as f:
                states = f['trajectory']; times = f['times']
                assert str(f['source_id'].item()) == str(seeds.iloc[j].source_id)
            np.testing.assert_array_equal(times, run.TIMES)
            np.testing.assert_array_equal(states[0], seeds.iloc[j][run.NAMES].to_numpy(float))
            assert states.shape == (501,6) and np.isfinite(states).all()
            errors = []
            for q in np.array_split(states[:, :3], 32):
                a, b, c = axis.evaluate(q)[1], lower.evaluate(q)[1], higher.evaluate(q)[1]
                errors.extend((np.linalg.norm(c-b,axis=1)/np.linalg.norm(a+c,axis=1)).tolist())
            item.update(maximum_extra_force_fraction=float(max(errors)), field_pass=max(errors)<.01,
                        radius_min_kpc=float(np.linalg.norm(states[:, :3],axis=1).min()),
                        radius_max_kpc=float(np.linalg.norm(states[:, :3],axis=1).max()),
                        path_sha256=run.digest(path))
            if row['numerical_pass'] and item['field_pass']:
                accepted.append(states); accepted_ids.append(seeds.iloc[j].source_id)
        checks.append(item)
        print('Field', j+1, '/', len(seeds), item['field_pass'], flush=True)
    run.save('field-check.json', dict(rows=checks, threshold_fraction=.01,
             input_hashes={str(p.relative_to(ROOT)):run.digest(p) for p in field_paths},
             passing_both=len(accepted), holdouts_opened=False))
    assert accepted, 'No passing added paths; retain failure record.'
    new = np.stack(accepted, axis=1)
    old_path = CACHE/'training-orbit-launches/full-field-orbits.npz'
    old_check = json.loads((HERE.parent/'training-orbit-launches/field-check.json').read_text())
    assert run.digest(old_path) == old_check['input_hashes'][str(old_path.relative_to(ROOT))]
    old_integrations = json.loads((HERE.parent/'training-orbit-launches/integration-results.json').read_text())
    assert all(r['numerical_checks_pass'] for r in old_integrations['rows'])
    assert all(r['field_gate_pass'] for r in old_check['rows'])
    with np.load(old_path) as f:
        old = f['trajectories']; old_ids = f['source_id']
        np.testing.assert_array_equal(f['times'], run.TIMES)
    combined = np.concatenate([old,new],axis=1)
    np.savez_compressed(OUT/'combined-passing-orbits.npz', times=run.TIMES, trajectories=combined,
                        source_id=np.r_[old_ids,np.array(accepted_ids,dtype=old_ids.dtype)])
    d = run.sample()
    # Exclude even failed new launch stars to keep the evaluation set identical.
    exclude = set(old_ids.tolist()) | set(seeds.source_id.tolist())
    d = d[(d.seed_eligibility=='passes_position_screen') & ~d.source_id.isin(exclude)].reset_index(drop=True)
    assert len(d)==27642-len(seeds) and not d.source_id.isin(exclude).any()
    x, v = d[run.NAMES[:3]].to_numpy(), d[run.NAMES[3:]].to_numpy()
    outputs = {}
    windows = [('all_positive',slice(1,None)),('half_sampling',slice(2,None,2)),('late_half',slice(251,None))]
    for name, window in windows:
        before = cov.coverage(x,v,old[window].reshape(-1,6))
        after = cov.coverage(x,v,combined[window].reshape(-1,6))
        assert np.all(after[0]<=before[0]) and np.all(after[1]>=before[1])
        outputs[name] = (before,after)
        print('Coverage',name, 'old',int((before[0][:,1]<=50).sum()),
              'expanded',int((after[0][:,1]<=50).sum()),'of',len(d),flush=True)
    base = outputs['all_positive']
    archive = pd.read_parquet(CACHE/'orbit-population-coverage/training-coverage.parquet')
    archive = d[['source_id']].merge(archive,on='source_id',validate='one_to_one')
    for j,r in enumerate(cov.RADII):
        np.testing.assert_array_equal(base[0][0][:,j],archive[f'min_velocity_gap_within_{r}_kpc'].to_numpy())
        np.testing.assert_array_equal(base[0][1][:,j],archive[f'samples_within_{r}_kpc'].to_numpy())
    for window in ['half_sampling','late_half']:
        for model in [0,1]:
            assert np.all(outputs[window][model][0]>=base[model][0])
            assert np.all(outputs[window][model][1]<=base[model][1])
    flat = combined[1:].reshape(-1,6)
    for i in np.linspace(0,len(d)-1,12,dtype=int):
        dx=np.linalg.norm(flat[:,:3]-x[i],axis=1); dv=np.linalg.norm(flat[:,3:]-v[i],axis=1)
        for j,r in enumerate(cov.RADII):
            mask=dx<=r; expected=dv[mask].min() if mask.any() else np.inf
            assert base[1][0][i,j]==expected and base[1][1][i,j]==mask.sum()
    cells=[]
    for lo,hi in [(.5,3.5),(3.5,5),(5,9)]:
        for zl,zh in [(0,.2),(.2,.5),(.5,1.5)]:
            mask=(d.mean_R_kpc>=lo)&((d.mean_R_kpc<=hi) if hi==9 else (d.mean_R_kpc<hi))
            mask&=(abs(d.mean_z_kpc)>=zl)&((abs(d.mean_z_kpc)<=zh) if zh==1.5 else (abs(d.mean_z_kpc)<zh))
            for sign in ['both','below','above']:
                selected=mask.copy()
                if sign=='below': selected&=d.mean_z_kpc<0
                if sign=='above': selected&=d.mean_z_kpc>=0
                row=dict(R_range=[lo,hi],absolute_z_range=[zl,zh],side=sign,stars=int(selected.sum()))
                for name, pair in outputs.items():
                    row[name]={model:int((vals[0][selected,1]<=50).sum()) for model,vals in zip(['old','expanded'],pair)}
                cells.append(row)
    all_mask=np.ones(len(d),bool)
    result=dict(scope='Training finite-orbit coverage only; no likelihood, gravity fit or holdout score.',
        target_stars=len(d), old_orbits=old.shape[1], proposed_new_orbits=len(seeds),
        passing_new_orbits=len(accepted), expanded_orbits=combined.shape[1],
        all_launch_stars_excluded=True, launch_samples_excluded=True,
        windows={name:{label:cov.summary(all_mask,*vals) for label,vals in zip(['old','expanded'],pair)} for name,pair in outputs.items()},
        cells=cells, old_coverage_exactly_reproduced=True, brute_force_controls=12,
        cutoffs_are_diagnostic_not_error_bars=True, holdouts_opened=False,
        selected_population_likelihood=False, orbital_stationarity_established=False,
        matched_ordinary_matter_library=False,
        hashes={str(p.relative_to(ROOT)):run.digest(p) for p in [HERE/'protocol.md',HERE/'run.py',Path(__file__),
           HERE/'selection.json',HERE/'integration.json',HERE/'field-check.json',old_path,OUT/'combined-passing-orbits.npz']})
    run.save('results.json',result)
    audit=d[['source_id']].copy()
    for name,pair in outputs.items():
        for label,vals in zip(['old','expanded'],pair):audit[f'{name}_{label}_covered']=vals[0][:,1]<=50
    audit.to_parquet(OUT/'target-coverage.parquet',index=False)
    print('Assessment complete',flush=True)

if __name__=='__main__': main()
