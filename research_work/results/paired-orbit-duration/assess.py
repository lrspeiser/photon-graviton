"""Audit paired long paths and equal-state-count coverage diagnostics."""
from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd
import astropy.units as u

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('paired_duration_driver',HERE/'run.py')
run=importlib.util.module_from_spec(spec);spec.loader.exec_module(run)
ROOT,CACHE,OUT=run.ROOT,run.CACHE,run.OUT

def main():
    info={kind:json.loads((HERE/f'{kind}-integration.json').read_text()) for kind in ['ordinary','full']}
    trajectories={};checks=[];hashes={}
    seed_table=pd.read_parquet(run.SOURCE)
    unit_myr=(1*u.kpc/(u.km/u.s)).to_value(u.Myr)
    phase_scales=[]
    for index in run.INDICES:
        seed=seed_table.iloc[index]
        R2=seed.x_kpc**2+seed.y_kpc**2
        relative=(seed.x_kpc*seed.vy_kms-seed.y_kpc*seed.vx_kms)/R2-run.orbit.OMEGA
        phase_scales.append(dict(seed_index=index,radius_kpc=float(np.sqrt(R2)),
            instantaneous_relative_angular_rate_kms_per_kpc=float(relative),
            constant_rate_phase_scale_Myr=float(2*np.pi/abs(relative)*unit_myr) if relative else None))
    for kind,data in info.items():
        assert data['completed'] and data['indices']==run.INDICES
        for p,sha in data['source_hashes'].items():assert run.digest(ROOT/p)==sha
        hashes[str((HERE/f'{kind}-integration.json').relative_to(ROOT))]=run.digest(HERE/f'{kind}-integration.json')
        for row in data['rows']:
            if row['complete_trajectory']:
                path=OUT/f"{kind}-{row['seed_index']}.npz"
                assert run.digest(path)==row['trajectory_sha256']
                with np.load(path) as f:
                    y=f['trajectory'];np.testing.assert_array_equal(f['times'],run.TIMES)
                    assert str(f['source_id'].item())==row['source_id']
                assert np.isfinite(y).all() and y.shape==(2001,6)
                np.testing.assert_array_equal(y[0],seed_table.iloc[row['seed_index']][run.exp.NAMES].to_numpy(float))
                trajectories[kind,row['seed_index']]=y
                hashes[str(path.relative_to(ROOT))]=run.digest(path)
    full=run.orbit.full
    fields=[full.CACHE/name for name in ['axisymmetric-reference-80-L128-N512.npz','fine.npz','finer.npz']]
    axis,lo,hi=[full.load_field(p) for p in fields]
    for p in fields:hashes[str(p.relative_to(ROOT))]=run.digest(p)
    for index in run.INDICES:
        row=next(x for x in info['full']['rows'] if x['seed_index']==index)
        check=dict(seed_index=index,field_pass=False,short_path_pass=False)
        if row['complete_trajectory']:
            y=trajectories['full',index]
            errors=[]
            for q in np.array_split(y[:,:3],128):
                a,b,c=axis.evaluate(q)[1],lo.evaluate(q)[1],hi.evaluate(q)[1]
                errors.extend((np.linalg.norm(c-b,axis=1)/np.linalg.norm(a+c,axis=1)).tolist())
            check.update(maximum_extra_force_difference=float(max(errors)),field_pass=max(errors)<.01)
            shortpath=CACHE/f'orbit-library-expansion/orbit-{index}.npz'
            shortaudit=json.loads((HERE.parent/'orbit-library-expansion/field-check.json').read_text())
            expected=next(x for x in shortaudit['rows'] if x['seed_index']==index)['path_sha256']
            assert run.digest(shortpath)==expected
            hashes[str(shortpath.relative_to(ROOT))]=expected
            with np.load(shortpath) as f:short=f['trajectory']
            dx=float(np.linalg.norm(y[:501,:3]-short[:,:3],axis=1).max())
            dv=float(np.linalg.norm(y[:501,3:]-short[:,3:],axis=1).max())
            check.update(short_path_difference_kpc=dx,short_path_difference_kms=dv,short_path_pass=dx<1e-4 and dv<.01)
            # Independent occupancy counts verify the inherited bin convention.
            pos=y[1:501,:3]
            cols=[np.hypot(pos[:,0],pos[:,1]),pos[:,2],np.arctan2(pos[:,1],pos[:,0])]
            bins=[run.old.R_EDGES,run.old.Z_EDGES,run.old.PHI_EDGES]
            idx=[np.digitize(v,b)-1 for v,b in zip(cols,bins)]
            measured=np.bincount(np.ravel_multi_index(idx,(8,8,8)),minlength=512)/len(pos)
            np.testing.assert_array_equal(measured,run.old.occupation(y[1:501]))
        checks.append(check)
        print('Field and short path',index,check['field_pass'],check['short_path_pass'],flush=True)
    pairs=[];excluded=[]
    for index in run.INDICES:
        statuses={kind:next(x for x in data['rows'] if x['seed_index']==index)['numerical_pass'] for kind,data in info.items()}
        extra=next(x for x in checks if x['seed_index']==index)
        good=all(statuses.values()) and extra['field_pass'] and extra['short_path_pass']
        (pairs if good else excluded).append(index)
    run.save('checks.json',dict(rows=checks,common_passing_pairs=pairs,excluded_pairs=excluded,input_hashes=hashes,
        ordinary_field_uncertainty='Same previously audited components; no new ordinary-field resolution certification.',holdouts_opened=False))
    assert pairs, 'No common passing pairs; checks retained, coverage not evaluated.'
    selection=json.loads((HERE.parent/'orbit-library-expansion/selection.json').read_text())
    for p,sha in selection['input_hashes'].items():
        assert run.digest(ROOT/p)==sha
        hashes[p]=sha
    d=run.exp.sample()
    launches=pd.read_parquet(CACHE/'orbit-library-expansion/launches.parquet')
    oldlaunch=pd.read_parquet(CACHE/'training-orbit-launches/launches.parquet')
    oldlaunchpath=CACHE/'training-orbit-launches/launches.parquet'
    hashes[str(oldlaunchpath.relative_to(ROOT))]=run.digest(oldlaunchpath)
    exclude=set(launches.source_id)|set(oldlaunch.source_id)
    d=d[(d.seed_eligibility=='passes_position_screen') & ~d.source_id.isin(exclude)].reset_index(drop=True)
    assert len(d)==27606 and (d.holdout_role=='training').all()
    assert not d.source_id.isin(exclude).any()
    x=d[run.exp.NAMES[:3]].to_numpy();v=d[run.exp.NAMES[3:]].to_numpy()
    windows=dict(early_equal250=np.arange(1,251),middle_equal250=np.arange(876,1126),late_equal250=np.arange(1751,2001),
        first_quarter500=np.arange(1,501),full_history_equal500=np.arange(4,2001,4),full_history2000=np.arange(1,2001))
    tables=[];gap_records={};cell_rows=[]
    for kind in ['ordinary','full']:
        states=np.stack([trajectories[kind,j] for j in pairs],axis=1)
        for name,indices in windows.items():
            flat=states[indices].reshape(-1,6)
            gaps,counts=run.exp.cov.coverage(x,v,flat)
            covered=gaps[:,1]<=50
            gap_records[kind,name]=gaps
            tables.append(dict(kind=kind,window=name,orbits=len(pairs),saved_times_per_orbit=len(indices),
                first_time_kpc_per_kms=float(run.TIMES[indices[0]]),last_time_kpc_per_kms=float(run.TIMES[indices[-1]]),
                target_stars=len(d),covered=int(covered.sum()),fraction=float(covered.mean()),spatially_covered=int((counts[:,1]>0).sum())))
            for i in np.linspace(0,len(d)-1,6,dtype=int):
                dx=np.linalg.norm(flat[:,:3]-x[i],axis=1);dv=np.linalg.norm(flat[:,3:]-v[i],axis=1)
                inside=dx<=.5;expected=dv[inside].min() if inside.any() else np.inf
                assert gaps[i,1]==expected and counts[i,1]==inside.sum()
            for loR,hiR in [(.5,3.5),(3.5,5),(5,9)]:
                for loz,hiz in [(0,.2),(.2,.5),(.5,1.5)]:
                    mask=(d.mean_R_kpc>=loR)&((d.mean_R_kpc<=hiR) if hiR==9 else (d.mean_R_kpc<hiR))
                    mask&=(abs(d.mean_z_kpc)>=loz)&((abs(d.mean_z_kpc)<=hiz) if hiz==1.5 else (abs(d.mean_z_kpc)<hiz))
                    for side in ['both','below','above']:
                        sel=mask.copy()
                        if side=='below':sel&=d.mean_z_kpc<0
                        if side=='above':sel&=d.mean_z_kpc>=0
                        cell_rows.append(dict(kind=kind,window=name,R_range=[loR,hiR],absolute_z_range=[loz,hiz],side=side,
                            stars=int(sel.sum()),covered=int((covered & sel).sum())))
            print('Coverage',kind,name,int(covered.sum()),'/',len(d),flush=True)
        for name in windows:
            assert np.all(gap_records[kind,'full_history2000']<=gap_records[kind,name])
    audit=d[['source_id']].copy()
    for (kind,name),gaps in gap_records.items():audit[kind+'_'+name]=gaps[:,1]<=50
    audit.to_parquet(OUT/'coverage.parquet',index=False)
    result=dict(scope='Common-passing training duration/coverage probes; not gravity ranking or a selected likelihood.',
        target_stars=len(d),selected_indices=run.INDICES,common_passing_pairs=pairs,excluded_pairs=excluded,
        tables=tables,cells=cell_rows,coverage_thresholds=dict(position_kpc=.5,velocity_kms=50),
        fixed_state_count_controls=True,all_108_launch_stars_excluded=True,holdouts_opened=False,
        orbit_population_stationarity_proved=False,selected_by_candidate_training_gaps=True,
        numerical_subset_is_not_unbiased=True,source_hashes=hashes,
        initial_phase_scales=phase_scales,
        phase_scale_scope='Supplementary kinematic diagnosis added after integrations started; not an orbit selection rule or measured return period.',
        code_hashes={str(p.relative_to(ROOT)):run.digest(p) for p in [HERE/'run.py',Path(__file__),HERE/'protocol.md']},
        coverage_sha256=run.digest(OUT/'coverage.parquet'))
    for p,sha in hashes.items():assert run.digest(ROOT/p)==sha
    run.save('results.json',result)

if __name__=='__main__':main()
