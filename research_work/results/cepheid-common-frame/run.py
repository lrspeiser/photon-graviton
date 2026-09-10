"""Declared training-only Cepheid reduction and frozen-field diagnostic."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactocentric, CartesianDifferential

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
P = json.loads((HERE/'protocol.json').read_text())
CACHE = ROOT/'research_work/data-cache/cepheid-stars'
for name, expected in P['frozen_files'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, name

def coords(d, scale=1., solar=None, delta=None):
    solar = solar or P['frame']
    delta = delta or {}
    sun_distance = solar.get('galcen_distance_kpc', np.hypot(solar['R_sun_kpc'], solar['z_sun_kpc']))
    frame = Galactocentric(galcen_distance=sun_distance*u.kpc,
        z_sun=solar['z_sun_kpc']*u.kpc,
        galcen_v_sun=CartesianDifferential(np.array(solar['v_sun_astropy_xyz_kms'])*u.km/u.s))
    c = SkyCoord(ra=d.ra.to_numpy()*u.deg, dec=d.dec.to_numpy()*u.deg,
        distance=d.distance_kpc.to_numpy()*scale*u.kpc,
        pm_ra_cosdec=(d.pmra.to_numpy()+delta.get('pmra',0))*u.mas/u.yr,
        pm_dec=(d.pmdec.to_numpy()+delta.get('pmdec',0))*u.mas/u.yr,
        radial_velocity=(d.radial_velocity.to_numpy()+delta.get('rv',0))*u.km/u.s).transform_to(frame)
    x,y,z = c.cartesian.xyz.to_value(u.kpc)
    vx,vy,vz = c.velocity.d_xyz.to_value(u.km/u.s)
    R = np.hypot(x,y)
    return np.c_[R, np.degrees(np.arctan2(y,-x)), z,
                 (x*vx+y*vy)/R, (y*vx-x*vy)/R, vz]

def variances(d, scale=1.):
    base=coords(d,scale)
    jac=np.stack([coords(d,scale,delta={k:1.})[:,3:]-base[:,3:] for k in ['pmra','pmdec','rv']],axis=2)
    distance_derivative=(coords(d,scale*1.0001)[:,3:]-coords(d,scale*.9999)[:,3:])/.0002
    sig=d[['pmra_error','pmdec_error','radial_velocity_error']].to_numpy()
    var=np.sum(jac**2*sig[:,None,:]**2,axis=2)
    var+=2*jac[:,:,0]*jac[:,:,1]*(sig[:,0]*sig[:,1]*d.pmra_pmdec_corr.to_numpy())[:,None]
    var+=(P['calibration']['distance_error_scenario_fraction']*distance_derivative)**2
    assert np.all(np.isfinite(var)) and np.all(var>=0)
    return var

def summarize_bins(c, var, membership):
    rows=[]
    for b in range(12):
        use=membership==b
        row=dict(bin=b,R_low_kpc=6+b,R_high_kpc=7+b,n=int(use.sum()))
        if use.sum()>=P['bins']['minimum_count_for_proxy']:
            R=float(c[use,0].mean())
            r2=float(np.mean(c[use,3]**2-var[use,0]))
            p2=float(np.mean(c[use,4]**2-var[use,1]))
            factor=1-R/P['velocity_proxy']['R_density_kpc']-2*R/P['velocity_proxy']['R_sigma_kpc']
            proxy2=p2-r2*factor
            raw2=float(np.mean(c[use,4]**2)-np.mean(c[use,3]**2)*factor)
            valid=(r2>=0 and p2>=0 and proxy2>0)
            row.update(R_mean_kpc=R,mean_vphi_kms=float(c[use,4].mean()),
                corrected_radial_second_moment=r2,corrected_azimuthal_second_moment=p2,
                jeans_proxy_kms=float(np.sqrt(proxy2)) if valid else None,
                uncorrected_jeans_proxy_kms=float(np.sqrt(raw2)) if raw2>0 else None,
                median_sigma_vphi_kms=float(np.median(np.sqrt(var[use,1]))),
                moment_valid=bool(valid))
        rows.append(row)
    return rows

def main():
    d=pd.read_parquet(CACHE/'gaia-dr3-dcep-parent-with-flags.parquet')
    bucket=d.source_id.map(lambda sid:int(hashlib.sha256(f'cepheid-common-frame-v1:{sid}'.encode()).hexdigest()[:8],16)%100)
    d['role']=np.where(bucket<60,'train',np.where(bucket<80,'validation','test'))
    calibrated=(d.mode_best_classification.eq('FUNDAMENTAL') & d.pf.gt(0)) | (d.mode_best_classification.eq('FIRST_OVERTONE') & d.p1_o.gt(0))
    eligible=d.flag_measurement_candidate & calibrated
    d[['source_id','role']].assign(measurement_and_mode_eligible=eligible).to_parquet(CACHE/'cepheid-common-frame-split.parquet',index=False)
    role_counts=d.role.value_counts().to_dict()
    eligible_counts=d.loc[eligible,'role'].value_counts().to_dict()
    # No held-out coordinates, velocities, residuals or outcome summaries below.
    t=d.loc[eligible & d.role.eq('train')].copy()
    fundamental=t.mode_best_classification.eq('FUNDAMENTAL')
    period=np.where(fundamental,t.pf,t.p1_o)
    cal=P['calibration']
    aa=np.where(fundamental,cal['fundamental']['intercept'],cal['first_overtone']['intercept'])
    bb=np.where(fundamental,cal['fundamental']['slope'],cal['first_overtone']['slope'])
    t['apparent_w']=t.int_average_g-cal['wesenheit_color_coefficient']*(t.int_average_bp-t.int_average_rp)
    t['absolute_w']=aa+bb*np.log10(period)
    t['distance_kpc']=10**((t.apparent_w-t.absolute_w-10)/5)
    assert np.all(np.isfinite(t.distance_kpc)) and np.all(t.distance_kpc>0)
    # Distance-modulus inversion sanity check for every processed star.
    assert np.max(abs((5*np.log10(t.distance_kpc)+10)-(t.apparent_w-t.absolute_w)))<1e-12
    c=coords(t)
    selected=(c[:,0]>=6)&(c[:,0]<=18)&(abs(c[:,1])<=30)&(abs(c[:,2])<=.5)&(abs(c[:,5])<=100)
    t=t.loc[selected].copy();c=c[selected]
    membership=np.clip(np.floor(c[:,0]-6).astype(int),0,11)
    var=variances(t)
    rows=summarize_bins(c,var,membership)
    for j,name in enumerate(['R_kpc','phi_deg','z_kpc','vR_kms','vphi_kms','vz_kms']):t[name]=c[:,j]
    t['bin']=membership
    t.to_parquet(CACHE/'cepheid-common-frame-train-selected.parquet',index=False)
    scenarios={str(scale):summarize_bins(coords(t,scale),variances(t,scale),membership) for scale in [.93,1.07]}
    frame_diagnostics={}
    for name,solar in {
        'Feng':dict(R_sun_kpc=8.275,z_sun_kpc=.025,v_sun_astropy_xyz_kms=[11.1,250.2,7.9]),
        'Eilers':dict(R_sun_kpc=8.122,z_sun_kpc=.025,v_sun_astropy_xyz_kms=[11.1,245.8,7.8])}.items():
        other=coords(t,solar=solar)
        frame_diagnostics[name]=dict(median_delta_R_kpc=float(np.median(other[:,0]-c[:,0])),
            median_delta_vphi_kms=float(np.median(other[:,4]-c[:,4])),
            delta_vphi_percentiles_kms=np.percentile(other[:,4]-c[:,4],[0,16,50,84,100]).tolist(),
            note='Exact coordinate transform of same selected training stars; not a circular-speed correction.')
    spec=importlib.util.spec_from_file_location('frozen_completion',HERE.parent/'conservative-field-completion/run.py')
    field=importlib.util.module_from_spec(spec);spec.loader.exec_module(field)
    outputs={}
    for refined in [False,True]:
        print('Building frozen fields, refined:',refined,flush=True)
        baryons=field.Baryons(refined);extra=field.Completion(baryons,refined)
        vb2=-c[:,0]*field.force(baryons,c[:,0],0)[:,0]
        vt2=vb2-c[:,0]*field.force(extra,c[:,0],0)[:,0]
        assert np.all(vb2>0) and np.all(vt2>0)
        pred={}
        for label,value in [('ordinary',vb2),('completion',vt2)]:
            pred[label]=[float(np.sqrt(value[membership==r['bin']].mean())) if r.get('moment_valid') else None for r in rows]
        outputs['refined' if refined else 'coarse']=pred
    scores={}
    for label,pred in outputs['refined'].items():
        residual=np.array([value-r['jeans_proxy_kms'] for value,r in zip(pred,rows) if value is not None])
        scores[label]=dict(bins=len(residual),rms_kms=float(np.sqrt(np.mean(residual**2))),bias_kms=float(residual.mean()))
        for row,value in zip(rows,pred):row[label+'_vc_kms']=value
    refinement={label:max(abs(a-b) for a,b in zip(outputs['coarse'][label],outputs['refined'][label]) if a is not None) for label in outputs['refined']}
    summary=dict(parent_role_counts=role_counts,measurement_mode_eligible_role_counts=eligible_counts,
        selected_training_stars=len(t),selected_modes=t.mode_best_classification.value_counts().to_dict(),
        selected_distance_percentiles_kpc=np.percentile(t.distance_kpc,[0,16,50,84,100]).tolist(),
        exact_coordinate_frame_diagnostics=frame_diagnostics,training_proxy_scores=scores,
        max_coarse_refined_bin_change_kms=refinement,
        held_out_outcomes_processed=False,new_physics_parameters_fitted=False,
        formula_status='Known empirical distance calibration, known coordinate transformations and approximate Jeans equation; frozen hypothetical/empirical additional field',
        limitations=['Random star split does not isolate spatial or survey systematics.', 'Published aggregate rotation curves already exposed.', '7 percent distance error is a sensitivity assumption; no complete covariance, extinction, metallicity or population likelihood.', 'Axisymmetry, density scale length and dispersion scale length imposed; vertical mixed moment omitted.', 'No exact original sample replication.', 'Training proxy discrepancies are not formal significances.'])
    for name,obj in [('results',summary),('training-bins',rows),('distance-sensitivity',scenarios)]:
        (HERE/f'{name}.json').write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
