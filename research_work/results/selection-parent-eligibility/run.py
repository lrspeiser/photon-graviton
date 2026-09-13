"""Audit a documented candidate eligibility rule against observed pilot targets."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache'
parent_path = CACHE / 'selection-fields/apogee2Object_300+00.fits'
main_path = CACHE / 'selection-observed-field/observed-main-metadata.parquet'
field_path = CACHE / 'selection-fields/allField.fits'
with fits.open(parent_path) as hd:
    t = hd[1].data
    cols = ['APOGEE_ID','RA','DEC','J','H','K','J_ERR','H_ERR','K_ERR','AK_TARG','AK_TARG_METHOD',
            'IRAC_4_5_ERR','WISE_4_5_ERR','TARG_4_5_ERR','TMASS_PHQUAL','TMASS_RDFLG',
            'TMASS_CCFLG','TMASS_GAL_CONTAM','TMASS_EXTKEY','TMASS_PROX']
    values = {}
    for c in cols:
        a = np.array(t[c])
        values[c] = np.char.strip(a.astype(str)) if a.dtype.kind in 'SU' else a.astype(a.dtype.newbyteorder('='))
    p = pd.DataFrame(values)
main = pd.read_parquet(main_path)
assert p.APOGEE_ID.is_unique and main.APOGEE_ID.is_unique
observed = p.APOGEE_ID.isin(main.APOGEE_ID)
assert observed.sum() == len(main)
field = fits.getdata(field_path,1)
field = field[np.char.strip(np.asarray(field['FIELD_NAME']).astype(str)) == '300+00']
assert len(field) == 1
center = SkyCoord(float(field['RA'][0])*u.deg, float(field['DEC'][0])*u.deg)
p['separation_deg'] = center.separation(SkyCoord(p.RA.to_numpy()*u.deg,p.DEC.to_numpy()*u.deg)).deg
# Known extinction convention used in the referenced APOGEE reader: AJ=2.5 AK.
p['jk0'] = p.J - p.K - 1.5*p.AK_TARG
finite = np.isfinite(p[['J','H','K','AK_TARG']]).all(axis=1) & p.AK_TARG.gt(-50)
cuts = {
    'finite_photometry_extinction': finite,
    'radius_at_most_0p8_deg': p.separation_deg.le(0.800000011920929),
    'JHK_errors_0_to_0p1': p[['J_ERR','H_ERR','K_ERR']].ge(0).all(axis=1) & p[['J_ERR','H_ERR','K_ERR']].le(.1).all(axis=1),
    'photometric_flags_A_or_B': p.TMASS_PHQUAL.str.fullmatch('[AB]{3}'),
    'read_flags_1_or_2': p.TMASS_RDFLG.str.fullmatch('[12]{3}'),
    'confusion_flags_zero': p.TMASS_CCFLG.eq('000'),
    'galaxy_contamination_zero': p.TMASS_GAL_CONTAM.eq('0'),
    'extkey_negative_sentinel': p.TMASS_EXTKEY.lt(0),
    'nearest_neighbor_at_least_6arcsec': p.TMASS_PROX.ge(6),
    'selected_midIR_error_0_to_0p1': p.TARG_4_5_ERR.between(0,.1),
    'H_10_to_13p3': p.H.between(10,float(np.float32(13.3))),
    'dereddened_color_at_least_0p5': p.jk0.ge(.5),
}
survives = pd.Series(True,index=p.index)
audit = []
for name,cut in cuts.items():
    cut = cut.fillna(False)
    survives &= cut
    audit.append(dict(cut=name, parent_individual_pass=int(cut.sum()),parent_cumulative_pass=int(survives.sum()),
                      observed_individual_fail=int((observed & ~cut).sum()),observed_cumulative_pass=int((observed & survives).sum())))
failures = []
for i in p.index[observed]:
    failed = [name for name,cut in cuts.items() if not cut.loc[i]]
    if failed:
        failures.append(dict(APOGEE_ID=p.loc[i,'APOGEE_ID'], failed=failed))
j = main[['APOGEE_ID','MIN_H','MAX_H','MIN_JK','MAX_JK']].merge(p[['APOGEE_ID','jk0']],on='APOGEE_ID',validate='one_to_one')
color_ok = j.jk0.ge(j.MIN_JK-1e-6) & j.jk0.le(j.MAX_JK+1e-6)
bins = []
without_midir = pd.Series(True,index=p.index)
for name,cut in cuts.items():
    if name != 'selected_midIR_error_0_to_0p1':
        without_midir &= cut
boundary_ids = main.loc[main.H.isin([float(np.float32(12.2)),float(np.float32(12.8)),float(np.float32(13.3))]),['APOGEE_ID','H','MIN_H','MAX_H']].to_dict('records')
for low,high in [(10,12.2),(12.2,12.8),(12.8,13.3)]:
    # Explicit endpoint convention; stored observed bounds are audited separately.
    hbin = p.H.gt(float(np.float32(low))) & p.H.le(float(np.float32(high)))
    for cmin,cmax in [(.5,.8),(.8,float('inf'))]:
        chosen = survives & hbin & p.jk0.ge(cmin) & p.jk0.lt(cmax)
        bins.append(dict(H_min=low,H_max=high,JK_min=cmin,JK_max=None if np.isinf(cmax) else cmax,
                         candidate_parent_count=int(chosen.sum()),observed_passing_count=int((chosen&observed).sum()),
                         parent_without_midIR_cut=int((without_midir&hbin&p.jk0.ge(cmin)&p.jk0.lt(cmax)).sum()),
                         observed_without_midIR_cut=int((observed&without_midir&hbin&p.jk0.ge(cmin)&p.jk0.lt(cmax)).sum())))
result = dict(scope='Candidate parent eligibility audit; weights not authorized by this diagnostic',
    input_sha256={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [parent_path,main_path,field_path]},
    quality_reference='https://arxiv.org/html/1708.00155#S4.T2',
    color_convention_reference='https://github.com/jobovy/apogee/blob/main/apogee/tools/read.py',
    candidate_parent_total=int(survives.sum()), observed_total=len(main), observed_pass=int((observed&survives).sum()),
    cuts=audit, failures=failures, bins=bins,
    H_endpoint_convention='lower exclusive, upper inclusive; observational bounds are kept separately',
    observed_H_boundary_cases=boundary_ids,
    sensitivity_without_midIR_cut=dict(parent=int(without_midir.sum()),observed=int((without_midir&observed).sum()),
                                     adopted=False),
    reconstructed_color_bound_disagreements=j.loc[~color_ok].to_dict('records'),
    extinction_methods=p.loc[observed,'AK_TARG_METHOD'].value_counts().to_dict(),
    missing_geometry='Central/off-axis camera masks and fiber allocation not reconstructed',
    heldout_kinematics_read=False, selection_weights_produced=False)
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['failures','input_sha256']},indent=2))
