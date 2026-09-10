"""Training-only normalization sensitivity; not a new SPS posterior."""
from pathlib import Path
import json,hashlib,shutil,argparse
import numpy as np
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
H=Path(__file__).resolve().parent;ROOT=H.parents[2]
parser=argparse.ArgumentParser();parser.add_argument('--updated-profile',action='store_true');args=parser.parse_args()
suffix='-updated-profile' if args.updated_profile else ''
C=ROOT/'research_work/data-cache/lens-stellar-populations'
p=Table.read(C/'table3.dat',format='ascii.cds',readme=str(C/'ReadMe')).to_pandas()
m=Table.read(C/'table4.dat',format='ascii.cds',readme=str(C/'ReadMe')).to_pandas()
assert len(p)==len(m)==85 and p.SDSS.is_unique and m.SDSS.is_unique
catalog=p.merge(m,on='SDSS',validate='one_to_one')
old=json.loads((H.parent/'lensing-data-readiness/lens-observations-and-image-models.json').read_text());old={r['Name']:r for r in old}
geometry=json.loads((H.parent/'lensing-data-readiness/conditional-geometry.json').read_text());geometry={r['Name']:r for r in geometry}
pilotfile=H.parent/('lens-training-pilot/predictions-refined'+suffix+'.json')
pilot=json.loads(pilotfile.read_text())
pilot=[r for r in pilot if r['seeing_fwhm_arcsec']==1.5 and r['cutoff_over_Re']==20]
names=sorted(set(r['Name'] for r in pilot));assert len(names)==33
catalog=catalog[catalog.SDSS.isin(names)].copy();assert len(catalog)==33
assert all(old[n]['role']=='training' for n in names)
# Reconstruct the publication's normalization ONLY, not the fictional universe.
reference=FlatLambdaCDM(H0=70,Om0=.3,Tcmb0=0)
rows=[];photo=[]
for _,r in catalog.iterrows():
    n=r.SDSS;z=old[n]['zFG'];D=geometry[n]['conditional_Dl_Mpc'];Lref=reference.luminosity_distance(z).value
    photo.append({k:r[k] for k in ['SDSS','Bmag','Vmag','f_Vmag','Imag','Hmag','Re(B)','Re(V)','Re(I)','MType','f_MType']})
    for imf,col,err in [('Chabrier','logMc','e_logMc'),('Salpeter','logMs','e_logMs')]:
        if not np.isfinite(r[col]):continue
        for branch,S in [('energy_loss_only',1.),('energy_loss_and_event_stretch',1+z)]:
            factor=D*D*(1+z)*S/Lref**2
            logmass=float(r[col])+np.log10(factor)
            for model in ['baryons','empirical_companion']:
                dyn=next(x for x in pilot if x['Name']==n and x['model']==model)
                rows.append(dict(Name=n,role='training',imf=imf,propagation_branch=branch,model=model,
                    published_log10_stellar_mass=float(r[col]),published_log10_mass_error=float(r[err]),normalization_factor=factor,
                    conditional_log10_stellar_mass=logmass,dispersion_inferred_mass_Msun=dyn['mass_from_sigma_Msun'],
                    log10_dynamical_to_conditional_stellar_mass=float(np.log10(dyn['mass_from_sigma_Msun'])-logmass),
                    old_Re_arcsec=old[n]['Re'],new_I_Re_arcsec=None if not np.isfinite(r['Re(I)']) else float(r['Re(I)'])))
summary={'training_systems':33,'extended_catalog_rows':85,'normalization_only':True,'scores':{},
    'band_coverage_training':{k:int(catalog[k].notna().sum()) for k in ['Bmag','Vmag','Imag','Hmag']},
    'missing_mass_training':int(catalog.logMc.isna().sum()),
    'photometric_error_column_present':False,
    'reference_cosmology_use':'Reconstruct published luminosity-distance normalization for a sensitivity conversion only; not assumed as fictional geometry.',
    'limitations':['Mass-to-light ratio, age/metallicity/dust/IMF and star-formation priors held fixed; not a re-fit of stellar populations.',
                  'Published uncertainties are not recalculated and omit these changes.',
                  'Earlier dynamical masses use the older size/profile approximation; no radius silently replaced.',
                  'No lensing or stellar holdout score opened; no parameter fitted to these comparisons.']}
for model in ['baryons','empirical_companion']:
    for imf in ['Chabrier','Salpeter']:
        for branch in ['energy_loss_only','energy_loss_and_event_stretch']:
            s=[x for x in rows if (x['model'],x['imf'],x['propagation_branch'])==(model,imf,branch)]
            delta=np.array([x['log10_dynamical_to_conditional_stellar_mass'] for x in s])
            summary['scores'][f'{model}/{imf}/{branch}']={'n':len(s),'median_mass_ratio':float(10**np.median(delta)),
                'median_log10_ratio':float(np.median(delta)),'log10_ratio_16_84':np.quantile(delta,[.16,.84]).tolist()}
summary['source_sha256']={n:hashlib.sha256((C/n).read_bytes()).hexdigest() for n in ['ReadMe','table3.dat','table4.dat']}
summary['origin']='https://cdsarc.cds.unistra.fr/ftp/J/ApJ/705/1099/'
summary['retrieved']='2026-09-10'
summary['updated_I_profile_used']=args.updated_profile
summary['pilot_sha256']=hashlib.sha256(pilotfile.read_bytes()).hexdigest()
if args.updated_profile:summary['limitations'][2]='Dynamical masses recomputed with updated I-band radius in the same spherical Hernquist approximation; no direct surface-brightness deprojection or population refit.'
for n,obj in [('results.json',summary),('normalization-sensitivity.json',rows)]:
    name=Path(n)
    (H/(name.stem+suffix+name.suffix)).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
# Pandas serialization retains missing photometry as null, not fabricated flux.
catalog[['SDSS','Bmag','Vmag','f_Vmag','Imag','Hmag','Re(B)','Re(V)','Re(I)','MType','f_MType']].to_json(H/'training-photometry.json',orient='records',indent=2)
raw=H/'source-tables';raw.mkdir(exist_ok=True)
for n in ['ReadMe','table3.dat','table4.dat']:shutil.copyfile(C/n,raw/n)
# Exact branch scaling check at fixed object, IMF and dynamical model.
for n in names:
    s=[x for x in rows if x['Name']==n and x['imf']=='Chabrier' and x['model']=='baryons']
    if len(s)==2:assert abs(s[1]['normalization_factor']/s[0]['normalization_factor']-(1+old[n]['zFG']))<1e-12
print(json.dumps(summary,indent=2))
