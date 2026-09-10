"""Role-filtered population normalization; no velocity or lens residual fitting."""
from pathlib import Path
import json,hashlib,argparse
import numpy as np
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
H=Path(__file__).resolve().parent;R=H.parent
parser=argparse.ArgumentParser();parser.add_argument('--role',choices=['validation','test'],default='validation');args=parser.parse_args()
if args.role=='test':H=R/'wave-lens-test'
protocol=json.loads((H/'protocol.json').read_text())
source=R/'lens-photometric-audit/source-tables'
p=Table.read(source/'table3.dat',format='ascii.cds',readme=str(source/'ReadMe')).to_pandas()
m=Table.read(source/'table4.dat',format='ascii.cds',readme=str(source/'ReadMe')).to_pandas()
catalog=p.merge(m,on='SDSS',validate='one_to_one').set_index('SDSS')
obsfile=R/'lensing-data-readiness/lens-observations-and-image-models.json'
geofile=R/'lensing-data-readiness/conditional-geometry.json'
obs=[r for r in json.loads(obsfile.read_text()) if r['role']==args.role]
geo={r['Name']:r for r in json.loads(geofile.read_text()) if r['role']==args.role}
reference=FlatLambdaCDM(H0=70,Om0=.3,Tcmb0=0)
rows=[];excluded=[];distance_errors=[]
for r in obs:
    name=r['Name'];reason=None
    if r['Mph']!='E':reason='Not early-type E in original catalog'
    elif not r['spectroscopic_dispersion_available']:reason='No spectroscopic dispersion'
    elif name not in catalog.index:reason='Missing population catalog row'
    elif not np.isfinite(catalog.loc[name,'logMs']):reason='Missing Salpeter mass'
    elif not np.isfinite(catalog.loc[name,'Re(I)']) or catalog.loc[name,'Re(I)']<=0:reason='Missing positive I-band radius'
    if reason:
        excluded.append({'Name':name,'reason':reason});continue
    row=catalog.loc[name];z=r['zFG'];D=geo[name]['conditional_Dl_Mpc']
    ref=reference.luminosity_distance(z).value
    independent=(1+z)*299792.458/70*quad(lambda s:1/np.sqrt(.3*(1+s)**3+.7),0,z,epsabs=1e-12)[0]
    distance_errors.append(abs(independent/ref-1))
    factor=D**2*(1+z)**2/ref**2
    rows.append({'Name':name,'role':args.role,'imf':'Salpeter','model':'baryons',
        'propagation_branch':'energy_loss_and_event_stretch','published_log10_stellar_mass':float(row['logMs']),
        'published_log10_mass_error':float(row['e_logMs']),'normalization_factor':factor,
        'conditional_log10_stellar_mass':float(row['logMs']+np.log10(factor)),
        'new_I_Re_arcsec':float(row['Re(I)'])})
assert len(rows)>0 and len(rows)+len(excluded)==len(obs) and max(distance_errors)<1e-10
if args.role=='validation':assert len(rows)==7 and len(excluded)==2
summary={'eligible':len(rows),'eligible_names':[r['Name'] for r in rows],'excluded':excluded,
    'reference_normalization_max_relative_check_error':max(distance_errors),
    'reference_cosmology_role':'Undo published luminosity normalization only, not the fictional model',
    'limitations':'Population M/L and age/dust/metallicity/IMF priors retained; not a new stellar-population posterior',
    'input_sha256':{str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in
        [source/'table3.dat',source/'table4.dat',source/'ReadMe',obsfile,geofile,H/'protocol.json']}}
for file,data in [('mass-inputs.json',rows),('preparation.json',summary)]:
    (H/file).write_text(json.dumps(data,indent=2)+'\n',newline='\n')
print(json.dumps(summary,indent=2))
