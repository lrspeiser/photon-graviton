"""Fit shared scaling on Sun/Earth, then predict Moon and other bodies."""
from pathlib import Path
import json,numpy as np
HERE=Path(__file__).resolve().parent
solar=json.loads((HERE/'solar-photon-budget-results.json').read_text())
planets=json.loads((HERE/'planet-well-energy-results.json').read_text())['rows']
Ms=1.9884e30;Rs=6.957e8;Ls=solar['nominal_solar_luminosity_W'];c=299792458.
alpha=solar['alpha_per_Mpc']/3.085677581491367e22
Pcs=Ls*(-np.expm1(-alpha*Rs))
data=[]
for r in planets:
    mu=r['mass_kg']/Ms;rad=r['radius_m']/Rs;lum=r['intercepted_solar_power_W']/Ls
    pc=r['intercepted_solar_power_W']*(-np.expm1(-alpha*r['radius_m']))/Pcs
    data.append(dict(body=r['body'],mu=mu,r=rad,l=lum,pc=pc,mass=r['mass_kg']))
e=next(r for r in data if r['body']=='Earth');results=[]
families=[('constant',0,0,0),('size',0,1,0),('well_depth',1,-1,0),('surface_acceleration',1,-2,0),('density',1,-3,0),('surface_radiative_flux',0,-2,1)]
for name,a,b,d in families:
    exponent=0. if name=='constant' else np.log(e['mu']/e['pc'])/(a*np.log(e['mu'])+b*np.log(e['r'])+d*np.log(e['l']))
    denominator=1-a*exponent
    assert abs(denominator)>1e-10
    rows=[]
    for v in data:
        mu=np.exp((np.log(v['pc'])+exponent*(b*np.log(v['r'])+d*np.log(v['l'])))/denominator)
        # Verify the implicit equation with the PREDICTED mass, not the observed target.
        rhs=v['pc']*(mu**a*v['r']**b*v['l']**d)**exponent
        assert abs(rhs/mu-1)<1e-10
        rows.append(dict(body=v['body'],role='calibration' if v['body']=='Earth' else 'exposed_prediction',observed_mass_kg=v['mass'],predicted_mass_kg=mu*Ms,predicted_over_observed=mu/v['mu']))
    if name!='constant':assert abs(next(r for r in rows if r['body']=='Earth')['predicted_over_observed']-1)<1e-10
    test=[r for r in rows if r['body']!='Earth']
    results.append(dict(family=name,exponent=exponent,predictor_powers=dict(mass=a,radius=b,luminosity=d),implicit_exponent_denominator=denominator,
        transfer_log10_RMS=float(np.sqrt(np.mean([np.log10(r['predicted_over_observed'])**2 for r in test]))),
        within_factor_two=sum(bool(.5<=r['predicted_over_observed']<=2) for r in test),prediction_count=len(test),rows=rows))
out=dict(status='Exploratory full-mass scaling; source/capture geometry stipulated; no independent support or energy-history solution',
         solar_K_seconds=Ms*c*c/Pcs,calibration=['Sun','Earth'],models=results)
(HERE/'capture-law-family-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
for r in results:
    print(r['family'],'p=',r['exponent'],'logRMS=',r['transfer_log10_RMS'],'factor2=',r['within_factor_two'])
    print(' '.join(v['body']+'='+format(v['predicted_over_observed'],'.4g') for v in r['rows'] if v['body']!='Earth'))

