"""Compatibility cost of exact lens constraints with constant orbital anisotropy."""
from pathlib import Path
import hashlib,json
import numpy as np
P=Path(__file__).resolve().parent;L=P.parent/'isotropic-galaxy-transfer'
branches=['original','reference','local','recycling']
out=dict(scope='Catalogue lens angle consumed to fix stellar mass; inner-star constant anisotropy fit; conditional exact-constraint cost, not calibrated significance or lens prediction',rows=[],summary=[])
paths=[L/'lensing.py']
budget=json.loads((P/'lens-bending-budget-results.json').read_text())
budget_rows={(r['branch'],r['Name'],r['population']):r for r in budget['rows']}
for branch in branches:
    freepath=L/('capacity-'+branch+'-optics-results.json')
    exactpath=L/('capacity-'+branch+'-exact-lens-optics-results.json')
    free=json.loads(freepath.read_text());exact=json.loads(exactpath.read_text())
    paths += [freepath,exactpath]
    old={(r['Name'],r['model']):r for r in free['rows']}
    assert len(exact['rows'])==12 and exact['C0_multiplier']==free['C0_multiplier']
    for r in exact['rows']:
        f=old[r['Name'],r['model']];pop=r['retention_mapping']['population']
        for key in ['observed_stellar_vrms','geometry','retention_mapping','lens_catalog_arcsec']:
            assert r[key]==f[key]
        assert abs(r['lens_fractional_residual'])<1e-7
        assert r['inner_chi2']>=f['inner_chi2']-1e-5
        assert r['optimizer_successes']>0 and not r['orbit_boundary']
        assert 1e7<r['mass_Msun']<1e14
        if branch in ['original','local']:
            required=budget_rows[branch,r['Name'],pop]['stellar_mass_multiplier_at_fixed_companions']
            assert abs(r['mass_Msun']/f['mass_Msun']-required)<1e-7
        out['rows'].append(dict(branch=branch,Name=r['Name'],population=pop,
            free_inner_chi2=f['inner_chi2'],exact_inner_chi2=r['inner_chi2'],delta_inner_chi2=r['inner_chi2']-f['inner_chi2'],
            free_beta=f['beta'],exact_beta=r['beta'],stellar_mass_multiplier=r['mass_Msun']/f['mass_Msun'],
            free_outer_standardized_residual=f['outer_conditional_standardized_residual'],exact_outer_standardized_residual=r['outer_conditional_standardized_residual'],
            exact_angle_fractional_residual=r['lens_fractional_residual']))
    for pop in ['Chabrier','Salpeter']:
        rr=[r for r in out['rows'] if r['branch']==branch and r['population']==pop]
        out['summary'].append(dict(branch=branch,population=pop,
            free_inner_chi2_sum=sum(r['free_inner_chi2'] for r in rr),exact_inner_chi2_sum=sum(r['exact_inner_chi2'] for r in rr),
            delta_inner_chi2_sum=sum(r['delta_inner_chi2'] for r in rr),
            free_outer_residual_square_sum=sum(r['free_outer_standardized_residual']**2 for r in rr),exact_outer_residual_square_sum=sum(r['exact_outer_standardized_residual']**2 for r in rr),
            systems_with_delta_inner_chi2_above_50=[r['Name'] for r in rr if r['delta_inner_chi2']>50],
            exact_beta_range=[min(r['exact_beta'] for r in rr),max(r['exact_beta'] for r in rr)]))
paths += [P/'lens-bending-budget-results.json',P/'capacity-lensing-results.json',P.parent/'slacs-outer-bin-check/protocol.json',P.parent/'slacs-resolved-input-audit/results.json']
out['input_sha256']={str(f.relative_to(P.parents[2])).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}
(P/'exact-lens-orbits-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['summary'],indent=2))
