"""Propagate the inner-kinematic mass posterior into the unchanged lens law."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy.special import logsumexp
from scipy.integrate import cumulative_trapezoid

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent/'slacs-component-refit'))
from model import ComponentModel, ARCSEC

paths = [HERE/'protocol.json', HERE.parent/'slacs-outer-predictive/protocol.json',
         HERE.parent/'slacs-resolved-input-audit/results.json',
         HERE.parent/'slacs-light-profile-audit/results.json',
         HERE.parent/'slacs-motion-lensing-pilot/results.json',
         HERE.parent/'lensing-data-readiness/conditional-geometry.json']
protocol, parproto, data, photo, base, geometry = [json.loads(f.read_text(encoding='utf-8')) for f in paths]
profiles = {r['Name']: r for r in photo['rows']}
geo = {r['Name']: r for r in geometry}
pilot = {r['Name']: r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20}
par = base['parameters']
A, p, astar = par['A'], par['p'], par['a_star_m_s2']*3.085677581491367e19/1e6
rows = []
for item in data['systems']:
    name = item['Name']
    if float(item['release_use_flag']) != 1 or profiles.get(name, {}).get('status') != 'conditional_components_available':
        continue
    Dl = geo[name]['conditional_Dl_Mpc']*1000
    ratio = geo[name]['conditional_Dls_over_Ds']
    a = pilot[name]['scale_a_kpc']
    edges = np.r_[item['inner_arcsec'], item['outer_arcsec'][-1]]*Dl/ARCSEC
    psf = item['psf_fwhm_arcsec']*Dl/ARCSEC/np.sqrt(8*np.log(2))
    comps = [dict(R=c['R_arcsec']*Dl/ARCSEC, n=c['n'], amp=c['amp_at_R'], bn=c['bn']) for c in profiles[name]['components']]
    model = ComponentModel(a, edges, psf, A, p, astar, 20, comps)
    y = np.array(item['vrms_kms'])[:-1]
    cov = np.array(item['covariance_kms_squared'])[:-1, :-1]
    chol = np.linalg.cholesky(cov)
    nm, nb = parproto['grids'][-1]
    lm = np.linspace(*np.log(parproto['mass_Msun_bounds']), nm)
    beta = np.linspace(*parproto['beta_bounds'], nb)
    m = np.exp(lm)/1e11
    coeff = np.array([model.coefficients(b) for b in beta])
    predictions = []
    for label, extra in [('baryons', False), ('empirical_extra', True)]:
        ll = np.empty((nb, nm))
        for j, (cb, cc) in enumerate(coeff):
            pred = np.sqrt(m[:, None]*cb[None, :] + (m[:, None]**p*cc[None, :] if extra else 0))
            white = np.linalg.solve(chol, (y-pred[:, :-1]).T)
            ll[j] = -.5*np.sum(white*white, axis=0)
        for prior in parproto['priors']:
            resolutions = []
            for step in [2, 1]:
                x, b = lm[::step], beta[::step]
                logp = ll[::step, ::step].copy()
                if prior == 'uniform_mass_and_beta':
                    logp += x[None, :]
                logp -= np.max(logp)
                marginal = np.trapezoid(np.exp(logp), b, axis=0)
                cdf = cumulative_trapezoid(marginal, x, initial=0)
                cdf /= cdf[-1]
                qlog = np.interp(protocol['mass_quantiles'], cdf, x)
                masses = np.exp(qlog)
                angles = [model.angle(float(v), Dl, ratio, extra) for v in masses]
                assert all(v is not None and np.isfinite(v) for v in angles)
                assert np.all(np.diff(angles) > 0)
                resolutions.append(dict(mass_nodes=len(x), beta_nodes=len(b),
                                        mass_quantiles_Msun=masses.tolist(), angle_quantiles_arcsec=angles))
            coarse, fine = resolutions
            dm = float(np.max(np.abs(np.log(coarse['mass_quantiles_Msun'])-np.log(fine['mass_quantiles_Msun']))))
            da = float(np.max(np.abs(np.array(coarse['angle_quantiles_arcsec'])-np.array(fine['angle_quantiles_arcsec']))))
            predictions.append(dict(model=label, prior=prior, resolutions=resolutions,
                                    max_log_mass_quantile_change=dm, max_angle_quantile_change_arcsec=da,
                                    numerical_pass=bool(dm <= protocol['maximum_coarse_fine_log_mass_quantile_change'] and
                                                        da <= protocol['maximum_coarse_fine_angle_quantile_change_arcsec'])))
    rows.append(dict(Name=name, catalog_SIE_arcsec=pilot[name]['catalog_SIE_arcsec'], predictions=predictions))
    print(json.dumps(rows[-1]), flush=True)
    out = dict(scope=protocol['scope'], quantile_levels=protocol['mass_quantiles'], rows=rows,
               complete=len(rows)==6,
               hashes={str(f.relative_to(ROOT)): hashlib.sha256(f.read_bytes()).hexdigest() for f in paths+[HERE/'run.py', HERE.parent/'slacs-component-refit/model.py', HERE.parent/'slacs-resolved-fit/model.py']})
    (HERE/'results.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
