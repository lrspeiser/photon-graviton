"""Reproducible, scoped electromagnetic audit. No expansion distance fit is used.

All observational scores are diagnostics on previously exposed data or published
summaries. Synthetic tests and prescribed fields are explicitly separated.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import re

import numpy as np
from scipy.constants import c, h, k, astronomical_unit
from scipy.integrate import quad, solve_ivp
from scipy.optimize import least_squares, minimize_scalar, brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
P = json.loads((HERE / 'protocol.json').read_text())
OBS = json.loads((HERE / 'observations.json').read_text())
MPC = 1e6 * astronomical_unit * 648000 / math.pi
LY = c * 31557600
ALPHA = P['alpha_per_mpc']
INPUT_PATHS = [
    'redshift_paper/all_164_groups.csv',
    'research_work/results/redshift-priority/maser-comparison.csv',
    'temporal_candidate_audit/data/spectral_aging.json',
    'time_first_principles/data/blondin2008.txt',
    'time_revision/data/firas.txt',
    'temporal_candidate_audit/data/published_measurements.json',
    'research_work/results/direct-conversion/results.json',
    'research_work/results/timing-foundation/input-audit.json',
    'research_work/results/timing-population/report.md',
    'research_work/results/atomic-line-response/report.md',
    'research_work/results/time-only-geometry/report.md',
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n',
                             encoding='utf-8', newline='\n')


def write_rows(name, rows):
    assert rows
    with (HERE / name).open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_rows(rel):
    with (ROOT / rel).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def power_log_stretch(p, A, ratio):
    """Exact energy-dependent drift; ratio is emitted E/E_ref.

    E_f/E_i = [1+p*A*ratio**p]**(-1/p). p=-1 is fixed absolute
    loss. If that model exhausts a photon, its redshift is undefined.
    """
    ratio = np.asarray(ratio, dtype=float)
    if p == 0:
        return np.zeros_like(ratio) + A
    term = p*A*ratio**p
    if np.any(term <= -1):
        raise ValueError('photon exhausted or outside positive-energy domain')
    return np.log1p(term)/p


def redshift_distances():
    out = []
    for row in read_rows(INPUT_PATHS[0]):
        D = float(row['catalog_distance_mpc'])
        z = float(row['observed_cmb_z'])
        prediction = math.expm1(ALPHA * D)
        out.append(dict(sample='SBF_group', object=row['group_pgc'], split=row['split'],
                        D_mpc=D, observed_z=z, predicted_z=prediction,
                        residual_c_dz_km_s=c/1000*(prediction-z)))
    for row in read_rows(INPUT_PATHS[1]):
        D = float(row['distance_mpc'])
        z = float(row['observed_cmb_z'])
        prediction = math.expm1(ALPHA * D)
        assert abs(prediction-float(row['exponential_z'])) < 1e-13
        out.append(dict(sample='maser_galaxy', object=row['name'], split='exposed_external',
                        D_mpc=D, observed_z=z, predicted_z=prediction,
                        residual_c_dz_km_s=c/1000*(prediction-z)))
    write_rows('redshift-predictions.csv', out)
    stats = []
    for sample, split in [('SBF_group', 'train'), ('SBF_group', 'validation'),
                          ('SBF_group', 'test'), ('maser_galaxy', 'exposed_external')]:
        errors = np.array([r['residual_c_dz_km_s'] for r in out
                           if r['sample'] == sample and r['split'] == split])
        stats.append(dict(sample=sample, split=split, n=len(errors),
                          rms_km_s=float(np.sqrt(np.mean(errors**2))),
                          mean_km_s=float(errors.mean())))
    old = json.loads((ROOT / INPUT_PATHS[6]).read_text())
    assert abs(stats[2]['rms_km_s']-old['statistics']['test']['conversion_path']['rms_km_s']) < 1e-9
    return dict(statistics=stats, interpretation='Fixed exposed-data prediction check; motion and covariance unresolved; no new parity claim')


def timing():
    raw = json.loads((ROOT / INPUT_PATHS[2]).read_text())
    # Independent parsing of the published table prevents confusing aging rate
    # with duration stretch or accidentally taking a theoretical table column.
    text = (ROOT / INPUT_PATHS[3]).read_text(encoding='utf-8')
    parsed = {}
    for line in text.splitlines():
        match = re.match(r'^\s*(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+\((\d+\.\d+)\)', line)
        if match and match[1] in {r[0] for r in raw}:
            parsed[match[1]] = [float(match[2]), float(match[4]), float(match[5])]
    assert len(parsed) == len(raw) == 35
    assert all(np.allclose(parsed[r[0]], r[1:], rtol=0, atol=1e-12) for r in raw)
    z, age, error = np.array([r[1:] for r in raw]).T
    out = [dict(object=r[0], z=r[1], observed_aging_rate=r[2], sigma=r[3],
                stationary_prediction=1, stretched_prediction=1/(1+r[1])) for r in raw]
    write_rows('spectral-aging-predictions.csv', out)
    results = []
    for label, selected in [('all', np.ones(35, dtype=bool)), ('low', z<.2), ('high', z>.2)]:
        zz, yy, ee = z[selected], age[selected], error[selected]
        for extra in [0., .1, .2]:
            sigma = np.sqrt(ee**2 + (extra*yy)**2)
            def score(b):
                return float(np.sum(((yy-(1+zz)**(-b))/sigma)**2))
            fit = minimize_scalar(score, bounds=(-2, 8), method='bounded', options={'xatol':1e-11})
            assert fit.success
            lo = brentq(lambda b:score(b)-score(fit.x)-1, -10, fit.x)
            hi = brentq(lambda b:score(b)-score(fit.x)-1, fit.x, 20)
            results.append(dict(sample=label, n=int(selected.sum()), extra_fractional_error=extra,
                                chi2_b0=score(0), chi2_b1=score(1), b_fit=float(fit.x),
                                chi2_fit=score(fit.x), formal_delta_chi2_one_interval=[lo, hi]))
    # The rounded published table cannot reproduce unrounded paper scores exactly.
    assert abs(results[0]['chi2_b0']-150.3) < 1
    assert abs(results[0]['chi2_b1']-27.0) < 1
    # Intrinsic source evolution exponent e and propagation exponent b enter as
    # their sum: age rate=(1+z)^[-(b+e)]. Data alone cannot separate them.
    np.testing.assert_allclose((1+z)**(-1), (1+z)**(-.3)*(1+z)**(-.7))
    return dict(scores=results, caveat='Published spectral-age summaries with source-template assumptions; no raw-spectrum refit, shared-template covariance, selection calibration or global significance. Extra-error rows are sensitivity examples only.',
                source_evolution_degeneracy='Only b_propagation + e_intrinsic is identifiable in this chosen family; b=0,e=1 gives the same timing predictions as b=1,e=0, but requires independent source physics.',
                DES_status='Author width summaries not independently refitted; the existing observer-time inference pilot failed its numerical gate and is not a validation result.')


def radio():
    r = OBS['radio_methanol']
    low, high = r['independent_centroid_groups']
    nu = np.array([low['frequency_GHz'], high['frequency_GHz']])
    z0 = low['z']
    A0 = math.log1p(z0)
    rows = []
    def anchored(p):
        # Set E_ref to the lower line; the common A is calibrated once per p.
        A = A0 if p == 0 else math.expm1(p*A0)/p
        zz = np.expm1(power_log_stretch(p, A, nu/nu[0]))
        # Conventional small differential shift reporting units, even for bad
        # models where the mismatch becomes large; not a physical gas velocity.
        offset = c/1000 * (zz[1]-zz[0])/(1+z0)
        return A, zz, offset
    for p in [-1, 0, 1, 2, 3]:
        A, zz, offset = anchored(p)
        rows.append(dict(p=p, calibrated_A=A, predicted_low_z=float(zz[0]),
                         predicted_high_z=float(zz[1]), observed_high_z=high['z'],
                         predicted_offset_km_s=float(offset),
                         offset_residual_in_quoted_sigma=float((offset-r['high_minus_low_velocity_offset_km_s'])/r['offset_sigma_km_s'])))
    fitted = brentq(lambda p:anchored(p)[2]-r['high_minus_low_velocity_offset_km_s'], -.1, .1, xtol=1e-15)
    formal = [brentq(lambda p:anchored(p)[2]-(r['high_minus_low_velocity_offset_km_s']+sign*r['offset_sigma_km_s']), -.1, .1, xtol=1e-15) for sign in [-1, 1]]
    fwhm = r['FWHM_km_s']
    sigma_fraction = fwhm / math.sqrt(8*math.log(2)) / (c/1000)
    eps_max = math.log1p(sigma_fraction**2)/A0
    noise = []
    for epsilon in [1e-3, 1e-6, 1e-9, 1e-12]:
        width = c/1000 * math.sqrt(math.expm1(A0*epsilon))
        noise.append(dict(fraction_per_jump=epsilon, mean_number=A0/epsilon,
                          predicted_rms_fractional_width_as_km_s=width))
    write_rows('radio-chromaticity.csv', rows)
    return dict(models=rows, fitted_p=float(fitted), formal_one_sigma_p_interval=formal,
                achromatic_quoted_sigma_residual=rows[1]['offset_residual_in_quoted_sigma'],
                p_caveat='Two-centroid conditional comparison with the lower line anchored; covariance and source systematics not fully supplied. Not a global p measurement or proof of photon conversion. A not independently inferred from distance.',
                jump_noise=noise, conditional_max_fraction_per_independent_jump=eps_max,
                linewidth_caveat='All observed Gaussian width generously assigned to independent multiplicative jumps; ordinary gas broadening would consume some allowance. Bound does not apply to coherent deterministic frequency conversion.',
                excluded='12.179 GHz profile mismatch; second 48 GHz centroid tied, not counted independently')


def planck(nu, temperature):
    return 2*h*nu**3/c**2/np.expm1(h*nu/(k*temperature))/1e-20


def thermal():
    data = np.loadtxt(ROOT / INPUT_PATHS[4])
    assert data.shape == (43, 5) and np.all(data[:, 3] > 0)
    nu, observed, sigma, gal = data[:,0]*100*c, data[:,2], data[:,3], data[:,4]
    reference = planck(nu, 2.725)
    wg = gal/sigma
    def fit_amplitude(amplitude):
        def predicted(T):
            return 1000*(amplitude*planck(nu,T)-reference)
        def profile(T, details=False):
            diff = (observed-predicted(T))/sigma
            g = float(wg@diff/(wg@wg))
            resid = diff-g*wg
            return (float(resid@resid), g, predicted(T)+g*gal) if details else float(resid@resid)
        fit = minimize_scalar(profile, bounds=(.5,4), method='bounded', options={'xatol':1e-12})
        assert fit.success
        chi2,g,pred = profile(fit.x,True)
        check = least_squares(lambda x:(predicted(x[0])+x[1]*gal-observed)/sigma,
                              [fit.x,g],bounds=([.5,-np.inf],[4,np.inf]),
                              ftol=1e-12,xtol=1e-12,gtol=1e-12)
        assert check.success and abs(check.fun@check.fun-chi2)<1e-5*max(1,chi2)
        return dict(amplitude=amplitude,T_color_K=float(fit.x),galaxy_coefficient=g,chi2_diagonal=chi2), pred
    q100 = math.exp(-ALPHA*100e6*LY/MPC)
    profiles, channels = [], []
    for q in [1., q100, .99, .5, 1/2.2]:
        for label, survival in [('number_preserved',1.),('extra_photon_removal',q**3)]:
            fit,pred = fit_amplitude(survival*q**-3)
            profiles.append(dict(q=q, model=label, photon_survival=survival, final_energy_fraction=q*survival,**fit))
            for i in range(len(nu)):
                channels.append(dict(q=q,model=label,frequency_Hz=nu[i],observed_residual_kJy_sr=observed[i],
                                     sigma_kJy_sr=sigma[i],predicted_residual_kJy_sr=pred[i]))
    write_rows('firas-predictions.csv', channels)
    # Dimensionless independent Planck integrals: initial temperature=1.
    moments = []
    for q in [.5,.9,.99]:
        number = quad(lambda x:x*x/math.expm1(x),1e-8,100,epsabs=1e-11)[0]
        energy = quad(lambda x:x**3/math.expm1(x),1e-8,100,epsabs=1e-11)[0]
        newN = quad(lambda x:x*x/q**3/math.expm1(x/q),1e-8,100*q,epsabs=1e-11)[0]
        newU = quad(lambda x:x**3/q**3/math.expm1(x/q),1e-8,100*q,epsabs=1e-11)[0]
        assert abs(newN/number-1)<1e-9 and abs(newU/energy-q)<1e-9
        # Concurrent removal at 3 times the fractional energy-loss rate:
        # dU/dA=-4U, conversion Q=U, removal sink=3U.
        lost = 1-q**4
        assert abs(q**4+lost/4+3*lost/4-1)<1e-14
        moments.append(dict(q=q,number_ratio=newN/number,energy_ratio=newU/energy,
                            concurrent_partial_transfer=lost/4,concurrent_removal_energy=3*lost/4))
    return dict(n=43,profiles=profiles,moment_checks=moments,
                caveat='Initially Planckian homogeneous fixed-volume bath only; no age/origin fitted, no full FIRAS covariance or calibration likelihood. Removal fixes this normalization by construction while changing photon count, flux and energy sink requirements.')


def wave_timing():
    def trace(mode,sign,launch,tol):
        def field(x,t):
            w = math.sin(math.pi*x)**2
            f,ft = (1.,0.) if mode=='static' else (1-math.exp(-.2*t),.2*math.exp(-.2*t))
            return 1+sign*.2*w*f, sign*.2*w*ft
        def ode(x,y):
            q,qt = field(x,y[0])
            return [1/q, -qt/q**2*y[1], qt/q**2]
        sol=solve_ivp(ode,(0,1),[launch,1,0],method='DOP853',rtol=tol,atol=tol*.01,max_step=.02)
        assert sol.success
        arrival,J,log_energy=sol.y[:,-1]
        return arrival,J,math.exp(log_energy)
    rows=[]
    for mode in ['static','evolving']:
        for sign in [-1,1]:
            for launch in [0.,.2,1.]:
                arrival,S,E=trace(mode,sign,launch,2e-12)
                refine=trace(mode,sign,launch,2e-10)
                assert abs(S-refine[1])<1e-9 and abs(S*E-1)<1e-10
                fd=(trace(mode,sign,launch+1e-5,2e-12)[0]-arrival)/1e-5
                assert abs(fd-S)<1e-6
                if mode=='static':assert abs(S-1)<1e-12
                if mode=='evolving' and sign<0:assert S>1 and E<1
                if mode=='evolving' and sign>0:assert S<1 and E>1
                for interval in [1e-5,.001,.01,.1]:
                    later=trace(mode,sign,launch+interval,2e-12)[0]
                    rows.append(dict(mode=mode,void_sign=sign,emission_time=launch,interval=interval,
                                     arrival_time=arrival,instantaneous_stretch=S,energy_survival=E,
                                     finite_event_stretch=(later-arrival)/interval,
                                     finite_relative_distortion=(later-arrival)/interval/S-1))
    write_rows('prescribed-time-field.csv',rows)
    return dict(n=len(rows),field='q=1+sign*0.2*sin(pi*x)^2*f(t), c=L=1, q=1 at both endpoints',
                max_static_error=max(abs(r['finite_event_stretch']-1) for r in rows if r['mode']=='static'),
                evolving_stretch_range=[min(r['instantaneous_stretch'] for r in rows if r['mode']=='evolving'),max(r['instantaneous_stretch'] for r in rows if r['mode']=='evolving')],
                largest_finite_event_distortion=max(abs(r['finite_relative_distortion']) for r in rows),
                status='Kinematic construction, not a sourced field, observed void or completed companion theory. Positive sign in this particular evolving family blueshifts; static profiles only delay.')


def grid_and_checks():
    bands=dict(radio=1e8,microwave=1e11,infrared=3e13,visible=6e14,ultraviolet=3e15,X_ray=2.4e17,gamma_ray=2.4e23)
    distances={'zero':0.,'1_mm':.001,'1_m':1.,'1_au':astronomical_unit,'100_au':100*astronomical_unit,
               '1_ly':LY,'1000_ly':1000*LY,'1_Mly':1e6*LY,'100_Mly':1e8*LY,'1000_Mly':1e9*LY}
    rows=[]
    for name,metres in distances.items():
        A=ALPHA*metres/MPC
        S=math.exp(A)
        loss=-math.expm1(-A)
        for band,nu in bands.items():
            rows.append(dict(band=band,distance=name,D_mpc=metres/MPC,emitted_Hz=nu,observed_Hz=nu/S,
                             emitted_eV=h*nu/1.602176634e-19,redshift=math.expm1(A),
                             energy_fraction_to_companions=loss,stationary_event_stretch=1,
                             proposed_event_stretch=S,proposed_flux_relative_to_inverse_square=math.exp(-2*A),
                             extra_removal_flux_relative_to_inverse_square=math.exp(-5*A),
                             status='synthetic_prediction_not_observed'))
    write_rows('seven-band-predictions.csv',rows)
    ode_errors=[]
    for p in [-1,0,1,2,3]:
        for e in [.5,1,2]:
            A=.01
            exact=e*math.exp(-float(power_log_stretch(p,A,e)))
            sol=solve_ivp(lambda x,y:[-y[0]**(p+1)],(0,A),[e],method='DOP853',rtol=1e-12,atol=1e-14)
            assert sol.success
            ode_errors.append(abs(sol.y[0,-1]/exact-1))
    assert max(ode_errors)<1e-10
    # Independently integrated spectral flux Jacobian, arbitrary test spectrum.
    flux_errors=[]
    for S in [1.,1.2,2.2]:
        source=quad(lambda nu:nu**2*math.exp(-nu),0,100)[0]
        observed=quad(lambda nu:(S*nu)**2*math.exp(-S*nu)/S,0,100/S)[0]
        flux_errors.append(abs(observed/source-S**-2))
    assert max(flux_errors)<1e-12
    # Explicit transfer then capture ledger, no permanent void storage implied.
    ledger_errors=[]
    for A in [0,.01,math.log(2),1.]:
        for K in [0,.5,3.]:
            photon=math.exp(-A)
            free=-math.expm1(-A)*math.exp(-K)
            deposited=-math.expm1(-A)*-math.expm1(-K)
            ledger_errors.append(abs(photon+free+deposited-1))
    assert max(ledger_errors)<1e-14
    beta=-300000/c
    # Purely illustrative approaching source: motion may beat positive transfer.
    example=math.exp(ALPHA)*math.sqrt((1+beta)/(1-beta))-1
    assert example<0
    return dict(n_predictions=len(rows),maximum_drift_ODE_relative_error=max(ode_errors),
                maximum_flux_integral_error=max(flux_errors),maximum_ledger_error=max(ledger_errors),
                illustrative_1Mpc_approaching_300km_s_total_z=example,
                conclusion='Positive path transfer is compatible with total blueshift when independent motion/endpoint factors dominate; do not fit arbitrary velocities to erase residuals.')


def local_and_multimessenger():
    old=json.loads((ROOT/INPUT_PATHS[5]).read_text())
    gamma_year=ALPHA*c/MPC*31557600
    clock=old['alpha_clock']
    bound=(abs(clock['alpha_dot_per_year'])+1.96*clock['sigma_per_year'])/gamma_year
    g=OBS['laboratory_gradient']
    sigma=math.hypot(g['stat_error_per_cm'],g['sys_error_per_cm'])
    clock_resid=(g['ordinary_prediction_per_cm']-g['frequency_gradient_per_cm'])/sigma
    # An extra clock-response coefficient is distinct from the path-loss rate.
    # No derived interaction demands dln(alpha_fs)/dt=gamma_year.
    speeds=[]
    for D in old['GW170817']['distance_sensitivity_Mpc']:
        for frac in [0.,1e-15,1e-12]:
            dt=D*MPC/c*frac
            speeds.append(dict(D_mpc=D,fractional_speed_mismatch_linearized=frac,propagation_delay_s=dt))
    drift=OBS['redshift_drift']
    return dict(loss_rate_per_year=gamma_year,
                optical_clock_coupling=dict(assumption='dln(alpha_fs)/dt = beta_clock * gamma; optional, not implied by universal time',
                                           conservative_gaussian_beta_abs_bound=bound,
                                           beta_zero_residual_sigma=-clock['alpha_dot_per_year']/clock['sigma_per_year'],
                                           beta_one_residual_sigma=(gamma_year-clock['alpha_dot_per_year'])/clock['sigma_per_year']),
                lab_clock_summary=dict(ordinary_prediction_residual_sigma=clock_resid,
                                       path_conversion_over_1cm=math.expm1(ALPHA*.01/MPC),
                                       caveat=g['limitations']),
                zero_drift_summary_residual_sigma=-drift['velocity_drift_m_s_year']/drift['sigma_m_s_year'],
                drift_caveat='A stationary conversion law predicts zero transfer drift for fixed path; no event stretch. A time-field branch has to derive its own drift, not inherit this zero.',
                multimessenger=speeds,
                observed_GW_gamma_gap=old['GW170817'],
                GRB090510=OBS['GRB090510'],
                multimessenger_caveat='Observed gaps combine emission lag and propagation. Constant equal speeds predict no propagation difference, not zero observed gap. A new companion need not be the measured GR tensor wave; equating them requires an action and source calculation.')


def exact_accumulating_time():
    """Post-initial-run explicit repair, tested under a separate amendment.

    Prescribe q=1/[1+c*kappa(s)*(t-t_ref)] and kappa=alpha*2sin²(pi*s/L).
    Then d t/d s=1/c+kappa(s)*(t-t_ref), so the arrival derivative is
    exp(integral kappa ds). This is a kinematic lapse, not a sourced metric.
    """
    import mpmath as mp
    mp.mp.dps=50
    rows=[]
    worst=0.
    for L_mpc in [1., 100e6*LY/MPC, 1000e6*LY/MPC]:
        A=ALPHA*L_mpc
        for launch_u in [0., .1, 1., 10.]:
            def ode(x,y):
                aw=A*2*math.sin(math.pi*x)**2
                return [1+aw*y[0], aw*y[1], -aw]
            sol=solve_ivp(ode,(0,1),[launch_u,1.,0.],method='DOP853',
                          rtol=1e-12,atol=1e-14,max_step=.02,dense_output=True)
            assert sol.success
            arrival,J,logE=sol.y[:,-1]
            S=math.exp(A)
            # Integrating-factor solution, independent quadrature.
            I=quad(lambda x:math.exp(-A*(x-math.sin(2*math.pi*x)/(2*math.pi))),0,1,epsabs=1e-13)[0]
            exact_arrival=S*(launch_u+I)
            worst=max(worst,abs(arrival-exact_arrival),abs(J/S-1),abs(math.exp(logE)*S-1))
            assert worst<1e-9
            seconds=L_mpc*MPC/c
            # High-precision finite interval check avoids subtracting nearly
            # equal ~billion-year arrivals with float64 and claiming precision.
            mA=mp.mpf(str(A)); mlaunch=mp.mpf(str(launch_u))
            delta=mp.mpf(10*86400)/mp.mpf(str(seconds))
            mI=mp.quad(lambda x:mp.exp(-mA*(x-mp.sin(2*mp.pi*x)/(2*mp.pi))),[0,1])
            t1=mp.exp(mA)*(mlaunch+mI)
            t2=mp.exp(mA)*(mlaunch+delta+mI)
            finite=float((t2-t1)/delta)
            assert abs(finite/S-1)<1e-14
            grid=np.linspace(0,1,301)
            q=1/(1+2*A*np.sin(np.pi*grid)**2*sol.sol(grid)[0])
            rows.append(dict(L_mpc=L_mpc,launch_time_in_L_over_c=launch_u,stretch=S,
                             finite_ten_day_arrival_interval_days=10*finite,
                             energy_survival=math.exp(logE),minimum_q_on_ray=float(min(q)),
                             same_metric_companion_survival=math.exp(logE),
                             fixed_path_redshift_drift=0.))
    write_rows('exact-time-candidate.csv',rows)
    gamma=ALPHA*c/MPC*31557600
    return dict(n=len(rows),maximum_ODE_vs_integrating_factor_error=worst,
                time_origin_domain_past_years_for_peak_kappa=1/(2*gamma),
                formula='q(s,t)=1/[1+c*kappa(s)*(t-t_ref)]; q=1 at endpoints where kappa=0',
                ten_day_interval_at_100Mly_days=rows[4]['finite_ten_day_arrival_interval_days'],
                construction='Known integrating-factor and lapse geometry, proposed inverse-affine time dependence; originality unverified. Chosen to enforce desired stretch, not derived from sources.',
                failures=['Same-metric massless companions also redshift: no-loss transport does not follow.',
                          'Prescribed pure-lapse, fixed-flat-space geometry lacks the required positive-energy GR source; modified gravity or spatial geometry must be derived.',
                          'Spatial time gradients exert forces on matter in a universal coupling.',
                          'Positive denominator has a finite past domain; cannot extrapolate this profile to an arbitrarily old universe.',
                          'No measured environmental kappa map, source evolution, actual photon-companion energy exchange, stable storage or full background solution.'],
                amendment_sha256=digest(HERE/'time-candidate-amendment.md'))


def main():
    manifest={rel:digest(ROOT/rel) for rel in INPUT_PATHS}
    manifest['protocol.json']=digest(HERE/'protocol.json')
    manifest['observations.json']=digest(HERE/'observations.json')
    assert manifest[INPUT_PATHS[0]]=='8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0'
    assert manifest[INPUT_PATHS[4]]=='df793c3dca09ebfa7dbc5aa0ec1951daa8884431bc30eff28a710d7516cf50fa'
    result=dict(status='Scoped joint audit; no completed all-observation theory',
                alpha_per_mpc=ALPHA,redshift=redshift_distances(),timing=timing(),radio=radio(),
                thermal=thermal(),time_field=wave_timing(),synthetic=grid_and_checks(),
                local_and_multimessenger=local_and_multimessenger(),exact_time_candidate=exact_accumulating_time())
    for rel in INPUT_PATHS:
        assert digest(ROOT/rel)==manifest[rel]
    save('input-manifest.json',dict(inputs=manifest,script_sha256=digest(Path(__file__)),
                                  source_status='Earlier sources reused, newly transcribed summaries already inspected before protocol; not blind'))
    save('results.json',result)
    print(json.dumps(dict(redshift=result['redshift']['statistics'],timing=result['timing']['scores'][:3],
                          radio=result['radio']['models'],radio_p=result['radio']['fitted_p'],
                          jump_limit=result['radio']['conditional_max_fraction_per_independent_jump'],
                          thermal=result['thermal']['profiles'],time_field=result['time_field'],
                          synthetic=result['synthetic'],local=result['local_and_multimessenger'],
                          exact_time_candidate=result['exact_time_candidate']),indent=2))


if __name__=='__main__':
    main()
