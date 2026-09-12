"""Real-filter stress test of wavelength-dependent source durations."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
base=ROOT/'research_work/generated/des-photometric-calibration/selected/filters/DES/DES-SN3YR_DECam'
manifest=json.loads((HERE.parent/'des-photometric-calibration/dependencies.json').read_text(encoding='utf-8'))
for v in manifest['selected']:
    if '/filters/' in '/'+v['archive_path']:assert hashlib.sha256((ROOT/v['cache_path']).read_bytes()).hexdigest()==v['sha256']
rows=[]
for band in 'griz':
    lam,T=np.loadtxt(base/f'DECam_{band}.dat',unpack=True)
    for b in [0.,1.]:
        for q in [-.5,0.,.5,1.]:
            widths={mode:[] for mode in ['fixed_observer_filter','matched_emission_filter']}
            for S in [1.,1.1,1.3,1.6,2.2]:
                for mode in widths:
                    observed_lam=lam if mode=='fixed_observer_filter' else S*lam
                    emitted_lam=observed_lam/S
                    assert emitted_lam.min()>1000 and emitted_lam.max()<20000
                    def count(t):
                        emitted_width=30*(emitted_lam/5500)**q
                        source=(emitted_lam/5500)**-2*np.exp(-np.log(2)*(2*t/S**b/emitted_width)**2)
                        # Relative photon counts: h*c and fixed geometric area cancel in ratios.
                        return float(np.trapezoid(T*observed_lam*source/S**(b+2),observed_lam))
                    peak=count(0);half=brentq(lambda t:count(t)/peak-.5,0,1000,xtol=1e-11)
                    width=2*half;widths[mode].append(width)
                    exponent=b-q if mode=='fixed_observer_filter' else b
                    expected=widths[mode][0]*S**exponent
                    error=abs(width/expected-1);assert error<1e-9
                    rows.append({'band':band,'b':b,'q':q,'S':S,'mode':mode,'observed_fwhm_days':width,'expected_exponent':exponent,'relative_width_scaling_error':error,'relative_count_peak':peak})
result={'scope':'Artificial wavelength-dependent source integrated through verified DES responses. No actual SN flux or physical source model.', 'source_model':'L_lambda proportional to (lambda_e/5500A)^-2 times exp[-ln(2)*(2*t_e/(30day*(lambda_e/5500A)^q))^2] over 1000-20000A, zero outside. Every tested passband maps strictly inside that range.', 'postulates':'Spectral stretch S; event stretch S^b; number conservation, fixed distance. F_lambda,o=S^-(b+2)*L_lambda,e(lambda_o/S,t_o/S^b) up to constant geometric dilution.', 'rows':rows,'maximum_relative_scaling_error':max(r['relative_width_scaling_error'] for r in rows),'observations_tested':False}
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'calculations':len(rows),'maximum_error':result['maximum_relative_scaling_error']}))
