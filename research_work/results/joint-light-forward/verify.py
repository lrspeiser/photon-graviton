"""Independent limiting-case checks for the combined forward calculator."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from transport import predict,load_filters,AB0_FNU,C_ANGSTROM_S,MPC_CM
HERE=Path(__file__).resolve().parent
filters=load_filters();bands=np.array(list('griz'));distance=1e-5
LNU=4*np.pi*(distance*MPC_CM)**2*AB0_FNU

def constant(t,lam):return np.broadcast_to(np.where((lam>=1000)&(lam<=20000),LNU*C_ANGSTROM_S/lam**2,0.),np.broadcast_shapes(t.shape,lam.shape))
base=predict(np.zeros(4),bands,distance,constant,alpha_per_mpc=0,filters=filters)
expected=np.array([1e11*10**(-.4*filters[b][3]) for b in bands])
normalization=float(np.max(abs(base['native_fluxcal']/expected-1)));assert normalization<1e-12
far=predict(np.zeros(4),bands,2*distance,constant,alpha_per_mpc=0,filters=filters)
assert np.allclose(far['native_fluxcal'],base['native_fluxcal']/4,rtol=1e-12)

def pulse(t,lam):return constant(t,lam)*np.exp(-np.log(2)*(2*t/30)**2)
rows=[]
for b in [0.,1.]:
    for S in [1.,1.2,2.2]:
        alpha=np.log(S)/100
        for band in bands:
            def value(time):return predict([time],[band],100,pulse,alpha_per_mpc=alpha,stretch_exponent=b,filters=filters)['native_fluxcal'][0]
            peak=value(0);half=brentq(lambda t:value(t)/peak-.5,0,200)
            reference=predict([0],[band],100,pulse,alpha_per_mpc=0,stretch_exponent=b,filters=filters)['native_fluxcal'][0]
            assert abs(2*half/(30*S**b)-1)<1e-10
            assert abs(peak/reference*S**b-1)<1e-12
            rows.append({'b':b,'S':S,'band':str(band),'fwhm_days':2*half,'peak_ratio':float(peak/reference)})
example=predict([0]*4,bands,30.660139,constant,filters=filters)
for key in ['photon_rate_per_cm2_s','ab_equivalent_fluxcal','native_fluxcal']:example[key]=example[key].tolist()
out={'scope':'Combined forward-model numerical checks with artificial flat-fnu sources. No SN observation fitted.','AB0_at_10pc_max_relative_error':normalization,'inverse_square_check':True,'pulse_checks':rows,'illustrative_100_million_light_year_source':example,'cause_derived':False,'source_luminosity_is_artificial':True}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'normalization_error':normalization,'pulse_cases':len(rows),'example_redshift':example['predicted_redshift']}))
