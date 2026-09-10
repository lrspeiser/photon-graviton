"""Established SBF extinction sensitivity; symbolic checks, no catalog edits."""
import json
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
m, color, E, RV, RI, a, b, pivot, dE, da, dEcal = s.symbols(
    'm color E R_V R_I a b pivot delta_E delta_a delta_E_cal', real=True)
mu = m-RI*E-(a+b*(color-(RV-RI)*E-pivot))
K = b*(RV-RI)-RI
derivative = s.simplify(s.diff(mu,E))
change = s.simplify(mu.subs({E:E+dE,a:a+da}, simultaneous=True)-mu)
assert s.simplify(derivative-K) == 0
assert s.simplify(change-(K*dE-da)) == 0
recalibrated = s.simplify(change.subs(da,K*dEcal))
assert s.simplify(recalibrated-K*(dE-dEcal)) == 0
assert s.simplify(recalibrated.subs(dE,dEcal)) == 0
out = {'status':'symbolic established-calibration sensitivity; no proposed gravity law or applied distance correction',
       'formula_provenance':'established SBF/color-extinction structure; algebraic restatement of known sensitivity, not original research',
       'assumptions':['fixed passband extinction coefficients','fixed linear color-calibration slope',
                      'same sensitivity K for target and calibrator in cancellation example',
                      'fixed observed fluctuation magnitude and observed color',
                      'no change in calibrator adopted distance in zero-point recalibration example'],
       'mu':str(mu),'d_mu_d_E':str(derivative),'finite_change_mu':str(change),
       'after_zero_point_recalibration':str(recalibrated),
       'common_shift_cancels':True,'symbolic_identities_verified':4}
(HERE/'sbf-extinction-response.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
print(json.dumps(out,indent=2))
