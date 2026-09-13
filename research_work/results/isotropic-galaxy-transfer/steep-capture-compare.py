"""Compare frozen capture predictions without selecting new parameters."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
old=json.loads((HERE/'predictions.json').read_text())
new=json.loads((HERE/'steep-capture-predictions.json').read_text())
lookup={(r['model'],r['galaxy']):r for r in old}
rows=[]
for r in new:
    before=lookup[r['model'],r['galaxy']]
    assert before['R_kpc']==r['R_kpc'] and before['observed_kms']==r['observed_kms']
    y=np.array(r['observed_kms']);p=np.array(r['predicted_kms']);p0=np.array(before['predicted_kms'])
    rows.append(dict(galaxy=r['galaxy'],model=r['model'],split=r['split'],
                     previous_RMSE_kms=float(np.sqrt(np.mean((p0-y)**2))),
                     revised_RMSE_kms=float(np.sqrt(np.mean((p-y)**2))),
                     previous_log_MSE=float(np.mean(np.log10(p0/y)**2)),
                     revised_log_MSE=float(np.mean(np.log10(p/y)**2))))
summary=[]
for model in ['transparent_control','attenuated']:
    for split in ['train','validation','test']:
        selected=[r for r in rows if r['model']==model and r['split']==split]
        summary.append(dict(model=model,split=split,n=len(selected),
                            galaxies_with_lower_kms_error=sum(r['revised_RMSE_kms']<r['previous_RMSE_kms'] for r in selected),
                            galaxies_with_lower_log_error=sum(r['revised_log_MSE']<r['previous_log_MSE'] for r in selected)))
checks=[]
for B,t in [(.2,-2.),(1.,0.),(1.,3.),(4.,-7.),(4.,7.)]:
    I2=t/(2*B**2*(B**2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3)
    exact=t/(4*B**2*(B**2+t*t)**2)+3*I2/(4*B**2)
    numeric=quad(lambda u:(B*B+u*u)**-3,-np.inf,t,epsabs=1e-12,epsrel=1e-11)[0]
    error=abs(exact/numeric-1)
    assert error<1e-8
    checks.append(dict(B=B,t=t,analytic=exact,quadrature=numeric,relative_error=error))
out=dict(summary=summary,rows=rows,upstream_integral_checks=checks)
(HERE/'steep-capture-comparison.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
