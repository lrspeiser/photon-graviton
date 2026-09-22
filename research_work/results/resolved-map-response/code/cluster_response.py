"""Apply a LOCAL source-unit susceptibility to observed X-COP density shells.

Not a cluster gravity fit: density reference is an uncalibrated common physical
parameter. The thermal dependence of exchange remains unspecified. No imported
hydrostatic mass, NFW profile or lensing reconstruction enters this calculation.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

def coefficients(h):
    kp=.5*np.exp(2*h);km=.5*np.exp(4*h*h)
    C_per_J=kp/(.5*km+.05*kp+.5*.05)
    dCdh=C_per_J*.5*((2-8*h)*km+2*.05)/(.5*km+.05*kp+.5*.05)
    return kp,km,C_per_J,dCdh

def main():
    ap=argparse.ArgumentParser(__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    if a.output.exists():raise FileExistsError('Preserve old results')
    data=json.loads(a.input.read_text());crit=brentq(lambda h:(8*h-2)*.5*np.exp(4*h*h)-.1,.25,.5);q=crit/(1-crit)
    h=np.linspace(0,1,1001);kp,km,c,d=coefficients(h);eps=1e-6;check=(coefficients(h[1:-1]+eps)[2]-coefficients(h[1:-1]-eps)[2])/(2*eps)
    rows=[]
    for name,cl in data['clusters'].items():
        den=cl['density'];r=.5*(np.array(den['r_in_kpc'])+np.array(den['r_out_kpc']));n=np.array(den['ne_cm3']);settings=[]
        for n0 in [1e-4,1e-3,1e-2]:
            hh=n/(n+n0);_,_,cc,der=coefficients(hh);cross=[]
            for j in range(len(n)-1):
                if (n[j]-q*n0)*(n[j+1]-q*n0)<0:
                    t=(np.log(q*n0)-np.log(n[j]))/(np.log(n[j+1])-np.log(n[j]));cross.append(float(np.exp(np.log(r[j])+t*(np.log(r[j+1])-np.log(r[j])))))
            settings.append(dict(reference_density_cm3=n0,threshold_density_cm3=q*n0,positive_shells=int((der>0).sum()),negative_shells=int((der<0).sum()),log_interpolated_crossing_radii_kpc=cross,min_C_per_source=float(cc.min()),max_C_per_source=float(cc.max()),C_per_source=cc.tolist(),log_density_derivative=(der*hh*(1-hh)).tolist()))
        rows.append(dict(cluster=name,density_shells=len(r),midpoint_kpc=r.tolist(),ne_cm3=n.tolist(),settings=settings))
    out=dict(scope=__doc__,input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),h_turnover=crit,n_over_n0_turnover=q,max_local_C_per_J=float(coefficients(crit)[2]),zero_gas_C_per_J=float(coefficients(0.)[2]),maximum_local_gain_over_zero=float(coefficients(crit)[2]/coefficients(0.)[2]-1),derivative_absolute_error=float(np.max(abs(check-d[1:-1]))),clusters=rows)
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('Turnover',crit,q,'max_gain',out['maximum_local_gain_over_zero'])
    for r in rows:print(r['cluster'],[(s['reference_density_cm3'],[round(x) for x in s['log_interpolated_crossing_radii_kpc']]) for s in r['settings']])
if __name__=='__main__':main()
