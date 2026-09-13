"""Necessary local photon supply for an optional homogeneous completion."""
from pathlib import Path
import json,hashlib
import sympy as sp
import numpy as np
HERE=Path(__file__).resolve().parent
p=HERE/'optical-source-results.json';source=json.loads(p.read_text())
f,q=sp.symbols('f q')
H=-(1-f)*sp.log(1-f)*sp.sqrt(1+f/(1+q*f))
series=sp.series(H,f,0,5).removeO().expand()
# R/alpha^2 = -6*h3 - 12*h4*f + O(f^2).
# u is proportional to R*(1-f)^2 and df/dz=1 at the observer.
r0=-6*series.coeff(f,3);r1=-12*series.coeff(f,4)
factor=sp.simplify(2-r1/r0)
value=float(factor.subs(q,source['q']))
u0=source['observer']['all_companion_u_J_m3'];required=value*u0
c=299792458.;kB=1.380649e-23;h=6.62607015e-34;Mpc=3.085677581491367e22
T=2.72548;dT=.00057
a=8*np.pi**5*kB**4/(15*h**3*c**3)
cmb=a*T**4
benchmarks=[]
for I in [45.,100.,170.]:
    ebl=4*np.pi*I*1e-9/c
    benchmarks.append(dict(UV_mm_intensity_nW_m2_sr=I,UV_mm_u_J_m3=ebl,
                           microwave_plus_UV_mm_u_J_m3=cmb+ebl,
                           required_to_benchmark_ratio=required/(cmb+ebl)))
out=dict(status='Conditional homogeneous temporal completion; not a general exclusion',
         source_input_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
         analytic_required_u_gamma_over_u_companion=str(factor),
         required_u_gamma_over_u_companion=value,
         required_photon_u_lower_bound_J_m3=required,
         required_growth_lower_source_W_m3=source['alpha_per_Mpc']*c/Mpc*required,
         required_bolometric_intensity_nW_m2_sr=c*required/(4*np.pi)*1e9,
         microwave_temperature_K=T,microwave_temperature_uncertainty_K=dT,
         microwave_u_J_m3=cmb,microwave_u_temperature_uncertainty_J_m3=4*cmb*dT/T,
         benchmarks=benchmarks,
         source_references=['https://arxiv.org/abs/0911.1955','https://ned.ipac.caltech.edu/level5/March03/Dwek2/Dwek6.html'])
(HERE/'optical-growth-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
