from pathlib import Path
import json
from model import *
# Scalar virial theorem for an isotropic self-gravitating Hernquist sphere.
cb,_=aperture_coefficients(3.,np.inf,0.,0.,.46,20000.,20.)
exact=G*1e11/(18*3);error=abs(cb/exact-1);assert error<1e-4,error
# Representative finite aperture plus Gaussian seeing: refine both integrations.
c=aperture_coefficients(3.,4.,1.5,.24,.46,22000.,20.)
f=aperture_coefficients(3.,4.,1.5,.24,.46,22000.,20.,n=8001,order=256)
change=max(abs(np.array(c)/f-1));assert change<.001,change
# Independent analytic projected Hernquist mass checks the Einstein root.
mass=3e11;a=3.;Dl=1e6;ratio=.6
angle_value=angle(mass,a,Dl,ratio,0.,.46,22000.,20.)
b=angle_value/ARCSEC*Dl;u=b/a
Q=np.arccosh(1/u)/np.sqrt(1-u*u) if u<1 else np.arccos(1/u)/np.sqrt(u*u-1)
projected=mass*u*u*(Q-1)/(1-u*u)
expected=ratio*4*G*projected/(C*C*b)*ARCSEC
lens_error=abs(angle_value/expected-1);assert lens_error<1e-7
Path(__file__).with_name('verification.json').write_text(json.dumps(dict(hernquist_lens_root_relative_error=lens_error,global_virial_relative_error=error,aperture_coefficients_relative_refinement=float(change),scope='Numerical checks only, not galaxy model validity'),indent=2)+'\n',encoding='utf-8',newline='\n')
print(error,change)
