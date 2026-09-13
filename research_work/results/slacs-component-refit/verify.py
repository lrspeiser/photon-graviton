from pathlib import Path
import json,numpy as np
from model import deproject,ComponentModel,G,C,ARCSEC
from scipy.special import gamma,gammainc
r=np.geomspace(.01,3,100);Re=1.2;bn=.6725;amp=2.;k=bn/Re**2
nu,total=deproject(r,[dict(R=Re,n=.5,amp=amp,bn=bn)])
exact=amp*np.exp(bn)*np.sqrt(k/np.pi)*np.exp(-k*r*r)
err=float(np.max(abs(nu/exact-1)));assert err<1e-9
r=np.geomspace(1e-5,1e5,4001);comps=[dict(R=2.,n=4.,amp=2.,bn=7.6697),dict(R=8.,n=1.,amp=.1,bn=1.6721)]
nu,total=deproject(r,comps);ref,_=deproject(r,comps,512)
keep=nu>nu.max()*1e-10;change=float(np.max(abs(nu[keep]/ref[keep]-1)));assert change<1e-5
mass=4*np.pi*np.trapezoid(r*r*nu,r);masserr=float(abs(mass/total-1));assert masserr<1e-4
m=ComponentModel(3.,np.array([0.,1.,2.]),.8,.24,.46,22000.,20.,comps)
theta=m.angle(3e11,1e6,.6,False);b=theta/ARCSEC*1e6
lum=np.array([2*np.pi*c['amp']*c['R']**2*c['n']*np.exp(c['bn'])*gamma(2*c['n'])/c['bn']**(2*c['n']) for c in comps])
frac=sum(L*gammainc(2*c['n'],c['bn']*(b/c['R'])**(1/c['n'])) for L,c in zip(lum,comps))/lum.sum()
expected=.6*4*G*3e11*frac/(C*C*b)*ARCSEC
lens_error=float(abs(theta/expected-1));assert lens_error<1e-4
Path(__file__).with_name('verification.json').write_text(json.dumps(dict(projected_sersic_lens_relative_error=lens_error,gaussian_abel_relative_error=err,deprojection_order_relative_change=change,spherical_mass_relative_error=masserr,scope='Known analytic Gaussian Abel inverse and two-component quadrature checks'),indent=2)+'\n',encoding='utf-8',newline='\n')
print(err,change,masserr)
