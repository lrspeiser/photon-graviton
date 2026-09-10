"""Check rotating canonical dynamics against an exact harmonic orbit."""
import numpy as np
from run import integrate,rotate,OMEGA,save

class Harmonic:
    frequency=123.
    def evaluate(self,x):
        return .5*self.frequency**2*np.sum(x*x,axis=1),-self.frequency**2*x

field=Harmonic()
y0=np.array([[1.,.2,.3,10.,160.,20.],[2.,-.4,.1,-30.,100.,-10.]])
times=np.linspace(0,.25,501)
t=times[:,None,None];w=field.frequency
inertial_x=y0[None,:,:3]*np.cos(w*t)+y0[None,:,3:]*np.sin(w*t)/w
inertial_p=-w*y0[None,:,:3]*np.sin(w*t)+y0[None,:,3:]*np.cos(w*t)
exact_x=rotate(inertial_x,-OMEGA*times[:,None])
exact_p=rotate(inertial_p,-OMEGA*times[:,None])
out,_=integrate(field,y0,2e-11,times)
dx=float(np.max(np.linalg.norm(out[:,:,:3]-exact_x,axis=2)))
dv=float(np.max(np.linalg.norm(out[:,:,3:]-exact_p,axis=2)))
assert dx<1e-7 and dv<1e-5
save('analytic-control.json',dict(known_control='Isotropic harmonic oscillator, exact inertial solution rotated analytically',
                                frequency_kms_per_kpc=w,bar_omega_kms_per_kpc=OMEGA,
                                max_position_error_kpc=dx,max_velocity_error_kms=dv,
                                astronomical_validation=False))
print('Exact harmonic rotating-orbit control passed:',dx,'kpc;',dv,'km/s')
