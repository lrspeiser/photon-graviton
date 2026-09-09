"""Synthetic connected, periodic 1D Hamiltonian field + ray calculation.
Dimensionless units: c=box length=K=1. This is not an astronomical fit.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
P=Path(__file__).resolve().parent;O=P/'results'

def setup(N):
 dx=1/N;x=np.arange(N)*dx;eps=1e-4;a=10.
 # Smooth dense region around the periodic boundary; a void in the middle.
 density=.5*(1-np.tanh((np.minimum(x,1-x)-.16)/.025))
 b2=100**2*density
 lap=(np.roll(np.eye(N),1,axis=1)+np.roll(np.eye(N),-1,axis=1)-2*np.eye(N))/dx**2
 L=np.block([[-lap+np.diag(a*a+b2),-a*a*np.sqrt(eps)*np.eye(N)],[-a*a*np.sqrt(eps)*np.eye(N),-lap+eps*a*a*np.eye(N)]])
 # ψ propagation speed is also c; mass-weighting gives the same Laplacian.
 rate=.007731496595524618
 vel=np.r_[rate*a*a/(a*a+b2),np.full(N,rate/np.sqrt(eps))]
 return dx,L,vel,x,density

def weights(pos,N):
 u=(pos%1)*N;i=int(np.floor(u));t=u-i
 ids=np.array([i-1,i,i+1,i+2])%N
 w=np.array([(1-t)**3,3*t**3-6*t*t+4,-3*t**3+3*t*t+3*t+1,t**3])/6
 dw=np.array([-3*(1-t)**2,9*t*t-12*t,-9*t*t+6*t+3,3*t*t])/6*N
 return ids,w,dw

def test_field(N,launch=0.):
 dx,L,v,x,density=setup(N);w2,U=eigh(L);w=np.sqrt(np.maximum(w2,0));vv=U.T@v
 def field(t):
  q=U@(vv*t*np.sinc(w*t/np.pi));dq=U@(vv*np.cos(w*t));return q,dq
 def ray(t,y):
  q,dq=field(t);ids,ww,dw=weights(y[0],N);chi=ww@q[ids];return [np.exp(-chi),-(ww@dq[ids])]
 def finish(t,y):return y[0]-.9
 finish.terminal=True;finish.direction=1
 ans=[]
 for start in [launch,launch+1e-5]:
  sol=solve_ivp(ray,[start,start+3.],[.1,0],events=finish,rtol=1e-10,atol=1e-12,max_step=.005)
  assert sol.t_events[0].size
  ans.append({'arrival':float(sol.t[-1]),'log_energy_ratio':float(sol.y[1,-1])})
 S=np.exp(-ans[0]['log_energy_ratio']);pulse=(ans[1]['arrival']-ans[0]['arrival'])/1e-5
 return {'N':N,'launch_time':launch,'frequency_stretch':float(S),'pulse_stretch':float(pulse),'pulse_vs_frequency_relative_difference':float(pulse/S-1),'arrival_time':ans[0]['arrival'],'initial_field_energy':float(.5*dx*(v@v)),'minimum_frequency_squared':float(w2[0])}

def coupled(N,energy_fraction):
 dx,L,v,x,density=setup(N);E0=.5*dx*(v@v);p0=energy_fraction*E0
 y0=np.r_[np.zeros(2*N),v,.1,p0]
 def rhs(t,y):
  q=y[:2*N];vel=y[2*N:4*N];pos,mom=y[-2:];ids,w,dw=weights(pos,N);chi=w@q[ids];speed=np.exp(-chi);energy=mom*speed
  # Conservative field-photon source from derivative of H_gamma=p exp(-chi).
  force=-L@q;np.add.at(force,ids,energy*w/dx)
  return np.r_[vel,force,speed,energy*(dw@q[ids])]
 def event(t,y):return y[-2]-.9
 event.terminal=True;event.direction=1
 sol=solve_ivp(rhs,[0,1.3],y0,method='DOP853',events=event,rtol=2e-9,atol=2e-11,max_step=.002)
 assert sol.success and sol.t_events[0].size
 totals=[];ph=[]
 for y in sol.y.T:
  q=y[:2*N];vel=y[2*N:4*N];ids,w,dw=weights(y[-2],N);ep=y[-1]*np.exp(-(w@q[ids]));ef=.5*dx*(vel@vel+q@L@q);totals.append(ep+ef);ph.append(ep)
 return {'N':N,'initial_photon_to_field_energy':energy_fraction,'frequency_stretch':float(ph[0]/ph[-1]),'photon_energy_change':float(ph[-1]-ph[0]),'total_relative_energy_error':float(np.max(np.abs(np.array(totals)/totals[0]-1))),'arrival':float(sol.t[-1]),'solver_evaluations':sol.nfev}
if __name__=='__main__':
 out={'test_ray':[test_field(n) for n in [81,161,321]],'launch_epochs':[test_field(161,t) for t in [10.,20.,40.]],'backreaction':[coupled(n,f) for n in [81,161] for f in [.001,.1]],'definition':'Periodic 1D box, smooth fixed-density slabs, positive gradient and quadratic field energies, two fields with epsilon=1e-4; cubic B-spline interpolation and its adjoint deposition. Local particle Hamiltonian Hgamma=p exp(-chi).','limitations':['Synthetic geometry and dimensionless energy ratios, not data calibration.','Photon packet is a 1D ray surrogate, not a complete electromagnetic or gravitational-wave action.','Matter density is fixed; field plus ray energy conserved but matter motion and gravity omitted.','No material clock law or absolute cosmological energy scale is derived.','Pulse-stretch equality tested only in the negligible-backreaction limit.']}
 (O/'connected_field_grid_tied_packet.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
