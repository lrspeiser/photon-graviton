"""Stage3: finite generated field, frozen-snapshot Hamilton rays and material force."""
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.interpolate import RectBivariateSpline
from common import save
from spatial import simulate

IMPACTS=[-2.,-1.,-.5,.5,1.,2.]

def rays(model,refined=False,zero=False):
    phi=model.phi*0 if zero else model.phi
    spline=RectBivariateSpline(model.axis,model.axis,phi,kx=3,ky=3,s=0)
    g=model.g;rows=[]
    def value(x,y,dx=0,dy=0):return float(spline.ev(x,y,dx=dx,dy=dy))
    def rhs(x,z):
        y,px,py,t=z
        norm=np.hypot(px,py);a=np.exp(-g*value(x,y))
        return [py/px,g*norm*norm/px*value(x,y,dx=1),
                g*norm*norm/px*value(x,y,dy=1),norm/(a*px)]
    for impact in IMPACTS:
        initial=[impact,np.exp(g*value(-8,impact)),0.,0.]
        sol=solve_ivp(rhs,[-8,8],initial,method='DOP853',
                      rtol=2e-11 if refined else 2e-9,
                      atol=2e-13 if refined else 2e-11,
                      max_step=.025 if refined else .1,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        x=np.linspace(-8,8,641);path=sol.sol(x);y,px,py,t=path
        h=np.exp(-g*spline.ev(x,y))*np.hypot(px,py)
        born=g*simpson(spline.ev(x,np.full_like(x,impact),dy=1),x=x)
        rows.append(dict(impact=impact,angle=np.arctan2(py[-1],px[-1]),
                         born_angle=born,ray_energy_error=np.max(abs(h-1)),
                         travel_time=t[-1],exit_y=y[-1],x=x[::8],y=y[::8]))
    return rows

def bend_metrics(rows):
    angles=np.array([r['angle'] for r in rows]);born=np.array([r['born_angle'] for r in rows])
    scale=max(np.max(abs(angles)),1e-30)
    return dict(max_bend=scale,mirror_relative=np.max(abs(angles+angles[::-1]))/scale,
                born_relative=np.max(abs(angles-born))/scale,
                all_attractive=bool(np.all(angles*np.array(IMPACTS)<0)),
                signed_inward=-angles*np.sign(IMPACTS))

def run(out):
    cases={};rays_by_case={}
    for n,b,energy,label in [(128,2,1.,'primary-128'),(192,2,1.,'primary-192'),
                            (256,2,1.,'primary-256'),(192,0,1.,'protected-clock'),
                            (192,1,1.,'optical-clock'),(192,2,0.,'source-off')]:
        result,model=simulate(n=n,dim=2,b=b,energy=energy,wavelength=3.)
        rr=rays(model);result['rays']=rr;result['bend']=bend_metrics(rr)
        radius=np.sqrt(sum(c*c for c in model.coords))
        result['central_fraction_of_receiving']=float(
            np.sum(model.receiving_density()[radius<3])*model.dv/max(model.sectors()[1],1e-30))
        # Same-position probes isolate composition dependence from position effects.
        points=np.array([[0.,1.],[0.,1.]])
        ph,grad,_=model.interpolation(model.phi,points)
        actions=np.array([.0005,.001]);omega=np.exp(-model.b*model.g*ph)
        acceleration=model.b*model.g*(actions*omega)[:,None]*grad
        result['same_position_probe_accelerations']=acceleration
        result['force_ratio_if_nonzero']=2. if np.max(abs(acceleration))>1e-20 else None
        if n==256:
            refined=rays(model,refined=True);straight=rays(model,zero=True)
            scale=result['bend']['max_bend']
            result['ray_refinement_relative']=max(abs(a['angle']-b['angle']) for a,b in zip(rr,refined))/scale
            result['zero_field_max_bend']=max(abs(a['angle']) for a in straight)
            result['refined_rays']=refined;result['zero_field_rays']=straight
            np.savez_compressed(out/'primary-final-field.npz',axis=model.axis,A=model.A,D=model.D,
                                phi=model.phi,Pi=model.Pi,X=model.X,P=model.P)
        save(out/(label+'.json'),result)
        cases[label]={k:v for k,v in result.items() if k not in
            ['snapshots','time','energy_history','momentum_history','central_receiving_energy']}
        rays_by_case[label]=rr
        print('2D',label,'energy',result['energy_relative_error'],
              'max bend',result['bend']['max_bend'],flush=True)
    coarse=cases['primary-192'];fine=cases['primary-256']
    scale=fine['bend']['max_bend']
    spatial=max(abs(a['angle']-b['angle']) for a,b in zip(coarse['rays'],fine['rays']))/scale
    weak_applicable=fine['config']['g']*fine['max_final_phi']<.02
    gates=dict(energy=fine['energy_relative_error']<1e-3,spatial_bending=spatial<.05,
               ray_refinement=fine['ray_refinement_relative']<1e-4,
               born_weak=(fine['bend']['born_relative']<.02 if weak_applicable else True),
               mirror=fine['bend']['mirror_relative']<.01,
               straight=fine['zero_field_max_bend']==0)
    physics=dict(nonzero_bending=scale>1e-10,attractive_bending=fine['bend']['all_attractive'],
                 net_conversion=fine['receiving_gain']>.001,
                 photon_exclusive=cases['source-off']['max_final_phi']<1e-10,
                 universal_matter_force=False,three_dimensional_gravity=False)
    return dict(stage='two_dimensional',gates=gates,numerical_pass=all(gates.values()),
                mechanism_gates=physics,mechanism_pass=all(physics.values()),
                spatial_bend_relative_change=spatial,born_weak_applicable=weak_applicable,
                cases=cases,
                note='Frozen finite periodic 2D snapshots, not evolving or isolated 3D astrophysical rays. Energy is fully counted; no input receiving field or halo.')
