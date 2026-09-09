"""Tree-level joint energy/angle kernel for a subluminal scalar companion."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'slow-companion'

def coordinates(y,v):
    loss=2*v*y/(1+v)
    x=1-loss
    # Algebraically stable form of (1+x*x-loss*loss/v**2)/(2*x).
    one_minus_mu=2*(1-v)*y*y/((1+v)*x)
    return x,loss,1-one_minus_mu

def moments(v):
    def integrand(y,which):
        x,loss,mu=coordinates(y,v)
        base=x*x*(1+mu*mu)/2
        return base*{'count':1.,'loss':loss,'loss2':loss*loss,'turn':1-mu,'angle2':np.arccos(np.clip(mu,-1,1))**2}[which]
    raw={key:quad(lambda y:integrand(y,key),0,1,epsabs=1e-15,epsrel=1e-10)[0] for key in ['count','loss','loss2','turn','angle2']}
    width=2*v/(1+v)
    return {'v_over_photon_speed':v,'minimum_surviving_energy_fraction':(1-v)/(1+v),
            'mean_fractional_loss_per_event':raw['loss']/raw['count'],
            'local_loss_second_to_first_moment':raw['loss2']/raw['loss'],
            'mean_one_minus_cos_theta':raw['turn']/raw['count'],
            'rms_deflection_degrees':np.sqrt(raw['angle2']/raw['count'])*180/np.pi,
            'rate_over_g_squared_E_cubed':width*raw['count']/(16*np.pi*v*v),
            'fractional_loss_rate_over_g_squared_E_cubed':width*raw['loss']/(16*np.pi*v*v)}

def main():
    rows=[moments(v) for v in [1e-6,1e-4,.01,.1,.5,.9,.99,.999999]]
    assert abs(rows[0]['mean_fractional_loss_per_event']/1e-6-10/11)<1e-5
    assert abs(rows[0]['mean_one_minus_cos_theta']-46/77)<1e-5
    assert abs(rows[0]['local_loss_second_to_first_moment']/1e-6-46/35)<1e-5
    assert abs(rows[-1]['mean_fractional_loss_per_event']-.25)<1e-5
    assert abs(rows[0]['fractional_loss_rate_over_g_squared_E_cubed']-1/(12*np.pi))<1e-6
    assert abs(rows[-1]['fractional_loss_rate_over_g_squared_E_cubed']-1/(192*np.pi))<1e-6
    rng=np.random.default_rng(40108);errors=[];ward=[];polarization=[];jacobian=[]
    for _ in range(300):
        v=float(rng.uniform(.001,.999));y=float(rng.uniform(.001,.999));E=float(rng.uniform(.1,10))
        x,loss,mu=coordinates(y,v);sin=np.sqrt(max(0,1-mu*mu))
        kin=np.array([0.,0.,E]);kout=E*x*np.array([sin,0.,mu]);q=kin-kout
        errors.append(abs(v*np.linalg.norm(q)-E*loss)/E)
        # Full electric-field vertex vanishes when either polarization is gauge.
        electric=lambda energy,momentum,e0,ev:energy*np.asarray(ev)-momentum*e0
        ward.append(float(np.linalg.norm(electric(E,kin,E,kin))))
        ei=[np.array([1.,0.,0.]),np.array([0.,1.,0.])]
        eo=[np.array([mu,0.,-sin]),np.array([0.,1.,0.])]
        pol=sum(np.dot(a,b)**2 for a in ei for b in eo)/2
        polarization.append(abs(pol-(1+mu*mu)/2))
        dm=1e-6
        w=lambda z:v*E*np.sqrt(1+x*x-2*x*z)
        # Near endpoint roots the finite-difference scale must be reduced.
        if 1+x*x-2*x*(mu+dm)>1e-4:
            numeric=(w(mu+dm)-w(mu-dm))/(2*dm)
            analytic=-v*v*E*x/loss
            jacobian.append(abs(numeric/analytic-1))
    assert max(errors)<1e-12 and max(ward)<1e-12 and max(polarization)<1e-12
    assert max(jacobian)<1e-4
    # Analytic phase-space normalization: integrate in x and compare y transform.
    phase=[]
    for row in rows[2:7]:
        v=row['v_over_photon_speed'];xmin=(1-v)/(1+v)
        def kernel(x):
            mu=(1+x*x-(1-x)**2/v**2)/(2*x)
            return x*x*(1+mu*mu)/2
        count=quad(kernel,xmin,1,epsabs=1e-12)[0]/(16*np.pi*v*v)
        phase.append(abs(count/row['rate_over_g_squared_E_cubed']-1))
    assert max(phase)<1e-10
    detailed_balance=[]
    for energy,ratio,temp in [(1.,.2,.5),(2.,.7,1.),(.1,.99,.2)]:
        f=lambda e:1/np.expm1(e/temp)
        forward=f(energy)*(1+f(energy*ratio))*(1+f(energy*(1-ratio)))
        backward=f(energy*ratio)*f(energy*(1-ratio))*(1+f(energy))
        detailed_balance.append(abs(forward/backward-1))
    assert max(detailed_balance)<1e-12
    result={'scope':'Preferred-frame scalar/electric dimension-five EFT, tree-level spontaneous gamma -> gamma + scalar, initially empty companion modes, all outgoing directions. No cosmological fit, physical medium or graviton identity is established.',
            'action':'L=(1+g phi) E_field^2/2-B^2/2+(dot phi)^2/2-v^2|grad phi|^2/2',
            'domain':'0<v<1, 1+g phi>0, photon and companion momenta well below EFT cutoff, perturbatively small rate/energy, preferred frame fixed; full matter/gravity completion absent',
            'kernel':'dGamma/dx=g^2 E^3 x^2(1+mu(x)^2)/(32 pi v^2), xmin<=x<=1',
            'scenarios':rows,'fractional_loss_ratio_at_double_photon_energy_same_g_v':8,
            'checks':{'random_energy_momentum_max_relative_error':max(errors),'gauge_vertex_max':max(ward),
                      'polarization_sum_max_error':max(polarization),'phase_space_coordinate_max_relative_error':max(phase),
                      'delta_jacobian_max_relative_error':max(jacobian),
                      'bosonic_thermal_detailed_balance_max_error':max(detailed_balance)},
            'verdict':'This specific scale-free coupling does not produce achromatic fractional energy loss. Slow companions make small energy packets but broad photon deflections; speed approaching photon speed narrows angles but leaves mean loss approaching one quarter.',
            'remaining':['Full interaction-derived environment and cutoff','Matter clocks, background stability and radiative corrections','Stimulated/inverse processes for a populated companion bath','Capture, support, gravity and source history']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'slow-companion-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
