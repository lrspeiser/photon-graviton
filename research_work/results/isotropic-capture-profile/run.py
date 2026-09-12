from pathlib import Path
import csv,json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid,quad,solve_ivp

HERE=Path(__file__).resolve().parent
R=10.;beta=1e-6;age=40.
results=[];profiles=[]
for count,order in [(321,48),(641,96)]:
    r=np.linspace(0,R,count);mu,wm=leggauss(order);z,w=leggauss(order);u=(z+1)/2;wu=w/2
    rr=r[:,None,None];mm=mu[None,:,None]
    sb=rr*mm+np.sqrt(R*R-rr*rr+rr*rr*mm*mm);s=sb*u[None,None,:]
    radii=np.sqrt(np.maximum(0,rr*rr+s*s-2*rr*mm*s))
    optical_base=sb[:,:,0]*np.sum(wu/(1+np.sqrt(1+radii*radii)),axis=2)
    for k0 in [.01,.1,1.,10.]:
        kappa=k0/(1+np.sqrt(1+r*r))
        q=4*np.pi*kappa*np.sum(np.exp(-k0*optical_base)*wm[None,:]/2,axis=1)
        enc=4*np.pi*cumulative_trapezoid(q*r*r,r,initial=0)
        integral=4*np.pi*cumulative_trapezoid(q*r,r,initial=0)
        psi=np.divide(enc,r,out=np.zeros_like(r),where=r>0)+integral[-1]-integral
        incident=4*np.pi*np.pi*R*R
        # Independent all-direction projected-chord integration.
        impact=R*u;half=np.sqrt(R*R-impact*impact)
        chord_r=np.sqrt(impact[:,None]**2+(half[:,None]*u[None,:])**2)
        tau=2*k0*half*np.sum(wu/(1+np.sqrt(1+chord_r*chord_r)),axis=1)
        absorbed=8*np.pi*np.pi*R*np.sum(wu*impact*(-np.expm1(-tau)))
        error=abs(enc[-1]/absorbed-1)
        if count==641:assert error<.002
        assert 0<=absorbed/incident<=1
        norm=beta/(age*psi[0]);mass=norm*age*enc[-1]
        fractions={str(a):float(np.interp(a,r,enc)/enc[-1]) for a in [1.,5.,9.]}
        snapshots=[]
        def piecewise(fn,edges,order):
            nodes,weights=leggauss(order);left=edges[:-1,None];width=np.diff(edges)[:,None]
            points=left+width*(nodes+1)/2
            return float(np.sum(fn(points)*weights*width/2))
        for b in [1.,5.,9.]:
            enclosed=float(np.interp(b,r,enc))*norm*age
            # Cylindrical mass from spherical shells, with a split at b.
            inside_edges=np.unique(np.r_[0,r[(r>0)&(r<b)],b])
            outside_edges=np.sqrt(np.unique(np.r_[b,r[r>b]])**2-b*b)
            inside=piecewise(lambda s:4*np.pi*s*s*np.interp(s,r,q),inside_edges,4)
            def projected_integrand(u):
                s=np.sqrt(b*b+u*u)
                return 4*np.pi*np.interp(s,r,q)*u*(s-u)
            outside=piecewise(projected_integrand,outside_edges,8)
            outside_check=piecewise(projected_integrand,outside_edges,16)
            projected_error=abs(outside-outside_check)/max(1e-100,inside+outside)
            assert projected_error<1e-9
            projected=norm*age*(inside+outside)
            snapshots.append(dict(radius=b,vc_km_s=float(299792.458*np.sqrt(enclosed/b)),deflection_arcsec=float(4*projected/b*180/np.pi*3600),projected_quadrature_relative_error=float(projected_error)))
        def psi_at(x):
            a=abs(x);return float(np.interp(a,r,psi)) if a<=R else float(enc[-1]/a)
        def phi(x,t):return -norm*(age+t)*psi_at(x)
        def ray(te):
            def rhs(x,y):
                n=np.exp(-2*phi(x,y[0]));pt=-norm*psi_at(x)
                return [n,-2*n*pt]
            sol=solve_ivp(rhs,[-20,20],[te,0.],rtol=1e-11,atol=1e-13,max_step=.025)
            assert sol.success;return sol.y[:,-1]
        to,logJ=ray(0.);S=np.exp(logJ+phi(20,to)-phi(-20,0))
        h=.001;early=ray(-h)[0];late=ray(h)[0]
        din=quad(lambda t:np.exp(phi(-20,t)),-h,h,epsabs=1e-14)[0]
        dout=quad(lambda t:np.exp(phi(20,t)),early,late,epsabs=1e-14)[0]
        cerr=abs(dout/din/S-1);assert cerr<1e-7
        row=dict(radial_nodes=count,quadrature_order=order,k0=k0,diameter_optical_depth=float(2*k0*(np.arcsinh(R)-R/(1+np.sqrt(1+R*R)))) ,
            captured_fraction=float(absorbed/incident),power_agreement=float(error),enclosed_fractions=fractions,outer_half_radius_fraction=1-fractions['5.0'],q9_over_q1=float(np.interp(9,r,q)/np.interp(1,r,q)),
            chosen_initial_depth=beta,chosen_age=age,required_intensity=float(norm),initial_effective_mass=float(mass),mass_growth_rate=float(norm*enc[-1]),snapshots=snapshots,
            measured_z=float(S-1),clock_error=float(cerr),shared_messenger_delay=0.)
        results.append(row)
        if count==641:
            profiles.extend(dict(k0=k0,radius=float(a),q_for_unit_intensity=float(b),enclosed_fraction=float(c/enc[-1]),potential_shape=float(d/psi[0])) for a,b,c,d in zip(r,q,enc,psi))
ref=[]
for k0 in [.01,.1,1.,10.]:
    lo=next(v for v in results if v['k0']==k0 and v['radial_nodes']==321);hi=next(v for v in results if v['k0']==k0 and v['radial_nodes']==641)
    err=max(abs(lo['enclosed_fractions'][k]-hi['enclosed_fractions'][k]) for k in hi['enclosed_fractions']);assert err<.003
    ref.append(dict(k0=k0,enclosed_fraction_difference=err,redshift_difference=abs(lo['measured_z']-hi['measured_z'])))
plummer={str(a):float((a**3/(1+a*a)**1.5)/(R**3/(1+R*R)**1.5)) for a in [1.,5.,9.]}
with (HERE/'profiles.csv').open('w',encoding='utf8',newline='') as fp:
    wr=csv.DictWriter(fp,fieldnames=list(profiles[0]),lineterminator='\n');wr.writeheader();wr.writerows(profiles)
(HERE/'results.json').write_text(json.dumps(dict(cases=results,refinement=ref,plummer_enclosed_fractions=plummer,scope='Isotropic absorption with fixed baryonic opacity, in-place retention and conditional weak-field response; no stellar fit or actual intensity calibration.'),indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps([v for v in results if v['radial_nodes']==641],indent=2))
