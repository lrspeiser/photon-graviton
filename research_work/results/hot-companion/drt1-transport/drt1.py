#!/usr/bin/env python3
"""DRT-1: reproducible hot-companion transport and canonical-coupling audit.

This is NOT the repository astronomy suite. Modules explicitly distinguish
exact ballistic transport, one conservative discrete collision model, and a
conditional canonical completion of the main-branch static field functional.
Python >=3.10; dependencies: numpy, scipy. No network or external data.
"""
from __future__ import annotations
import argparse, json, math, platform, time
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp, quad
from scipy.special import roots_legendre, exp1, erfc

HEAD = 'ee5a5f77d46371d05b5960fdf47d34cc3157ebb9'
G = 6.67430e-11
C = 299792458.0
PC = 3.085677581491367e16
YR = 365.25*86400
U_SPEED = 169443.36869
A = 6.297889390049439e-11
ELL = A*U_SPEED/2  # algebraically consistent adopted normalization


def finite(obj):
    if isinstance(obj, np.ndarray): return obj.tolist()
    if isinstance(obj, (np.integer, np.floating)): return obj.item()
    if isinstance(obj, dict): return {k:finite(v) for k,v in obj.items()}
    if isinstance(obj, (tuple,list)): return [finite(v) for v in obj]
    return obj


def moment_bridge():
    """Independent luminosity and gravity-kernel sums for prescribed hot emitters."""
    rng=np.random.default_rng(260926)
    x=rng.normal(size=(137,3))*4
    mass=rng.uniform(.05,2,size=137)
    heat=rng.uniform(0,12,size=137)
    probes=rng.normal(size=(31,3))*.6
    # Dimensionless constants chosen independently to test factors, not fit.
    ell, u, grav=0.731, 1.27, .183
    errU=[]; errF=[]
    for p in probes:
        r=p-x; d=np.linalg.norm(r,axis=1); n=r/d[:,None]
        P=ell*mass*heat
        density=np.sum(P/(4*np.pi*u*d*d))
        flux=np.sum((P/(4*np.pi*d*d))[:,None]*n,axis=0)
        S=grav*np.sum(heat*mass/d**2)
        gh=-grav*np.sum((heat*mass/d**2)[:,None]*n,axis=0)
        errU.append(abs((4*np.pi*grav*u/ell)*density-S)/S)
        errF.append(np.linalg.norm(-(4*np.pi*grav/ell)*flux-gh)/max(np.linalg.norm(gh),1e-30))
    assert max(errU)<1e-13 and max(errF)<1e-13
    return dict(emitters=137,probes=31,max_scalar_relative_error=max(errU),max_vector_relative_error=max(errF),
                relation='S_hot=4*pi*G*u*U_hot/ell; g_hot=-4*pi*G*F_hot/ell',
                scope='Stationary, transparent, incoherent hot luminosity ell*k*rho; cold coherence not replaced.')


def pair_collision(f, strength, dt):
    """Exact x+ + x- <-> y+ + y- reaction. Energy and momentum conserved.
    f axis 0 is +x,-x,+y,-y. Assumes equal packet energy at one frequency.
    strength has reciprocal density-time units; arbitrary hypothesis, not fitted.
    """
    E=f.sum(axis=0); px=f[0]-f[1]; py=f[2]-f[3]
    old=f[0]+f[1]
    eq=np.divide(E*E+px*px-py*py,2*E,out=np.zeros_like(E),where=E>0)
    p=eq+(old-eq)*np.exp(-strength*E*dt)
    out=np.stack(((p+px)/2,(p-px)/2,(E-p+py)/2,(E-p-py)/2))
    # Do not clip: conservation checks should expose any numerical problem.
    assert np.min(out)>-2e-13
    return out


def stream(f):
    out=np.empty_like(f)
    out[0]=np.roll(f[0],1,axis=0);out[1]=np.roll(f[1],-1,axis=0)
    out[2]=np.roll(f[2],1,axis=1);out[3]=np.roll(f[3],-1,axis=1)
    return out


def crossing_test(outdir, n=160, run_local=True):
    """Counterpropagating packets, exact lattice streaming, reversible elastic reaction."""
    length=20.; dx=length/n; dt=dx; steps=round(6/dt)
    z=(np.arange(n)-n/2)*dx; X,Y=np.meshgrid(z,z,indexing='ij')
    sigma=.8
    f0=np.zeros((4,n,n));f0[0]=np.exp(-((X+4)**2+Y**2)/(2*sigma**2))
    f0[1]=np.exp(-((X-4)**2+Y**2)/(2*sigma**2))
    rows=[]
    for kappa in [0.,1.,10.]:
        f=f0.copy(); E0=f.sum()*dx*dx
        p0=2*np.array([(f[0]-f[1]).sum(),(f[2]-f[3]).sum()])*dx*dx
        maxE=0.;maxP=0.
        for it in range(steps):
            f=pair_collision(f,kappa,dt/2)
            f=stream(f)
            f=pair_collision(f,kappa,dt/2)
            E=f.sum()*dx*dx
            pp=2*np.array([(f[0]-f[1]).sum(),(f[2]-f[3]).sum()])*dx*dx
            maxE=max(maxE,abs(E/E0-1));maxP=max(maxP,np.linalg.norm(pp-p0)/E0)
        frac=(f[2]+f[3]).sum()/f.sum()
        rows.append(dict(collision_coefficient=kappa,transverse_energy_fraction=frac,
                         relative_energy_drift=maxE,relative_momentum_residual=maxP,
                         final_time=steps*dt,min_intensity=float(f.min()),grid_n=n))
        np.savez_compressed(outdir/f'crossing_n{n}_k{kappa:g}.npz',x=z,energy=f.sum(axis=0),directions=f)
    # Single outgoing beam cannot self-scatter through this reaction.
    one=f0.copy();one[1]=0
    nochange=np.max(abs(pair_collision(one,1e5,100)-one))
    # Local independent ODE check on random nonequilibrium cells.
    rng=np.random.default_rng(764); localerr=0
    for _ in range(24 if run_local else 0):
        ff=rng.uniform(.01,2,4)
        def rhs(t,v):
            r=1.7*(v[0]*v[1]-v[2]*v[3]);return np.array([-r,-r,r,r])
        sol=solve_ivp(rhs,(0,.8),ff,rtol=2e-11,atol=2e-13,method='DOP853')
        pred=pair_collision(ff[:,None],1.7,.8)[:,0]
        localerr=max(localerr,np.max(abs(pred-sol.y[:,-1])))
    assert max(r['relative_energy_drift'] for r in rows)<1e-12
    assert max(r['relative_momentum_residual'] for r in rows)<1e-12
    assert localerr<1e-9 and nochange<1e-12
    return dict(rows=rows,local_ode_cases=24 if run_local else 0,max_local_ode_discrepancy=localerr,
                single_beam_change=nochange,
                scope='Four angular bins, 2D geometry, one explicitly reversible 2-to-2 reaction; not general 3D opacity.')


def hot_shell():
    """3D exact shell quadrature. Total luminosity 4*pi, speed=radius=1."""
    mu,w=roots_legendre(512)
    rows=[]
    for r in [0.,.1,.3,.5,.7,.9]:
        d2=1+r*r-2*r*mu
        # Luminosity per source solid angle=1; integral azimuth=2*pi.
        energy=.5*np.sum(w/d2)
        flux=.5*np.sum(w*(r-mu)/d2**1.5)
        radial_stress=.5*np.sum(w*(r-mu)**2/d2**2)
        analytic=1 if r==0 else np.arctanh(r)/r
        rows.append(dict(radius=r,energy_relative_to_center=energy,
                         analytic_energy=analytic,outward_flux=flux,
                         radial_angular_second_moment_fraction=radial_stress/energy))
    # Switch shell off at t=0. At position r, retarded rays remain for t<distance.
    switched=[]
    for r in [0.,.5]:
        for t in [0.,.4,.75,1.,1.25,1.6]:
            if r==0: fraction=1. if t<1 else 0.
            else:
                upper=min(1.,(1+r*r-t*t)/(2*r))
                if upper<=-1: fraction=0.
                else:
                    val=quad(lambda m: .5/(1+r*r-2*r*m),-1,upper,epsabs=1e-12)[0]
                    fraction=val/(np.arctanh(r)/r)
            switched.append(dict(radius=r,time=t,remaining_fraction=fraction))
    assert max(abs(x['outward_flux']) for x in rows)<1e-10
    assert max(abs(x['energy_relative_to_center']/x['analytic_energy']-1) for x in rows)<1e-11
    return dict(steady=rows,after_switch_off=switched,
                scope='Transparent shell; zero flux is crossing intensity, not a trapped reservoir.')


def slab_detailed_balance():
    """Exact two-direction conservative slab, u=1, x in [0,1].
    Incoming f_+(0)=L, f_-(1)=R; uniform volumetric production q, split equally.
    Collision is kappa(E/2-f). J'=q; E'=-kappa*J.
    """
    x=np.linspace(0,1,1001); rows=[]
    for left,right,q,name in [(1.,1.,0.,'external_isotropic'),(1.,0.,0.,'external_one_sided'),
                              (1.,1.,1.,'external_plus_internal_source')]:
        for k in [0.,1.,3.,10.,30.]:
            j0=(2*(left-right)-q*(1+k/2))/(2+k)
            J=j0+q*x
            E=2*left-j0-k*(j0*x+q*x*x/2)
            fp=(E+J)/2;fm=(E-J)/2
            boundary_err=max(abs(fp[0]-left),abs(fm[-1]-right))
            balance=(fp[-1]+fm[0])-(left+right+q)
            # Quadratic narrow-band carrier: p/E=2/u; medium force density=2*k*J for u=1.
            medium_impulse_rate=2*k*(j0+q/2)
            momentum_boundary=2*(left-right-fp[-1]+fm[0])
            rows.append(dict(case=name,opacity=k,mean_energy=2*left-j0-k*(j0/2+q/6),
                             center_energy=float(E[500]),max_energy=float(E.max()),
                             transmission=float(fp[-1]),reflection=float(fm[0]),
                             energy_balance_residual=balance,
                             momentum_balance_residual=momentum_boundary-medium_impulse_rate,
                             boundary_residual=boundary_err))
            assert fp.min()>-1e-12 and fm.min()>-1e-12
    # Independent discrete-ordinate steady boundary linear solves; convergence.
    errors=[]
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import spsolve
    k=10.;L=R=1.;q=1.
    j0=(2*(L-R)-q*(1+k/2))/(2+k)
    for N in [64,128,256]:
        h=1/N; mat=lil_matrix((2*N,2*N));b=np.zeros(2*N)
        # centers; upwind streaming. Vacuum-boundary consistency first order.
        for i in range(N):
            mat[i,i]=1/h+k/2;mat[i,N+i]=-k/2;b[i]=q/2
            if i>0:mat[i,i-1]=-1/h
            else:b[i]+=L/h
            mat[N+i,N+i]=1/h+k/2;mat[N+i,i]=-k/2;b[N+i]=q/2
            if i<N-1:mat[N+i,N+i+1]=-1/h
            else:b[N+i]+=R/h
        f=spsolve(mat.tocsr(),b);xx=(np.arange(N)+.5)*h
        Ean=2*L-j0-k*(j0*xx+q*xx*xx/2)
        errors.append(dict(N=N,max_energy_error=float(np.max(abs(f[:N]+f[N:]-Ean)))))
    assert max(abs(z['energy_balance_residual']) for z in rows)<1e-12
    assert max(abs(z['momentum_balance_residual']) for z in rows)<1e-12
    return dict(rows=rows,independent_grid_convergence=errors,
                scope='Stationary reciprocal elastic slab; externally lit gas cannot accumulate an isotropic bath merely by increasing opacity. Source positions matter.')


def moving_source():
    """Ballistic ray ages. Two explicit emission-velocity conventions compared.
    No force or collision; source history is prescribed. Momentum/energy requirements
    of Galilean-inherited quadratic packets are audited separately, not hidden.
    """
    rng=np.random.default_rng(528);N=24000
    age=rng.uniform(0,2.,N)
    n=rng.normal(size=(N,3));n/=np.linalg.norm(n,axis=1)[:,None]
    V=np.array([1.3,0,0]); present=np.array([0.,0,0.])
    emitted=-age[:,None]*V
    fixed=emitted+age[:,None]*n
    inherited=emitted+age[:,None]*(V+n)
    # Algebraic exact identities rather than Monte Carlo centroid accuracy.
    centered_error=np.max(abs(inherited-age[:,None]*n))
    max_front_error=np.max(abs(np.linalg.norm(inherited,axis=1)-age))
    fixed_relative_centroid=np.mean(fixed,axis=0)
    inherit_relative_centroid=np.mean(inherited,axis=0)
    # After stop, an old packet's cloud center coasts with old source velocity.
    t_after=.5
    old_birth=-age[:,None]*V
    old_after=old_birth+(age+t_after)[:,None]*(V+n)
    old_center=V*t_after
    old_error=np.max(abs(old_after-old_center-(age+t_after)[:,None]*n))
    # Exact antithetic direction quadrature for energy/momentum budget E=v^2/2,m=1.
    n2=np.concatenate([n,-n],axis=0)
    rates=[]
    for beta in [0.,.5,1.,1.3,3.,6.]:
        vv=n2+np.array([beta,0,0]);ee=.5*np.sum(vv*vv,axis=1)
        enhancement=ee.mean()/.5
        rates.append(dict(V_over_u=beta,lab_energy_over_rest=enhancement,
                          bulk_kinetic_share=beta*beta/(1+beta*beta),
                          mean_packet_momentum_x=float(vv[:,0].mean())))
    return dict(packets=N,max_comoving_identity_error=centered_error,max_causal_front_error=max_front_error,
                fixed_medium_centroid=fixed_relative_centroid,source_inherited_centroid=inherit_relative_centroid,
                stopped_source_old_cloud_center=old_center,max_old_cloud_identity_error=old_error,
                energy_cost_of_quadratic_inheritance=rates,
                radius_200Myr_kpc=U_SPEED*(200e6*YR)/(1000*PC),
                scope='Inheriting source velocity makes comoving memory kinematically exact, but lab energy and momentum must also be inherited; not a completed source ledger.')


def dwarf_and_field_coupling():
    """Dwarf: superposition survives transport; the nonlinear force response still changes."""
    mu,w=roots_legendre(800)
    rows=[]
    gi=1e-4; a=1.
    isolated=gi+np.sqrt(a*gi)
    for ratio in [0.,1.,3.,10.,100.,1000.]:
        ge=ratio*gi
        mag=np.sqrt(gi*gi+ge*ge+2*gi*ge*mu)
        radial=gi+np.sqrt(a/mag)*(gi+ge*mu)
        # Uniform external vector projection integrates to zero.
        actual=.5*np.dot(w,radial)
        approx=gi*(1+5/6*np.sqrt(a/ge)) if ge>0 else isolated
        rows.append(dict(external_over_internal=ratio,mean_internal_force=actual,
                         force_relative_to_isolated=actual/isolated,linearized=approx))
    # Explicit linear transport leaves a tagged dwarf intensity unchanged.
    rng=np.random.default_rng(18);selfI=rng.uniform(.2,1,(57,12));host=rng.uniform(.1,2,(57,12))
    # Floating cancellation avoided as a false exact assertion for high host amplitude.
    superposition=[]
    for scale in [0,1,100,10000]:
        recovered=(selfI+scale*host)-scale*host
        superposition.append(dict(host_scale=scale,relative_intensity_error=float(np.max(abs(recovered-selfI))/selfI.max())))
    # External tidal differential is distinct and nonzero even with uniform field removed.
    tide=[dict(tidal_gradient=t,differential_acceleration_at_unit_separation=t) for t in [0.,.001,.01,.1]]
    return dict(force_rows=rows,tagged_transport_superposition=superposition,tidal_controls=tide,
                scope='Transparent linear transport does not erase dwarf light; retaining that light does not remove nonlinear external-field response. Angular force diagnostic, not a solved dwarf galaxy.')


def canonical_completion():
    """Test one specified first-order kinetic completion, NOT all field theories.
    Transport normalization maps S_ex=4pi G*u/ell*(U-|F|/u).
    Main-branch static field energy is g^3/(12pi G*a)-S_ex*g/(4pi G).
    Adding this as interaction energy to a canonical free hot occupation U gives
    H_chi=(1-b)U+b*|M|, M=sum(energy*n), b=u*g/ell=2g/a.
    At the monochromatic quadratic carrier shell, M=(u/2)*P.
    Occupation variations preserve that shell; no derivative at |M|=0 is assumed.
    """
    rows=[]
    rng=np.random.default_rng(81)
    p=rng.normal(size=(100,3));p/=np.linalg.norm(p,axis=1)[:,None]
    p=np.concatenate([p,-p],axis=0)  # exactly balanced counterflow
    U=np.linalg.norm(p,axis=1).sum();P=np.linalg.norm(p.sum(axis=0))
    for ga in [0.,.1,.25,.49,.5,1.,3.,10.]:
        b=2*ga
        H=(1-b)*U+b*P
        # Add equal opposite occupations preserving total momentum exactly.
        dH=2*(1-b)
        rows.append(dict(g_over_a=ga,b=b,hot_energy_coefficient=H/U,
                         energy_change_per_added_counterpropagating_pair=dH))
    # Derived threshold and homogeneous infinitesimal-counterpair audit.
    threshold=A/2
    # Positive saturation is tested only as diagnostic, not proposed as derived repair.
    sat=[]
    from scipy.optimize import brentq
    for Q,S in [(0.,.01),(0.,1.),(.01,1.),(1.,1.),(10.,10.)]:
        g=math.sqrt(Q+S)
        gs=brentq(lambda z:z*z-Q-S*np.exp(-2*z),0,math.sqrt(Q+S)+1)
        sat.append(dict(Q=Q,S=S,original_force=g,saturated_force=gs,
                        residual_without_source_Q0=(Q==0 and gs>0)))
    joint=[]
    from scipy.optimize import minimize_scalar
    for y in [.01,.0625,.25,.5625,1.,4.,25.]:
        opt=minimize_scalar(lambda x:y+(2/3)*x**3-2*x*y,bounds=(0,2*np.sqrt(y)+1),method='bounded',options={'xatol':1e-13})
        exact=y-(4/3)*y**1.5
        joint.append(dict(density_y=y,optimal_g_over_a=opt.x,analytic_g_over_a=np.sqrt(y),total_energy=opt.fun,analytic_total_energy=exact,total_energy_over_free_energy=exact/y))
        assert abs(opt.fun-exact)<1e-10*max(1,abs(exact))
    return dict(counterstream_rows=rows,joint_field_energy_minimization=joint,threshold_m_s2=threshold,
                conditional_hamiltonian='H_chi=(1-2g/a)*U+(2g/a)*|M|; M=(u/2)*P on the central quadratic carrier shell',
                diagnostic_bounded_alternative=sat,
                scope='This fails the minimal canonical kinetic+static-energy identification. A second-order wave action, constrained medium, or interaction not equal to this energy requires a separate derivation.')


def energy_mass_budget():
    """Transparent spherical emission; benchmark ordinary energy/c^2 coupling."""
    rows=[]
    for r_kpc in [1,30,100,1000]:
        r=r_kpc*1000*PC
        rows.append(dict(radius_kpc=r_kpc,companion_energy_mass_over_baryon_mass=ELL*r/(U_SPEED*C*C),
                         crossing_years=r/U_SPEED/YR))
    age=13e9*YR
    return dict(rows=rows,all_cold_output_retained_13Gyr_fraction=ELL*age/C**2,
                scope='One illustrative standard mass-energy normalization; no independent invisible halo, no extra fitted gravity multiplier. Active stress or modified coupling cannot silently be inferred.')



def spectrum_age_join():
    """Conditional join of SGM-3's emitted power spectrum and quadratic carrier.
    z=omega*L/u; b=pi/16 as in prior barrier-spectrum test.
    u_g(z)=u*sqrt(2z) for the prior D=u*L/2 choice.
    For an age T, only z>=(r/u/T)^2/2 has arrived.
    This is not a claim that the internal gate's spectrum is already the outgoing one.
    """
    b=np.pi/16
    rows=[]
    for r_kpc in [30.,100.,1000.]:
        for age_Gyr in [1.,5.,13.]:
            r=r_kpc*1000*PC;T=age_Gyr*1e9*YR
            zmin=.5*(r/U_SPEED/T)**2
            root=zmin/b
            exact=exp1(root)/np.sqrt(2*np.pi*b)
            # Independent log-frequency quadrature.
            numerical=quad(lambda y:np.exp(-np.exp(y)/b)/np.sqrt(2*np.pi*b),
                           np.log(zmin),np.log(max(zmin+100*b,100*b)),epsabs=1e-12,epsrel=2e-11)[0]
            fraction=erfc(np.sqrt(root))
            rows.append(dict(radius_kpc=r_kpc,age_Gyr=age_Gyr,minimum_arriving_frequency=zmin,
                             density_relative_to_monochromatic_u=exact,
                             luminosity_fraction_arrived=fraction,
                             independent_integral_relative_error=abs(numerical-exact)/max(exact,1e-290)))
    high=erfc(600000/U_SPEED/np.sqrt(2*b))
    return dict(rows=rows,spectrum_b=b,mean_emitted_power_weighted_group_speed_km_s=U_SPEED*np.sqrt(2*b/np.pi)/1000,
                power_fraction_group_speed_above_600_km_s=high,
                stationary_IR_integrand='w(z)/u_g(z) proportional to 1/z; logarithmic divergence without an infrared cutoff or finite age',
                scope='Conditional combination of two separate SGM-3 choices; finite age regularizes density but changes the one-speed interpretation.')


def gas_own_emission():
    """Equal-size slab diagnostic. Cold gas still emits ell per mass.
    Background illumination is the same in both slabs and cancels in the excess.
    Galaxies: q_s=(1+k)*M_s, transparent. Gas: q_g=M_g, opacity tau.
    Source-generated integrated energy = q*(1/2+tau/12).
    """
    rows=[]
    for ratio in [10.,20.,50.]:
        for k in [30.,100.,150.]:
            critical=max(0.,6*((1+k)/ratio-1))
            rows.append(dict(gas_to_stellar_mass=ratio,stellar_heat_weight=k,
                             gas_opacity_where_source_energy_overtakes=critical))
    example=[]
    for opacity in [0.,10.,24.3,30.]:
        star=101*.5;gas=20*(.5+opacity/12)
        example.append(dict(gas_opacity=opacity,stellar_source_energy=star,gas_source_energy=gas,gas_over_stellar=gas/star))
    return dict(thresholds=rows,example=example,
                scope='Equal geometries and stipulated opacities, not measured cluster opacities or a lensing prediction; gas is not falsely treated as nonemitting.')


def recoil_microcheck():
    """Independent exact elastic encounter, quadratic carrier represented by m=1.
    Unlike infinite-mass gray scattering, this test retains finite target recoil.
    It is NOT an executed 3D transport simulation through moving gas.
    """
    rng=np.random.default_rng(9007);rows=[]
    for mass in [1.,10.,1000.,1e6]:
        n=10000
        v=rng.normal(size=(n,3));v/=np.linalg.norm(v,axis=1)[:,None]
        vg=rng.normal(scale=.05,size=(n,3))
        center=(v+mass*vg)/(1+mass)
        rel=np.linalg.norm(v-vg,axis=1)
        ns=rng.normal(size=(n,3));ns/=np.linalg.norm(ns,axis=1)[:,None]
        vp=center+mass/(1+mass)*rel[:,None]*ns
        vgp=center-1/(1+mass)*rel[:,None]*ns
        oldE=.5*np.sum(v*v,axis=1)+.5*mass*np.sum(vg*vg,axis=1)
        newE=.5*np.sum(vp*vp,axis=1)+.5*mass*np.sum(vgp*vgp,axis=1)
        oldP=v+mass*vg;newP=vp+mass*vgp
        errE=np.max(abs(newE-oldE)/oldE)
        errP=np.max(np.linalg.norm(newP-oldP,axis=1)/np.maximum(np.linalg.norm(oldP,axis=1),1e-30))
        rows.append(dict(target_to_packet_effective_mass=mass,encounters=n,max_relative_energy_error=errE,max_relative_momentum_error=errP))
        assert errE<1e-12 and errP<1e-12
    return dict(rows=rows,scope='Exact two-body kinetic-energy and momentum check; finite recoil generally changes carrier speed and frequency.')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=Path('results'))
    args=parser.parse_args();out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    start=time.time()
    tests=[('moment_bridge',moment_bridge),('crossing',lambda:crossing_test(out)),('crossing_refined',lambda:crossing_test(out,320,False)),('crossing_finest',lambda:crossing_test(out,640,False)),('shell',hot_shell),
           ('gas_countercheck',slab_detailed_balance),('moving_source',moving_source),
           ('dwarf',dwarf_and_field_coupling),('canonical_completion',canonical_completion),
           ('energy_mass_budget',energy_mass_budget),('spectrum_age_join',spectrum_age_join),
           ('gas_own_emission',gas_own_emission),('recoil_microcheck',recoil_microcheck)]
    result={'provenance':{'repo':'lrspeiser/photon-graviton','head_read':HEAD,
                         'u_m_s':U_SPEED,'a_m_s2':A,'ell_W_kg':ELL,
                         'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'interpretation':'Transport and conditional coupling tests; no astronomical refit; assertions are numerical checks, not physical acceptance.'}
    for name,fn in tests:
        t=time.time();print('RUN',name,flush=True);result[name]=fn();print('OK',name,round(time.time()-t,3),flush=True)
    result['wall_seconds']=time.time()-start
    (out/'results.json').write_text(json.dumps(finite(result),indent=2)+'\n')
    print('WROTE',out/'results.json',flush=True)

if __name__=='__main__':main()
