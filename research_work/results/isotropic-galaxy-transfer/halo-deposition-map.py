"""Translate fitted NFW configurations into finite-aperture deposition targets."""
import json,hashlib
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.special import gamma,gammainc
P=Path(__file__).resolve().parent
def read(n):return json.loads((P/n).read_text())
halos=[r for r in read('nfw-geometry-results.json')['rows'] if r['model']=='NFW_exact_lens']
cp=read('third-radiation-retention-results.json')['models']['attenuated'];assert cp['q']==1/3
companions=read('third-retention-optics-results.json')['rows'];mu,w=leggauss(96)
rows=[]
for dep in companions:
    a=dep['capture_scale_kpc'];eta=dep['retention_mapping']['eta'];C0=2*eta*cp['C_Msun_kpc3'];k=cp['k0_per_kpc']
    x=np.geomspace(1e-8,1e7,24001);r=a*x;t=x[:,None]*mu;B2=1+x[:,None]**2*(1-mu**2);B=np.sqrt(B2)
    tau=k*a*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3));J=.5*np.sum(np.exp(-np.maximum(tau,0))*w,axis=1)
    rho=C0*J/(1+x*x)**2;mass=4*np.pi*(rho[0]*r[0]**3/3+cumulative_trapezoid(r*r*rho,r,initial=0));M=PchipInterpolator(np.log(r),mass)
    T=np.pi*k*a/2;sigma=np.pi*a*a*(T**(2/3)*gamma(1/3)*gammainc(1/3,T)-1+np.exp(-T));total=C0/k*sigma
    assert abs(mass[-1]/total-1)<2e-4
    for h in halos:
        if h['Name']!=dep['Name']:continue
        A=h['halo_A_Msun'];rs=h['rs_kpc'];Re=rs/h['rs_over_Re'];rE=dep['lens_catalog_arcsec']/206264.80624709636*h['Dl_Mpc']*1000
        density=A/(4*np.pi*rs**3)
        def hm(q):
            u=q/rs;return A*(np.log1p(u)-u/(1+u))
        samples=[]
        for label,R in [('0.5Re',.5*Re),('Re',Re),('2Re',2*Re),('5Re',5*Re),('10kpc',10.),('Einstein_sphere',rE)]:
            target=hm(R);current=float(M(np.log(R)))
            samples.append(dict(aperture=label,radius_kpc=R,NFW_enclosed_Msun=target,original_deposit_enclosed_Msun=current,NFW_to_original_enclosed=target/current,required_energy_J=target*1.98847e30*299792458**2,total_inventory_to_target=total/target))
        # Minimum mass moved to reproduce density inside a finite diagnostic aperture.
        # This permits arbitrary redistribution outside and does not specify dynamics.
        R=5*Re;rr=np.geomspace(max(1e-8*a,1e-9*R),R,12001)
        targetrho=A/(4*np.pi*rs**3)/(rr/rs*(1+rr/rs)**2)
        currentrho=np.exp(PchipInterpolator(np.log(r),np.log(rho))(np.log(rr)))
        difference=4*np.pi*rr**2*(targetrho-currentrho)
        deficit=float(np.trapezoid(np.maximum(difference,0),rr));excess=float(np.trapezoid(np.maximum(-difference,0),rr))
        row=dict(Name=h['Name'],geometry=h['geometry'],population=dep['retention_mapping']['population'],rs_kpc=rs,rho_s_Msun_kpc3=density,A_Msun=A,stellar_mass_Msun=h['stellar_mass_Msun'],beta=h['beta'],beta_boundary=h['beta_boundary'],scale_boundary=h['scale_boundary'],fraction_boundary=h['fraction_boundary'],stellar_chi2=h['stellar_chi2'],Re_kpc=Re,original_capture_scale_kpc=a,eta=eta,original_deposit_total_Msun=total,total_integral_relative_error=abs(mass[-1]/total-1),samples=samples,finite_target_radius_kpc=R,finite_target_mass_Msun=hm(R),inventory_sufficient_inside_5Re=bool(total>=hm(R)),minimum_moved_mass_inside_5Re_Msun=max(deficit,excess),minimum_moved_fraction_of_original_total=max(deficit,excess)/total,added_mass_inside_5Re_Msun=deficit,removed_mass_inside_5Re_Msun=excess)
        rows.append(row)
out=dict(scope='Required finite-aperture halo profiles; original fixed companion inventory, not measured photon supply',standard_comparison='Standard-geometry target is compared with the original physical companion profile without changing its inventory or scale. This is a hypothetical redistribution target, not a recomputed standard-geometry companion prediction.',rows=rows,input_sha256={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['nfw-geometry-results.json','third-retention-optics-results.json','third-radiation-retention-results.json']})
(P/'halo-deposition-map-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
lines=['# Fitted halo configurations and companion deposition targets','',
       'These are the exact numerical configurations of our restricted best fits, not uniquely measured halos. Some have poor stellar fits, extreme scale radii or parameter boundaries. They are not the published best-fit configurations from a full lens-image analysis.','',
       '## Halo definition','',
       'Known NFW density: rho_h(r)=rho_s/[x(1+x)^2], x=r/r_s. Known enclosed mass: M_h(<r)=A[ln(1+x)-x/(1+x)], A=4*pi*rho_s*r_s^3. Both r_s and rho_s specify the entire spherical fitting profile. r_s is its slope-transition radius, not the halo edge. The untruncated NFW total mass diverges logarithmically; A is a normalization, not total mass. A physical finite halo requires an outer truncation not measured here.','',
       '## Configurations under the standard benchmark geometry','',
       '| Galaxy | r_s kpc | rho_s Msun/kpc^3 | Halo mass inside 10 kpc, Msun | Stellar mass, Msun | Limits |','|---|---:|---:|---:|---:|---|']
for q in rows:
    if q['population']!='Chabrier' or q['geometry']!='standard_flat_FLRW':continue
    m=next(s for s in q['samples'] if s['aperture']=='10kpc');flags=', '.join(k for k in ['beta','scale','fraction'] if q[k+'_boundary']) or 'no bound hit'
    lines.append(f"| {q['Name']} | {q['rs_kpc']:.5g} | {q['rho_s_Msun_kpc3']:.5g} | {m['NFW_enclosed_Msun']:.5g} | {q['stellar_mass_Msun']:.5g} | {flags} |")
lines += ['', 'The table uses the standard geometry because it gave better NFW fits in the preceding comparison. This geometry is an alternative benchmark only. Our own geometry has separately fitted halo parameters in the JSON and should not be mixed with these numbers. J1402 still fits poorly; J1621 has a scale at its upper bound; J1204/J1402 have very compact fitted scale radii. No confidence intervals or realistic halo-population priors were obtained.','',
          '## Could our fixed deposit inventory occupy those regions?','',
          'Under our assumed ordinary gravitational coupling, rho_dep=u_dep/c^2. A deposited component matching the same density and relevant stresses would have the same gravitational effect regardless of its origin. That is a conditional equivalence, not evidence that capture produces that distribution.','',
          '| Galaxy | Original total deposit inventory, Msun | Target halo mass inside 5 Re, Msun | Inventory / target | Minimum inventory fraction moved to reproduce density inside 5 Re |','|---|---:|---:|---:|---:|']
for q in rows:
    if q['population']!='Chabrier' or q['geometry']!='standard_flat_FLRW':continue
    move=f"{100*q['minimum_moved_fraction_of_original_total']:.3f}%" if q['inventory_sufficient_inside_5Re'] else 'Insufficient fixed inventory'
    lines.append(f"| {q['Name']} | {q['original_deposit_total_Msun']:.5g} | {q['finite_target_mass_Msun']:.5g} | {q['original_deposit_total_Msun']/q['finite_target_mass_Msun']:.3f} | {move} |")
lines += ['', 'This is a finite-aperture bookkeeping test, not a formation result. Five effective radii is a declared diagnostic region, not a measured halo boundary; portions extend beyond stellar constraints. The original profile is evaluated at the target physical radii without changing its existing inventory or scale. That comparison does not recalculate the one-third law under standard geometry. Both population cases and both halo geometries are retained, with masses at 0.5, 1, 2, 5 Re, 10 kpc, and the sphere at the Einstein radius. A sphere is not the projected cylinder measured by lensing.','',
          'Minimum moved mass is max(integrated positive density deficit, integrated positive density excess) inside 5 Re, assuming unrestricted rearrangement outside, a nonnegative reservoir, and enough total inventory. It is not an efficiency, migration probability or measured amount of movement. It applies only when the reported inventory suffices. Matching this finite radial interval does not reproduce the full line-of-sight lensing signal; mass beyond it and stellar mass must also be consistent.','',
          '## What the one-third law would need to produce','',
          'Our fixed law remains eta=X^(1/3)/(1+X^(1/3)), with rho0(r)=2*C_half*eta*J(r)/[1+(r/a)^2]^2. The Hill function and attenuation mathematics are known; their companion interpretation is proposed. The required redistribution is Delta_rho(r)=rho_h(r)-rho0(r). Positive Delta_rho means deposit more locally; negative means relocate some elsewhere. Cumulative target/original ratios in the JSON locate the mismatch without changing eta.','',
          'Equivalently a proposed effective capture-history model would need u_dep(r)=c^2*rho_h(r). If its local retained input power is q_ret(r), lossless accumulation requires integral q_ret(r,t)dt=c^2*rho_h(r), with transport and any work included. We have not measured that input or history. The fitted C_half inventory is not an independently established photon budget; passing its inventory test cannot establish sufficient real photon energy.','',
          'For a spherically symmetric redistribution flux j_E, energy continuity gives partial_t u+1/r^2 partial_r(r^2 j_E)=sources-sinks. Integrated over redistribution with no net local sources or sinks, 4*pi*r^2*integral j_E dt=-c^2*Delta_M(<r). This is known continuity mathematics: a positive required enclosed-mass change demands net inward energy transport. It does not supply a force law or a timescale for that transport.','',
          'Deposits also require a viable support/stress prescription. Placing positive energy at the fitted radii is not enough to ensure it stays there or produces the assumed lensing response. The model still lacks a demonstrated capture-and-migration law yielding these NFW-like targets. Merely naming the halo companion energy would relabel the inferred component without explaining it.','',
          '## Verification','',
          'The angular-capture density is integrated directly and compared with its analytic total capture-area expression, with relative tolerance 2e-4. Known spherical mass and E=Mc^2 conversions are used. SHA-256 hashes identify all three input results. No fit was rerun or silently selected on this inventory calculation. The prior geometry and stellar fit limitations continue to apply. All six goals remain open.']
(P/'halo-deposition-map-report.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines[7:31]))
