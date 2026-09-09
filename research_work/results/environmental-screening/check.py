"""Independent restoring-field profile and wave-interface checks."""
from pathlib import Path
import json,hashlib,os
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.constants import c,hbar,e,parsec

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'environmental-screening'


def transmission(B,w):
    if w==1:return 1/(1+B*B)
    if w<1:return 1/(1+np.sinh(2*B*np.sqrt(1-w*w))**2/(4*w*w*(1-w*w)))
    return 1/(1+np.sin(2*B*np.sqrt(w*w-1))**2/(4*w*w*(w*w-1)))


def boundary_match(B,w):
    # Set m=v=1, thickness d=2B. Unknowns r, A, B_inside, t.
    d=2*B;k=w
    if w==1:
        mat=np.array([[1,-1,0,0],[-1j*k,0,-1,0],[0,1,d,-1],[0,0,1,-1j*k]],dtype=complex)
    else:
        q=np.sqrt(complex(w*w-1));ep=np.exp(1j*q*d);em=np.exp(-1j*q*d)
        mat=np.array([[1,-1,-1,0],[-1j*k,-1j*q,1j*q,0],
                      [0,ep,em,-1],[0,1j*q*ep,-1j*q*em,-1j*k]],dtype=complex)
    rhs=np.array([-1,-1j*k,0,0],dtype=complex)
    solution=np.linalg.solve(mat,rhs)
    return float(abs(solution[0])**2),float(abs(solution[3])**2)


def main():
    protocol=json.loads((HERE/'protocol.json').read_text());checks=protocol['independent_checks'];waves=[];profiles=[];energy=[]
    for B in checks['dimensionless_barrier_strengths']:
        for w in checks['frequency_over_mass']:
            reflected,transmitted=boundary_match(B,w);exact=transmission(B,w)
            assert abs(transmitted/exact-1)<checks['relative_transmission_tolerance']
            assert abs(reflected+transmitted-1)<checks['flux_balance_tolerance']
            waves.append({'strength':B,'frequency_over_mass':w,'transmission':transmitted,'reflection':reflected,'relative_formula_error':abs(transmitted/exact-1),'flux_balance_error':abs(reflected+transmitted-1)})
    for B in checks['static_ode_strengths']:
        sol=solve_ivp(lambda x,y:[y[1],B*B*y[0]],[0,1],[1.,0.],method='DOP853',rtol=1e-12,atol=1e-13,max_step=.02)
        assert sol.success
        suppression=1/sol.y[0,-1];exact=1/np.cosh(B)
        assert abs(suppression/exact-1)<checks['relative_static_ode_tolerance']
        profiles.append({'strength':B,'numerical_center_suppression':float(suppression),'relative_error':float(abs(suppression/exact-1))})
        # a=v=K=1, m=B; A linear in time with chosen value and slope.
        A=.7;Adot=.03
        f=lambda x:np.cosh(B*x)/np.cosh(B)
        fx=lambda x:B*np.sinh(B*x)/np.cosh(B)
        growth=A*Adot*quad(lambda x:fx(x)**2+B*B*f(x)**2,-1,1,epsabs=1e-12)[0]
        inward=2*A*Adot*B*np.tanh(B)
        assert abs(growth-inward)<1e-10
        energy.append({'strength':B,'interior_energy_growth':float(growth),'incoming_boundary_power':float(inward)})
    illustrations=[];cfg=protocol['illustrations'];a=cfg['half_width_kpc']*1000*parsec;v=cfg['wave_speed_over_c']*c
    for f in cfg['center_suppressions']:
        B=float(np.arccosh(1/f));mass_frequency=B*v/a
        illustrations.append({'center_suppression':f,'required_strength':B,'illustrative_mass_angular_frequency_per_second':mass_frequency,
                              'illustrative_cutoff_energy_eV':hbar*mass_frequency/e,
                              'atomic_frequency_ratio_for_unit_boundary_field':(1+f)**-2,
                              'transmission_by_frequency_ratio':{str(w):float(transmission(B,w)) for w in cfg['wave_frequency_over_mass']}})
    source=ROOT/'research_work/results/conversion-first/results.json';prior=json.loads(source.read_text())
    old_rate=prior['fitted']['c_alpha_km_s_per_mpc']*1000/(1e6*parsec)
    result={'scope':protocol['scope'],'checks_pass':True,'matched_wave_cases':waves,'profile_cases':profiles,'boundary_energy_checks':energy,
            'illustrations':illustrations,'historical_effective_rate_per_second_for_scale_comparison':old_rate,
            'high_frequency_transmission_lower_bound_at_omega_over_m_2':48/49,
            'absorption_in_this_model':0,'density_coupling_or_actual_clock_bound_derived':False,'capture_or_gravity_derived':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'check.py',source]}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'environmental-screening-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'wave_cases':len(waves),'profile_cases':len(profiles),'max_transmission_error':max(x['relative_formula_error'] for x in waves),'illustrations':illustrations,'historical_rate':old_rate},indent=2))


if __name__=='__main__':main()
