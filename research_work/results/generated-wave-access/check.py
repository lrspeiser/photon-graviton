"""Generate companion waves with the photon Hamiltonian and measure slab access."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = Path(os.environ.get('PHOTON_GRAVITON_RESULTS', ROOT/'research_work/generated'))/'generated-wave-access'


def simulate(cells, sigma, mass, protocol):
    cfg = protocol['parameters']
    xmin, xmax = cfg['domain']
    dx = (xmax-xmin)/cells
    # Cell centers put slab and monitor boundaries on faces at both resolutions.
    x = xmin + (np.arange(cells)+0.5)*dx
    K, v = cfg['field_inertia'], cfg['field_speed']
    z = (x-cfg['conversion_center'])/cfg['conversion_half_width']
    g = np.where(abs(z)<1, np.cos(np.pi*z/2)**2, 0.0)
    slab = (x>cfg['slab_edges'][0]) & (x<cfg['slab_edges'][1])
    m2 = mass**2*slab
    left, right = [int(round((edge-xmin)/dx))-1 for edge in cfg['monitor_edges']]
    region = slice(left+1, right+1)
    energy0 = cfg['photon_initial_energy']
    # Additional states: incident, reflected, transmitted, incoming from far left,
    # source work in monitored region, integrated net force from fixed profiles/grid.
    initial = np.zeros(2*cells+8)
    initial[2*cells:2*cells+2] = [cfg['photon_initial_x'], energy0]

    def quantities(y):
        phi, u = y[:cells], y[cells:2*cells]
        X, P = y[2*cells:2*cells+2]
        W = np.exp(-0.5*((x-X)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        weight = W*g
        n = 1 + dx*np.dot(weight, phi)
        n_X = dx*np.dot(weight*(x-X)/sigma**2, phi)
        source = (P/n**2)*weight
        gradient = (np.roll(phi,-1)-phi)/dx
        average_u = (u+np.roll(u,-1))/2
        # In a massless region these are directional wave-energy flux magnitudes.
        # Their difference is exactly the discrete energy flux at a cell face.
        go_left = K*v/4*(average_u+v*gradient)**2
        go_right = K*v/4*(average_u-v*gradient)**2
        return phi, u, P, n, n_X, source, gradient, go_left, go_right

    def rhs(t, y):
        phi,u,P,n,n_X,source,gradient,L,R = quantities(y)
        acceleration = v*v*(np.roll(phi,-1)-2*phi+np.roll(phi,1))/dx**2-m2*phi+source/K
        central_gradient = (np.roll(phi,-1)-np.roll(phi,1))/(2*dx)
        photon_force = P*n_X/n**2
        # Exact rate of the chosen discrete total momentum. The static profiles
        # and spatial grid can exchange momentum; their energy source is zero.
        external_force = photon_force-K*dx*np.dot(acceleration,central_gradient)
        rates = [L[right],R[right],L[left],R[left],dx*np.dot(source[region],u[region]),external_force]
        return np.concatenate([u, acceleration, [1/n,photon_force], rates])

    times = np.linspace(0,cfg['end_time'],cfg['samples'])
    sol = solve_ivp(rhs, [0,cfg['end_time']], initial, t_eval=times,
                    method='DOP853',rtol=2e-9,atol=2e-12,max_step=0.75*dx/v)
    assert sol.success, sol.message
    total, stored, optical, mom_errors, boundary_energies, incident_amplitude = [], [], [], [], [], []
    for y in sol.y.T:
        phi,u,P,n,n_X,source,gradient,L,R = quantities(y)
        density = K/2*(u*u+m2*phi*phi) + K*v*v/4*(gradient**2+np.roll(gradient,1)**2)
        total.append(dx*density.sum()+P/n)
        stored.append(dx*density[region].sum())
        optical.append(n)
        momentum = P-K*dx*np.dot(u,(np.roll(phi,-1)-np.roll(phi,1))/(2*dx))
        mom_errors.append(abs(momentum-energy0-y[-1]))
        boundary_energies.append(dx*(density[:8].sum()+density[-8:].sum()))
        incident_amplitude.append((u[right]+u[right+1])/2+v*gradient[right])
    final=sol.y[:,-1]
    incident, reflected, transmitted, left_incoming, work, force_impulse = final[-6:]
    residual=np.asarray(stored)-(sol.y[-6]-sol.y[-5]-sol.y[-4]+sol.y[-3]+sol.y[-2])
    fractions = {'reflected':float(reflected/incident),'transmitted':float(transmitted/incident),
                 'temporarily_inside':float(stored[-1]/incident),'incoming_from_left':float(left_incoming/incident)}
    checks = protocol['checks']
    energy_drift=max(abs(np.asarray(total)-energy0))/energy0
    balance=max(abs(residual))/incident
    assert energy_drift<checks['relative_total_energy_tolerance'],energy_drift
    assert balance<checks['local_balance_over_incident_tolerance'],balance
    assert max(mom_errors)/energy0<checks['relative_total_energy_tolerance']
    assert max(boundary_energies)/energy0<1e-12, 'Periodic boundary contamination'
    assert min(optical)>checks['minimum_effective_optical_factor']
    assert abs(work)/incident<checks['maximum_direct_source_work_over_incident']
    samples=[]
    for target in [15,25,35,50]:
        j=int(np.argmin(abs(times-target)))
        samples.append({'time':float(times[j]),'incident_energy':float(sol.y[-6,j]),
                        'reflected_energy':float(sol.y[-5,j]),'transmitted_energy':float(sol.y[-4,j]),
                        'energy_inside_monitors':float(stored[j])})
    row = {'cells':cells,'dx':dx,'smoothing_width':sigma,'restoring_frequency':mass,
            'slab_static_center_suppression':float(1/np.cosh(mass/v)),
            'photon_initial_energy':energy0,'photon_final_energy':float(final[2*cells+1]/optical[-1]),
            'photon_energy_lost':float(energy0-final[2*cells+1]/optical[-1]),
            'incident_energy':float(incident),'reflected_energy':float(reflected),'transmitted_energy':float(transmitted),
            'energy_inside_monitors_at_end':float(stored[-1]),'fractions_of_incident':fractions,
            'maximum_relative_total_energy_drift':float(energy_drift),'maximum_local_balance_error_over_incident':float(balance),
            'maximum_momentum_balance_error_over_initial_photon_energy':float(max(mom_errors)/energy0),
            'maximum_boundary_strip_energy_over_initial_photon_energy':float(max(boundary_energies)/energy0),
            'fixed_profile_and_grid_impulse':float(force_impulse),'direct_source_work_inside_monitors':float(work),
            'minimum_effective_optical_factor':float(min(optical)),'final_photon_x':float(final[2*cells]),
            'time_samples':samples,'rhs_evaluations':sol.nfev}
    return row, np.asarray(incident_amplitude)


def spectral_prediction(amplitude, dt, mass, v, half_width):
    # Zero padding refines quadrature of the finite pulse's Fourier transform;
    # it is not extra temporal data. Keep the DC component (do not demean).
    count=32768
    omega=2*np.pi*np.fft.rfftfreq(count,dt)
    power=abs(np.fft.rfft(amplitude,n=count))**2
    power[1:-1]*=2
    weight=power/power.sum()
    if mass==0:
        transmission=np.ones_like(omega)
    else:
        B=mass*half_width/v; w=omega/mass
        transmission=np.zeros_like(w)
        low=(w>0)&(w<1); high=w>1; threshold=w==1
        transmission[low]=1/(1+np.sinh(2*B*np.sqrt(1-w[low]**2))**2/(4*w[low]**2*(1-w[low]**2)))
        transmission[high]=1/(1+np.sin(2*B*np.sqrt(w[high]**2-1))**2/(4*w[high]**2*(w[high]**2-1)))
        transmission[threshold]=1/(1+B*B)
    return {'predicted_transmitted_fraction':float(np.dot(weight,transmission)),
            'incident_energy_fraction_below_cutoff':float(weight[omega<mass].sum()),
            'median_incident_angular_frequency':float(omega[np.searchsorted(np.cumsum(weight),.5)]),
            'frequency_below_which_99_percent_of_energy_lies':float(omega[np.searchsorted(np.cumsum(weight),.99)]),
            'frequency_quadrature_spacing':float(omega[1])}


def main():
    protocol=json.loads((HERE/'protocol.json').read_text())
    grid=protocol['grid']; runs=[]; source_signals={}
    for sigma in grid['packet_smoothing_width']:
        for mass in grid['restoring_frequency']:
            for cells in grid['cells']:
                row,signal=simulate(cells,sigma,mass,protocol);runs.append(row)
                if mass==0: source_signals[(cells,sigma)]=signal
                print(f"N={cells} sigma={sigma} m={mass}: {row['fractions_of_incident']}",flush=True)
    convergence=[]
    for sigma in grid['packet_smoothing_width']:
        for mass in grid['restoring_frequency']:
            a,b=[r for r in runs if r['smoothing_width']==sigma and r['restoring_frequency']==mass]
            delta=max(abs(a['fractions_of_incident'][key]-b['fractions_of_incident'][key]) for key in a['fractions_of_incident'])
            incident_change=abs(a['incident_energy']/b['incident_energy']-1)
            convergence.append({'smoothing_width':sigma,'restoring_frequency':mass,
                                'maximum_absolute_fraction_change':delta,'relative_incident_energy_change':incident_change})
            assert delta<protocol['checks']['max_absolute_fraction_change_on_refinement'],convergence[-1]
            assert incident_change<protocol['checks']['max_relative_incident_energy_change_on_refinement'],convergence[-1]
    spectral_checks=[]
    cfg=protocol['parameters']
    for row in runs:
        spectrum=spectral_prediction(source_signals[(row['cells'],row['smoothing_width'])],
                                     cfg['end_time']/(cfg['samples']-1),row['restoring_frequency'],
                                     cfg['field_speed'],(cfg['slab_edges'][1]-cfg['slab_edges'][0])/2)
        actual=row['fractions_of_incident']['transmitted']
        spectrum['relative_transmission_difference']=abs(actual/spectrum['predicted_transmitted_fraction']-1)
        row['independent_source_spectrum_comparison']=spectrum
        spectral_checks.append(spectrum['relative_transmission_difference'])
    spectral_intervals=[]
    for sigma in grid['packet_smoothing_width']:
        for mass in grid['restoring_frequency']:
            a,b=[r for r in runs if r['smoothing_width']==sigma and r['restoring_frequency']==mass]
            transmission=b['fractions_of_incident']['transmitted']
            remaining=b['fractions_of_incident']['temporarily_inside']
            mesh_change=abs(a['fractions_of_incident']['transmitted']-transmission)
            predicted=b['independent_source_spectrum_comparison']['predicted_transmitted_fraction']
            # The frequency-domain prediction includes the entire outgoing tail;
            # at finite time some energy is still inside. Mesh change is an
            # empirical error estimate, not a mathematically certified bound.
            lower=transmission-mesh_change-1e-12
            upper=transmission+remaining+mesh_change+1e-12
            assert lower<=predicted<=upper,(sigma,mass,lower,predicted,upper)
            spectral_intervals.append({'smoothing_width':sigma,'restoring_frequency':mass,
                                       'predicted_eventual_transmission':predicted,
                                       'finite_window_plus_mesh_comparison_interval':[lower,upper],
                                       'agreement_with_finite_time_and_mesh_allowance':True})
    result={'scope':protocol['scope'],'checks_pass':True,'runs':runs,'grid_convergence':convergence,
            'source_spectrum_checks':spectral_intervals,
            'absorption_or_permanent_deposit_channel':False,'density_profile_and_gravity_dynamically_derived':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'check.py']}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'generated-wave-access-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'checks_pass':True,'grid_convergence':convergence,'max_spectral_transmission_difference':max(spectral_checks)},indent=2))


if __name__=='__main__':
    main()
