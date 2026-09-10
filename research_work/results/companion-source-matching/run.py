"""Match the frozen lens-source mass to production and cold-state kinematics.

Known relativistic identities; conditional same-species assumptions, no fit.
"""
from pathlib import Path
import hashlib
import json
import mpmath as mp

mp.mp.dps=110
HERE=Path(__file__).resolve().parent
CAL=HERE.parent/'wave-lens-training/parameter-grid/selected-calibration.json'
cal=json.loads(CAL.read_text())
mass=mp.mpf(str(cal['field_mass_eV_c2'])) # numerical value of rest energy in eV
c=mp.mpf('299792458')
h=mp.mpf('6.62607015e-34')
electronvolt=mp.mpf('1.602176634e-19')
energies=[('1 GHz photon',h*mp.mpf('1e9')/electronvolt),
          ('1 meV photon',mp.mpf('.001')),('1 eV photon',mp.mpf(1)),
          ('1 keV photon',mp.mpf(1000)),('1 MeV photon',mp.mpf('1e6'))]

def value(x):
    return float(x)

production=[]; residuals=[]
for name,energy in energies:
    angle=2*mp.asin(mass/(2*energy))
    boost=2*energy/mass
    # Stable near-light-speed subtraction; binary64 would simply return zero.
    speed_deficit=1/(boost**2*(1+mp.sqrt(1-1/boost**2)))
    invariant=2*energy**2*(1-mp.cos(angle))
    residuals.append(abs(invariant/mass**2-1))
    assert residuals[-1]<mp.mpf('1e-45')
    production.append(dict(input=name,each_photon_energy_eV=value(energy),
        required_full_opening_angle_rad=value(angle),
        produced_scalar_energy_eV=value(2*energy),
        produced_scalar_lorentz_factor=value(boost),
        one_minus_scalar_speed_over_c=value(speed_deficit),
        classification='Vacuum two-photon to one on-shell scalar; kinematic condition, not rate'))

storage=[]
for speed_kms in [100,200,300]:
    beta=mp.mpf(speed_kms)*1000/c
    gamma=1/mp.sqrt(1-beta**2)
    cold_energy=gamma*mass
    for label,energy in energies:
        retained=cold_energy/energy
        count=energy/cold_energy
        assert abs(retained*count-1)<mp.mpf('1e-100')
        storage.append(dict(incoming_companion_energy_eV=value(energy),
            energy_scale_label=label.replace('photon','energy scale'),
            illustrative_final_speed_kms=speed_kms,
            cold_energy_per_quantum_eV=value(cold_energy),
            kinetic_energy_over_rest_energy=value(gamma-1),
            energy_fraction_retained_in_same_one_cold_quantum=value(retained),
            cold_quanta_for_all_incoming_energy=value(count)))

# Alternative cold production: two opposing photons with energy mc^2/2 each.
frequency=mass*electronvolt/(2*h)
result=dict(classification='Conditional same-field consistency check; no observational scores or photon-supply budget',
    calibration_path=str(CAL.relative_to(HERE.parents[2])),
    calibration_sha256=hashlib.sha256(CAL.read_bytes()).hexdigest(),
    frozen_mass_eV_c2=value(mass),production=production,storage=storage,
    head_on_cold_production=dict(each_photon_energy_eV=value(mass/2),
        each_photon_frequency_hz=value(frequency),
        wavelength_m=value(c/frequency)),
    verification=dict(arithmetic_decimal_precision=mp.mp.dps,
        max_independent_invariant_relative_error=float(max(residuals)),
        all_energy_partition_reciprocity_checks_pass=True),
    assumptions=['Locally vacuum photons and usual special-relativistic dispersion',
                 'One on-shell scalar in two-photon production diagnostic',
                 'Same scalar rest mass as frozen nonrelativistic lens-source fit',
                 'Cold-state counting is an occupation/particle diagnostic, not proof of particle-number conservation'])
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))
