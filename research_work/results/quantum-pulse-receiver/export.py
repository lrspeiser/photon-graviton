"""Energy/momentum-preserving receiver-to-null-mode isometry; no emission rate."""
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
unit=.05
n=np.arange(96,105)
energy=2*n*unit
U=50.
er=U+unit*np.arange(201)
r=np.exp(-((er-55)/.5)**2/4)*np.exp(1j*er)
r/=np.linalg.norm(r)
a=np.exp(-((energy-10)/.22)**2/4)
a/=np.linalg.norm(a)
size=len(er)+2*int(n.max())
er_final=U+unit*np.arange(size)
companion_energy=er_final-U
joint=np.zeros((len(n),len(n),size),complex)
for i,shift1 in enumerate(n):
    for j,shift2 in enumerate(n):
        shift=shift1+shift2
        joint[i,j,shift:shift+len(r)]=a[i]*a[j]*np.exp(4j*energy[j])*r
flat=joint.reshape(len(n)**2,size)
rho_before=flat@flat.conj().T
# Isometry relabels r -> |ground U> |massless energy r-U>.
# Ground factor is fixed. All old distinct receiver labels remain orthogonal.
exported=flat.copy()
rho_after=exported@exported.conj().T
export_mean=float(np.sum(abs(exported)**2*companion_energy[None,:]))
seed=float(abs(r)**2@(er-U))
photon_loss=float(abs(a)**2@energy)  # two photons each lose half
energy_error=float(np.max(abs(er_final-(U+companion_energy))))
momentum_error=float(np.max(abs((er_final-U)-companion_energy)))

# If full export occurs after photon 1, receiver 2 starts at sharp energy U.
# Different photon-2 losses lead to orthogonal exported companion energies.
second_reset_rho=np.diag(abs(a)**2).astype(complex)
coherent_joint=np.zeros((len(n),len(er)+int(n.max())),complex)
for i,shift in enumerate(n):
    coherent_joint[i,shift:shift+len(r)]=a[i]*r
second_fresh_coherent_rho=coherent_joint@coherent_joint.conj().T
def offdiagonal_l1(rho):
    return float(np.sum(abs(rho-np.diag(np.diag(rho)))))
assert energy_error<1e-12 and momentum_error<1e-12
assert np.max(abs(rho_before-rho_after))<1e-12
assert abs(export_mean-seed-photon_loss)<1e-12
result=dict(scope='Asymptotic export map; local interaction and emission lifetime absent',
    receiver_ground_energy=U,initial_receiver_excess_energy=seed,
    photon_supplied_energy=photon_loss,total_exported_energy=export_mean,
    basis_energy_error=energy_error,basis_momentum_error=momentum_error,
    photon_density_matrix_change_after_unconditional_export=float(np.max(abs(rho_before-rho_after))),
    second_photon_offdiagonal_l1_with_coherent_receiver=offdiagonal_l1(second_fresh_coherent_rho),
    second_photon_offdiagonal_l1_after_ground_reset=offdiagonal_l1(second_reset_rho),
    ground_reset_restores_initial_coherence=False,
    null_companion_identity_is_graviton=False,
    rate_distance_timing_or_galaxy_data_fitted=False)
(HERE/'export-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
