"""Energy and momentum ledger for receiver-assisted threshold production (RB-1).

Per galaxy and for the fixed retained inventory M_ret. Energies are in units of
M_ret c^2 unless stated. Nothing here is fitted. The heavy-receiver limit is
used throughout; receiver recoil is bounded by q/2 of the absorbed energy with
q = E'/(M_R c^2), negligible for any eV-scale quantum on nucleon receivers.
"""
import numpy as np
from kinematics import C_KMS, REFERENCE_CHI, bound_energy_fraction, drag_coefficient

GYR = 0.9777922216807891            # Gyr per kpc/(km/s)


def ledger(pot, r_sites, receiver_mass, J, zeta, disk_fraction, M_ret, push_limit=.1):
    """Energy, momentum and angular-momentum bookkeeping.

    receiver_mass: baryonic receiver mass per site (Msun); J, zeta: incident
    mean intensity and net outward flux fraction at the sites; disk_fraction:
    rotating-disk share of the receivers at each site (the rest carry no net
    rotation). Absorption is proportional to receiver_mass*J.
    """
    r, m = np.asarray(r_sites, float), np.asarray(receiver_mass, float)
    W = m*J
    ve, vc = pot.vesc(r), pot.vcirc(r)
    fE = bound_energy_fraction(ve)                        # bound rest energy / absorbed energy
    efficiency = float(np.sum(W*fE)/np.sum(W))
    multiplier = 1/efficiency
    bound = W*ve**3
    bound /= bound.sum()
    kinetic = float(np.sum(bound*.3*(ve/C_KMS)**2))       # <v^2> = 3/5 v_esc^2 in the ball
    orbital = float(np.sum(bound*.4*pot.phi(r))/C_KMS**2)  # <E> = (2/5) Phi(r0)
    kappa = drag_coefficient(REFERENCE_CHI)
    Jdisk = float(np.sum(m*disk_fraction*r*vc))           # receivers' net angular momentum
    lever = float(np.sum(W*disk_fraction*r*vc)/np.sum(W))  # absorption-weighted r v_phi
    loss_ideal = 4/3*M_ret*lever/Jdisk                    # gray drag, every quantum retained
    loss_declared = kappa*multiplier*M_ret*lever/Jdisk
    Krec = float(np.sum(m*vc*vc)/2)
    drain_declared = kappa*multiplier*M_ret*float(np.sum(W*vc*vc)/np.sum(W))/Krec
    g = vc*vc/r
    # a_push = |zeta| (absorbed power per receiver mass)/c; power per mass = E_abs J/(T sum W).
    T_min = push_limit**-1*np.abs(zeta)*J*multiplier*M_ret*C_KMS/(g*W.sum())*GYR
    order = np.argsort(r)
    cw = np.cumsum(W[order])/W.sum()
    core = order[(cw > .005) & (cw < .995)]
    return dict(
        retained_rest=1., absorbed=multiplier, absorbed_ideal=1.,
        bound_efficiency=efficiency, supply_multiplier=multiplier,
        bound_kinetic_at_injection=kinetic, bound_orbital_energy=orbital,
        escaping_rest_plus_kinetic=multiplier - 1 - kinetic,
        receiver_recoil_bound='<= q/2 of absorbed energy, q=E\'/(M_R c^2)',
        drag_coefficient_declared=kappa, drag_coefficient_gray=4/3,
        receiver_angular_momentum=Jdisk, absorption_weighted_r_vphi_kpc_kms=lever,
        receiver_angular_momentum_loss_ideal=loss_ideal,
        receiver_angular_momentum_loss_declared=loss_declared,
        receiver_kinetic_drain_declared=drain_declared,
        push_history_gyr_median=float(np.median(T_min[core])) if len(core) else float('nan'),
        push_history_gyr_max=float(np.max(T_min[core])) if len(core) else float('nan'),
        mean_abs_zeta=float(np.sum(W*np.abs(zeta))/W.sum()))
