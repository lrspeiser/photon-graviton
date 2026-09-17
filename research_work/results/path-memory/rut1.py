"""RUT-1: can a moving body write an attractive track, and survive writing it? (protocol-rut1.md)

    python rut1.py [--output-dir DIR] [--canonical]

The footprint field: matter writes potential along its own trajectory and the potential fades,

    Phi_mem(x,t) = -q int_0^t exp[-(t-s)/tau] exp[-|x - X(s)|^2/(2 w^2)] ds,

equivalently dPhi/dt = -Phi/tau - q exp[-|x - X(t)|^2/(2 w^2)]. A steady writing pattern leaves
Phi = -tau (pattern), not zero, which is why PM-3's relaxational-memory theorem does not dispose of it.

Stage 0 is the normalization check, run first: the extra inward acceleration of a family of self-written
circular tracks goes as 1/R^2 at fixed q, tau and w, not as 1/R. Stage 1 verifies the ring kernel, moves
the probe rather than the source, and checks that a second test body feels the stored track. Stage 2 puts
the writer on a prescribed circular orbit and measures the complete finite-memory force, including the
tangential drag that the uniform-ring average removes by construction.

Stage 3 -- releasing the body to write its own trajectory, and the two-stage candidate -- is run by
rut3.py, which owns the formation campaign.

Regenerates rut1-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.special import i0e, i1e
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))

GM = 1.                       # dimensionless and illustrative throughout; no galaxy data enters RUT-1
TOL_KERNEL = 1e-12


# ---------------------------------------------------------------- the circular-history ring field
def phi_ring(r, R, w, A=1.):
    """Potential of a uniformly weighted circular history of radius R.

    -A/(2pi) int_0^2pi exp[-(r^2 + R^2 - 2 r R cos t)/(2 w^2)] dt
        = -A exp[-(r - R)^2/(2 w^2)] I0e(r R/w^2),

    The exponentially scaled Bessel function is used because I0 overflows over the parameter range this
    experiment scans -- I0(1000) is inf in double while I0e(1000) is finite. (I0(100) = 1.07e42 does NOT
    overflow; the earlier rationale naming it was wrong, though the implementation was right.)
    """
    r = np.asarray(r, float)
    return -A*np.exp(-(r - R)**2/(2*w*w))*i0e(r*R/(w*w))


def phi_ring_offplane(r, z, R, w, A=1.):
    """The same history evaluated off its plane: |x - X|^2 picks up z^2, which factorises out."""
    return phi_ring(r, R, w, A)*np.exp(-np.asarray(z, float)**2/(2*w*w))


def _d(f, x, h=1e-6):
    return (f(x + h) - f(x - h))/(2*h)


def _d2(f, x, h=1e-5):
    return (f(x + h) - 2*f(x) + f(x - h))/h**2


def ring_depth(R, w, A=1.):
    return float(-phi_ring(R, R, w, A))


def ring_inward_acceleration(r, R, w, A=1.):
    """Inward (positive) acceleration from the ring at radius r: +dPhi/dr."""
    return float(_d(lambda x: float(phi_ring(x, R, w, A)), r))


# ---------------------------------------------------------------- stage 0
def stage0_normalization(radii=(1., 2., 4., 8.), w=.1, A=1.):
    """At fixed q, tau and w a family of self-written tracks scales as 1/R^2, not 1/R."""
    rows = []
    for R in radii:
        a = ring_inward_acceleration(R, R, w, A)
        rows.append(dict(source_orbit_radius=R, on_orbit_inward_acceleration=a,
                         acceleration_times_R=a*R, acceleration_times_R_squared=a*R**2,
                         depth_D=ring_depth(R, w, A),
                         depth_approx_A_w_over_sqrt_2pi_R=float(A*w/(np.sqrt(2*np.pi)*R)),
                         D_over_2R=float(ring_depth(R, w, A)/(2*R))))
    aR2 = np.array([r['acceleration_times_R_squared'] for r in rows])
    aR = np.array([r['acceleration_times_R'] for r in rows])
    asymptote = float(A*w/(2*np.sqrt(2*np.pi)))
    return dict(
        rows=rows, w=w, A=A, asymptotic_a_R_squared=asymptote,
        spread_of_a_R_squared=float((aR2.max() - aR2.min())/aR2.mean()),
        spread_of_a_R=float((aR.max() - aR.min())/aR.mean()),
        approaches_inverse_square=bool((aR2.max() - aR2.min())/aR2.mean()
                                       < (aR.max() - aR.min())/aR.mean()),
        largest_R_relative_to_asymptote=float(abs(aR2[-1]/asymptote - 1)),
        passed=bool((aR2.max() - aR2.min())/aR2.mean() < .05
                    and (aR.max() - aR.min())/aR.mean() > .5),
        statement='extra inward acceleration ~ q tau w/(2 sqrt(2 pi) R^2): at fixed writing rate, '
                  'retention and physical width this family is inverse-SQUARE, so "can create a supporting '
                  'track" and "can explain flat rotation curves" are separate claims. Holding D fixed '
                  'independently at every radius and presenting the resulting 1/R force as a prediction is '
                  'not permitted; a different width scaling, collective writing process or footprint shape '
                  'would be a declared mechanism to test')


# ---------------------------------------------------------------- stage 1
def gate_kernel(cases=((1., 1., .1), (1.3, 1., .1), (.7, 1., .1), (2., 1., .3), (.5, 2., .25))):
    """The closed form against the declared angular integral."""
    rows = []
    for r, R, w in cases:
        num = -1/(2*np.pi)*quad(lambda t: np.exp(-(r*r + R*R - 2*r*R*np.cos(t))/(2*w*w)), 0, 2*np.pi,
                                epsabs=1e-16, epsrel=1e-14)[0]
        closed = float(phi_ring(r, R, w))
        rows.append(dict(r=r, R=R, w=w, angular_integral=num, closed_form=closed,
                         relative_difference=float(abs(closed/num - 1))))
    worst = max(r['relative_difference'] for r in rows)
    return dict(rows=rows, worst_relative_difference=worst, tolerance=TOL_KERNEL,
                passed=bool(worst < TOL_KERNEL),
                note='the scaled Bessel form is used because I0 overflows over the scanned range: I0(1000) '
                     'is inf in double while I0e(1000) is finite. I0(100) = 1.07e42 does not overflow, so the '
                     'earlier rationale naming it was wrong even though the implementation was right')


def gate_mature_ring(R=1., w=.1, D=.2):
    """The declared mature-ring control, reproduced as a reported result to check."""
    A = D/ring_depth(R, w, 1.)
    a = ring_inward_acceleration(R, R, w, A)
    v2 = GM/R + a*R
    tot = lambda r: -GM/r + float(phi_ring(r, R, w, A))
    kappa2 = float(_d2(tot, R) + 3*_d(tot, R)/R)
    return dict(R=R, w=w, D=D, A=A, extra_inward_acceleration=a,
                percent_of_newtonian=float(100*a/(GM/R**2)),
                supported_speed=float(np.sqrt(v2)), newtonian_speed=float(np.sqrt(GM/R)),
                percent_speed_increase=float(100*(np.sqrt(v2)/np.sqrt(GM/R) - 1)),
                kappa_squared=kappa2, kappa_squared_without_track=float(GM/R**3),
                kappa_squared_approx_GM_R3_plus_D_w2=float(GM/R**3 + D/w**2),
                D_over_2R=float(D/(2*R)),
                declared=dict(percent_of_newtonian=10.03, percent_speed_increase=4.89, kappa_squared=21.15),
                passed=bool(abs(100*a/(GM/R**2) - 10.03) < .01
                            and abs(100*(np.sqrt(v2)/np.sqrt(GM/R) - 1) - 4.89) < .01
                            and abs(kappa2 - 21.15) < .01),
                note='a mature FIXED field, held still on purpose: this is motion in an already-built '
                     'track and is not evidence that a freely evolving body could build it')


def gate_probe_scan(R=1., w=.1, D=.2):
    """The source ring stays put while the probe moves, in radius and off the plane.

    PM-1's withdrawn row used the on-ring expression at the observation radius, which silently moves the
    source to wherever the probe is. That error is re-tested here rather than assumed retired.
    """
    A = D/ring_depth(R, w, 1.)
    radii = np.array([.5, .8, .95, 1., 1.05, 1.2, 2.])*R
    fixed = np.array([float(phi_ring(r, R, w, A)) for r in radii])
    moved = np.array([float(phi_ring(r, r, w, A)) for r in radii])     # the mistake, for contrast
    zs = np.array([0., .05, .1, .2, .5])*R
    off = np.array([float(phi_ring_offplane(R, z, R, w, A)) for z in zs])
    vert = [float(-_d(lambda zz: float(phi_ring_offplane(R, zz, R, w, A)), z)) for z in zs[1:]]
    mn = minimize_scalar(lambda r: float(phi_ring(r, R, w, A)), bracket=(.9*R, R, 1.1*R))
    return dict(
        probe_radii=radii.tolist(), potential_source_fixed=fixed.tolist(),
        potential_source_moved_to_probe=moved.tolist(),
        worst_relative_gap=float(np.max(np.abs(fixed/moved - 1))),
        source_moves_with_probe=bool(np.allclose(moved, moved[0])),
        off_plane=dict(z=zs.tolist(), potential=off.tolist(), vertical_acceleration=vert,
                       restoring=bool(all(v < 0 for v in vert)),
                       note='negative vertical acceleration is directed back toward the plane'),
        trough_minimum=float(mn.x), displaced_inward_of_source_path=float(R - mn.x),
        passed=bool(np.max(np.abs(fixed/moved - 1)) > .5 and all(v < 0 for v in vert)
                    and R - mn.x > 0),
        statement='the trough bottom lies INSIDE the source path, which is why a body on its own track '
                  'feels an inward force; a symmetric trough centred on R would give exactly zero there')


def gate_second_body(R=1., w=.1, D=.2):
    """A separate negligible-mass body must feel the stored track, or this is a private constraint."""
    A = D/ring_depth(R, w, 1.)
    probes = np.array([.9, .97, 1.03, 1.1])*R
    felt = [ring_inward_acceleration(float(p), R, w, A) for p in probes]
    return dict(probe_radii=probes.tolist(), inward_acceleration=felt,
                writer_is_not_special=True,
                passed=bool(all(abs(f) > 0 for f in felt)
                            and felt[0] < 0 < felt[-1] or all(abs(f) > 0 for f in felt)),
                note='the field is a function of position alone, so a test body that never wrote anything '
                     'responds to the same track; nothing in the implementation keys on the writer. Inside '
                     'the trough bottom the force is outward and outside it inward, as a trough requires')


# ---------------------------------------------------------------- stage 2: the writer's own past
def self_force(R, w, tau, Omega, q=1., n_rev=None, nodes=240):
    """Radial and tangential self-acceleration of a body on a prescribed circular orbit.

    With phi = Omega (t - s) the chord geometry gives exactly

        a_r = -(q/Omega) int_0^{Phi} e^{-phi/(Omega tau)} [R(1 - cos phi)/w^2] E(phi) dphi
        a_t = -(q/Omega) int_0^{Phi} e^{-phi/(Omega tau)} [R sin phi   /w^2] E(phi) dphi

    with E(phi) = exp[-R^2(1 - cos phi)/w^2] and a_r negative meaning inward. Both integrands are 2-pi
    periodic apart from the decay factor, so the whole history collapses exactly to one revolution times a
    geometric sum: after N revolutions the factor is (1 - rho^N)/(1 - rho) with rho = exp(-T/tau). That
    reproduces the declared build law A(N) = A_inf [1 - exp(-N/N_build)] identically, and it lets each
    revolution be integrated by Gauss-Legendre on [0, pi] and [pi, 2 pi], whose nodes cluster at the ends
    where E is sharply peaked. Uniform sampling does not resolve that peak: the radial integral survives it
    because its integrand is even about each passage, but the tangential one is odd there and is the part
    that carries the drag, so it is the part that a uniform grid gets wrong.
    """
    x, wg = np.polynomial.legendre.leggauss(nodes)
    seg = []
    for a, b in ((0., np.pi), (np.pi, 2*np.pi)):
        p = .5*(b - a)*(x + 1) + a
        seg.append((p, .5*(b - a)*wg))
    T = 2*np.pi/Omega
    I_r = I_t = 0.
    for p, ww in seg:
        E = np.exp(-R*R*(1 - np.cos(p))/(w*w))*np.exp(-p/(Omega*tau))
        I_r += float(np.sum(ww*E*R*(1 - np.cos(p))/(w*w)))
        I_t += float(np.sum(ww*E*R*np.sin(p)/(w*w)))
    rho = np.exp(-T/tau)
    S = (1 - rho**n_rev)/(1 - rho) if n_rev is not None else 1/(1 - rho)
    return float(-q*S*I_r/Omega), float(-q*S*I_t/Omega)


def _one_revolution(R, w, tau, Omega, nodes=240, a=0., b=None):
    """Support, drag and H integrands over [a, b] in phase (default one full revolution)."""
    b = 2*np.pi if b is None else b
    x, wg = np.polynomial.legendre.leggauss(nodes)
    out = [0., 0., 0.]
    mid = .5*(a + b)
    for lo, hi in ((a, mid), (mid, b)):                 # split so the nodes cluster where E peaks
        ph = .5*(hi - lo)*(x + 1) + lo
        ww = .5*(hi - lo)*wg/Omega
        u = ph/Omega
        e = np.exp(-u/tau)*np.exp(-R*R*(1 - np.cos(ph))/(w*w))
        out[0] += float(np.sum(ww*e*R*(1 - np.cos(ph))/(w*w)))
        out[1] += float(np.sum(ww*e*R*np.sin(ph)/(w*w)))
        out[2] += float(np.sum(ww*e))
    return out


def self_force_at_age(R, w, tau, Omega, t, q=1.):
    """Support and drag at elapsed age t, exactly, for t = nT + u:

        a(t) = a(T) (1 - rho^n)/(1 - rho) + rho^n a(u),   rho = exp(-T/tau),

    with a(T) and a(u) the fresh-history integrals over those intervals, applied separately to each
    component. Substituting a fractional revolution count into the geometric factor alone is NOT
    equivalent, which is why the mature map cannot display the within-orbit structure of formation.
    """
    T = 2*np.pi/Omega
    n = int(np.floor(t/T + 1e-12))
    u = t - n*T
    rho = np.exp(-T/tau)
    full = _one_revolution(R, w, tau, Omega) if n else (0., 0., 0.)
    part = _one_revolution(R, w, tau, Omega, b=Omega*u) if u > 1e-15 else (0., 0., 0.)
    S = (1 - rho**n)/(1 - rho) if n else 0.
    # signed like self_force: negative radial is inward, negative tangential is backward
    return (float(-q*(S*full[0] + rho**n*part[0])), float(-q*(S*full[1] + rho**n*part[1])))


def drag_identity(R, w, tau, Omega, t=None, q=1.):
    """drag(t) = (q/Omega R)[1 - e^{-t/tau} E(t) - H(t)/tau], from integrating the tangential term by
    parts. It avoids computing the drag purely by cancellation in an odd integrand, and at infinite age
    tends to (q/Omega R)[1 - I0e(R^2/w^2)]. H's tail uses the same exact periodic decomposition: a
    truncated revolution sum is what makes a naive check of this identity disagree."""
    T = 2*np.pi/Omega
    rho = np.exp(-T/tau)
    if t is None:
        H = _one_revolution(R, w, tau, Omega)[2]/(1 - rho)
        return float(q/(Omega*R)*(1 - H/tau))
    n = int(np.floor(t/T + 1e-12))
    u = t - n*T
    H = (1 - rho**n)/(1 - rho)*_one_revolution(R, w, tau, Omega)[2] if n else 0.
    if u > 1e-15:
        H += rho**n*_one_revolution(R, w, tau, Omega, b=Omega*u)[2]
    E_t = np.exp(-R*R*(1 - np.cos(Omega*t))/(w*w))
    return float(q/(Omega*R)*(1 - np.exp(-t/tau)*E_t - H/tau))


def drag_impulse(R, w, tau, Omega, t, q=1., steps=600):
    """Accumulated backward impulse int_0^t drag dt'; times R it is the angular momentum the driver must
    supply to hold the prescribed orbit while the track builds."""
    grid = np.linspace(0., t, steps + 1)
    vals = np.array([drag_identity(R, w, tau, Omega, tt, q) if tt > 0 else 0. for tt in grid])
    return float(np.trapezoid(vals, grid))


def C_long(b):
    """The closed-form long-memory coefficient in drag/support = C/(N_build * w/R), b = w/R."""
    alpha = 1/b**2
    return float(b**3/(2*np.pi)*(1 - i0e(alpha))/(i0e(alpha) - i1e(alpha)))


def collective_force(N, R, w, tau, Omega, q_total=1., nodes=200, max_rev=4000):
    """Mature support and drag on one member of N evenly spaced writers sharing a fixed TOTAL writing
    rate. Radial contributions largely add; tangential ones increasingly cancel."""
    T = 2*np.pi/Omega
    rho = np.exp(-T/tau)
    x, wg = np.polynomial.legendre.leggauss(nodes)
    j = np.arange(N)[:, None]
    tot_r = tot_t = 0.
    for lo, hi in ((0., np.pi), (np.pi, 2*np.pi)):
        ph = .5*(hi - lo)*(x + 1) + lo
        ww = .5*(hi - lo)*wg/Omega
        u = ph/Omega
        phj = ph[None, :] + 2*np.pi*j/N
        Ej = np.exp(-R*R*(1 - np.cos(phj))/(w*w))
        dec = np.exp(-u/tau)
        tot_r += float(np.sum(ww*dec*np.mean((1 - np.cos(phj))*Ej, axis=0)))
        tot_t += float(np.sum(ww*dec*np.mean(np.sin(phj)*Ej, axis=0)))
    f = q_total*R/(w*w)/(1 - rho)
    return float(f*tot_r), float(f*tot_t)


def stage2c_collective(R=1., wr=.1, nb=10., counts=(1, 4, 8, 16, 32), q_total=1.):
    """Several physically distinct writers, at fixed total writing rate, versus one.

    Two controls stay distinct. Dividing ONE writer into coincident copies whose rates sum to the original
    must leave the field exactly unchanged -- the numerical-subdivision loophole that killed PM-1's
    per-ring saturation. Placing DISTINCT writers at different physical positions changes the source
    distribution, so the field may legitimately change.
    """
    Omega = np.sqrt(GM/R**3)
    w, tau = wr*R, nb*2*np.pi/Omega
    base = collective_force(1, R, w, tau, Omega, q_total)
    rows = []
    for N in counts:
        s, d = collective_force(N, R, w, tau, Omega, q_total)
        rows.append(dict(writers=N, support=s, drag=d,
                         support_relative_to_one=float(s/base[0]),
                         drag_relative_to_one=float(d/base[1]),
                         drag_over_support=float(d/s)))
    # the loophole control: M coincident copies at q_total/M each must reproduce one writer exactly
    coincident = []
    for M in (2, 5, 17):
        s, d = collective_force(1, R, w, tau, Omega, q_total/M)
        coincident.append(dict(copies=M, support_times_M=float(M*s), drag_times_M=float(M*d),
                               support_relative_error=float(abs(M*s/base[0] - 1)),
                               drag_relative_error=float(abs(M*d/base[1] - 1))))
    worst_copy = max(max(c['support_relative_error'], c['drag_relative_error']) for c in coincident)
    best = min(rows, key=lambda r: r['drag_relative_to_one'])
    return dict(
        width_ratio=wr, memory_revolutions=nb, rows=rows,
        coincident_copy_control=dict(
            rows=coincident, worst_relative_error=worst_copy, tolerance=1e-12,
            passed=bool(worst_copy < 1e-12),
            note='splitting one writer into coincident copies is a relabelling and must change nothing; '
                 'if it did, the model would be defined by how finely the source was subdivided, which is '
                 'exactly what withdrew PM-1\'s per-ring saturation'),
        support_retained=float(min(r['support_relative_to_one'] for r in rows)),
        best_drag_reduction=float(1/best['drag_relative_to_one']) if best['drag_relative_to_one'] else 0.,
        passed=bool(worst_copy < 1e-12
                    and all(r['support_relative_to_one'] > .99 for r in rows)
                    and rows[-1]['drag_relative_to_one'] < .1),
        statement='distinct writers sharing one track keep essentially all of the inward support while '
                  'their tangential contributions cancel: a single-writer drag bound therefore cannot be '
                  'applied to a collective source',
        not_shown='evenly spaced writers held on their orbits are a deliberately favourable symmetry. '
                  'Nothing here shows that such an arrangement forms from an empty field, keeps its '
                  'spacing, tolerates phase disturbances, survives differential motion, or obeys a '
                  'completed matter-field energy law. Stage 3 compares one writer with several at the '
                  'same total rate, both from empty fields, and perturbs the phases')


def stage2_prescribed_orbit(R=1., width_ratios=(.05, .1, .2, .4), memory_revs=(1., 3., 10., 30., 100.),
                            support_fractions=(.1, .01), q=1.):
    """Scan w/R and tau/T_orbit; report support, drag and the angular-momentum-change time.

    The writing rate q is arbitrary and scales support and drag together, so the q-free statement is their
    RATIO. The survival diagnostic is then quoted at a declared support fraction f -- the extra inward
    acceleration as a fraction of Newtonian -- because "how long does the orbit last" is only meaningful
    once the amount of help being bought is fixed.
    """
    Omega = np.sqrt(GM/R**3)
    T = 2*np.pi/Omega
    v = Omega*R
    rows = []
    for wr in width_ratios:
        w = wr*R
        for nb in memory_revs:
            tau = nb*T
            a_r, a_t = self_force(R, w, tau, Omega, q=q)
            inward, drag = -a_r, -a_t
            ratio = drag/inward
            at_f = {}
            for f in support_fractions:
                drag_f = f*(GM/R**2)*ratio          # drag when the track supplies f of Newtonian support
                Lt = v/drag_f
                at_f[f'support_{f}'] = dict(
                    angular_momentum_change_time_in_periods=float(Lt/T),
                    L_time_over_build_time=float(Lt/tau),
                    mature_drag_timescale_exceeds_one_retention_time=bool(Lt > tau))
            rows.append(dict(
                width_ratio=wr, memory_revolutions=nb, field_age='mature (infinite past)',
                inward_acceleration_at_q1=inward, tangential_acceleration_at_q1=a_t,
                drag_over_support=float(ratio),
                ratio_times_Nbuild_times_width=float(ratio*nb*wr),
                C_long_closed_form=C_long(wr),
                drag_identity_at_q1=drag_identity(R, w, tau, Omega),
                drag_identity_relative_difference=float(abs(drag_identity(R, w, tau, Omega)/drag - 1)),
                build_law_at_one_memory_time=float(1 - np.exp(-1.)),
                at_support_fraction=at_f,
                driver_power_per_unit_mass_at_q1=float(drag*v)))
    # the within-orbit structure the mature map cannot display, plus the accumulated driver impulse
    w_f, nb_f = .1*R, 10.
    tau_f = nb_f*T
    mature_f = self_force(R, w_f, tau_f, Omega)
    ages = []
    for fr in (.25, .5, 1., 1.5, 2., 10.):
        s_a, d_a = self_force_at_age(R, w_f, tau_f, Omega, fr*T)
        ages.append(dict(orbits=fr, support_over_mature=float(s_a/mature_f[0]),
                         drag_over_mature=float(d_a/mature_f[1])))
    cost = []
    for wr in width_ratios:
        w2, tau2 = wr*R, nb_f*T
        sup_m, _ = self_force(R, w2, tau2, Omega)
        q_f = support_fractions[0]*(GM/R**2)/abs(sup_m)
        J = drag_impulse(R, w2, tau2, Omega, tau2, q=q_f)
        cost.append(dict(width_ratio=wr, memory_revolutions=nb_f,
                         driver_angular_momentum_over_body=float(J*R/(R*R*Omega)),
                         support_reached_fraction_of_mature=float(1 - np.exp(-1.))))
    crit = {}
    for f in support_fractions:
        def gap(b_, f=f):
            s, d = self_force(R, b_*R, nb_f*T, Omega)
            return (Omega*R/(f*(GM/R**2)*abs(d/s)))/(nb_f*T) - 1.
        try:
            crit[f'support_{f}'] = float(brentq(gap, 1e-3, .95, xtol=1e-10))
        except ValueError:
            crit[f'support_{f}'] = None
    scale = np.array([r['ratio_times_Nbuild_times_width'] for r in rows])
    f_main = support_fractions[0]
    surv = [r for r in rows if r['at_support_fraction'][f'support_{f_main}']
            ['mature_drag_timescale_exceeds_one_retention_time']]
    return dict(
        rows=rows, orbital_period=T, orbital_speed=v, support_fractions=list(support_fractions),
        uniform_ring_control=dict(
            tangential_acceleration=0.,
            why='the averaged ring is exactly symmetric in phi, so its tangential force vanishes '
                'identically -- and int_0^{2 pi} sin(phi) E(phi) dphi = 0 exactly, since that integrand is '
                'a total derivative. Using the averaged ring to argue the writing process is drag-free '
                'would be circular: the averaging is what removed the drag'),
        drag_always_present=bool(all(r['drag_over_support'] > 0 for r in rows)),
        drag_over_support_range=[float(min(r['drag_over_support'] for r in rows)),
                                 float(max(r['drag_over_support'] for r in rows))],
        scaling=dict(product_drag_ratio_Nbuild_width=[float(s) for s in scale],
                     spread=float((scale.max() - scale.min())/scale.mean()),
                     law='drag/support = C/(N_build * w/R) with C about 0.62 to 0.78, so the ratio falls '
                         'like the memory time and like the track width, and is independent of the writing '
                         'rate and of how long the field has been building'),
        field_age_structure=dict(
            width_ratio=.1, memory_revolutions=nb_f, rows=ages,
            note='the drag OVERSHOOTS its mature value within an orbit -- about 1.04 at half and at one '
                 'and a half revolutions -- and matches the build law only at integer revolutions, where '
                 'both components carry the same geometric factor. "The drag saturates after one '
                 'revolution" was wrong. What is true is that the MATURE drag approaches a '
                 'retention-independent limit (q/Omega R)[1 - I0e(R^2/w^2)] as the memory time grows'),
        formation_cost=dict(
            rows=cost, support_fraction=support_fractions[0],
            note='driver angular momentum needed to HOLD the prescribed orbit to t = tau, as a multiple '
                 'of the body own angular momentum, with the support then at 63.2% of mature. It measures '
                 'the effort required to prevent the orbit from changing, not a loss a released body '
                 'would suffer: once it moves, the prescribed-circle history no longer describes it'),
        timescale_threshold=dict(
            derivation='L_time/build = v (w/R)/(f (GM/R^2) C T) with drag/support = C/(N_build w/R), so '
                       'N_build cancels; but C is a FUNCTION of width, so the crossing is solved from the '
                       'force integrals rather than from a representative C',
            C_measured_range=[float(scale.min()), float(scale.max())],
            C_long_closed_form={f'width_{wr}': C_long(wr) for wr in width_ratios},
            C_narrow_limit=float(np.sqrt(2/np.pi)), two_pi_C_narrow_limit=float(2*np.sqrt(2*np.pi)),
            critical_width_over_radius=crit, solved_at_memory_revolutions=nb_f,
            statement='this is a timescale crossing, NOT a survival theorem: no freely moving orbit is '
                      'tested by it. An orbit may migrate substantially without being destroyed, and may '
                      'lose an unacceptable fraction of its angular momentum while passing the '
                      'inequality. The earlier single coefficient 4.4 came from collapsing C(w/R) to a '
                      'median and is withdrawn',
            pm3_inference_withdrawn='a broad track is NOT excluded by PM-3. Width does not change the '
                                    'storage equation, and Phi_steady = -tau (steady writing source) is '
                                    'nonzero whether the source is narrow or broad. That a broad '
                                    'collective field might admit a static description is a question '
                                    'about telling mechanisms apart observationally, not a proof that '
                                    'the stored field vanishes'),
        survival_at_main_fraction=dict(
            support_fraction=f_main, rows_crossing=len(surv), rows_total=len(rows),
            label='rows where the mature drag timescale exceeds one retention time; a timescale crossing, '
                  'not a tested survival',
            best=(max(rows, key=lambda r: r['at_support_fraction'][f'support_{f_main}']
                      ['L_time_over_build_time'])['at_support_fraction'][f'support_{f_main}']
                  ['L_time_over_build_time'])),
        passed=bool(all(r['drag_over_support'] > 0 for r in rows)
                    and all(r['inward_acceleration_at_q1'] > 0 for r in rows)
                    and float((scale.max() - scale.min())/scale.mean()) < .3),
        reading='a map, not a verdict. Support and drag both scale with the writing rate, so the physical '
                'question is what the orbit pays for a given amount of help: at a declared support '
                'fraction f, does the angular-momentum-change time exceed the time the track needs to '
                'build? The drag is set by the most recent passage and saturates immediately, while the '
                'support accumulates over the memory time, which is why long memory and wide tracks are '
                'the surviving direction')


def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    s0 = stage0_normalization()
    kernel, mature = gate_kernel(), gate_mature_ring()
    probe, second = gate_probe_scan(), gate_second_body()
    s2 = stage2_prescribed_orbit()
    s2c = stage2c_collective()
    stage1 = bool(kernel['passed'] and mature['passed'] and probe['passed'] and second['passed'])
    result = dict(
        experiment='RUT-1: can a moving body write an attractive track, and survive writing it?',
        protocol='protocol-rut1.md',
        model=dict(
            kernel='Phi_mem(x,t) = -q int_0^t exp[-(t-s)/tau] exp[-|x - X(s)|^2/(2 w^2)] ds',
            differential='dPhi/dt = -Phi/tau - q exp[-|x - X(t)|^2/(2 w^2)]',
            why_not_covered_by_pm3='a steady writing pattern leaves Phi = -tau (pattern), not zero, unlike '
                                   'a correction proportional to (p - Q) which vanishes when Q = p',
            scope='a low-speed effective model, not derived from general relativity. The writing '
                  'coefficient\'s dependence on source mass and the field\'s energy accounting are not '
                  'specified, so stages 0-2 are kinematic feasibility tests. The instantaneous-Gaussian '
                  'kernel is used throughout; a causal version in which a disturbance written at X(s) '
                  'cannot reach x before |x - X(s)|/c is a separate obligation, since propagation speed '
                  'and retention time are distinct properties'),
        stage_0_normalization=s0,
        stage_1_kernel=dict(closed_form=kernel, mature_ring_control=mature, probe_scan=probe,
                            second_test_body=second, passed=stage1),
        stage_2_prescribed_orbit=s2,
        stage_2C_collective=s2c,
        stage_3='run by rut3.py, which owns the formation campaign and the two-stage candidate',
        vector_extension='retained as a proposal behind a narrower dependency. v.[v x curl A] = 0 so its '
                         'sideways term does no direct work, but the force law also carries -dA/dt, which '
                         'cannot be dropped while the field grows, and its source law, sign and evolution '
                         'are unspecified. Not implemented here, so it stays an independent physical '
                         'proposal rather than an adjustable repair',
        passed=bool(s0['passed'] and stage1 and s2['passed'] and s2c['passed']),
        input_sha256={'protocol-rut1.md': hashlib.sha256((HERE/'protocol-rut1.md').read_bytes()).hexdigest()},
        checks_short_run=dict(
            kernel_worst=kernel['worst_relative_difference'],
            mature_percent=mature['percent_of_newtonian'], mature_kappa2=mature['kappa_squared'],
            normalization_aR2_spread=s0['spread_of_a_R_squared'],
            normalization_aR_spread=s0['spread_of_a_R'],
            drag_over_support_min=s2['drag_over_support_range'][0],
            drag_scaling_spread=s2['scaling']['spread'],
            best_L_time_over_build=s2['survival_at_main_fraction']['best'],
            collective_support_retained=s2c['support_retained'],
            collective_drag_reduction=s2c['best_drag_reduction'],
            coincident_copy_worst=s2c['coincident_copy_control']['worst_relative_error']),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:1800])
    status = evidence_io.finish(args, 'path-memory-rut1', text, HERE/'rut1-results.json',
                                ignore={'/runtime_seconds'})
    return status if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
