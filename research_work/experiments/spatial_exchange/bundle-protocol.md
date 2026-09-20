# SE-B: evolving light bundles and massive test bodies

Declared 20 September 2026 before execution. Use the pinned SE-2 background
equations and SE-P test-body Hamiltonian. Integrate background and probes at the
same four RK4 stage times. The probes have infinitesimal statistical weight:
their mass/momentum normalization sets trajectories but adds no finite energy
or recoil to the background ledger. This is a declared test-body limit.

Initial light travels along +x from x=-1.5, centered at y=1,z=0. Nine rays:
center and +/-y,+/-z offsets for each spacing .05 and .025. Canonical p=(1,0,0),
m=0. Detector plane x=1.5. Interpolate the first crossing within a step linearly;
record position, coordinate time, momentum, Hamiltonian velocity and coordinate
energy. Flat reference flight time is 3. Coordinate energy ratios are not yet
observed frequency shifts; an observer/clock model is needed.

Six massive test bodies start equally spaced on radius2 in the xy plane with
tangential p/m=.02 and m=1. Report positions, velocities and radial displacement
at T4. This short interval is not an orbital-stability or flat-rotation test.
Spherical probe radius .9 by default. Stop before any probe kernel touches the
boundary. Background L16,T4,radius .9,excitation .001 per source as in SE-2.

Seven runs: n32 dt.02 emission-only (theta0,chi0,mix0), X emission chi200,
mixed emission theta pi/4 chi200, Y emission theta pi/2 chi200, then Y at
dt.01, Y at n40 dt.01, and Y with probe radius1.2. Keep all other parameters
fixed. Initial fields remain empty. No frozen lens or external wave reservoir.

Use detector transverse differences to estimate the 2x2 parallel-ray transport
Jacobian M at each bundle spacing. Report trace distortion 1-tr(M)/2, symmetric
trace-free components, antisymmetric rotation and 1/abs(det M). These are
transport diagnostics in the specified model coordinates, not cosmological
convergence/shear or astronomical magnification until source/observer geometry
is derived. Report central velocity deflection atan2(v_y,v_x), arrival offset
t-3 and transverse detector displacement. Do not substitute these for complete
image positions, magnifications and observed time delays.

Numerical gates: all background SE-2 energy/cone/edge/clearance gates; all rays
cross by T4; probe cone excess<1e-10 at recorded steps. Two bundle spacings must
agree in M to max absolute .001. Where a matching prior SE-1/SE-2 state exists,
background final-state difference must be <1e-12 (unaltered by probes).
Y time/space comparisons require central bend agreement within 1%/5%, using
max(abs(finer bend),1e-8), arrival offset within1e-4/.001, M within .001/.005.
Report probe-radius differences without accepting a point-ray interpretation.

Archive every background endpoint, all probe-step states/velocities/energies,
all detector events and failed gates. Independently reconstruct crossings and
bundle diagnostics before interpretation. Source refinement continues separately;
these calculations cannot establish long-time stability or observational success.
Hamiltonian optics and bundle derivatives are established mathematics, credited;
the shared response law is our candidate assumption, not an old-theory target.
