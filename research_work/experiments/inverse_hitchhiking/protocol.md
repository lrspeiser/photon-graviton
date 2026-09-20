# IH-1: observation-constrained coherent hitchhiking

Declared before IH-1 fits. Work under the owner's hypothesis that a matter-fed
companion swirl and temporary photon coupling are the right physical direction.
Find the response that this mechanism must deliver. A calibrated response is
not yet a derivation of that mechanism. No dark matter, expanding background,
distance adjustment, random photon release, or claims of unprecedented theory.

## Equations

Use a deterministic coupling fraction f in [0,1], not a sampled Bernoulli state:
ell df/ds=f_eq-f, f_eq=chi/(1+chi), a_H=a_sat*f_eq.
The observational primary is the local equilibrium limit ell -> 0. Finite
memory is left explicit and unfitted; exact relaxation is a numerical control.
For identical initial data this closure has no stochastic angular spread.

For measured ordinary-matter acceleration g_b:
v_pred^2=r*(g_b+a_H).
Inverse diagnostic f_required=(v_observed^2/r-g_b)/a_sat, with chi_required=
f_required/(1-f_required) only when 0<f_required<1. Keep inadmissible negative
and >=1 values, not clip them into apparent successes. Pointwise inversion
copies the observations and is not an independent prediction.

Define x=ln(g_b/1e-10 m s^-2), y=ln(r/10 kpc),
z=ln(M_source/1e11 Msun), w=tanh(ln(M_source/1e12 Msun)/2).
Five families:
L: ln chi=c+q*x
R: ln chi=c+q*x+s*y
M: ln chi=c+q*x+m*z
RM: ln chi=c+q*x+s*y+m*z
T: ln chi=c+q*x+s*y+m*z+d*w
Four fixed capacities a_sat={3e-9,1e-8,3e-8,1e-7} m/s^2.
These are 20 response variants, not 20 novel theories. q bounds [-2,3],
s,m [-2,2], d [-10,10], c [-25,5].
Use three starts (the q=0.5 scale guess, and two seeded perturbations, seed
190921). Preserve every solver status and boundary hit. Run galaxy-only and
joint fits for each variant: 40 fits / 120 optimization attempts.
No inferred dark-matter mass enters the model.

## Inputs and splits

SPARC: original rotation-model ZIP and catalog, 149 existing selected galaxies
and existing 89/29/31 split. Select names/splits from the frozen selection;
read observed velocities and ordinary component velocities from original data.
M/L_disk=0.5, M/L_bulge=0.7. Source-mass proxy is
(0.5*L_3.6+1.33*M_HI), in solar units; its simplified bulge treatment is explicit.
This proxy is an input to a response law, never an added gravitating halo.
Use all positive-radius rows for those galaxies and report exact count.

X-COP: archived density, gas mass, measured stars when available and the frozen
median stellar/gas fraction from the seven measured systems otherwise.
No hydrostatic inferred total mass or NFW product is used, evaluated or fitted.
Reimplement the established pressure reduction from cl2_sources.py without its
excluded model comparisons. Integrate dP_e/dr=-mu*m_p*n_e*g (mu=0.6) inward;
fit one nonnegative boundary pressure per cluster as a nuisance. Report this
also for transfer clusters: their shape is tested conditional on that boundary.
Grid 600 radii; selected-law verification 1200.
Sort 12 cluster names; index modulo 4 equal 0 or 1 -> training, 2 -> validation,
3 -> test (6/3/3). Published radius/density calibrations stay fixed.

All these observations were exposed in earlier project work. They are withheld
from this fit's parameter estimation, not genuinely new or historically blind
data. No new observations or source geometry are invented.

## Objective and selection

Galaxy term = equal-galaxy mean squared velocity error / (20 km/s)^2.
Cluster term = chi_squared per pressure point / 10, using quoted diagonal
errors and profiled boundary pressure. Joint loss is the sum, with both terms
reported separately. Choose one joint winner by validation sum, tie by smaller
parameter count then identifier. Choose a galaxy-only winner by galaxy
validation error. Freeze choices before scoring test targets or Coma.
An operational comparison target is galaxy validation RMSE <=22 km/s and
cluster validation chi_squared/point <=10; it is not a universal physical
acceptance rule, and it is not changed after seeing scores.

## Light from the same response

Postulate in static weak-field spherical geometry:
alpha(b)=4/c_light^2 int_0^infinity [g_b(r)+eta*a_H(r)]*b/r dz.
Primary eta=1. This is stipulated common optical coupling, not a derivation
from literal quantum graviton attachment. A deterministic potential description
Phi'_eff=g_b+eta*a_H and H_gamma=c_light*|p|*exp(2 Phi_eff/c_light^2)
makes stationary test-ray energy conserved and gives the displayed first-order
bending. The test-particle Hamiltonian is |p|^2/(2m)+m Phi_eff with eta=1.
Full source, binding, field and recoil energies remain owed; these Hamiltonians
do not magically fund the field.

Use Coma's existing ordinary gas model (core296 kpc, beta0.75, truncation3 Mpc)
with the two previously declared gas/stellar normalizations. Common ordinary
mass stays fixed outside3 Mpc. Apply the fitted coupling to radius9 Mpc,
then continue its acceleration as r^-2 for a finite asymptotic potential.
Declare reach sensitivities3 and30 Mpc, not selected on shear.
Report eta={0,0.25,1,4}; eta=0 is ordinary-matter optical reference.
Use alpha/b-alpha' as the geometry-free tangential shear shape.
Compare the six archived Kubo figure-reconstructed tangential shear bins with
one fitted nonnegative amplitude. It absorbs unknown source geometry; it is
NOT an absolute lensing prediction or a distance fit. Retain negative bins,
plotted errors, missing covariance and digitization limitations.
No X-COP cluster has matched shear inputs in this campaign.

## Verification and reporting

Before fits: analytic relaxation and occupancy bounds; pointwise inversion;
unit conversion; point-mass deflection/shear; independent finite difference
of projected deflection; synthetic recovery on a noiseless identified local
family. Verify input/source SHA256 and archive predictions, residuals, all
starts, selected law, inverse requirements, and boundary nuisance values.
After fits: selected pressure quadrature600/1200; Coma projection quadrature
128/256 and derivative step0.002/0.001. Relative numerical target 0.5 percent
on a maximum-profile scale; preserve failures.

Publish equations with fitted coefficients, all block scores, a diagram/plot
from actual predictions, and exactly which assumptions supply any match.
The square-root-acceleration limit resembles Milgrom1983 and is not claimed
as our invention. Logistic saturation, relaxation, potential mechanics and
line-of-sight lensing are credited standard mathematics. Data sources:
SPARC https://arxiv.org/abs/1606.09251 ;
X-COP https://arxiv.org/abs/1805.00042 ;
Kubo https://arxiv.org/abs/0709.0506 ;
Milgrom https://adsabs.harvard.edu/pdf/1983ApJ...270..365M .
