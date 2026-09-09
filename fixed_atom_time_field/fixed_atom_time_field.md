# A propagation field with fixed atomic standards

Research note, 8 September 2026. This is a proposed effective model for the fictional nonexpanding universe. It derives consequences from explicit postulates. It is not a derivation from established physics, and it does not yet prove exact protection of atomic behavior in a microscopic theory.

## The starting requirement

Take matter's background spatial geometry to be static. Atomic transition frequencies, particle masses and material ruler lengths are stipulated to be constant in matter's local units. For the homogeneous calculation, clocks measure t. Ordinary gravitational and kinematic clock comparisons can be added later; they are not part of this first flat-background calculation. No density screening is introduced.

Introduce a positive, dimensionless propagation field n(x,t), normalized to n=1 at the present observer. The field controls the relation between a traveling wave's spatial variation and its temporal oscillation:

\[
H_\gamma(\mathbf x,\mathbf p,t)=\frac{c_0|\mathbf p|}{n(\mathbf x,t)}.
\]

This photon Hamiltonian is a postulate. It supplies a specific meaning for the proposed time effect relative to atomic clocks: traveling light's frequency and group speed depend on n. It is not merely a relabeling of time, because matter standards are held fixed. A homogeneous n changes coordinate frequency without changing coordinate wave number, so “stretching” here initially means longer temporal periods, not a longer coordinate wavelength.

Hamilton's equations give

\[
\dot{\mathbf x}=\frac{c_0}{n}\hat{\mathbf p},\qquad
\dot{\mathbf p}=\frac{H_\gamma}{n}\nabla n,\qquad
\frac{d\ln\omega}{dt}=-\partial_t\ln n.
\]

The last derivative on the right is a partial derivative at the ray's position. In this Hamiltonian, a static spatial variation changes propagation speed and direction but does not generate cumulative frequency change between fixed matter clocks. A temporally evolving background is essential. Time-varying photonic media provide an established mathematical analogue for momentum conservation with frequency conversion [1,2]; they are not evidence that intergalactic vacuum behaves this way.

## A new dynamical route to the exponential

Give n a simple positive field energy through the Lagrangian density

\[
\mathcal L_n=\frac K2\left[(\partial_t n)^2-v_n^2|\nabla n|^2\right]-V(n),\qquad K>0.
\]

K sets the energy scale; vₙ is the propagation speed of field disturbances. Together with a photon ensemble described by the preceding Hamiltonian, this is a local classical field–ray model. It is not yet a quantized electromagnetic theory coupled consistently to atoms.

In a homogeneous region with constant V and negligible radiation backreaction, its equation is

\[
\ddot n=0,\qquad n(t)=1+\gamma t.
\]

The constant γ is initial field velocity. Choosing the positive direction gives redshift. Neither its sign nor its magnitude is predicted by the free field alone.

For a source at fixed geometric separation R, emitting at tₑ and received today at tᵣ=0,

\[
R=\int_{t_e}^{0}\frac{c_0\,dt}{1+\gamma t}
=\frac{c_0}{\gamma}\ln\frac1{n(t_e)}.
\]

Homogeneity conserves photon momentum, so νᵣ/νₑ=nₑ/nᵣ. With an unchanged atomic reference this yields

\[
\boxed{1+z=\frac{n_r}{n_e}=\exp\left(\frac{\gamma R}{c_0}\right)}.
\]

Thus the supplied law is obtained with γ/c₀=κ and κ=0.000077315 per million light-years. At today's normalization, γ=7.7315×10⁻¹¹ per year. The numerical coefficient remains fitted, rather than independently derived.

The new ingredient relative to the preceding note is **which field has the simple kinetic energy**. A canonical free χ=ln n gives n=e^(γt) and the linear distance law 1+z=1+γR/c₀. A canonical free n gives the exponential distance law. In χ variables, the latter kinetic term is K e^(2χ)(∂χ)²/2, still positive. Therefore reproducing an affine n does not universally require the unstable potential that arose when we insisted on canonical χ. These are distinct choices of dynamics; a field redefinition must transform the kinetic term as well. No principle yet uniquely selects one choice.

The affine history has n=0 at t=−1/γ, where the propagation Hamiltonian becomes singular. This is a boundary of the approximation, not a derived Big Bang or an established universe age. A global model needs a continuation or additional dynamics.

## Whole-signal stretching follows from the same rays

Hold R fixed and differentiate the light-travel integral between neighboring emission and reception events:

\[
\frac{dt_r}{dt_e}=\frac{n_r}{n_e}=1+z.
\]

Because atomic clocks measure t in this model, the result predicts longer arrival intervals while leaving the emitting and receiving clocks locally unchanged. For affine n it holds for finite intervals too:

\[
1+\gamma t_r=e^{\gamma R/c_0}(1+\gamma t_e).
\]

The script integrates two rays independently at dimensionless distance X=γR/c₀=0.1. Their measured pulse stretch is 1.10517091809084 versus the analytic 1.10517091807565, a relative numerical difference of 1.37×10⁻¹¹. This checks the ray mapping. Actual source light-curve dynamics are an additional modeling task.

## Energy exchange changes the history

For an isotropic homogeneous bath with conserved adiabatic mode occupations in fixed volume, radiation energy density is uγ=A/n. Its loss is

\[
\dot u_\gamma=-\frac{\dot n}{n}u_\gamma.
\]

Varying the combined field–radiation energy gives

\[
\boxed{K\ddot n+V'(n)=\frac{u_\gamma}{n}},\qquad
\boxed{\mathcal E=\frac K2\dot n^2+V(n)+\frac A n=\text{constant}}.
\]

For V=0 the radiation's energy loss increases the field's kinetic energy. A finite radiation bath therefore prevents n from remaining exactly affine. The exponential law is a negligible-backreaction limit, not an exact result for arbitrary radiation density. We have given an energy reservoir to a physical field; we have not established that time itself consumes energy.

Define η=uγ,0/(Kγ²), using present n=1 and ṅ=γ. The code integrates the conservative system backward along rays at fixed geometric distances. At X=7, where the unperturbed stretch would be 1096.633:

| η | Calculated stretch | Difference from exponential |
|---:|---:|---:|
| 0 | 1096.633 | Numerical agreement |
| 0.000001 | 1095.440 | −0.109% |
| 0.0001 | 986.314 | −10.06% |
| 0.001 | 459.162 | −58.13% |
| 0.1 | 0.241 | Blueshift at this distance |

Across the runs, maximum relative energy error was below 1.6×10⁻¹⁵. These are synthetic, dimensionless numerical experiments, not measurements or new fits to astronomical data.

The homogeneous bath gives an analytic constraint. Nonnegative field kinetic energy requires

\[
n\ge\frac{\eta}{\tfrac12+\eta},\qquad
\boxed{(1+z)_{max}=1+\frac1{2\eta}}.
\]

The past trajectory reaches a positive minimum of n, turns, and eventually produces decreasing redshift; this replaces the free model's singular boundary in this reduced system. It does not establish that the universe is cyclic or solve gravitational backreaction. To obtain even a stretch of 1100.7 requires η≤0.0004547, equivalently present field kinetic energy at least 1099.7 times this bath's present radiation energy. This is only a necessary condition; keeping the exponential accurate over that range requires a stronger condition. Those energy densities must ultimately be coupled to gravity, not treated as invisible by assumption.

## Atomic invariance: what is assumed and what remains to derive

The effective model deliberately leaves atomic standards fixed. It therefore accomplishes the requested ray-level derivation without shrinking rulers or changing the reference transition. However, it has not yet accomplished exact local atomic invariance from microscopic physics.

Real propagating photons and the electromagnetic interactions within atoms belong to the same underlying theory. Specifying a Hamiltonian for free photons and fixed atomic levels separately does not prove their consistency once emission, absorption, magnetic interactions and quantum corrections are included. A minimal change to Maxwell's magnetic term can leave the leading electrostatic hydrogen potential unchanged, but that alone does not protect all atomic observables. Scalar couplings to light and matter can have nontrivial spectroscopic consequences [3].

Thus the postulate “atoms remain unchanged” must become either an exact symmetry demonstrated in a full action or an approximation whose residual effects are calculated and compared to precision data. The present effective Hamiltonian does neither yet. It also selects a preferred background frame; a relativistic completion, matter dynamics and gravitational-wave propagation are outstanding.

This is not screening: n acts everywhere, including the Solar System. With fixed rulers and clocks the model predicts a physical local drift cγ=c₀/n, currently d ln cγ/dt=−γ. A fixed-size ideal cavity compared with an atomic clock would therefore drift at −7.73×10⁻¹¹ per year if the same nondispersive propagation law applies to cavity light and its material boundaries remain fixed. Real cavity dimensions and dispersion need calculation before interpreting measurements. Defining an SI metre through c cannot remove a dimensionless cavity/atom frequency-ratio prediction. This is a direct test of the fixed-atom branch, beyond accumulated intergalactic redshift.

For the exact prescribed affine history, sources at fixed R also have zero redshift drift: d z/dtᵣ=0. That distinguishes it from canonical rolling χ, for which d z/dtᵣ=γz at today's normalization. Moving sources and gravitational effects require separate terms. Radiation backreaction changes these predictions.

## Temperature and scope

For adiabatic, frequency-independent changes preserving occupations, an existing Planck distribution transforms to Tᵣ=Tₑ/(1+z) in the fixed atomic standards stipulated here. In the fixed-coordinate-volume model uγ∝n³T⁴∝1/n, consistently with the energy equation. Arbitrary radiation is not made thermal by this operation. No temperature normalization, thermalization mechanism, CMB angular structure, preferred-axis alignment or galaxy rotation law is derived in this note.

## Milestone outcome

We now have a positive-energy effective propagation-field model that derives the requested exponential and whole-signal stretching in its free-background limit. It also derives a conservative correction when radiation affects the field, exposing why distant extrapolation is not automatic. The next milestone is the microscopic protection of atomic standards, including the cavity/atom comparison, followed by coupling the field's energy to gravity. Exact atomic invariance, global stability, the fitted rate and observational viability remain open.

Sources:

1. Ortega-Gomez et al., *A tutorial on the conservation of momentum in photonic time-varying media* (2023). https://arxiv.org/abs/2301.03333
2. Moussa et al., *Observation of Temporal Reflections and Broadband Frequency Translations at Photonic Time-Interfaces* (2023 publication; 2022 preprint). Laboratory analogue only. https://arxiv.org/abs/2208.07236
3. van de Bruck, Mifsud & Nunes, *The variation of the fine-structure constant from disformal couplings* (2015). https://arxiv.org/abs/1510.00200

Reproduce the synthetic checks with Python, NumPy and SciPy using `calculate.py`; numerical results are in `results.json`. No new real-data likelihood or spacecraft residual fit was performed in this step.
