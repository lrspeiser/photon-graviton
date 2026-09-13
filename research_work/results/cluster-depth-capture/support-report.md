# Can particles sustain the depth-capture profile?

Status: conditional dynamical support test, not an observational fit. The previous goal turn made progress by projecting this profile into lensing and motion predictions. This turn tests whether one proposed form of stored companions can sustain that profile. All six objectives remain open.

## Question and assumptions

The injection profile says where incoming energy is captured. It does not establish where that energy remains. Here assume captured companions become nonrelativistic collisionless bound particles, with a spherical, steady distribution and no preferred velocity direction. Their density is required to retain the shape q(r) calculated from kappa=chi*(-Phi)^4. Ordinary matter supplies a fixed Plummer potential; deposit self-gravity is neglected. These are optional storage assumptions, not locked-in companion properties. The test does not address trapped waves, field stress, interacting matter, or an evolving population.

Known Eddington inversion determines the unique isotropic distribution f(E) for a specified density and potential under these boundary conditions. A negative f represents a negative particle population, so it invalidates that pairing of assumptions. This use of the method is standard; see [Sanchez Almeida, Plastino and Trujillo (2024)](https://arxiv.org/abs/2407.16519). We use the mathematical consistency criterion, not that paper's dark-matter potentials or cosmological assumptions.

In units GM=a=1, relative potential Psi=1/sqrt(1+r^2), binding energy E=Psi-v^2/2 and normalized density rho=q/A:

    rho(Psi)=4 pi sqrt(2) integral_0^Psi f(E) sqrt(Psi-E) dE,
    f(E)=1/(sqrt(8) pi^2) integral_0^E rho''(Psi)/sqrt(E-Psi) dPsi.

The boundary term vanishes because rho'(0)=0. These are known formulas. The optional capture law supplies rho, which is our model-specific input. Multiplying rho by a positive mass normalization cannot remove a negative f in a fixed potential.

## Results

| A | f at maximum binding E=1 | Negative values in 1001 sampled energies? |
|---:|---:|---|
| 0.1 | 0.393239 | none found |
| 1 | 0.037199 | none found |
| 1.5 | -0.041817 | yes |
| 1.833275 | -0.067856 | yes |
| 3 | -0.076005 | yes |
| 10 | 0.006069 | yes, at intermediate energies |

The sampled range is E=0.01 to 1. Positive samples alone are not proof of positivity at every energy, dynamical stability, or a formation mechanism. In the analytic optically thin limit, rho=Psi^4 and f=12 B(3,1/2) E^(5/2)/(sqrt(8) pi^2), which is positive throughout 0<E<=1. Here B denotes the beta function, not the surface-curvature coefficient in the lensing report.

The high-binding endpoint crosses zero near A=1.19134954 on the bracket tested. This is a diagnostic crossing, not a proven universal support boundary: at A=10 the endpoint is positive again, but the distribution remains negative elsewhere. The A=1.5 failure occurs before the 3D central-curvature transition at A=1.83327 and the projected transition at A=2.74533. Thus an acceptable-looking central density and lensing profile are insufficient to establish particle support.

## Verification

support.py preserves all results in support-results.json. Density splines use the analytic endpoint derivative, rho'(1)=-2 exp(-pi A/4)[-2+pi A/4+A^2/6]. The inversion removes the integrable square-root singularity by substituting Psi=E(1-t^2).

Doubling the potential grid from 2048 to 4096 intervals, angular quadrature from 256 to 512, and inversion quadrature from 256 to 512 changes sampled f by at most 2.8e-7 relative to its case's maximum absolute value. The endpoint crossing changes by 2.5e-8 in A. Reintegrating f, including all negative values without clipping, recovers the target density at five potentials per case to a relative discrepancy below 4e-7. An independent analytic thin-limit distribution agrees to 5.5e-6 at A=1e-6 over the sampled energies.

These checks support the numerical signs. Reconstruction of a signed distribution is an inversion check, not evidence that a negative distribution is physical.

## Consequence for the model

Retain weak capture as an isotropic particle-storage candidate. Do not fit the stronger frozen injection profiles as isotropic particle halos in this fixed well. To use stronger capture, calculate redistribution after capture, a velocity distribution with directional preferences, field or interaction support, or self-consistent gravity before interpreting the profile as a persistent halo. None of these alternatives is established by this test. The next physical calculation must connect capture to storage rather than independently choosing a favorable gravitational profile. No observational holdouts were opened.
