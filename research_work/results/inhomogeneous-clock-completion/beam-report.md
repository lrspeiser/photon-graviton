# Brightness and angular size from the same clock geometry

The retained evolving lapse has a calculable beam-area effect. Its source-to-observer geometric area is generally not 4 pi times coordinate distance squared. Thus the empirical Euclidean brightness comparisons cannot automatically be called predictions of this particular curved-spacetime candidate. This calculation supplies its missing narrow-beam relation; it does not fit a cosmic field or explain companion production.

## Derivation using established ray optics

Keep the known metric form ds^2=-dt^2/n(x,t)^2+dx^2+dy^2+dz^2 and the previously prescribed n=1+beta*t*sin^2(pi*x). Transverse translations are symmetries, so transverse photon momentum p_y is conserved. For the central ray from x=0 to x=D, let p(x) be its local spatial momentum magnitude. To first order in a small angle,

    dy/dx = p_y/p(x)
    I = integral_0^D dx/p(x).

The source opening angle is p_y/p_e. The observer opening angle for the reverse pencil is p_y/p_o. Fixed transverse spatial rods give

    D_G = p_e I,  D_A = p_o I
    S = p_e/p_o = 1+z
    D_G = S D_A.

D_G describes the observer-screen area per solid angle emitted at the source; D_A describes physical source size per angle seen by the observer. Their roles cannot be interchanged. The plane symmetry gives the same transverse factor in both screen directions; these are infinitesimal pencils around the central ray, not all-sky distances for a general galaxy distribution.

Conserved photons lose a factor S in photon energy and a factor S in arrival rate, so

    F = L/[4 pi D_G^2 S^2]
    D_L = S D_G = S^2 D_A.

This is the known Etherington distance-duality relation, recovered from the explicit ray map. It is not unique to expansion or to our hypothesis. Its geometric assumptions and relation to photon conservation are discussed by [More et al.](https://arxiv.org/abs/1612.08784). An interacting non-geodesic conversion theory needs its own optical transport and cannot automatically inherit this derivation.

The earlier flux expression L/[4 pi D^2 S^2] remains a conditional Euclidean-area benchmark. For this metric it applies only if the actual D_G equals D. The earlier specification already required a derived beam area in curved geometry; the present calculation evaluates it instead of assuming equality.

## Same redshift, different beam areas

For one normalized cell D=1, endpoint clocks have n=1 and the previously derived stretch is S=exp(beta/2), independent of emission time t_e. Normalize p_e=1. The existing central-ray equations give

    dt/dx=n
    d ln(omega)/dx=-partial_t n
    p(x)=n(x,t(x))*omega(x)
    D_G=integral dx/[n omega].

The last integral depends on when the ray traverses the field even though S does not. The following uses selected dimensionless field parameters, not measured cosmic epochs or a universal field fit:

| Redshift | t_e | D_G / D | D_A / D | Flux relative to the Euclidean S^-2 benchmark | Magnitude change |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 0 | 1.000919 | 0.909926 | 0.998165 | +0.001994 |
| 0.1 | 1 | 0.920649 | 0.836954 | 1.179809 | -0.179529 |
| 1.0 | 0 | 1.044066 | 0.522033 | 0.917368 | +0.093640 |
| 1.0 | 1 | 0.731704 | 0.365852 | 1.867796 | -0.678324 |

At the same z=1, the later crossing is brighter than the Euclidean benchmark while the earlier crossing is slightly dimmer. A choice of field epoch is therefore part of the prediction, not a harmless coordinate relabeling when n is kept in this fixed functional form. Changing the origin consistently would also change the profile's initial state. We have not imposed a universe age or fitted a separate epoch per object.

These rows vary beta to set each diagnostic redshift; they are not one globally specified sky realization. They cannot be applied as a universal correction to the 960 supernovae, nor interpreted as a demonstrated solution to their residuals.

## Direct check and physical consequence

beam.py integrates the central ray and I, then independently propagates slightly tilted rays forward and backward using the full Hamilton equations with conserved transverse momentum. The opening angle is 1e-5 radians. Those screen displacements reproduce D_G and D_A within 1e-7 relative error, including the flat-space control. Reciprocity follows exactly from the defined distances and is not a separate observational validation.

The same optical geometry now predicts brightness and apparent transverse size. This connects the first and third objectives: focusing cannot be adjusted for brightness while leaving lensing or angular-size predictions untouched. However, this plane-symmetric lapse is not a galaxy/cluster lens model, has no demonstrated dynamical source and still produces no separate companion current without an additional interaction.

**Decision:** retain the Euclidean observational fits as conditional benchmarks, not as full tests of the lapse candidate. A future metric candidate must specify one spatial/temporal field and propagate beam areas as well as frequencies and arrival intervals through it. A redshift relation alone does not determine all three observables for this family. No final holdout, new observational fit, or derived capture law was added; all six objectives remain open.
