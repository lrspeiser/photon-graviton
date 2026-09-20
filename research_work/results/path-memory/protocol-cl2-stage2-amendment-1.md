# CL-2 stage 2, amendment 1: gate G-N3's comparison target
Declared 19 September 2026 before the archived run, after the numerical development of the disk-force
routine showed that gate G-N3 cannot pass as written for any correct code.

G-N3 as declared: "a compact ring set reproduces the point-mass force to 1e-6 at ten times its radius".
A thin ring of radius a is not a point mass in its own plane: its in-plane force at R > a is
G M/R^2 [1 + (3/4)(a/R)^2 + ...] (established: the multipole expansion of the ring potential), so at
ten radii the difference is 7.5e-3, and the development run measured 7.6e-3 for the uniform annulus
0.9 to 1.1 kpc at 10 kpc. The gate compares with the wrong reference, not with a wrong force.

The amendment keeps G-N3 in the archive as declared, where it fails, and adds G-N3b: the same annulus
against the exact annulus force, the ring force integrated over the annulus by adaptive quadrature,
to 1e-6 at 10 kpc; and G-N3c: the point-mass comparison at one thousand radii, where the ring's
quadrupole is 7.5e-7, to 1e-5. Nothing else in the protocol changes; the development values are not
the archive.
