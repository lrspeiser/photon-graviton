# SE-Rot3: frozen-source vector stiffness and a conservative timestep bound

Declared 20 September 2026. Freeze phi=0, q on the source ring, orbital p=0,
and internal Q,P; vary vector A only. This is a forced linear vector subsystem,
not the complete dynamic-rotor stability problem. Its Hessian is

    L = -discrete_laplacian + omega^2
        + epsilon^2 sum_i D_i^* (|Q_i|^2 I-Q_i Q_i^T) D_i,
    D_i A = sum W_i curl_centered A.

Adjoints use the grid volume inner product, so D_i^* b=curl(W_i b/dV).
Every added quadratic form is nonnegative; L>=omega^2 I. A conservative upper
bound is omega^2+12/h^2 + (3 epsilon^2/h^2) sum_i |Q_i|^2 sum W_i^2/dV.
This follows from the centered curl Fourier norm <=sqrt(3)/h and Cauchy-Schwarz.
For this frozen subsystem only, RK4 on imaginary frequencies is stable if
dt sqrt(lambda_max)<=2 sqrt(2). Report half that bound using the conservative
upper bound, rather than silently adopting it for the complete nonlinear model.

Eight fixtures: n12 and n24 in L8,radius1.2, six radius.7 source positions,
Q_i=sqrt(.001) e_x, omega=.2, epsilon=0,2,10,100. Use scipy eigsh largest
eigenpair with tol1e-9, seeded start20260930. Require relative residual<1e-7,
eigenvalue within upper bound+1e-7 and >=omega^2-1e-7. Check operator symmetry
on a seeded pair, scaled error<1e-10. At epsilon0 verify the largest eigenvalue
against the exact discrete checkerboard eigenvalue omega^2+12/h^2, error<1e-7.

Record that this averaged derivative coupling is nonlocal across each finite
source kernel. The local vacuum principal cone is unchanged, but these frozen
spectral checks do not establish causal propagation through matter, nonlinear
stability, dynamic-rotor torque transfer or a safe full-system timestep. Preserve
all values and failures, with source hashes. Eigenvalue methods, quadratic-form
bounds and RK4 stability are established mathematics; the coupling is a candidate.
