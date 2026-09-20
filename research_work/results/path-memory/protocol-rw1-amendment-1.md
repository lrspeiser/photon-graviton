# RW-1 amendment 1: the shell average's series carries 40 terms
Declared 20 September 2026 after the synthetic gates and before any calculation on data.

## What the gates showed
The protocol's two-branch shell average used the binomial series in x = 2 r r'/(r^2 + r'^2 + w^2) with 20 terms
where x <= 0.2. Gate K5's second part compares the two branches where 0.15 < x < 0.3, which evaluates the series out
to x = 0.3; there the 20-term series misses by 7.1e-8 for p = 2 (its slowest convergence), above the declared 1e-8,
while the closed form agrees with adaptive quadrature to 4e-11. The production kernel (series only where x <= 0.2)
was unaffected: K5's first part, the two-branch kernel against adaptive quadrature at 525 points, passed at 2.1e-8
against 1e-6.

## The amendment
The series carries 40 terms; nothing else changes. With 40 terms the branches agree in the window to 2.5e-10 and
every synthetic gate K1 to K6 passes (K1 2.5e-11, K2 1.2e-10, K3 2.5e-11, K4 5.3e-12, K5 2.1e-8 and 2.5e-10,
K6 1.5e-4). The thresholds are unchanged.
